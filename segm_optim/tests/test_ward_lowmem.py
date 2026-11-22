#### eHabitat+ HFT model: test of the classification with little memory (ward_lowmem.py)
#### License: CC BY-SA 3.0
"""
Checks ward_lowmem.py, which classifies the segments into HFTs when a threshold has more than LOWMEM_SEGMENTS
segments, against the classification with the matrix of the distances:

- the weighted Ward of the nearest-neighbour chain gives the merges and heights of the usual algorithm (an O(n³)
  brute force with weights, and scipy without them), as an hclust object with the conventions of R (merge, increasing
  heights, order), also with ties and with a merge that rounding left below one of its groups;
- the cophenetic correlation and the silhouette widths (with the neighbour HFT) are those computed with the matrices;
- the command line (the files that moranvar_plots.R writes and reads) and its errors;
- with R: class_tree() of moranvar_plots.R gives the same tree (merge, heights, order, labels, cuts of the tree),
  cophenetic correlation and silhouettes (sil_widths) with and without little memory, with CLASS_METHOD 'ward_area',
  'ward' and the classifications of 'ward_area_kmeans', switches above the limit, stops if Python fails, and its
  LOWMEM_SEGMENTS is that of hft_state.py.

    cd segm_optim && python -s tests/test_ward_lowmem.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEGM = os.path.dirname(HERE)
sys.path.insert(0, SEGM)  # segm_optim: ward_lowmem.py, hft_state.py

import hft_state  # noqa: E402
import ward_lowmem as wl  # noqa: E402


def tree(x, w):
	mi, mj, mh = wl.nn_chain(np.ascontiguousarray(x, dtype=np.float64), np.asarray(w, dtype=np.float64))
	return wl.to_hclust(mi, mj, mh, len(w))


def groups_of_merge(merge, n):
	"""Rows of each group of an hclust merge matrix."""
	g = []
	for a, b in merge:
		g.append((g[a - 1] if a > 0 else [-a - 1]) + (g[b - 1] if b > 0 else [-b - 1]))
	return g


def cut(merge, n, k):
	"""Group 0..k-1 of each row after the first n - k merges (as cutree(k=)), numbered by their first row."""
	lab = np.arange(n)
	for rows in groups_of_merge(merge[:n - k], n):
		lab[rows] = min(rows)
	return np.unique(lab, return_inverse=True)[1]


def canon(lab):
	"""A partition numbered by the order of its first row, to compare partitions."""
	_, first, inv = np.unique(lab, return_index=True, return_inverse=True)
	rank = np.argsort(np.argsort(first))
	return rank[inv]


def ward_brute(x, w):
	"""Weighted Ward by brute force: merges in order (rows of the two groups) and heights
	sqrt(2 w_a w_b / (w_a + w_b) ||c_a - c_b||²)."""
	groups = [([i], x[i].astype(float), float(w[i])) for i in range(len(w))]
	merges, heights = [], []
	while len(groups) > 1:
		best = None
		for a in range(len(groups)):
			for b in range(a + 1, len(groups)):
				ra, ca, wa = groups[a]
				rb, cb, wb = groups[b]
				d = 2 * wa * wb / (wa + wb) * np.sum((ca - cb) ** 2)
				if best is None or d < best[0]:
					best = (d, a, b)
		d, a, b = best
		ra, ca, wa = groups[a]
		rb, cb, wb = groups[b]
		merges.append((ra, rb))
		heights.append(np.sqrt(d))
		groups[a] = (ra + rb, (wa * ca + wb * cb) / (wa + wb), wa + wb)
		del groups[b]
	return merges, np.array(heights)


def check_hclust(merge, height, order, n):
	"""An hclust object that R accepts: each row once, each group used once after it is made, increasing heights,
	the conventions of hcass2 and the order of the leaves of every group together, left before right."""
	assert merge.shape == (n - 1, 2) and height.shape == (n - 1,) and order.shape == (n,)
	assert np.all(np.diff(height) >= 0)
	used = [-int(v) for v in merge.ravel() if v < 0]
	assert sorted(used) == list(range(1, n + 1))
	steps = [int(v) for v in merge.ravel() if v > 0]
	assert sorted(steps) == list(range(1, n - 1))
	for t, (a, b) in enumerate(merge, start=1):
		assert a < t and b < t
		if a < 0 and b < 0:
			assert -a < -b  # two rows: the lower first
		elif a > 0 or b > 0:
			assert a < 0 or (b > 0 and a < b)  # a row before a group; two groups: the earlier first
	assert sorted(order.tolist()) == list(range(1, n + 1))
	pos = {int(v): i for i, v in enumerate(order)}
	for t, rows in enumerate(groups_of_merge(merge, n)):
		p = sorted(pos[r + 1] for r in rows)
		assert p[-1] - p[0] == len(p) - 1  # together in the dendrogram
		a = merge[t, 0]
		left = groups_of_merge(merge, n)[a - 1] if a > 0 else [-a - 1]
		assert max(pos[r + 1] for r in left) < min(pos[r + 1] for r in rows if r not in left)


def test_tree():
	rng = np.random.default_rng(1)
	# with weights: the merges and heights of the brute force
	for n, p in ((2, 1), (3, 2), (25, 3), (60, 5)):
		x = rng.normal(size=(n, p))
		w = rng.integers(1, 50, size=n).astype(float)
		merge, height, order = tree(x, w)
		check_hclust(merge, height, order, n)
		bm, bh = ward_brute(x, w)
		assert np.allclose(height, bh, rtol=1e-12, atol=0)
		got = [sorted(g) for g in groups_of_merge(merge, n)]
		assert got == [sorted(a + b) for a, b in bm]
	# without weights: those of scipy (hclust ward.D2 with members = 1)
	from scipy.cluster.hierarchy import fcluster, linkage
	x = rng.normal(size=(400, 6))
	merge, height, order = tree(x, np.ones(400))
	check_hclust(merge, height, order, 400)
	z = linkage(x, 'ward')
	assert np.allclose(height, z[:, 2], rtol=1e-12, atol=0)
	for k in (2, 3, 7, 30, 399):
		assert np.array_equal(canon(cut(merge, 400, k)), canon(fcluster(z, k, 'maxclust')))
	print('test_tree: OK')


def test_ties():
	"""Many equal distances (a grid with repeated points): still a valid tree with the heights of the brute force."""
	x = np.array([[i % 4, i // 4 % 3] for i in range(36)], dtype=float)
	w = np.ones(36)
	merge, height, order = tree(x, w)
	check_hclust(merge, height, order, 36)
	assert np.allclose(height, np.sort(ward_brute(x, w)[1]), rtol=1e-12, atol=1e-12)
	assert np.sum(height == 0) == 24  # the repeated points first
	print('test_ties: OK')


def test_to_hclust():
	"""Conventions of R and a merge left below one of its groups by rounding."""
	# made: (3, 4) at 1.0, (0, 2) at 2.0, (0, 3) at 1.0 - 1e-15 (below (3, 4): takes 1.0, after it), (0, 1) at 5.0
	mi = np.array([3, 0, 0, 0])
	mj = np.array([4, 2, 3, 1])
	mh = np.array([1.0, 2.0, 1.0 - 1e-15, 5.0])
	merge, height, order = wl.to_hclust(mi, mj, mh, 5)
	check_hclust(merge, height, order, 5)
	assert merge.tolist() == [[-4, -5], [-1, -3], [1, 2], [-2, 3]]
	assert height.tolist() == [1.0, 2.0, 2.0, 5.0]
	assert order.tolist() == [2, 4, 5, 1, 3]
	print('test_to_hclust: OK')


def cophenetic_brute(merge, height, n):
	c = np.zeros((n, n))
	for t, (a, b) in enumerate(merge):
		ga = groups_of_merge(merge, n)[a - 1] if a > 0 else [-a - 1]
		gb = groups_of_merge(merge, n)[b - 1] if b > 0 else [-b - 1]
		c[np.ix_(ga, gb)] = height[t]
		c[np.ix_(gb, ga)] = height[t]
	return c


def test_cophenetic():
	from scipy.spatial.distance import pdist, squareform
	rng = np.random.default_rng(2)
	for n, p in ((3, 2), (7, 1), (120, 4), (301, 3)):  # odd and even n (rows paired in _pair_sums)
		x = rng.normal(size=(n, p)) + rng.integers(0, 3, size=(n, 1)) * 3
		w = rng.integers(1, 30, size=n).astype(float)
		merge, height, order = tree(x, w)
		iu = np.triu_indices(n, 1)
		ref = np.corrcoef(cophenetic_brute(merge, height, n)[iu], squareform(pdist(x))[iu])[0, 1]
		assert abs(wl.cophenetic_cor(x, merge, height, order) - ref) < 1e-12, n
	assert np.isnan(wl.cophenetic_cor(np.zeros((2, 1)), np.array([[-1, -2]]), np.array([0.0]), np.array([1, 2])))
	print('test_cophenetic: OK')


def silhouette_brute(x, lab):
	"""sil_width of moranvar_plots.R: s and neighbour of each row."""
	d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(-1))
	k = lab.max() + 1
	nh = np.bincount(lab, minlength=k).astype(float)
	md = np.stack([d[:, lab == g].sum(1) for g in range(k)], 1) / nh
	own = md[np.arange(len(lab)), lab]
	a = own * nh[lab] / np.maximum(nh[lab] - 1, 1)
	md[np.arange(len(lab)), lab] = np.inf
	nb = np.argmin(md, 1)  # the first in a tie
	b = md[np.arange(len(lab)), nb]
	mx = np.maximum(a, b)
	s = np.where(mx > 0, (b - a) / np.where(mx > 0, mx, 1), 0)
	s[nh[lab] == 1] = 0
	return s, nb


def test_silhouettes():
	rng = np.random.default_rng(3)
	n = 157
	x = rng.normal(size=(n, 3))
	labs = [rng.integers(0, 2, n), rng.integers(0, 9, n), np.r_[np.zeros(n - 1, int), 1]]  # the last: a row alone
	x[5] = x[6]  # two rows at distance 0
	labs.append(np.r_[np.arange(2), np.zeros(n - 2, int)])
	g = np.asfortranarray(np.stack(labs, 1).astype(np.int32))
	s, nb = wl.silhouettes(x, g, (g.max(0) + 1).astype(np.int64))
	for c, lab in enumerate(labs):
		rs, rnb = silhouette_brute(x, lab)
		assert np.allclose(s[:, c], rs, rtol=0, atol=1e-13), c
		assert np.array_equal(nb[:, c], rnb), c
	# equal mean distances to two HFTs: the first one is the neighbour
	x = np.array([[0.0], [-1.0], [1.0], [0.0]])
	s, nb = wl.silhouettes(x, np.asfortranarray(np.array([[0], [1], [2], [0]], np.int32)), np.array([3]))
	assert nb[0, 0] == 1 and s[0, 0] == 1.0
	print('test_silhouettes: OK')


def run_cli(args):
	return subprocess.run([sys.executable, '-s', os.path.join(SEGM, 'ward_lowmem.py')] + [str(a) for a in args],
	                      capture_output=True, text=True)


def test_cli():
	"""The files of moranvar_plots.R: matrices by columns, little-endian."""
	rng = np.random.default_rng(4)
	n, p = 50, 3
	x = rng.normal(size=(n, p))
	w = rng.integers(1, 9, size=n).astype(float)
	with tempfile.TemporaryDirectory() as d:
		x.ravel(order='F').astype('<f8').tofile(os.path.join(d, 'x.f8'))
		w.astype('<f8').tofile(os.path.join(d, 'w.f8'))
		r = run_cli(['tree', d, n, p])
		assert r.returncode == 0, r.stderr
		merge = np.fromfile(os.path.join(d, 'merge.i4'), '<i4').reshape((n - 1, 2), order='F')
		height = np.fromfile(os.path.join(d, 'height.f8'), '<f8')
		order = np.fromfile(os.path.join(d, 'order.i4'), '<i4')
		em, eh, eo = tree(x, w)
		assert np.array_equal(merge, em) and np.array_equal(height, eh) and np.array_equal(order, eo)
		assert np.fromfile(os.path.join(d, 'coph.f8'), '<f8')[0] == wl.cophenetic_cor(x, em, eh, eo)
		g = np.stack([cut(em, n, 2), cut(em, n, 5)], 1).astype('<i4')
		g.ravel(order='F').tofile(os.path.join(d, 'groups.i4'))
		r = run_cli(['silhouette', d, n, p, 2])
		assert r.returncode == 0, r.stderr
		s = np.fromfile(os.path.join(d, 's.f8'), '<f8').reshape((n, 2), order='F')
		nb = np.fromfile(os.path.join(d, 'neighbour.i4'), '<i4').reshape((n, 2), order='F')
		es, enb = wl.silhouettes(x, np.asfortranarray(g), np.array([2, 5]))
		assert np.array_equal(s, es) and np.array_equal(nb, enb)
		# errors: values missing, NA, a segment without pixels, a classification with one group
		assert run_cli(['tree', d, n + 1, p]).returncode != 0
		w[3] = 0
		w.astype('<f8').tofile(os.path.join(d, 'w.f8'))
		assert 'each with pixels' in run_cli(['tree', d, n, p]).stderr
		x[2, 1] = np.nan
		x.ravel(order='F').astype('<f8').tofile(os.path.join(d, 'x.f8'))
		assert 'NA' in run_cli(['tree', d, n, p]).stderr
		x[2, 1] = 0
		x.ravel(order='F').astype('<f8').tofile(os.path.join(d, 'x.f8'))
		np.zeros(n, '<i4').tofile(os.path.join(d, 'groups.i4'))
		assert 'k >= 2' in run_cli(['silhouette', d, n, p, 1]).stderr
		assert run_cli(['other', d, n, p]).returncode != 0
	print('test_cli: OK')


def test_limit():
	"""The same limit in moranvar_plots.R and hft_state.py, ward_lowmem.py in the code of the plots step, and the
	memory with the matrices is that of the limit at most."""
	with open(os.path.join(SEGM, 'moranvar_plots.R'), encoding='utf-8') as f:
		m = re.search(r'^LOWMEM_SEGMENTS <- (\d+)$', f.read(), re.M)
	assert m and int(m.group(1)) == hft_state.LOWMEM_SEGMENTS
	assert 'ward_lowmem.py' in hft_state.CODE['plots']
	assert abs(hft_state.class_memory_gb(6000) - 2.52) < 0.01
	assert hft_state.class_memory_gb(55000) == hft_state.class_memory_gb(hft_state.LOWMEM_SEGMENTS) < 4
	print('test_limit: OK')


R_TEST = r"""
suppressMessages(library(vegan))
source('plot_utils.R')
for (e in parse('moranvar_plots.R', keep.source=FALSE)) {
	if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
	    ((is.call(e[[3]]) && identical(e[[3]][[1]], as.name('function'))) || identical(as.character(e[[2]]), 'LOWMEM_SEGMENTS')))
		eval(e, globalenv())
}
set.seed(7)
ns <- 400
nv <- 4
m <- matrix(rnorm(ns * nv), ns) + sample(0:3, ns, TRUE) * 2
v <- matrix(rexp(ns * nv, 5), ns)
hri <- data.frame(1, sprintf('1%06d', seq_len(ns)), m, v, sample(5:300, ns, TRUE), v)
names(hri) <- c('ecoregion', 'segm_id', paste0('v', 1:nv, 'pamean'), paste0('v', 1:nv, 'pavar'), 'sumpamask', paste0('v', 1:nv, 'pavar2'))
fail <- character(0)
check <- function(ok, what) if (!isTRUE(ok)) fail <<- c(fail, what)
for (method in c('ward_area', 'ward')) {
	a <- class_tree(hri, method, lowmem=FALSE)
	b <- class_tree(hri, method, lowmem=TRUE)
	check(!a$lowmem && b$lowmem && is.null(b$dist), paste(method, 'lowmem'))
	check(identical(a$tree$merge, b$tree$merge), paste(method, 'merge'))
	check(identical(a$tree$order, b$tree$order), paste(method, 'order'))
	check(isTRUE(all.equal(a$tree$height, b$tree$height, tolerance=1e-12)), paste(method, 'heights'))
	check(identical(a$tree$labels, b$tree$labels), paste(method, 'labels'))
	cuts <- function(tr) lapply((1:9) / 10, function(f) cutree(tr, h=f * max(tr$height)))
	cls <- cuts(a$tree)
	check(identical(cls, cuts(b$tree)), paste(method, 'cuts'))
	check(identical(lapply(2:30, function(k) cutree(a$tree, k=k)), lapply(2:30, function(k) cutree(b$tree, k=k))), paste(method, 'k'))
	check(abs(cor(cophenetic(a$tree), a$dist) - b$coph) < 1e-10, paste(method, 'cophenetic correlation'))
	# also one HFT, one HFT per segment (NA) and those of 'ward_area_kmeans' (not cuts of the tree)
	cls <- c(cls, list(rep(1L, ns), seq_len(ns)), lapply(2:6, function(k) classify_k(a, k, 'ward_area_kmeans')))
	sa <- sil_widths(a, as.matrix(dist(a$z)), cls)
	sb <- sil_widths(b, NULL, cls)
	for (i in seq_along(cls)) {
		check(identical(is.na(sa[[i]]$s), is.na(sb[[i]]$s)), paste(method, 'silhouette NA', i))
		check(max(c(0, abs(sa[[i]]$s - sb[[i]]$s)), na.rm=TRUE) < 1e-12, paste(method, 'silhouette', i))
		check(identical(sa[[i]]$neighbour, sb[[i]]$neighbour), paste(method, 'neighbour', i))
	}
	pdf(NULL)
	check(!inherits(try(plot(b$tree, hang=-1)), 'try-error') && !inherits(try(rect.hclust(b$tree, k=3)), 'try-error'), paste(method, 'plot'))
	dev.off()
}
Sys.setenv(EHAB_LOWMEM_SEGMENTS=as.character(ns - 1))
check(class_tree(hri, 'ward_area')$lowmem, 'little memory above the limit')
Sys.setenv(EHAB_LOWMEM_SEGMENTS=as.character(ns))
check(!class_tree(hri, 'ward_area')$lowmem, 'the matrices up to the limit')
Sys.unsetenv('EHAB_LOWMEM_SEGMENTS')
check(!class_tree(hri, 'ward_area')$lowmem && LOWMEM_SEGMENTS > ns, 'default limit')
Sys.setenv(EHAB_PYTHON=file.path(tempdir(), 'no_python'))
check(inherits(try(class_tree(hri, 'ward_area', lowmem=TRUE), silent=TRUE), 'try-error'), 'an error if Python fails')
if (length(fail)) {
	cat('FAILED:', fail, sep='\n')
	quit(status=1)
}
cat('R OK\n')
"""


def test_r():
	# the R of the same environment (that of conda next to its Python), or the one in the PATH
	rscript = shutil.which('Rscript', path=os.path.dirname(sys.executable)) or shutil.which('Rscript')
	if not rscript:
		print('test_r: skipped (Rscript not found)')
		return
	with tempfile.NamedTemporaryFile('w', suffix='.R', delete=False) as f:
		f.write(R_TEST)
	try:
		env = dict(os.environ, EHAB_PYTHON=sys.executable, PYTHONNOUSERSITE='1')
		env.pop('EHAB_LOWMEM_SEGMENTS', None)
		r = subprocess.run([rscript, f.name], cwd=SEGM, env=env, capture_output=True, text=True)
	finally:
		os.unlink(f.name)
	assert r.returncode == 0 and 'R OK' in r.stdout, r.stdout + r.stderr
	print('test_r: OK')


if __name__ == '__main__':
	test_tree()
	test_ties()
	test_to_hclust()
	test_cophenetic()
	test_silhouettes()
	test_cli()
	test_limit()
	test_r()
	print('All the tests passed')
