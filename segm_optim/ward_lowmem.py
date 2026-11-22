#### eHabitat+ HFT model: classification of many segments with little memory (Ward, cophenetic correlation, silhouettes)
#### License: CC BY-SA 3.0
"""
The classification of the segments into HFTs (moranvar_plots.R) builds with hclust the matrix of the distances between
every pair of segments, and so do the cophenetic correlation and the silhouette widths: about 70 bytes x segments²
(hft_state.CLASS_BYTES), 200 GB with 55,000 segments. Above hft_state.LOWMEM_SEGMENTS segments moranvar_plots.R
calls this script instead, which computes the same results with memory proportional to the number of segments, the
distances computed when needed (time still grows with segments², in parallel with numba):

- tree: Ward weighted by the number of pixels of each segment with the nearest-neighbour chain (Murtagh 1983,
  https://doi.org/10.1093/comjnl/26.4.354; Mullner 2011, https://arxiv.org/abs/1109.2378), which merges reciprocal nearest
  neighbours and gives the tree of the usual algorithm (Ward is reducible). The distance between two groups is that of
  hclust(method = 'ward.D2', members = w) on the distances sqrt(2 w_i w_j / (w_i + w_j)) ||x_i - x_j|| between the rows
  of x, computed from the centroids of the groups, so the merges and heights are those of class_tree() (up to
  rounding). Returned as the merge matrix, heights and order of an hclust object, and with the cophenetic correlation:
  the Pearson correlation between the cophenetic distances and the Euclidean distances between the rows of x, over
  every pair (cor(cophenetic(tree), dist(x)) of moranvar_plots.R).
- silhouette: the silhouette widths of the segments in several classifications at once, and the neighbour HFT of each
  segment (sil_width() of moranvar_plots.R), with the Euclidean distances between the rows of x.

moranvar_plots.R writes the inputs to a folder DIR and reads the outputs from it (little-endian binary files, matrices
by columns as R stores them):

    python ward_lowmem.py tree DIR N P          # x.f8 (N x P), w.f8 (N) -> merge.i4 (N-1 x 2), height.f8, order.i4, coph.f8
    python ward_lowmem.py silhouette DIR N P C  # x.f8 (N x P), groups.i4 (N x C, 0..k-1) -> s.f8 (N x C), neighbour.i4 (N x C)
"""

import os
import sys

import numpy as np
from numba import njit, prange


@njit(parallel=True, cache=True)
def nn_chain(x, w):
	"""Weighted Ward of the rows of X (N x P) with weights W (N) by the nearest-neighbour chain. Returns the N - 1
	merges in the order they are made: the two groups (each known by its lowest row; the merged group keeps the
	lowest) and the height sqrt(2 w_a w_b / (w_a + w_b) ||c_a - c_b||²), c: the weighted mean of the rows of the
	group. Ties: the first group in the list of active groups, but the previous group of the chain when it is as near
	(so that the chain ends in reciprocal neighbours)."""
	n, p = x.shape
	c = x.copy()
	wt = w.copy()
	act = np.arange(n)  # active groups
	pos = np.arange(n)  # position of each group in act
	nact = n
	chain = np.empty(n, np.int64)
	clen = 0
	mi = np.empty(n - 1, np.int64)
	mj = np.empty(n - 1, np.int64)
	mh = np.empty(n - 1, np.float64)
	buf = np.empty(n, np.float64)
	m = 0
	while nact > 1:
		if clen == 0:
			chain[0] = act[0]
			clen = 1
		a = chain[clen - 1]
		wa = wt[a]
		for t in prange(nact):
			b = act[t]
			s = 0.0
			for k in range(p):
				df = c[b, k] - c[a, k]
				s += df * df
			buf[t] = wa * wt[b] / (wa + wt[b]) * s
		best = np.inf
		bi = -1
		for t in range(nact):
			b = act[t]
			if b != a and buf[t] < best:
				best = buf[t]
				bi = b
		if clen >= 2:
			prev = chain[clen - 2]
			if buf[pos[prev]] <= best:
				bi = prev
				best = buf[pos[prev]]
		if clen >= 2 and bi == chain[clen - 2]:  # reciprocal nearest neighbours: merged
			i = min(a, bi)
			j = max(a, bi)
			mi[m] = i
			mj[m] = j
			mh[m] = np.sqrt(2.0 * best)
			m += 1
			wi = wt[i]
			wj = wt[j]
			for k in range(p):
				c[i, k] = (wi * c[i, k] + wj * c[j, k]) / (wi + wj)
			wt[i] = wi + wj
			pj = pos[j]
			last = act[nact - 1]
			act[pj] = last
			pos[last] = pj
			nact -= 1
			clen -= 2
		else:
			chain[clen] = bi
			clen += 1
	return mi, mj, mh


@njit(cache=True)
def to_hclust(mi, mj, mh, n):
	"""The merges of nn_chain as an hclust object of R: merge (N-1 x 2: -j is the row j, from 1, and a positive number
	the group made in that step), height (increasing) and order (the rows from left to right in the dendrogram), with
	the conventions of hclust (hcass2): two rows, the lower first; a row and a group, the row first; two groups, the
	earlier first; the order lists the left branch of each merge before the right one. The merges are sorted by
	height, those of the same height in the order they were made. A merge whose height came out lower than that of
	one of its groups by rounding takes the height of that group, so that heights increase along every branch."""
	h = mh.copy()
	step = np.full(n, -1, np.int64)  # merge (in the order made) that formed the group of each row that keeps a group
	for s in range(n - 1):
		for g in (mi[s], mj[s]):
			if step[g] >= 0 and h[step[g]] > h[s]:
				h[s] = h[step[g]]
		step[mi[s]] = s
		step[mj[s]] = -1
	o = np.argsort(h, kind='mergesort')  # stable: the merges of a group keep their order, so each one after its groups
	node = np.empty(n, np.int64)  # the group of each row that keeps a group, as in merge
	for g in range(n):
		node[g] = -(g + 1)
	merge = np.empty((n - 1, 2), np.int32)
	height = np.empty(n - 1, np.float64)
	for t in range(n - 1):
		s = o[t]
		a = node[mi[s]]
		b = node[mj[s]]
		if a > 0 and b < 0:
			a, b = b, a
		elif a > 0 and b > 0 and a > b:
			a, b = b, a
		merge[t, 0] = a
		merge[t, 1] = b
		height[t] = h[s]
		node[mi[s]] = t + 1
	# order: the leaves of the root from left to right (depth first, left branch first)
	order = np.empty(n, np.int32)
	stack = np.empty(n, np.int64)
	top = 0
	stack[0] = n - 1  # the last merge (1-based), the root
	top = 1
	k = 0
	while top > 0:
		top -= 1
		v = stack[top]
		if v < 0:
			order[k] = -v
			k += 1
		else:
			stack[top] = merge[v - 1, 1]
			stack[top + 1] = merge[v - 1, 0]
			top += 2
	return merge, height, order


@njit(cache=True)
def _row_pairs(xo, gap, p, mc, md):
	"""Sums over the pairs (p, q > p) of the rows of XO in the order of the dendrogram: of the distance d, of
	(c - MC)(d - MD) and of (d - MD)², c being the cophenetic distance, the highest GAP between p and q."""
	n, nv = xo.shape
	sd = 0.0
	scd = 0.0
	sdd = 0.0
	cmax = -np.inf
	for q in range(p + 1, n):
		if gap[q - 1] > cmax:
			cmax = gap[q - 1]
		s = 0.0
		for k in range(nv):
			df = xo[q, k] - xo[p, k]
			s += df * df
		d = np.sqrt(s)
		sd += d
		scd += (cmax - mc) * (d - md)
		sdd += (d - md) * (d - md)
	return sd, scd, sdd


@njit(parallel=True, cache=True)
def _pair_sums(xo, gap, mc, md):
	n = xo.shape[0]
	sd = 0.0
	scd = 0.0
	sdd = 0.0
	for t in prange((n + 1) // 2):  # rows t and n - 1 - t together: the same work in every iteration
		a, b, c = _row_pairs(xo, gap, t, mc, md)
		sd += a
		scd += b
		sdd += c
		u = n - 1 - t
		if u != t:
			a, b, c = _row_pairs(xo, gap, u, mc, md)
			sd += a
			scd += b
			sdd += c
	return sd, scd, sdd


@njit(cache=True)
def _sizes_gaps(merge, height, order, n):
	"""Rows of each group (size), and the cophenetic distance between the neighbouring leaves of the dendrogram: the
	height of the merge that joins them (gap[t], between the leaves t and t + 1 of order)."""
	size = np.empty(n - 1, np.int64)
	first = np.empty(n - 1, np.int64)  # position in order of the first row of each group
	pos = np.empty(n + 1, np.int64)
	for t in range(n):
		pos[order[t]] = t
	gap = np.empty(n - 1, np.float64)
	for t in range(n - 1):
		sa = 1
		fa = 0
		sb = 1
		a = merge[t, 0]
		b = merge[t, 1]
		if a < 0:
			fa = pos[-a]
		else:
			sa = size[a - 1]
			fa = first[a - 1]
		if b < 0:
			fb = pos[-b]
		else:
			sb = size[b - 1]
			fb = first[b - 1]
		size[t] = sa + sb
		first[t] = min(fa, fb)
		# the left group (a) ends just before the right one (b) in order
		gap[fa + sa - 1] = height[t]
	return size, gap


def cophenetic_cor(x, merge, height, order):
	"""Pearson correlation between the cophenetic distances of the tree and the Euclidean distances between the rows
	of X, over every pair of rows (cor(cophenetic(tree), dist(x)) in R), without the matrix of the distances: the
	cophenetic distance between two leaves is the highest merge between them in the order of the dendrogram."""
	n = x.shape[0]
	if n < 3:
		return np.nan
	size, gap = _sizes_gaps(merge, height, order, n)
	left = np.where(merge[:, 0] < 0, 1, size[np.maximum(merge[:, 0], 1) - 1])
	pairs = left * (size - left).astype(np.float64)  # pairs of rows whose cophenetic distance is each height
	npairs = n * (n - 1) / 2
	mc = float(np.sum(pairs * height) / npairs)
	scc = float(np.sum(pairs * (height - mc) ** 2))
	xo = np.ascontiguousarray(x[order - 1])
	sd = _pair_sums(xo, gap, mc, 0.0)[0]
	md = sd / npairs
	_, scd, sdd = _pair_sums(xo, gap, mc, md)
	if scc <= 0 or sdd <= 0:
		return np.nan
	return scd / np.sqrt(scc * sdd)


@njit(parallel=True, cache=True)
def silhouettes(x, groups, nlev):
	"""Silhouette widths of the rows of X (N x P) in C classifications, GROUPS (N x C, the group 0..k-1 of each row,
	k = NLEV[c] >= 2): as sil_width() of moranvar_plots.R, a = mean distance to the other rows of its group, b = the
	smallest mean distance to the rows of another group (its neighbour; ties: the first group), s = (b - a) / max(a, b),
	0 for a row alone in its group. Returns s and the neighbour (N x C)."""
	n, nv = x.shape
	nc = groups.shape[1]
	kmax = 1
	for c in range(nc):
		kmax = max(kmax, nlev[c])
	cnt = np.zeros((nc, kmax), np.float64)
	for c in range(nc):
		for j in range(n):
			cnt[c, groups[j, c]] += 1.0
	sil = np.empty((n, nc), np.float64)
	nb = np.empty((n, nc), np.int32)
	bs = 64
	nblk = (n + bs - 1) // bs
	for blk in prange(nblk):
		d = np.empty(n, np.float64)
		acc = np.empty(kmax, np.float64)
		for i in range(blk * bs, min(n, (blk + 1) * bs)):
			for j in range(n):
				s = 0.0
				for k in range(nv):
					df = x[j, k] - x[i, k]
					s += df * df
				d[j] = np.sqrt(s)
			for c in range(nc):
				k = nlev[c]
				acc[:k] = 0.0
				for j in range(n):
					acc[groups[j, c]] += d[j]
				g = groups[i, c]
				nh = cnt[c, g]
				a = acc[g] / nh * nh / max(nh - 1.0, 1.0)
				b = np.inf
				nbi = -1
				for h in range(k):
					if h != g and cnt[c, h] > 0:
						v = acc[h] / cnt[c, h]
						if v < b:
							b = v
							nbi = h
				mx = max(a, b)
				sil[i, c] = (b - a) / mx if mx > 0 else 0.0
				if nh == 1:
					sil[i, c] = 0.0
				nb[i, c] = nbi
	return sil, nb


def _read(path, dtype, shape):
	a = np.fromfile(path, dtype=dtype)
	if a.size != int(np.prod(shape)):
		sys.exit(f'ERROR: {path} has {a.size} values, not {int(np.prod(shape))}')
	return a.reshape(shape, order='F')


def _write(path, a, dtype):
	np.asarray(a, dtype=dtype).ravel(order='F').tofile(path)


def main(argv):
	if len(argv) < 4 or argv[0] not in ('tree', 'silhouette'):
		sys.exit(__doc__)
	cmd, folder = argv[0], argv[1]
	n, p = int(argv[2]), int(argv[3])
	x = np.ascontiguousarray(_read(os.path.join(folder, 'x.f8'), '<f8', (n, p)))
	if not np.all(np.isfinite(x)):
		sys.exit('ERROR: the values of the segments for the classification have NA or infinite values')
	if cmd == 'tree':
		w = _read(os.path.join(folder, 'w.f8'), '<f8', (n,)).copy()
		if n < 2 or not np.all(w > 0):
			sys.exit('ERROR: the classification needs 2 segments or more, each with pixels')
		mi, mj, mh = nn_chain(x, w)
		merge, height, order = to_hclust(mi, mj, mh, n)
		_write(os.path.join(folder, 'merge.i4'), merge, '<i4')
		_write(os.path.join(folder, 'height.f8'), height, '<f8')
		_write(os.path.join(folder, 'order.i4'), order, '<i4')
		_write(os.path.join(folder, 'coph.f8'), [cophenetic_cor(x, merge, height, order)], '<f8')
	else:
		nc = int(argv[4])
		groups = np.asfortranarray(_read(os.path.join(folder, 'groups.i4'), '<i4', (n, nc)))
		nlev = groups.max(axis=0) + 1 if n else np.zeros(nc, np.int64)
		if groups.min(initial=0) < 0 or np.any(nlev < 2):
			sys.exit('ERROR: every classification for the silhouettes needs groups 0..k-1 with k >= 2')
		sil, nb = silhouettes(x, groups, nlev.astype(np.int64))
		_write(os.path.join(folder, 's.f8'), sil, '<f8')
		_write(os.path.join(folder, 'neighbour.i4'), nb, '<i4')


if __name__ == '__main__':
	main(sys.argv[1:])
