"""Revisión independiente de compare_features.py: recalcula sus números con código propio.

    /opt/anaconda3/envs/ehabpy/bin/python -s review_check.py

1. Descriptores de los segmentos (media, mediana, SD, MAD, CV) con pandas groupby desde los píxeles leídos aquí.
2. Ward ponderado por píxeles con un algoritmo distinto (voraz por centroides: en cada paso une el par de grupos con
   menor aumento de la suma de cuadrados n_a n_b / (n_a + n_b) |c_a - c_b|^2) y con el código real de
   moranvar_plots.R (std_stats + class_tree, hclust ward.D2 con members) en R.
3. Medidas en el umbral 0.5 con k = 3 y 9 desde los píxeles: R2, JM con covarianza completa (np.cov), silueta
   (sklearn), R2 robusto, ARI de los píxeles (sklearn) y la estabilidad bootstrap (réplicas con los segmentos
   repetidos como puntos separados, no fundidos).
4. ARI entre umbrales de los píxeles (sklearn) con k = 9.
5. Mediana tipificada como std_stats() con STAT_AGG = 'median' (SD de los píxeles estimada con la mediana y la
   varianza de cada segmento) frente a la del script.
6. Bootstrap pareado mediana - media con 200 réplicas nuevas por umbral (¿supera el ruido la ventaja de la mediana?).
7. Diferencias mediana - media en todos los k (umbral 0.5) con las dos tipificaciones de la mediana.
8. ¿Difieren los mapas de media y mediana sobre todo por el TWI? (ARI sin el TWI y con mediana solo en el TWI).
Ficheros auxiliares en review_tmp/ (hri para R, particiones de R, réplicas bootstrap).
"""

import glob
import os
import subprocess
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd
from osgeo import gdal
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score, silhouette_samples

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, 'algorithm_comparison', 'data')
sys.argv = sys.argv[:1]
sys.path.insert(0, HERE)
import compare_features as cf  # solo para comparar con sus resultados

gdal.UseExceptions()
t0 = time.time()


def say(*a):
	print(*a, flush=True)


# --- 1. Píxeles ---
nr, nc = pd.read_csv(os.path.join(DATA, 'grid.csv')).iloc[0]
vars_ = pd.read_csv(os.path.join(DATA, 'vars.csv'))
P = np.empty((nr * nc, 7))
for i in range(7):
	f = sorted(glob.glob(os.path.join(ROOT, 'inVars', f'Guajares_{i + 1}_*.tif')))
	assert len(f) == 1
	ds = gdal.Open(f[0])
	b = ds.GetRasterBand(1)
	x = b.ReadAsArray().astype(float).ravel()
	nd = b.GetNoDataValue()
	bad = (x == nd) | ~np.isfinite(x)
	other = np.sort(np.unique(x[~bad]))[:2]
	say(f'{os.path.basename(f[0])[:40]}: {ds.RasterYSize}x{ds.RasterXSize} nodata={nd} n_nodata={bad.sum()} '
	    f'min={x[~bad].min():.4g} max={x[~bad].max():.4g} (lowest values {other})')
	x[bad] = np.nan
	P[:, i] = x
names = list(vars_['short'])


def segdata(t):
	seg = np.fromfile(os.path.join(DATA, f'seg_0{t}.bin'), dtype=np.int32)
	ok = (seg > 0) & np.all(np.isfinite(P), 1)
	Xo = P[ok]
	Z = (Xo - Xo.mean(0)) / Xo.std(0, ddof=1)
	g = seg[ok]
	return seg, ok, Xo, Z, g


T = 5
seg, ok, Xo, Z, g = segdata(T)
say(f'\nThreshold 0.{T}: {ok.sum()} pixels, {len(np.unique(g))} segments')
df = pd.DataFrame(Z, columns=names)
df['g'] = g
G = df.groupby('g')
mean = G[names].mean().to_numpy()
median = G[names].median().to_numpy()
sd = G[names].std(ddof=0).to_numpy()
mad = 1.4826 * (df[names] - G[names].transform('median')).abs().groupby(df.g).median().to_numpy()
do = pd.DataFrame(Xo, columns=names)
do['g'] = g
om, osd = do.groupby('g')[names].mean().to_numpy(), do.groupby('g')[names].std(ddof=0).to_numpy()
cv = osd / np.abs(om)
n = G.size().to_numpy().astype(float)

s = cf.Segments(T, P.reshape(nr, nc, 7), (nr, nc))
for nm, a, b in [('mean', mean, s.mean), ('median', median, s.median), ('sd', sd, s.sd), ('mad', mad, s.mad),
                 ('cv', cv, s.cv), ('n', n, s.n)]:
	say(f'  descriptor {nm}: max |mine - script| = {np.nanmax(np.abs(a - b)):.2e}')


# --- 2. Ward ponderado independiente ---

def ward_naive(X, w, ks):
	"""Ward de los píxeles partiendo de los segmentos: une el par con menor aumento de la suma de cuadrados."""
	lab = np.arange(len(X))
	C, W = X.astype(float).copy(), w.astype(float).copy()
	act = list(range(len(X)))
	out = {}
	while len(act) > 1:
		A = np.array(act)
		d2 = ((C[A][:, None] - C[A][None]) ** 2).sum(2)
		cost = np.outer(W[A], W[A]) / np.add.outer(W[A], W[A]) * d2
		np.fill_diagonal(cost, np.inf)
		a, b = np.unravel_index(np.argmin(cost), cost.shape)
		i, j = A[a], A[b]
		C[i] = (W[i] * C[i] + W[j] * C[j]) / (W[i] + W[j])
		W[i] += W[j]
		lab[lab == j] = i
		act.remove(j)
		if len(act) in ks:
			out[len(act)] = np.unique(lab, return_inverse=True)[1]
	return out


def same_partition(a, b):
	return adjusted_rand_score(a, b) == 1.0


KS = list(range(2, 21))
say('\nWeighted Ward: independent greedy centroid Ward vs the script (k = 2..20)')
script_parts = {}
for nm in ['mean', 'median', 'mean+sd', 'mean+cv', 'median+mad', 'mean+sd_half']:
	F = s.features(nm)[0]
	mine = ward_naive(F, s.n, KS)
	sc = cf.cut_tree(cf.ward(F, s.n), s.ns, KS)
	script_parts[nm] = sc
	eq = [k for k in KS if same_partition(mine[k], sc[k])]
	say(f'  {nm}: identical partitions at {len(eq)}/{len(KS)} k' + ('' if len(eq) == len(KS) else
	    f'; differ at {[k for k in KS if k not in eq]}'))

# R: std_stats() + class_tree() de moranvar_plots.R sobre un hri con media/mediana y varianza de cada segmento
tmp = os.path.join(HERE, 'review_tmp')
os.makedirs(tmp, exist_ok=True)
var_o = do.groupby('g')[names].var(ddof=0).to_numpy()
med_o = do.groupby('g')[names].median().to_numpy()
for stat, M in [('mean', om), ('median', med_o)]:
	hri = pd.DataFrame({'ecoregion': 5, 'segm_id': np.arange(len(n)) + 1})
	for j, v in enumerate(names):
		hri[v + 'pamean'] = M[:, j]
	for j, v in enumerate(names):
		hri[v + 'pavar'] = var_o[:, j]
	hri['sumpamask'] = n.astype(int)
	for j, v in enumerate(names):
		hri[v] = M[:, j] * n
	hri.to_csv(os.path.join(tmp, f'hri_{stat}.csv'), index=False)
rcode = r'''
ex <- parse(file.path('%s', 'segm_optim', 'moranvar_plots.R'))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name('<-')) && as.character(e[[2]]) %%in%% c('std_stats', 'class_tree')) eval(e)
for (stat in c('mean', 'median')) {
	hri <- read.csv(file.path('%s', paste0('hri_', stat, '.csv')))
	ct <- class_tree(hri, 'ward_area')
	out <- sapply(2:20, function(k) cutree(ct$tree, k=k))
	write.csv(out, file.path('%s', paste0('r_', stat, '.csv')), row.names=FALSE)
	write.csv(ct$z, file.path('%s', paste0('rz_', stat, '.csv')), row.names=FALSE)
}
''' % (ROOT, tmp, tmp, tmp)
with open(os.path.join(tmp, 'r_check.R'), 'w') as fh:
	fh.write(rcode)
r = subprocess.run(['/opt/anaconda3/envs/ehabpy/bin/Rscript', os.path.join(tmp, 'r_check.R')], capture_output=True,
                   text=True)
if r.returncode:
	say(r.stderr)
say('\nmoranvar_plots.R class_tree(ward_area) in R vs the script:')
for stat in ['mean', 'median']:
	R = pd.read_csv(os.path.join(tmp, f'r_{stat}.csv')).to_numpy()
	eq = [k for i, k in enumerate(KS) if same_partition(R[:, i], script_parts[stat][k])]
	ari_k = {k: adjusted_rand_score(np.repeat(R[:, i], n.astype(int)), np.repeat(script_parts[stat][k], n.astype(int)))
	         for i, k in enumerate(KS)}
	say(f'  {stat}: identical at {len(eq)}/{len(KS)} k; pixel ARI k=3 {ari_k[3]:.3f}, k=8 {ari_k[8]:.3f}, '
	    f'k=9 {ari_k[9]:.3f}; min over k {min(ari_k.values()):.3f}')
	if stat == 'median':
		Rz = pd.read_csv(os.path.join(tmp, f'rz_{stat}.csv')).to_numpy()
		ratio = np.nanmedian((s.median - np.average(s.median, axis=0, weights=s.n)) /
		                     np.where(np.abs(Rz) > 1e-9, Rz, np.nan), axis=0)
		say('  std_stats(median) divides by s; script divides by the pixel SD. Ratio s / pixel SD per variable: ' +
		    ', '.join(f'{v} {x:.3f}' for v, x in zip(names, ratio)))
		Rcl = {k: R[:, i] for i, k in enumerate(KS)}


# --- 3. Medidas desde los píxeles ---

def metrics_pix(cl, Zp, gp, nseg_w, D):
	"""R2, JM medio y mínimo con np.cov de los píxeles de cada HFT, silueta (sklearn) por área y R2 robusto."""
	L = cl[gp]
	k = L.max() + 1
	gm = Zp.mean(0)
	r2 = np.mean([np.sum(np.bincount(L, minlength=k) * (np.bincount(L, Zp[:, j], k) / np.bincount(L, minlength=k)
	                                                   - gm[j]) ** 2) / np.sum((Zp[:, j] - gm[j]) ** 2) for j in range(7)])
	mus, covs = [], []
	for h in range(k):
		Y = Zp[L == h]
		mus.append(Y.mean(0))
		covs.append(np.cov(Y.T) + 1e-6 * np.eye(7))
	jm = []
	for a, b in combinations(range(k), 2):
		Cm = (covs[a] + covs[b]) / 2
		dm = mus[a] - mus[b]
		B = dm @ np.linalg.solve(Cm, dm) / 8 + 0.5 * np.log(np.linalg.det(Cm) / np.sqrt(np.linalg.det(covs[a]) *
		                                                                             np.linalg.det(covs[b])))
		jm.append(2 * (1 - np.exp(-B)))
	sil = np.sum(silhouette_samples(D, cl, metric='precomputed') * nseg_w) / nseg_w.sum()
	dz = pd.DataFrame(Zp)
	dz['L'] = L
	medh = dz.groupby('L').transform('median').to_numpy()
	rr = np.mean(1 - np.abs(Zp - medh).sum(0) / np.abs(Zp - np.median(Zp, 0)).sum(0))
	return {'r2': r2, 'jm_mean': np.mean(jm), 'jm_min': np.min(jm), 'sil': sil, 'r2_rob': rr}


tab = pd.read_csv(os.path.join(HERE, f'medidas_0{T}.csv'))
gi = np.searchsorted(np.unique(g), g)  # segmento 0..ns-1 de cada píxel
D = squareform(pdist(mean))
mypart = {nm: ward_naive(s.features(nm)[0], s.n, [3, 9]) for nm in ['mean', 'median', 'mean+sd', 'mean+sd_half']}
say(f'\nMetrics at threshold 0.{T} from the pixels (independent partitions), mine vs medidas_0{T}.csv:')
for nm in ['mean', 'median', 'mean+sd', 'mean+sd_half']:
	for k in [3, 9]:
		m = metrics_pix(mypart[nm][k], Z, gi, s.n, D)
		row = tab[(tab.set == nm) & (tab.k == k)].iloc[0]
		say(f'  {nm:13s} k={k}: ' + '  '.join(f'{x} {m[x]:.4f}/{row[x]:.4f}' for x in m))

say('\nPixel ARI vs mean (sklearn on expanded pixels), mine vs ari_ref:')
for nm in ['median', 'mean+sd', 'mean+sd_half']:
	for k in [3, 9]:
		a = adjusted_rand_score(mypart['mean'][k][gi], mypart[nm][k][gi])
		say(f'  {nm} k={k}: {a:.4f} / {tab[(tab.set == nm) & (tab.k == k)].ari_ref.iloc[0]:.4f}')

# Bootstrap: segmentos repetidos como puntos separados (no fundidos), mismas réplicas (semilla 1000 + t)
say('\nBootstrap ARI (k=3, 9), duplicates as separate points, mine vs medidas:')
for nm in ['mean', 'median']:
	F = s.features(nm)[0]
	full = ward_naive(F, s.n, [3, 9])
	rng = np.random.default_rng(1000 + T)
	acc = {3: [], 9: []}
	for _ in range(cf.NBOOT):
		idx = rng.integers(0, s.ns, s.ns)
		parts = ward_naive(F[idx], s.n[idx], [3, 9])
		inb = np.unique(idx)
		first = np.array([np.flatnonzero(idx == u)[0] for u in inb])
		for k in [3, 9]:
			pb = parts[k][first]
			acc[k].append(adjusted_rand_score(np.repeat(pb, s.n[inb].astype(int)),
			                                  np.repeat(full[k][inb], s.n[inb].astype(int))))
	for k in [3, 9]:
		say(f'  {nm} k={k}: {np.mean(acc[k]):.4f} / {tab[(tab.set == nm) & (tab.k == k)].boot.iloc[0]:.4f}')

# --- 4. ARI entre umbrales (k = 9) ---
say('\nCross-threshold pixel ARI at k=9 (mean of the 10 pairs), mine vs estabilidad_umbrales.csv:')
stab = pd.read_csv(os.path.join(HERE, 'estabilidad_umbrales.csv'))
maps = {}
for t in [3, 4, 5, 6, 7]:
	sg, okt, _, Zt, gt = segdata(t)
	st = cf.Segments(t, P.reshape(nr, nc, 7), (nr, nc))
	for nm in ['mean', 'median']:
		lab = ward_naive(st.features(nm)[0], st.n, [9])[9]
		full = np.full(nr * nc, -1)
		full[okt] = lab[np.searchsorted(np.unique(gt), gt)]
		maps[(t, nm)] = full
for nm in ['mean', 'median']:
	v = []
	for a, b in combinations([3, 4, 5, 6, 7], 2):
		c = (maps[(a, nm)] >= 0) & (maps[(b, nm)] >= 0)
		v.append(adjusted_rand_score(maps[(a, nm)][c], maps[(b, nm)][c]))
	say(f'  {nm}: {np.mean(v):.4f} / {stab[(stab.set == nm) & (stab.k == 9)].ari_thr.iloc[0]:.4f}')

# --- 5. Mediana tipificada como std_stats() del pipeline ---
say('\nMedian set with the pipeline std_stats(median) scaling (R class_tree) vs the script median set:')
for k in [3, 9]:
	m = metrics_pix(np.unique(Rcl[k], return_inverse=True)[1], Z, gi, s.n, D)
	row = tab[(tab.set == 'median') & (tab.k == k)].iloc[0]
	a_script = adjusted_rand_score(Rcl[k][gi], script_parts['median'][k][gi])
	a_mean = adjusted_rand_score(Rcl[k][gi], mypart['mean'][k][gi])
	say(f'  k={k}: ARI with script median {a_script:.3f}, with mean {a_mean:.3f} (script median vs mean '
	    f'{row.ari_ref:.3f}); ' + '  '.join(f'{x} {m[x]:.4f}/{row[x]:.4f}' for x in m))

# --- 6. ¿Es la ventaja de la mediana en estabilidad bootstrap mayor que el ruido? (más réplicas, pareadas) ---
NB2 = 200
say(f'\nPaired bootstrap median - mean, {NB2} resamples per threshold (ward of the script, verified above):')
rows = []
for t in [3, 4, 5, 6, 7]:
	st = cf.Segments(t, P.reshape(nr, nc, 7), (nr, nc))
	F = {nm: st.features(nm)[0] for nm in ['mean', 'median']}
	ref = {nm: cf.cut_tree(cf.ward(F[nm], st.n), st.ns, KS) for nm in F}
	rng = np.random.default_rng(777 + t)
	for b in range(NB2):
		mlt = np.bincount(rng.integers(0, st.ns, st.ns), minlength=st.ns)
		inb = np.flatnonzero(mlt > 0)
		for nm in F:
			parts = cf.cut_tree(cf.ward(F[nm][inb], st.n[inb] * mlt[inb]), len(inb), KS)
			for k in KS:
				rows.append((t, b, nm, k, cf.ari(parts[k], ref[nm][k][inb], st.n[inb])))
bt = pd.DataFrame(rows, columns=['threshold', 'rep', 'set', 'k', 'ari'])
w = bt.pivot_table(index=['threshold', 'rep', 'k'], columns='set', values='ari').reset_index()
w['d'] = w['median'] - w['mean']
cell = w.groupby(['threshold', 'k']).d.mean()
say(f'  per (threshold, k) mean difference: median {cell.median():+.4f}, share > 0 {np.mean(cell > 0):.2f} '
    f'(script with 30 resamples: +0.013, 0.73)')
for t, gt_ in w.groupby('threshold'):
	rep = gt_.groupby('rep').d.mean()  # una réplica = una unidad (promedio de los k)
	say(f'  threshold 0.{t}: mean over k and resamples {rep.mean():+.4f}, SE over resamples {rep.std(ddof=1) / np.sqrt(len(rep)):.4f}')
rep = w.groupby(['threshold', 'rep']).d.mean().groupby('rep').mean()
say(f'  all thresholds: {rep.mean():+.4f} +- {1.96 * rep.std(ddof=1) / np.sqrt(len(rep)):.4f} (95% CI half-width, resamples as units)')
bt.to_csv(os.path.join(HERE, 'review_tmp', 'boot_pareado.csv'), index=False)

# --- 7. Mediana con la tipificación de std_stats(median): diferencias con la media en todos los k (umbral 0.5) ---
say('\nThreshold 0.5, all k: median vs mean with the script scaling and with the pipeline std_stats(median) scaling')
R = pd.read_csv(os.path.join(tmp, 'r_median.csv')).to_numpy()
mm_ = []
for i, k in enumerate(KS):
	a = cf.measures(R[:, i] - 1, s, s.features('median')[0])
	b = tab[(tab.set == 'median') & (tab.k == k)].iloc[0]
	c = tab[(tab.set == 'mean') & (tab.k == k)].iloc[0]
	mm_.append({x: (a[x] - c[x], b[x] - c[x]) for x in ['r2', 'jm_mean', 'jm_min', 'sil', 'r2_rob']})
for x in ['r2', 'jm_mean', 'jm_min', 'sil', 'r2_rob']:
	p = np.array([r[x][0] for r in mm_])
	q = np.array([r[x][1] for r in mm_])
	say(f'  d_{x}: pipeline scaling median {np.median(p):+.4f} (better {np.mean(p > 1e-9):.2f}), '
	    f'script scaling {np.median(q):+.4f} (better {np.mean(q > 1e-9):.2f})')
say(f'\nDone in {time.time() - t0:.0f} s')

# --- 8. ¿Difieren los mapas de media y mediana sobre todo por el TWI? (umbral 0.5) ---
say('\nThreshold 0.5: pixel ARI mean vs median at the same k, all 7 variables / median only for TWI / without TWI')
mix = s.mean.copy()
mix[:, 6] = s.median[:, 6]
pm = {nm: cf.cut_tree(cf.ward(F_, s.n), s.ns, KS) for nm, F_ in
      [('mean', s.mean), ('median', s.median), ('mix', mix), ('mean6', s.mean[:, :6]), ('median6', s.median[:, :6])]}
for k in [3, 5, 9, 15]:
	say(f'  k={k}: all {cf.ari(pm["mean"][k], pm["median"][k], s.n):.3f}  TWI-only median '
	    f'{cf.ari(pm["mean"][k], pm["mix"][k], s.n):.3f}  without TWI {cf.ari(pm["mean6"][k], pm["median6"][k], s.n):.3f}')
allk = [(cf.ari(pm['mean'][k], pm['median'][k], s.n), cf.ari(pm['mean'][k], pm['mix'][k], s.n),
         cf.ari(pm['mean6'][k], pm['median6'][k], s.n)) for k in KS]
say('  median over k 2-20: all %.3f, TWI-only median %.3f, without TWI %.3f' % tuple(np.median(allk, 0)))
say(f'\nDone in {time.time() - t0:.0f} s')
