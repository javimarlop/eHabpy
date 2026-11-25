"""Compara los descriptores de los segmentos con que se clasifican en HFTs: media o mediana de cada variable y
si añadir su desviación típica (SD), MAD o coeficiente de variación (CV) dentro del segmento mejora la clasificación.

    /opt/anaconda3/envs/ehabpy/bin/python -s compare_features.py [UMBRALES] [NBOOT]

UMBRALES: índices separados por comas (por defecto 3,4,5,6,7); NBOOT: réplicas bootstrap de los segmentos para la
estabilidad (por defecto 30). Usa los datos de ../algorithm_comparison/data (segment_sums.R) y los píxeles de
../inVars/Guajares_<i>_*.tif (los 7 de la raíz de inVars, en el orden de vars.csv).

Todos los conjuntos usan el mismo algoritmo, el de CLASS_METHOD = 'ward_area': Ward ponderado por el nº de píxeles de
cada segmento. Se implementa (agnes() de compare_algorithms.py) con la fórmula de Lance-Williams con pesos: la
distancia inicial entre los segmentos i y j es 2 n_i n_j / (n_i + n_j) |x_i - x_j|^2 (el doble del aumento de la suma
de cuadrados dentro de los grupos al unirlos, contando sus píxeles) y al unir dos grupos se actualiza con sus pesos
acumulados. Es exactamente el Ward de todos los píxeles partiendo de los segmentos (como hclust(..., members = n)
en moranvar_plots.R); replicar los segmentos n veces daría el mismo árbol pero con ~10^6 puntos.

Conjuntos de descriptores de cada segmento (V = 7 variables):
  mean          media de cada variable (referencia: el método actual)
  median        mediana (STAT_AGG = 'median')
  mean+sd       media y SD de los píxeles del segmento
  median+mad    mediana y MAD (mediana de las desviaciones absolutas a la mediana, x 1.4826: SD si fuera normal)
  mean+cv       media y CV = SD / |media| en unidades originales
  mean+sd_half  como mean+sd con el bloque de SD escalado para pesar la mitad que el de medias (ver abajo)
  mean+sd_pix   (extra) media y SD en unidades de la SD de los píxeles, sin re-tipificar: la escala "natural"
Tipificación: los descriptores de posición (media, mediana) son los de std_stats(): cada variable tipificada con la
SD de los píxeles (dentro + entre segmentos). Los de dispersión (SD, MAD, CV) son una propiedad del segmento sin
unidades comparables a la posición, y se tipifican con su media y SD ponderadas por píxeles entre los segmentos:
cada uno tiene varianza 1 entre segmentos (los de posición algo menos, la parte de la varianza de los píxeles que
está entre segmentos, ~0.8-0.9), así que mean+sd da a la dispersión algo más de peso total que a la posición
(peso = suma de varianzas ponderadas del bloque, ver pesos.csv). mean+sd_half escala el bloque de SD para que su peso
sea 0.5 del de las medias; mean+sd_pix deja la SD en unidades de píxel (peso mucho menor, la SD varía poco).

Medidas (las de compare_algorithms.py, con los píxeles tipificados):
  r2       varianza de los píxeles entre HFTs / total, media de las 7 variables. Es exactamente lo que maximiza
           (de forma voraz) Ward sobre las medias ponderadas por área: favorece a 'mean' por construcción
  jm_mean  Jeffries-Matusita medio entre pares de HFTs (normales multivariantes con la covarianza completa de
           los píxeles: ya tiene en cuenta la dispersión dentro de los HFTs)
  jm_min   JM del par de HFTs menos separable
  sil      silueta media de los segmentos (distancias entre sus medias tipificadas, como compare_algorithms.py;
           favorece a 'mean'), ponderada por área
  n_eff    nº efectivo de HFTs, exp(entropía de las fracciones de área)
Medidas añadidas que no favorecen a 'mean' por construcción:
  r2_rob   R2 robusto de los píxeles con medianas: 1 - sum|x - mediana del HFT| / sum|x - mediana global|, media
           de las variables (favorece algo a 'median')
  sd_eta2  parte de la varianza (ponderada por área) de la SD de los segmentos que está entre HFTs, media de las
           variables: ¿agrupa la clasificación segmentos de heterogeneidad interna parecida? (favorece a +sd)
  sil_own  silueta en el propio espacio de descriptores de cada conjunto (no comparable entre espacios en sentido
           estricto; se da como referencia)
  ari_ref  acuerdo (índice de Rand ajustado por área) con 'mean' con el mismo nº de HFTs
Estabilidad:
  boot     ARI medio por área entre la partición de los segmentos de una réplica bootstrap de los segmentos (los
           mismos remuestreos para todos los conjuntos) y la de todos, en los segmentos de la réplica, mismo k
  ari_thr  ARI de los píxeles entre los mapas de HFTs de dos umbrales de segmentación con el mismo k, media de los
           pares de umbrales: un mapa que depende menos del umbral es más estable
  cortes   nº de HFTs de los cortes p10-p90 del árbol (fracción de su altura máxima, SIMIL = 3) y su variación
           entre las réplicas bootstrap

Salidas en esta carpeta: medidas_0k.csv (todas las medidas por conjunto y k), resumen.csv (diferencia mediana con
'mean' al mismo k y % de casos mejores), estabilidad_umbrales.csv, cortes.csv, cv_por_variable.csv, pesos.csv,
media_vs_mediana.csv, fig1_curvas_0k.png, fig2_resumen.png y fig3_cv_estabilidad.png (cada figura con su .csv).
"""

import csv
import glob
import os
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd
from osgeo import gdal
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

gdal.UseExceptions()
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, 'algorithm_comparison', 'data')
INVARS = os.path.join(ROOT, 'inVars')
KMAX = 20
KS = list(range(2, KMAX + 1))
THR = [int(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else [3, 4, 5, 6, 7]
NBOOT = int(sys.argv[2]) if len(sys.argv) > 2 else 30
NODATA = -9999
SETS = ['mean', 'median', 'mean+sd', 'median+mad', 'mean+cv', 'mean+sd_half', 'mean+sd_pix']
METRICS = ['r2', 'jm_mean', 'jm_min', 'sil', 'n_eff', 'r2_rob', 'sd_eta2', 'sil_own', 'boot']
BETTER = ['r2', 'jm_mean', 'jm_min', 'sil', 'r2_rob', 'sd_eta2', 'sil_own', 'boot']  # mayor = mejor
CUTS = [('p10', .1), ('p25', .25), ('p50', .5), ('p75', .75), ('p90', .9)]


# --- Píxeles y descriptores de los segmentos ---

def read_pixels():
	"""Valores de las 7 variables (filas x columnas x 7), NaN sin datos; comprueba tamaño y orden con data/."""
	nr, nc = pd.read_csv(os.path.join(DATA, 'grid.csv')).iloc[0]
	vars_ = pd.read_csv(os.path.join(DATA, 'vars.csv'))
	P = []
	for i, (v, short) in enumerate(zip(vars_['var'], vars_['short']), 1):
		f = glob.glob(os.path.join(INVARS, f'Guajares_{i}_*.tif'))
		assert len(f) == 1, f
		name = os.path.basename(f[0])[:-4].replace('-', '_')
		assert name == v, (name, v)  # mismo orden que vars.csv
		ds = gdal.Open(f[0])
		b = ds.GetRasterBand(1)
		assert (ds.RasterYSize, ds.RasterXSize) == (nr, nc), (f[0], ds.RasterYSize, ds.RasterXSize)
		assert b.GetNoDataValue() == NODATA
		x = b.ReadAsArray().astype(float)
		x[x == NODATA] = np.nan
		P.append(x)
	return np.stack(P, -1), list(vars_['short']), (nr, nc)


def group_median(x, g, ng):
	"""Mediana de x en cada grupo g (0..ng-1)."""
	o = np.lexsort((x, g))
	xs = x[o]
	cnt = np.bincount(g, minlength=ng)
	st = np.concatenate([[0], np.cumsum(cnt)[:-1]])
	lo, hi = st + (cnt - 1) // 2, st + cnt // 2
	med = np.full(ng, np.nan)
	ok = cnt > 0
	med[ok] = (xs[lo[ok]] + xs[hi[ok]]) / 2
	return med


def wstd(X, w):
	"""Media y SD ponderadas por w de cada columna de X."""
	mu = np.average(X, axis=0, weights=w)
	return mu, np.sqrt(np.average((X - mu) ** 2, axis=0, weights=w))


class Segments:
	"""Píxeles tipificados y descriptores de los segmentos de un umbral."""
	def __init__(self, t, P, shape):
		d = pd.read_csv(os.path.join(DATA, f'sums_0{t}.csv'), dtype={'segm_id': str})
		self.rows = np.flatnonzero(d['n'].to_numpy() > 0)
		d = d.iloc[self.rows]
		nv = self.nv = P.shape[-1]
		seg = np.fromfile(os.path.join(DATA, f'seg_0{t}.bin'), dtype=np.int32)
		assert seg.size == shape[0] * shape[1]
		X = P.reshape(-1, nv)
		ok = (seg > 0) & np.all(np.isfinite(X), 1)
		rowmap = np.full(seg.max() + 1, -1)
		rowmap[self.rows + 1] = np.arange(len(self.rows))
		self.pix = np.flatnonzero(ok)  # índice en la malla de cada píxel usado
		self.g = rowmap[seg[ok]]  # segmento (0..ns-1) de cada píxel usado
		assert self.g.min() >= 0
		Xo = X[ok]
		ns = self.ns = len(self.rows)
		self.Z = (Xo - Xo.mean(0)) / Xo.std(0, ddof=1)  # como scale() de segment_sums.R
		self.n = np.bincount(self.g, minlength=ns).astype(float)
		self.N = self.n.sum()
		# comprobación con las sumas de segment_sums.R (mismo orden de píxeles y variables)
		S = np.stack([np.bincount(self.g, self.Z[:, j], ns) for j in range(nv)], 1)
		assert np.allclose(self.n, d['n'].to_numpy()), 'n'
		Scsv = d[[f's{i}' for i in range(1, nv + 1)]].to_numpy()
		assert np.allclose(S, Scsv, atol=1e-6 * np.abs(Scsv).max()), 'sums'
		qc = [c for c in d.columns if c.startswith('q_')]
		self.pairs = np.array([[int(x) - 1 for x in c.split('_')[1:]] for c in qc])
		self.diag = np.array([qc.index(f'q_{i}_{i}') for i in range(1, nv + 1)])
		self.S, self.Q = S, d[qc].to_numpy()
		# descriptores en unidades de la SD de los píxeles (Z) y CV en unidades originales
		self.mean = S / self.n[:, None]
		self.sd = np.sqrt(np.maximum(np.stack([np.bincount(self.g, self.Z[:, j] ** 2, ns) for j in range(nv)], 1)
		                             / self.n[:, None] - self.mean ** 2, 0))
		self.median = np.stack([group_median(self.Z[:, j], self.g, ns) for j in range(nv)], 1)
		self.mad = 1.4826 * np.stack([group_median(np.abs(self.Z[:, j] - self.median[self.g, j]), self.g, ns)
		                              for j in range(nv)], 1)
		om = np.stack([np.bincount(self.g, Xo[:, j], ns) for j in range(nv)], 1) / self.n[:, None]
		osd = np.sqrt(np.maximum(np.stack([np.bincount(self.g, Xo[:, j] ** 2, ns) for j in range(nv)], 1)
		                         / self.n[:, None] - om ** 2, 0))
		self.omean, self.osd = om, osd
		self.pmin, self.pmax = Xo.min(0), Xo.max(0)
		self.cv = osd / np.abs(om)
		self.D = squareform(pdist(self.mean))  # distancias entre medias tipificadas (silueta de compare_algorithms.py)
		# para r2_rob: valores de cada variable ordenados y desviación absoluta a la mediana global
		self.order = [np.argsort(self.Z[:, j], kind='stable') for j in range(nv)]
		self.zsort = [self.Z[o, j] for j, o in enumerate(self.order)]
		self.gabs = np.array([np.abs(self.Z[:, j] - np.median(self.Z[:, j])).sum() for j in range(nv)])

	def features(self, name):
		"""Matriz de descriptores de un conjunto y peso de cada bloque (suma de varianzas ponderadas)."""
		loc = self.median if name.startswith('median') else self.mean
		if '+' not in name:
			return loc, (wstd(loc, self.n)[1] ** 2).sum(), 0.0
		disp = {'sd': self.sd, 'mad': self.mad, 'cv': self.cv, 'sd_half': self.sd, 'sd_pix': self.sd}[name.split('+')[1]]
		mu, sdv = wstd(disp, self.n)
		if name.endswith('sd_pix'):
			disp = disp - mu
		else:
			disp = (disp - mu) / np.where(sdv > 0, sdv, 1)
		wl = (wstd(loc, self.n)[1] ** 2).sum()
		if name.endswith('sd_half'):
			disp = disp * np.sqrt(0.5 * wl / (wstd(disp, self.n)[1] ** 2).sum())
		return np.hstack([loc, disp]), wl, (wstd(disp, self.n)[1] ** 2).sum()


# --- Ward ponderado, cortes y medidas (de compare_algorithms.py) ---

def ward(X, w):
	"""Uniones (i, j, altura) del Ward ponderado por w (Lance-Williams con pesos; ver el docstring)."""
	n = len(X)
	D = squareform(pdist(X))
	D = 2 * np.outer(w, w) / np.add.outer(w, w) * D ** 2
	W = w.astype(float).copy()
	np.fill_diagonal(D, np.inf)
	merges = []
	for _ in range(n - 1):
		i, j = sorted(np.unravel_index(np.argmin(D), D.shape))
		h, di, dj = D[i, j], D[i].copy(), D[j].copy()
		new = ((W[i] + W) * di + (W[j] + W) * dj - W * h) / (W[i] + W[j] + W)
		D[i], D[:, i] = new, new
		D[j], D[:, j] = np.inf, np.inf
		D[i, i] = np.inf
		W[i] += W[j]
		merges.append((i, j, np.sqrt(h)))
	return merges


def cut_tree(merges, n, ks):
	lab, out = np.arange(n), {}
	for m, (i, j, _) in enumerate(merges):
		lab[lab == j] = i
		if n - m - 1 in ks:
			out[n - m - 1] = np.unique(lab, return_inverse=True)[1]
	return out


def cut_height(merges, n, f):
	hmax = max(h for _, _, h in merges)
	return n - sum(h <= f * hmax for _, _, h in merges)


def rowsum(x, cl, k):
	out = np.zeros((k,) + x.shape[1:])
	np.add.at(out, cl, x)
	return out


def silhouette(D, cl):
	"""Anchura de silueta de cada punto (0 en un grupo de un solo punto, como sklearn y cluster::silhouette)."""
	k, n = cl.max() + 1, len(cl)
	if k < 2 or k >= n:
		return np.zeros(n)
	nh = np.bincount(cl, minlength=k)
	md = rowsum(D, cl, k) / nh[:, None]
	idx = np.arange(n)
	a = md[cl, idx] * nh[cl] / np.maximum(nh[cl] - 1, 1)
	md[cl, idx] = np.inf
	b = md.min(0)
	s = np.where(np.maximum(a, b) > 0, (b - a) / np.maximum(a, b), 0)
	s[nh[cl] == 1] = 0
	return s


def ari(x, y, w):
	_, x = np.unique(x, return_inverse=True)
	_, y = np.unique(y, return_inverse=True)
	T = np.zeros((x.max() + 1, y.max() + 1))
	np.add.at(T, (x, y), w)
	cb = lambda v: np.sum(v * (v - 1) / 2)
	e = cb(T.sum(1)) * cb(T.sum(0)) / cb(np.array([w.sum()]))
	return (cb(T) - e) / ((cb(T.sum(1)) + cb(T.sum(0))) / 2 - e)


def ari_pix(x, y):
	"""ARI entre dos etiquetados de los mismos píxeles (enteros >= 0)."""
	kx, ky = x.max() + 1, y.max() + 1
	T = np.bincount(x.astype(np.int64) * ky + y, minlength=kx * ky).reshape(kx, ky).astype(float)
	cb = lambda v: np.sum(v * (v - 1) / 2)
	e = cb(T.sum(1)) * cb(T.sum(0)) / cb(np.array([T.sum()]))
	return (cb(T) - e) / ((cb(T.sum(1)) + cb(T.sum(0))) / 2 - e)


def measures(cl, s, F):
	"""Medidas de la partición cl de los segmentos s (F: descriptores del conjunto, para sil_own)."""
	_, cl = np.unique(cl, return_inverse=True)
	k = cl.max() + 1
	Nh, Sh, Qh = rowsum(s.n, cl, k), rowsum(s.S, cl, k), rowsum(s.Q, cl, k)
	St = s.S.sum(0)
	tot = s.Q[:, s.diag].sum(0) - St ** 2 / s.N
	r2 = np.mean(((Sh ** 2 / Nh[:, None]).sum(0) - St ** 2 / s.N) / tot)
	mu = Sh / Nh[:, None]
	C = np.zeros((k, s.nv, s.nv))
	C[:, s.pairs[:, 0], s.pairs[:, 1]] = Qh
	C[:, s.pairs[:, 1], s.pairs[:, 0]] = Qh
	cov = (C - Nh[:, None, None] * mu[:, :, None] * mu[:, None, :]) / np.maximum(Nh - 1, 1)[:, None, None]
	cov += 1e-6 * np.eye(s.nv)
	a, b = np.array(list(combinations(range(k), 2))).T
	Cm = (cov[a] + cov[b]) / 2
	dm = mu[a] - mu[b]
	ld = np.linalg.slogdet(cov)[1]
	B = np.einsum('pi,pi->p', dm, np.linalg.solve(Cm, dm[:, :, None])[:, :, 0]) / 8 + \
		0.5 * (np.linalg.slogdet(Cm)[1] - 0.5 * (ld[a] + ld[b]))
	jm = 2 * (1 - np.exp(-B))
	p = Nh / s.N
	sil = silhouette(s.D, cl)
	sil_own = silhouette(squareform(pdist(F)), cl)
	# sd_eta2: parte entre HFTs de la varianza ponderada de la SD de los segmentos
	msd = rowsum(s.n[:, None] * s.sd, cl, k) / Nh[:, None]
	gsd = np.average(s.sd, axis=0, weights=s.n)
	tsd = (s.n[:, None] * (s.sd - gsd) ** 2).sum(0)
	sd_eta2 = np.mean((Nh[:, None] * (msd - gsd) ** 2).sum(0) / tsd)
	# r2_rob: medianas de los píxeles de cada HFT (desde el orden de cada variable, ordenación estable por HFT)
	L = cl[s.g]
	cnt = np.bincount(L, minlength=k)
	st = np.concatenate([[0], np.cumsum(cnt)[:-1]])
	lo, hi = st + (cnt - 1) // 2, st + cnt // 2
	rr = []
	for j in range(s.nv):
		Ls = L[s.order[j]]
		xs = s.zsort[j][np.argsort(Ls, kind='stable')]
		med = (xs[lo] + xs[hi]) / 2
		rr.append(1 - np.abs(s.Z[:, j] - med[L]).sum() / s.gabs[j])
	return {'hfts': k, 'r2': r2, 'jm_mean': jm.mean(), 'jm_min': jm.min(), 'sil': np.sum(sil * s.n) / s.N,
	        'n_eff': np.exp(-np.sum(p * np.log(p))), 'r2_rob': np.mean(rr), 'sd_eta2': sd_eta2,
	        'sil_own': np.sum(sil_own * s.n) / s.N}


def bootstrap(F, s, ref, nboot, seed):
	"""ARI medio (por k) entre la partición de cada réplica bootstrap de los segmentos y la de todos (ref), y nº de
	HFTs de los cortes de cada réplica. Las réplicas dependen solo de la semilla: las mismas para cada conjunto."""
	rng = np.random.default_rng(seed)
	A = {k: [] for k in KS}
	cuts = {p: [] for p, _ in CUTS}
	for _ in range(nboot):
		m = np.bincount(rng.integers(0, s.ns, s.ns), minlength=s.ns)
		inb = np.flatnonzero(m > 0)
		tree = ward(F[inb], s.n[inb] * m[inb])
		parts = cut_tree(tree, len(inb), [k for k in KS if k <= len(inb)])
		for k, p in parts.items():
			A[k].append(ari(p, ref[k][inb], s.n[inb]))
		for p, f in CUTS:
			cuts[p].append(cut_height(tree, len(inb), f))
	return {k: np.mean(v) for k, v in A.items()}, cuts


# --- Comparación ---

def cv_report(s, names, t):
	"""Cómo se comporta el CV de cada variable entre los segmentos."""
	rows = []
	for j, v in enumerate(names):
		m, cv = s.omean[:, j], s.cv[:, j]
		w = s.n / s.N
		pos = np.sum(w[m > 0])
		rows.append({'threshold': t, 'var': v, 'pix_min': s.pmin[j], 'pix_max': s.pmax[j], 'seg_mean_min': m.min(),
		             'seg_mean_max': m.max(), 'area_mean_pos': pos, 'min_abs_mean': np.abs(m).min(),
		             'cv_p50': np.median(cv), 'cv_p95': np.quantile(cv, .95), 'cv_max': cv.max(),
		             'cv_max_over_p50': cv.max() / np.median(cv),
		             'cv_sd_w_over_iqr': wstd(cv[:, None], s.n)[1][0] / (np.quantile(cv, .75) - np.quantile(cv, .25)),
		             'segments_cv_gt_1': int(np.sum(cv > 1))})
	return rows


def run(t, P, shape, names):
	s = Segments(t, P, shape)
	print(f'\n=== Threshold 0.{t}: {s.ns} segments, {s.N:.0f} pixels ===', flush=True)
	rows, parts, wts, cuts, mm = [], {}, [], [], []
	# media frente a mediana: correlación ponderada de los segmentos y asimetría de los píxeles
	for j, v in enumerate(names):
		x, y = s.mean[:, j], s.median[:, j]
		c = np.cov(x, y, aweights=s.n)
		skew = np.bincount(s.g, (s.Z[:, j] - s.mean[s.g, j]) ** 3, s.ns) / s.n / np.maximum(s.sd[:, j], 1e-12) ** 3
		mm.append({'threshold': t, 'var': v, 'r_mean_median': c[0, 1] / np.sqrt(c[0, 0] * c[1, 1]),
		           'mean_abs_diff_pixsd': np.average(np.abs(x - y), weights=s.n),
		           'between_sd_mean': np.sqrt(c[0, 0]), 'between_sd_median': np.sqrt(c[1, 1]),
		           'within_skew_median': np.median(skew), 'within_skew_p90abs': np.quantile(np.abs(skew), .9)})
	for name in SETS:
		F, wl, wd = s.features(name)
		wts.append({'threshold': t, 'set': name, 'w_loc': wl, 'w_disp': wd, 'ratio': wd / wl})
		tree = ward(F, s.n)
		pk = cut_tree(tree, s.ns, KS)
		boot, bc = bootstrap(F, s, pk, NBOOT, seed=1000 + t)
		c = {'threshold': t, 'set': name}
		for p, f in CUTS:
			c[p] = cut_height(tree, s.ns, f)
			c[p + '_boot_med'] = np.median(bc[p])
			c[p + '_boot_iqr'] = np.subtract(*np.quantile(bc[p], [.75, .25]))
			c[p + '_boot_same'] = np.mean(np.array(bc[p]) == c[p])
		cuts.append(c)
		for k in KS:
			parts[(name, k)] = pk[k]
			rows.append({'threshold': t, 'set': name, 'k': k, **measures(pk[k], s, F), 'boot': boot[k]})
		print(f'{name}: cuts p25/p50 = {c["p25"]}/{c["p50"]}, disp/loc weight {wd / wl:.2f}', flush=True)
	tab = pd.DataFrame(rows)
	tab['ari_ref'] = [ari(parts[(r.set, r.k)], parts[('mean', r.k)], s.n) for r in tab.itertuples()]
	cd = pd.DataFrame(cuts)
	# ARI de cada conjunto en su corte p25/p50 con 'mean' en el suyo
	for p in ['p25', 'p50']:
		kref = int(cd.loc[cd.set == 'mean', p].iloc[0])
		refp = cut_tree(ward(s.features('mean')[0], s.n), s.ns, [kref])[kref]
		cd[p + '_ari_mean'] = [ari(cut_tree(ward(s.features(n)[0], s.n), s.ns, [int(kk)])[int(kk)], refp, s.n)
		                       for n, kk in zip(cd.set, cd[p])]
	tab.to_csv(os.path.join(HERE, f'medidas_0{t}.csv'), index=False)
	return s, tab, parts, cd, pd.DataFrame(wts), pd.DataFrame(mm), cv_report(s, names, t)


def summary(tab):
	ref = tab[tab.set == 'mean'].set_index(['threshold', 'k'])[METRICS]
	x = tab[tab.set != 'mean'].copy()
	r = ref.loc[list(zip(x.threshold, x.k))].to_numpy()
	for i, m in enumerate(METRICS):
		x['d_' + m] = x[m].to_numpy() - r[:, i]
	out = []
	for name, g in x.groupby('set', sort=False):
		row = {'set': name, 'cases': len(g), 'ari_ref': g.ari_ref.median()}
		for m in METRICS:
			row['d_' + m] = g['d_' + m].median()
		for m in BETTER:
			row['better_' + m] = np.mean(g['d_' + m] > 1e-9)
		out.append(row)
	return x, pd.DataFrame(out)


def thr_stability(keep):
	"""ARI de los píxeles entre los mapas de HFTs de cada par de umbrales con el mismo k."""
	rows = []
	ts = sorted(keep)
	for name in SETS:
		for k in KS:
			v = []
			for a, b in combinations(ts, 2):
				sa, sb = keep[a][0], keep[b][0]
				common, ia, ib = np.intersect1d(sa.pix, sb.pix, assume_unique=True, return_indices=True)
				v.append(ari_pix(keep[a][2][(name, k)][sa.g[ia]], keep[b][2][(name, k)][sb.g[ib]]))
			rows.append({'set': name, 'k': k, 'ari_thr': np.mean(v), 'ari_thr_min': np.min(v)})
	return pd.DataFrame(rows)


# --- Figuras ---

INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#6250d6']
COL = dict(zip(SETS[1:], SLOTS), mean=INK)
MARK = dict(zip(SETS, ['o', 's', '^', 'v', 'D', 'P', 'X']))
LABEL = {'mean': 'media (actual)', 'median': 'mediana', 'mean+sd': 'media + SD', 'median+mad': 'mediana + MAD',
         'mean+cv': 'media + CV', 'mean+sd_half': 'media + SD (peso 0.5)', 'mean+sd_pix': 'media + SD sin re-tipificar'}
TITLE = {'r2': 'R² de los píxeles (sesgado hacia la media)', 'jm_mean': 'JM medio entre pares de HFTs',
         'jm_min': 'JM mínimo (par menos separable)', 'sil': 'Silueta en el espacio de medias (por área)',
         'n_eff': 'Nº efectivo de HFTs', 'r2_rob': 'R² robusto (medianas, desv. absolutas)',
         'sd_eta2': 'Parte entre HFTs de la varianza de la SD de los segmentos',
         'sil_own': 'Silueta en el espacio propio', 'boot': 'Estabilidad bootstrap (ARI)'}
SHORT = {'r2': 'Δ R² píxeles\n(sesgado a la media)', 'jm_mean': 'Δ JM medio', 'jm_min': 'Δ JM mínimo',
         'sil': 'Δ silueta\n(espacio de medias)', 'r2_rob': 'Δ R² robusto\n(medianas)', 'sd_eta2': 'Δ parte entre HFTs\nde la SD (sesgado a +SD)',
         'boot': 'Δ estabilidad\nbootstrap', 'ari_ref': 'ARI con la media\n(mismo k)'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#c3c2b7', 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.labelcolor': INK2, 'text.color': INK, 'axes.titleweight': 'bold', 'axes.titlesize': 9.5,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.6, 'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb',
                     'savefig.facecolor': '#fcfcfb', 'legend.frameon': False})


def save_csv(df, png):
	"""CSV con los datos dibujados junto a cada figura (separado por espacios, colores entre comillas)."""
	df.to_csv(png[:-4] + '.csv', sep=' ', index=False, quoting=csv.QUOTE_NONNUMERIC)


def fig_curves(tab, t, cd):
	fig, axs = plt.subplots(3, 3, figsize=(13, 11))
	keep = []
	for ax, m in zip(axs.flat, METRICS):
		for name in SETS:
			g = tab[tab.set == name].sort_values('k')
			ax.plot(g.k, g[m], color=COL[name], lw=2.4 if name == 'mean' else 1.4, marker=MARK[name], ms=3.5,
			        label=LABEL[name], zorder=3 if name == 'mean' else 2)
			keep.append(pd.DataFrame({'metric': m, 'set': name, 'k': g.k, 'value': g[m], 'colour': COL[name]}))
		for p in ['p25', 'p50']:
			k = int(cd.loc[cd.set == 'mean', p].iloc[0])
			ax.axvline(k, color=MUTED, lw=0.8, ls='--', zorder=0)
			ax.text(k, 1, f' {p}', transform=ax.get_xaxis_transform(), color=MUTED, fontsize=7.5, va='top')
		ax.set_title(TITLE[m], loc='left')
		ax.set_xticks(range(2, KMAX + 1, 2))
		ax.set_xlabel('Número de HFTs')
	axs[0, 0].legend(fontsize=8, loc='lower right')
	fig.suptitle(f'Descriptores de los segmentos: medidas de los HFTs (umbral 0.{t}, Ward ponderado por área)',
	             x=0.05, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.05, 0.952, 'Línea negra: el método actual (media). Discontinuas: cortes p25 y p50 de su árbol',
	         color=INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.95))
	png = os.path.join(HERE, f'fig1_curvas_0{t}.png')
	fig.savefig(png, dpi=140)
	plt.close(fig)
	save_csv(pd.concat(keep), png)


def fig_summary(dif, ts):
	order = SETS[1:]
	ms = ['r2', 'jm_mean', 'jm_min', 'sil', 'r2_rob', 'sd_eta2', 'boot', 'ari_ref']
	fig, axs = plt.subplots(1, len(ms), figsize=(16, 4.6), sharey=True)
	keep = []
	for ax, m in zip(axs, ms):
		col = m if m == 'ari_ref' else 'd_' + m
		for i, name in enumerate(order):
			g = dif[dif.set == name]
			per = g.groupby('threshold')[col].median()
			med = g[col].median()
			ax.plot([per.min(), per.max()], [i, i], color=GRID, lw=3, solid_capstyle='round', zorder=1)
			ax.scatter(per.values, [i] * len(per), s=12, color=MUTED, zorder=2, linewidths=0)
			c = SLOTS[0] if (m == 'ari_ref' or med > 0) else SLOTS[1]
			ax.scatter([med], [i], s=46, color=c, zorder=3, edgecolors='#fcfcfb', linewidths=1.5)
			keep.append({'metric': m, 'set': name, 'median_all': med, 'colour': c,
			             **{f'median_0{k}': v for k, v in per.items()}})
		if m != 'ari_ref':
			ax.axvline(0, color=INK2, lw=0.9)
		ax.set_title(SHORT[m], loc='left', fontsize=8.5)
		ax.grid(axis='y', visible=False)
	axs[0].set_yticks(range(len(order)), [LABEL[o] for o in order])
	axs[0].invert_yaxis()
	fig.suptitle('Diferencia con la media (método actual) con el mismo nº de HFTs (2 a %d HFTs, umbrales %s)'
	             % (KMAX, ', '.join(f'0.{x}' for x in ts)), x=0.01, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.01, 0.9, 'Punto grande: mediana de todos los casos (azul, mayor que la media; naranja, menor). '
	         'Puntos grises y barra: mediana de cada umbral', color=INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.89))
	png = os.path.join(HERE, 'fig2_resumen.png')
	fig.savefig(png, dpi=140)
	plt.close(fig)
	save_csv(pd.DataFrame(keep), png)


def fig_cv_stab(s, names, tab_all, stab, t):
	fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
	keep = []
	ax = axs[0]
	cv = np.where(s.cv > 0, s.cv, np.nan)
	ax.boxplot([cv[:, j][np.isfinite(cv[:, j])] for j in range(s.nv)], orientation='horizontal', whis=(0, 100), widths=0.6,
	           medianprops={'color': SLOTS[1]})
	ax.set_xscale('log')
	ax.set_yticks(range(1, s.nv + 1), names)
	ax.invert_yaxis()
	ax.set_title(f'CV de los segmentos por variable (umbral 0.{t})', loc='left')
	ax.set_xlabel('CV = SD / |media| (escala log; bigotes: mín. y máx.)')
	for j, v in enumerate(names):
		keep.append(pd.DataFrame({'panel': 'cv', 'var': v, 'x': s.cv[:, j], 'n': s.n, 'colour': SLOTS[1]}))
	ax = axs[1]
	g = tab_all.groupby(['set', 'k'])['boot'].median().reset_index()
	for name in SETS:
		h = g[g.set == name]
		ax.plot(h.k, h.boot, color=COL[name], marker=MARK[name], ms=3.5, lw=2.4 if name == 'mean' else 1.4,
		        label=LABEL[name])
		keep.append(pd.DataFrame({'panel': 'boot', 'var': name, 'x': h.k, 'n': h.boot, 'colour': COL[name]}))
	ax.set_title('Estabilidad bootstrap de los segmentos (mediana de umbrales)', loc='left')
	ax.set_xlabel('Número de HFTs')
	ax.set_ylabel('ARI con la partición de todos')
	ax.legend(fontsize=7.5)
	ax = axs[2]
	for name in SETS:
		h = stab[stab.set == name]
		ax.plot(h.k, h.ari_thr, color=COL[name], marker=MARK[name], ms=3.5, lw=2.4 if name == 'mean' else 1.4)
		keep.append(pd.DataFrame({'panel': 'ari_thr', 'var': name, 'x': h.k, 'n': h.ari_thr, 'colour': COL[name]}))
	ax.set_title('Estabilidad entre umbrales de segmentación', loc='left')
	ax.set_xlabel('Número de HFTs')
	ax.set_ylabel('ARI de los píxeles, media de pares de umbrales')
	for a in axs[1:]:
		a.set_xticks(range(2, KMAX + 1, 2))
	fig.tight_layout()
	png = os.path.join(HERE, 'fig3_cv_estabilidad.png')
	fig.savefig(png, dpi=140)
	plt.close(fig)
	save_csv(pd.concat(keep), png)


if __name__ == '__main__':
	t0 = time.time()
	P, names, shape = read_pixels()
	print('Pixels read:', shape, 'variables:', names)
	keep, tabs, cds, wts, mms, cvs = {}, [], [], [], [], []
	for t in THR:
		s, tab, parts, cd, wt, mm, cv = run(t, P, shape, names)
		keep[t] = (s, tab, parts)
		tabs.append(tab)
		cds.append(cd)
		wts.append(wt)
		mms.append(mm)
		cvs += cv
		fig_curves(tab, t, cd)
		print(f'{time.time() - t0:.0f} s', flush=True)
	tab = pd.concat(tabs)
	dif, summ = summary(tab)
	if len(keep) < 2:
		sys.exit('Cross-threshold stability needs at least two thresholds')
	stab = thr_stability(keep)
	st = stab.groupby('set').agg(ari_thr=('ari_thr', 'median'), ari_thr_min=('ari_thr_min', 'median')).reindex(SETS)
	ref = stab[stab.set == 'mean'].set_index('k').ari_thr
	stab['d_ari_thr'] = stab.ari_thr.to_numpy() - ref.loc[stab.k].to_numpy()
	summ['d_ari_thr'] = summ.set.map(stab.groupby('set').d_ari_thr.median())
	summ['better_ari_thr'] = summ.set.map(stab.groupby('set').d_ari_thr.apply(lambda v: np.mean(v > 1e-9)))
	summ.to_csv(os.path.join(HERE, 'resumen.csv'), index=False)
	stab.to_csv(os.path.join(HERE, 'estabilidad_umbrales.csv'), index=False)
	cd = pd.concat(cds)
	cd.to_csv(os.path.join(HERE, 'cortes.csv'), index=False)
	pd.concat(wts).to_csv(os.path.join(HERE, 'pesos.csv'), index=False)
	pd.concat(mms).to_csv(os.path.join(HERE, 'media_vs_mediana.csv'), index=False)
	cvd = pd.DataFrame(cvs)
	cvd.to_csv(os.path.join(HERE, 'cv_por_variable.csv'), index=False)
	pd.set_option('display.width', 250)
	pd.set_option('display.max_columns', 40)
	print(f'\nMedian difference with mean at the same k and share of cases better, thresholds {THR}, k 2-{KMAX}:')
	print(summ.round(3).to_string(index=False))
	print('\nCross-threshold pixel ARI (median over k):')
	print(st.round(3).to_string())
	print('\nCuts (HFTs) per set and threshold, bootstrap share giving the same k, ARI vs mean at its own cut:')
	print(cd[['threshold', 'set', 'p10', 'p25', 'p50', 'p75', 'p25_boot_same', 'p50_boot_same', 'p25_boot_iqr',
	          'p50_boot_iqr', 'p25_ari_mean', 'p50_ari_mean']].round(2).to_string(index=False))
	print('\nWeights (sum of area-weighted variances) of the blocks:')
	print(pd.concat(wts).groupby('set')[['w_loc', 'w_disp', 'ratio']].median().round(3).to_string())
	print('\nMean vs median of the segments (median over thresholds):')
	print(pd.concat(mms).groupby('var').median(numeric_only=True).drop(columns='threshold').round(3).to_string())
	print('\nCV per variable:')
	print(cvd.groupby('var').median(numeric_only=True).drop(columns='threshold').round(3).to_string())
	fig_summary(dif, THR)
	s5 = keep[5][0] if 5 in keep else keep[THR[0]][0]
	fig_cv_stab(s5, names, tab, stab, 5 if 5 in keep else THR[0])
	print(f'\nDone in {time.time() - t0:.0f} s. Outputs in', HERE)
