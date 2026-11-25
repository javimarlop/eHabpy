"""Compara algoritmos de clasificación no supervisada de los segmentos en HFTs por la separabilidad de las clases.

    Rscript segment_sums.R            # una vez: sumas de los píxeles de cada segmento en data/
    python compare_algorithms.py [UMBRALES] [UMBRAL_MAPAS]

UMBRALES: índices de los umbrales separados por comas (por defecto los de data/); UMBRAL_MAPAS: el de las
figuras de detalle (por defecto 5). Necesita numpy, pandas, scipy, scikit-learn, matplotlib y artlib
(Fuzzy ART, pip install artlib).

Todos los algoritmos clasifican los segmentos por su media de cada variable tipificada con los píxeles
(como CLASS_METHOD = 'ward_area') y, salvo los que no admiten pesos, ponderan cada segmento por su nº de
píxeles. Algoritmos de Tabla_Algoritmos.xlsx probados:

  Con el nº de HFTs fijado (2 a KMAX):
    ward_area         referencia: AGNES con el criterio de Ward ponderado por área (el método actual)
    ward_area_kmeans  ward_area + k-medias ponderadas desde sus HFTs (opción actual de CLASS_METHOD)
    agnes_average     AGNES, enlace medio (UPGMA) ponderado por área
    agnes_complete    AGNES, enlace completo
    kmeans            k-medias ponderadas, k-means++ con 50 inicios (mejor de 5 semillas por la inercia)
    bisecting_kmeans  k-medias bisectivas ponderadas: árbol divisivo, cortes anidados
    som               SOM por lotes ponderado con k neuronas (malla lo más cuadrada posible): cada neurona un HFT
    hsom              SOM jerárquico: SOM de ~nº de segmentos / 3 neuronas y Ward ponderado sobre las neuronas
  Eligen ellos el nº de HFTs (se barren sus parámetros):
    dbscan            DBSCAN ponderado (min_samples en píxeles); el ruido pasa al HFT de media más próxima
    hdbscan           HDBSCAN (sin pesos); ídem
    fuzzy_art         Fuzzy ART (artlib), segmentos presentados al azar en proporción a su área
    dual_vigilance    Dual Vigilance Fuzzy ART (artlib; la versión no distribuida de DDVFA)
    gng               Growing Neural Gas: los HFTs son las componentes conexas del grafo
No probados: ART1 (solo datos binarios) y DDVFA (solo en Julia, AdaptiveResonance.jl).

Medidas (las mismas que compare_classification.R, con los píxeles de ../inVars tipificados):
  r2       varianza de los píxeles entre HFTs / total, media de las variables
  jm_mean  distancia de Jeffries-Matusita media entre pares de HFTs (0-2, normales multivariantes)
  jm_min   JM del par de HFTs menos separable
  sil      anchura de silueta media de los segmentos (distancias entre medias tipificadas), ponderada por área
  n_eff    nº efectivo de HFTs, exp(entropía de las fracciones de área): k si todos tienen la misma área.
           Aislar HFTs diminutos y extremos sube jm_mean y jm_min aunque casi toda el área quede en un HFT:
           esas particiones tienen un n_eff mucho menor que k y un R2 bajo
Resumen (resumen.csv): diferencia mediana con ward_area con el mismo nº de HFTs, % de casos mejores en cada
medida y % de casos mejores en todas a la vez (all4: R2, JM medio, JM mínimo y silueta).
  ari_ref  acuerdo (índice de Rand ajustado, por área) con ward_area con el mismo nº de HFTs
"""

import os
import sys
import time
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN, HDBSCAN, KMeans
from sklearn.metrics import silhouette_samples
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from artlib import DualVigilanceART, FuzzyART

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
KMAX = 20
SEEDS = range(5)
THR = [int(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else \
	sorted(int(f[6]) for f in os.listdir(DATA) if f.startswith('sums_0'))
TMAP = int(sys.argv[2]) if len(sys.argv) > 2 else 5
METRICS = ['r2', 'jm_mean', 'jm_min', 'sil', 'n_eff']
SEPAR = METRICS[:4]  # medidas de separabilidad (mayor = mejor)


# --- Datos y medidas ---

class Segments:
	"""Sumas de los píxeles tipificados de cada segmento de un umbral (segment_sums.R)."""
	def __init__(self, t):
		d = pd.read_csv(os.path.join(DATA, f'sums_0{t}.csv'), dtype={'segm_id': str})
		self.rows = np.flatnonzero(d['n'].to_numpy() > 0)  # filas de hri_results con píxeles
		d = d.iloc[self.rows]
		self.nv = sum(c[0] == 's' and c[1:].isdigit() for c in d.columns)
		qc = [c for c in d.columns if c.startswith('q_')]
		self.pairs = np.array([[int(x) - 1 for x in c.split('_')[1:]] for c in qc])
		self.diag = np.array([qc.index(f'q_{i}_{i}') for i in range(1, self.nv + 1)])
		self.n = d['n'].to_numpy(float)
		self.S = d[[f's{i}' for i in range(1, self.nv + 1)]].to_numpy()
		self.Q = d[qc].to_numpy()
		self.M = self.S / self.n[:, None]  # media tipificada de cada segmento
		self.D = squareform(pdist(self.M))
		self.N = self.n.sum()


def rowsum(x, cl, k):
	out = np.zeros((k,) + x.shape[1:])
	np.add.at(out, cl, x)
	return out


def measures(cl, s):
	"""R2, JM medio y mínimo y silueta de la partición cl de los segmentos s."""
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
	sil = silhouette_samples(s.D, cl, metric='precomputed') if 1 < k < len(cl) else np.zeros(len(cl))
	p = Nh / s.N
	return {'hfts': k, 'r2': r2, 'jm_mean': jm.mean(), 'jm_min': jm.min(), 'sil': np.sum(sil * s.n) / s.N, 'n_eff': np.exp(-np.sum(p * np.log(p)))}


def ari(x, y, w):
	"""Índice de Rand ajustado entre las particiones x e y de los segmentos, contando sus píxeles (w)."""
	_, x = np.unique(x, return_inverse=True)
	_, y = np.unique(y, return_inverse=True)
	T = np.zeros((x.max() + 1, y.max() + 1))
	np.add.at(T, (x, y), w)
	cb = lambda v: np.sum(v * (v - 1) / 2)
	e = cb(T.sum(1)) * cb(T.sum(0)) / cb(np.array([w.sum()]))
	return (cb(T) - e) / ((cb(T.sum(1)) + cb(T.sum(0))) / 2 - e)


# --- Algoritmos con el nº de HFTs fijado ---

def agnes(X, w, method):
	"""Uniones (i, j, altura) del agrupamiento aglomerativo de los puntos X con pesos (members) w, como
	hclust(..., members = w): 'ward' (ward.D2 sobre las distancias de Ward entre segmentos, como
	moranvar_plots.R), 'average', 'complete' o 'single' (distancias euclídeas)."""
	n = len(X)
	D = squareform(pdist(X))
	if method == 'ward':
		D = 2 * np.outer(w, w) / np.add.outer(w, w) * D ** 2  # 2 x incremento de la suma de cuadrados
	W = w.astype(float).copy()
	np.fill_diagonal(D, np.inf)
	merges = []
	for _ in range(n - 1):
		i, j = sorted(np.unravel_index(np.argmin(D), D.shape))
		h, di, dj = D[i, j], D[i].copy(), D[j].copy()
		if method == 'ward':
			new = ((W[i] + W) * di + (W[j] + W) * dj - W * h) / (W[i] + W[j] + W)
		elif method == 'average':
			new = (W[i] * di + W[j] * dj) / (W[i] + W[j])
		else:
			new = np.maximum(di, dj) if method == 'complete' else np.minimum(di, dj)
		D[i], D[:, i] = new, new
		D[j], D[:, j] = np.inf, np.inf
		D[i, i] = np.inf
		W[i] += W[j]
		merges.append((i, j, np.sqrt(h) if method == 'ward' else h))
	return merges


def cut_tree(merges, n, ks):
	"""Particiones de un árbol de agnes() con k grupos, para cada k de ks."""
	lab, out = np.arange(n), {}
	for m, (i, j, _) in enumerate(merges):
		lab[lab == j] = i
		if n - m - 1 in ks:
			out[n - m - 1] = np.unique(lab, return_inverse=True)[1]
	return out


def cut_height(merges, n, f):
	"""Nº de grupos del corte del árbol a la fracción f de su altura máxima (los cortes de SIMIL = 3)."""
	hmax = max(h for _, _, h in merges)
	return n - sum(h <= f * hmax for _, _, h in merges)


def wkmeans(X, w, cl, iter=100):
	"""k-medias ponderadas desde la partición cl (como wkmeans de moranvar_plots.R)."""
	for _ in range(iter):
		lev = np.unique(cl)
		cen = np.array([np.average(X[cl == c], axis=0, weights=w[cl == c]) for c in lev])
		new = lev[np.argmin(((X[:, None, :] - cen[None]) ** 2).sum(2), axis=1)]
		if np.all(new == cl):
			break
		cl = new
	return cl


def kmeans(X, w, k):
	"""k-medias ponderadas: la de menor inercia de varias semillas, cada una con 50 inicios k-means++."""
	fits = [KMeans(k, n_init=50, random_state=sd).fit(X, sample_weight=w) for sd in SEEDS]
	return min(fits, key=lambda f: f.inertia_).labels_, [f.labels_ for f in fits]


def bisecting(X, w, ks, seed=0):
	"""k-medias bisectivas ponderadas: parte en dos (k-medias, 10 inicios) el grupo de mayor suma de
	cuadrados ponderada, hasta max(ks) grupos; las particiones de cada k están anidadas."""
	lab, out = np.zeros(len(X), int), {}
	for k in range(2, max(ks) + 1):
		sse = [np.sum(w[lab == c] * ((X[lab == c] - np.average(X[lab == c], axis=0, weights=w[lab == c])) ** 2).sum(1))
		       if np.sum(lab == c) > 1 else -1 for c in range(k - 1)]
		m = np.flatnonzero(lab == int(np.argmax(sse)))
		km = KMeans(2, n_init=10, random_state=seed).fit(X[m], sample_weight=w[m])
		lab[m[km.labels_ == 1]] = k - 1
		out[k] = lab.copy()
	return out


def som(X, w, rows, cols, s_end, iters=60):
	"""SOM por lotes (Kohonen) ponderado por w, malla rows x cols, vecindad gaussiana de radio de max(rows, cols)/2
	a s_end, iniciado sobre los dos primeros ejes de la PCA ponderada (determinista). Prototipos y neurona
	ganadora de cada punto."""
	g = np.array([(r, c) for r in range(rows) for c in range(cols)], float)
	G2 = ((g[:, None] - g[None]) ** 2).sum(2)
	mu = np.average(X, axis=0, weights=w)
	val, vec = np.linalg.eigh(np.cov((X - mu).T, aweights=w))
	span = lambda m: np.linspace(-1, 1, m) if m > 1 else np.zeros(1)
	P = mu + np.outer(span(cols)[g[:, 1].astype(int)], 2 * np.sqrt(val[-1]) * vec[:, -1]) + \
		np.outer(span(rows)[g[:, 0].astype(int)], 2 * np.sqrt(val[-2]) * vec[:, -2])
	s0 = max(max(rows, cols) / 2, s_end)
	for t in range(iters):
		sg = s0 * (s_end / s0) ** (t / (iters - 1))
		bmu = np.argmin(((X[:, None] - P[None]) ** 2).sum(2), axis=1)
		H = np.exp(-G2[bmu] / (2 * sg ** 2)) * w[:, None]
		hit = H.sum(0) > 0
		P[hit] = (H.T @ X)[hit] / H.sum(0)[hit, None]
	return P, np.argmin(((X[:, None] - P[None]) ** 2).sum(2), axis=1)


def grid_shape(k):
	r = max(d for d in range(1, int(np.sqrt(k)) + 1) if k % d == 0)
	return r, k // r


def som_k(X, w, k):
	"""SOM de k neuronas: cada neurona con segmentos es un HFT (radio final 0.3, casi k-medias)."""
	return np.unique(som(X, w, *grid_shape(k), s_end=0.3)[1], return_inverse=True)[1]


def hsom(X, w, ks):
	"""SOM jerárquico: SOM de unas nº de puntos / 3 neuronas (radio final 1) y Ward ponderado por el área de las
	neuronas sobre sus prototipos; cada segmento va al HFT de su neurona."""
	side = max(2, int(round(np.sqrt(len(X) / 3))))
	P, bmu = som(X, w, side, side, s_end=1.0)
	used, node = np.unique(bmu, return_inverse=True)
	wn = np.bincount(node, w)
	parts = cut_tree(agnes(P[used], wn, 'ward'), len(used), [k for k in ks if k <= len(used)])
	return {k: p[node] for k, p in parts.items()}


# --- Algoritmos que eligen el nº de HFTs ---

def assign_noise(X, w, lab):
	"""Los puntos sin grupo (-1) pasan al grupo de media ponderada más próxima."""
	lev = np.unique(lab[lab >= 0])
	if len(lev) == 0:
		return None
	cen = np.array([np.average(X[lab == c], axis=0, weights=w[lab == c]) for c in lev])
	lab = lab.copy()
	no = lab < 0
	lab[no] = lev[np.argmin(((X[no, None] - cen[None]) ** 2).sum(2), axis=1)]
	return np.unique(lab, return_inverse=True)[1]


def dbscan_runs(s):
	for frac in [0.005, 0.01, 0.02, 0.05, 0.1]:
		for q in np.linspace(0.02, 0.5, 25):
			eps = np.quantile(s.D[np.triu_indices(len(s.n), 1)], q)
			lab = DBSCAN(eps=eps, min_samples=int(frac * s.N), metric='precomputed').fit(s.D, sample_weight=s.n).labels_
			yield f'min_samples={frac:g}N eps=q{q:.2f}', lab, np.sum(s.n[lab < 0]) / s.N


def hdbscan_runs(s):
	for mcs in range(2, 31):
		for ms in [1, 3, 5]:
			# HDBSCAN cambia la diagonal de la matriz de distancias que recibe
			lab = HDBSCAN(min_cluster_size=mcs, min_samples=ms, metric='precomputed').fit(s.D.copy()).labels_
			yield f'min_cluster_size={mcs} min_samples={ms}', lab, np.sum(s.n[lab < 0]) / s.N


def art_runs(s, dual):
	for sd in SEEDS:
		order = np.random.default_rng(sd).choice(len(s.n), 3 * len(s.n), p=s.n / s.N)
		for rho in np.round(np.arange(0.30, 0.951, 0.025), 3):
			for lower in ([rho - 0.1, rho - 0.2] if dual else [None]):
				if lower is not None and lower <= 0:
					continue
				base = FuzzyART(rho=rho, alpha=1e-3, beta=1.0)
				art = DualVigilanceART(base, rho_lower_bound=lower) if dual else base
				Xp = art.prepare_data(s.M)
				art.fit(Xp[order])
				par = f'rho={rho:g}' + (f' rho_lower={lower:g}' if dual else '') + f' seed={sd}'
				yield par, art.predict(Xp), 0.0


def gng(X, w, max_nodes, age_max, seed, n_iter=20000, eps_b=0.05, eps_n=0.006, lam=100, alpha=0.5, dec=0.995):
	"""Growing Neural Gas (Fritzke, 1995) con los puntos presentados al azar en proporción a w. Los grupos son las
	componentes conexas del grafo final; cada punto va al de su nodo más próximo."""
	rng = np.random.default_rng(seed)
	sample = X[rng.choice(len(X), n_iter, p=w / w.sum())]
	P = np.zeros((max_nodes, X.shape[1]))
	P[:2] = X[rng.choice(len(X), 2, replace=False)]
	alive = np.zeros(max_nodes, bool)
	alive[:2] = True
	err = np.zeros(max_nodes)
	age = np.full((max_nodes, max_nodes), -1)  # -1: sin arista
	for t, x in enumerate(sample, 1):
		d2 = np.where(alive, ((P - x) ** 2).sum(1), np.inf)
		s1, s2 = np.argsort(d2)[:2]
		nb = age[s1] >= 0
		age[s1, nb] += 1
		age[nb, s1] += 1
		err[s1] += d2[s1]
		P[s1] += eps_b * (x - P[s1])
		P[nb] += eps_n * (x - P[nb])
		age[s1, s2] = age[s2, s1] = 0
		age[age > age_max] = -1
		alive &= (age >= 0).any(1)
		if t % lam == 0 and alive.sum() < max_nodes:
			q = np.argmax(np.where(alive, err, -1))
			nq = np.flatnonzero(age[q] >= 0)
			f = nq[np.argmax(err[nq])]
			r = np.flatnonzero(~alive)[0]
			P[r] = (P[q] + P[f]) / 2
			age[r] = age[:, r] = -1
			age[q, f] = age[f, q] = -1
			age[q, r] = age[r, q] = age[r, f] = age[f, r] = 0
			alive[r] = True
			err[q] *= alpha
			err[f] *= alpha
			err[r] = err[q]
		err *= dec
	idx = np.flatnonzero(alive)
	_, comp = connected_components(age[np.ix_(idx, idx)] >= 0, directed=False)
	node = idx[np.argmin(((X[:, None] - P[idx][None]) ** 2).sum(2), axis=1)]
	return comp[np.searchsorted(idx, node)]


def gng_runs(s):
	for mn in [10, 20, 40, 80]:
		for am in [20, 50, 100]:
			for sd in SEEDS[:3]:
				yield f'max_nodes={mn} age_max={am} seed={sd}', gng(s.M, s.n, mn, am, sd), 0.0


FIXED = ['ward_area', 'ward_area_kmeans', 'agnes_average', 'agnes_complete', 'kmeans', 'bisecting_kmeans', 'som', 'hsom']
AUTO = {'dbscan': dbscan_runs, 'hdbscan': hdbscan_runs, 'fuzzy_art': lambda s: art_runs(s, False),
        'dual_vigilance': lambda s: art_runs(s, True), 'gng': gng_runs}
LABEL = {'ward_area': 'Ward ponderado (actual, ward_area)', 'ward_area_kmeans': 'Ward ponderado + k-medias',
         'agnes_average': 'AGNES, enlace medio', 'agnes_complete': 'AGNES, enlace completo',
         'agnes_single': 'AGNES, enlace simple', 'kmeans': 'k-medias ponderadas',
         'bisecting_kmeans': 'k-medias bisectivas', 'som': 'SOM (una neurona por HFT)', 'hsom': 'SOM jerárquico (HSOM)',
         'dbscan': 'DBSCAN ponderado', 'hdbscan': 'HDBSCAN', 'fuzzy_art': 'Fuzzy ART',
         'dual_vigilance': 'Dual Vigilance Fuzzy ART', 'gng': 'Growing Neural Gas'}


# --- Comparación ---

def run(t):
	s = Segments(t)
	ns, ks = len(s.n), list(range(2, KMAX + 1))
	print(f'\n=== Threshold 0.{t}: {ns} segments, {s.N:.0f} pixels ===')
	tree = agnes(s.M, s.n, 'ward')
	ref = cut_tree(tree, ns, ks)
	cuts = {p: cut_height(tree, ns, f) for p, f in [('p10', .1), ('p25', .25), ('p50', .5), ('p75', .75), ('p90', .9)]}
	print('HFTs of the cuts of ward_area:', cuts)
	parts = {('ward_area', k): ref[k] for k in ks}
	t0 = time.time()
	for k in ks:
		parts[('ward_area_kmeans', k)] = wkmeans(s.M, s.n, ref[k])
		parts[('kmeans', k)], runs = kmeans(s.M, s.n, k)
		parts[('kmeans_seeds', k)] = runs
		parts[('som', k)] = som_k(s.M, s.n, k)
	for mth in ['average', 'complete', 'single']:
		for k, p in cut_tree(agnes(s.M, s.n, mth), ns, ks).items():
			parts[('agnes_' + mth, k)] = p
	for k, p in bisecting(s.M, s.n, ks).items():
		parts[('bisecting_kmeans', k)] = p
	for k, p in hsom(s.M, s.n, ks).items():
		parts[('hsom', k)] = p
	rows = []
	for (mth, k), p in parts.items():
		if mth == 'kmeans_seeds':
			continue
		rows.append({'threshold': t, 'method': mth, 'params': '', 'k': k, **measures(p, s),
		             'ari_ref': ari(p, ref[k], s.n), 'noise': 0.0})
	print(f'fixed-k methods: {time.time() - t0:.0f} s')
	for mth, gen in AUTO.items():
		t0 = time.time()
		for par, lab, noise in gen(s):
			lab = assign_noise(s.M, s.n, lab)
			if lab is None or len(np.unique(lab)) < 2:
				continue
			m = measures(lab, s)
			r = {'threshold': t, 'method': mth, 'params': par, 'k': m['hfts'], **m, 'noise': noise}
			r['ari_ref'] = ari(lab, ref[m['hfts']], s.n) if m['hfts'] in ref else np.nan
			rows.append(r)
		print(f'{mth}: {time.time() - t0:.0f} s')
	tab = pd.DataFrame(rows)
	tab.to_csv(os.path.join(HERE, f'medidas_0{t}.csv'), index=False)

	# Estabilidad de k-medias entre semillas (ARI medio por área entre las 5 soluciones de cada k)
	stab = [np.mean([ari(a, b, s.n) for a, b in combinations(parts[('kmeans_seeds', k)], 2)]) for k in ks]
	print('k-means: mean ARI between seeds, k = 2..%d: %s' % (KMAX, np.round(stab, 2)))
	return s, tab, parts, cuts


def compare(tab):
	"""Diferencia de cada medida con ward_area con el mismo nº de HFTs (k <= KMAX)."""
	ref = tab[tab.method == 'ward_area'].set_index(['threshold', 'k'])[METRICS]
	x = tab[(tab.method != 'ward_area') & (tab.k <= KMAX)].copy()
	r = ref.loc[list(zip(x.threshold, x.k))].to_numpy()
	for i, m in enumerate(METRICS):
		x['d_' + m] = x[m].to_numpy() - r[:, i]
	return x


def summary(dif):
	out = []
	for mth, g in dif.groupby('method', sort=False):
		row = {'method': mth, 'runs': len(g), 'k_range': f'{g.k.min()}-{g.k.max()}', 'ari_ref': g.ari_ref.median()}
		for m in METRICS:
			row['d_' + m] = g['d_' + m].median()
		for m in SEPAR:
			row['better_' + m] = np.mean(g['d_' + m] > 1e-9)
		row['better_all4'] = np.mean(np.all(g[['d_' + m for m in SEPAR]].to_numpy() > 1e-9, axis=1))
		out.append(row)
	return pd.DataFrame(out)


# --- Figuras ---

INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
COL = dict(zip(FIXED[1:], SLOTS), ward_area=INK)
COL.update(dict(zip(AUTO, SLOTS)))
MARK = dict(zip(FIXED, ['o', 's', '^', 'v', 'D', 'P', 'X', '*']))
MARK.update(dict(zip(AUTO, ['o', 's', '^', 'D', 'X'])))
TITLE = {'r2': 'R² de los píxeles (media de las variables)', 'jm_mean': 'JM medio entre pares de HFTs (máx. 2)',
         'jm_min': 'JM mínimo: el par de HFTs menos separable',
         'sil': 'Silueta media de los segmentos (por área)', 'n_eff': 'Nº efectivo de HFTs (por su área)'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#c3c2b7', 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.labelcolor': INK2, 'text.color': INK, 'axes.titleweight': 'bold', 'axes.titlesize': 10,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.6, 'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb',
                     'savefig.facecolor': '#fcfcfb', 'legend.frameon': False})


def mark_cuts(ax, cuts):
	for p in ['p25', 'p50']:
		ax.axvline(cuts[p], color=MUTED, lw=0.8, ls='--', zorder=0)
		ax.text(cuts[p], 1, f' {p}', transform=ax.get_xaxis_transform(), color=MUTED, fontsize=7.5, va='top')


def fig_fixed(tab, t, ns, cuts):
	fig, axs = plt.subplots(3, 2, figsize=(11, 11.5))
	for ax, m in zip(axs.flat, METRICS):
		for mth in FIXED:
			g = tab[tab.method == mth].sort_values('k')
			ax.plot(g.k, g[m], color=COL[mth], lw=2.4 if mth == 'ward_area' else 1.5, marker=MARK[mth], ms=4,
			        label=LABEL[mth], zorder=3 if mth == 'ward_area' else 2)
		mark_cuts(ax, cuts)
		ax.set_title(TITLE[m], loc='left')
		ax.set_xticks(range(2, KMAX + 1, 2))
		ax.set_xlabel('Número de HFTs')
	axs[2, 1].set_axis_off()
	axs[2, 1].legend(*axs[0, 0].get_legend_handles_labels(), loc='center left', fontsize=9.5, title='Algoritmo',
	                 alignment='left', title_fontproperties={'weight': 'bold'})
	fig.suptitle(f'Separabilidad de los HFTs con el nº de HFTs fijado (umbral 0.{t}, {ns} segmentos)',
	             x=0.06, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.06, 0.952, 'Línea negra: el método actual. Líneas discontinuas: cortes p25 y p50 de su árbol',
	         color=INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.95))
	fig.savefig(os.path.join(HERE, f'fig1_k_fijo_0{t}.png'), dpi=150)
	plt.close(fig)


def fig_auto(tab, t, ns, cuts):
	fig, axs = plt.subplots(3, 2, figsize=(11, 11.5))
	ref = tab[tab.method == 'ward_area'].sort_values('k')
	for ax, m in zip(axs.flat, METRICS):
		ax.plot(ref.k, ref[m], color=INK, lw=2.4, label=LABEL['ward_area'], zorder=3)
		for mth in AUTO:
			g = tab[(tab.method == mth) & (tab.k <= KMAX)]
			ax.scatter(g.k + (list(AUTO).index(mth) - 2) * 0.12, g[m], s=16, color=COL[mth], marker=MARK[mth],
			           alpha=0.55, linewidths=0, label=LABEL[mth], zorder=2)
		mark_cuts(ax, cuts)
		ax.set_title(TITLE[m], loc='left')
		ax.set_xticks(range(2, KMAX + 1, 2))
		ax.set_xlabel('Número de HFTs que da el algoritmo')
	axs[2, 1].set_axis_off()
	axs[2, 1].legend(*axs[0, 0].get_legend_handles_labels(), loc='center left', fontsize=9.5, markerscale=1.8,
	                 title='Algoritmo', alignment='left', title_fontproperties={'weight': 'bold'})
	fig.suptitle(f'Algoritmos que eligen el nº de HFTs: cada punto, una combinación de parámetros (umbral 0.{t}, {ns} segmentos)',
	             x=0.06, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.06, 0.952, 'Comparados con el método actual con el mismo nº de HFTs (línea negra)', color=INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.95))
	fig.savefig(os.path.join(HERE, f'fig2_k_automatico_0{t}.png'), dpi=150)
	plt.close(fig)


def fig_summary(dif):
	"""Diferencia mediana con ward_area con el mismo nº de HFTs, por método y umbral."""
	order = [m for m in FIXED[1:] + ['agnes_single'] + list(AUTO) if m in set(dif.method)]
	fig, axs = plt.subplots(1, 5, figsize=(14, 5.4), sharey=True, gridspec_kw={'width_ratios': [1, 1, 1, 1, 0.7]})
	thr = sorted(dif.threshold.unique())
	for ax, m in zip(axs, SEPAR):
		for i, mth in enumerate(order):
			g = dif[dif.method == mth]
			per = g.groupby('threshold')['d_' + m].median()
			ax.plot([per.min(), per.max()], [i, i], color=GRID, lw=3, solid_capstyle='round', zorder=1)
			ax.scatter(per.values, [i] * len(per), s=12, color=MUTED, zorder=2, linewidths=0)
			med = g['d_' + m].median()
			ax.scatter([med], [i], s=46, color=SLOTS[0] if med > 0 else SLOTS[1], zorder=3, edgecolors='#fcfcfb',
			           linewidths=1.5)
		ax.axvline(0, color=INK2, lw=0.9)
		ax.set_title(TITLE[m].split(':')[0].split(' (')[0], loc='left', fontsize=9.5)
		ax.grid(axis='y', visible=False)
	# casos en que el método mejora las cuatro medidas a la vez
	all4 = [100 * np.mean(np.all(dif[dif.method == o][['d_' + m for m in SEPAR]].to_numpy() > 1e-9, axis=1)) for o in order]
	axs[4].barh(range(len(order)), all4, height=0.55, color=SLOTS[0])
	for i, v in enumerate(all4):
		axs[4].text(v + 1, i, f'{v:.0f} %', va='center', fontsize=8, color=INK2)
	axs[4].set_xlim(0, 100)
	axs[4].set_title('Mejor en las 4 a la vez', loc='left', fontsize=9.5)
	axs[4].grid(axis='y', visible=False)
	axs[0].set_yticks(range(len(order)), [LABEL[o] for o in order])
	axs[0].invert_yaxis()
	fig.suptitle('Diferencia con el método actual con el mismo nº de HFTs (2 a %d HFTs, umbrales %s)'
	             % (KMAX, ', '.join(f'0.{x}' for x in thr)), x=0.02, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.02, 0.905, 'Punto grande: mediana de todos los casos (azul, mejor que el actual; naranja, peor). '
	         'Puntos grises y barra: mediana de cada umbral', color=INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.9))
	fig.savefig(os.path.join(HERE, 'fig3_resumen.png'), dpi=150)
	plt.close(fig)


def fig_maps(s, tab, parts, t, k, methods):
	"""Mapas de los HFTs de varios métodos con k HFTs; los HFTs de cada uno se numeran como los de ward_area
	con los que comparten más área, para compararlos."""
	nr, nc = pd.read_csv(os.path.join(DATA, 'grid.csv')).iloc[0]
	seg = np.fromfile(os.path.join(DATA, f'seg_0{t}.bin'), dtype=np.int32).reshape(nr, nc)
	rowmap = np.full(seg.max() + 1, -1)
	rowmap[s.rows + 1] = np.arange(len(s.rows))
	ok = np.argwhere(seg > 0)
	(r0, c0), (r1, c1) = ok.min(0), ok.max(0) + 1
	seg = seg[r0:r1, c0:c1]
	ref = parts[('ward_area', k)]
	cols = SLOTS + ['#86b6ef', '#f4a582', '#7fd8b4', '#f6cd73', '#9085e9', '#c3c2b7']
	cmap = ListedColormap(cols[:max(k, 2)] if k <= len(cols) else plt.get_cmap('tab20').colors[:k])
	fig, axs = plt.subplots(2, 2, figsize=(11.5, 9.2))
	for ax, mth in zip(axs.flat, methods):
		p = parts[(mth, k)]
		# renumeración: emparejamiento voraz por área compartida con los HFTs de ward_area
		T = np.zeros((p.max() + 1, ref.max() + 1))
		np.add.at(T, (p, ref), s.n)
		lab, free, used = np.full(p.max() + 1, -1), set(range(k)), set()
		for i, j in sorted(np.ndindex(T.shape), key=lambda ij: -T[ij]):
			if lab[i] < 0 and j not in used:
				lab[i] = j
				used.add(j)
		for i in np.flatnonzero(lab < 0):
			lab[i] = min(free - used)
			used.add(lab[i])
		img = np.where(seg > 0, lab[p[np.maximum(rowmap[seg], 0)]], -1).astype(float)
		img[(seg == 0) | (rowmap[seg] < 0)] = np.nan
		ax.imshow(img, cmap=cmap, vmin=-0.5, vmax=cmap.N - 0.5, interpolation='nearest')
		r = tab[(tab.method == mth) & (tab.k == k) & (tab.params == '')].iloc[0]
		ax.set_title(f'{LABEL[mth]}\nR² {r.r2:.2f} · JM medio {r.jm_mean:.2f}, mín. {r.jm_min:.2f} · silueta {r.sil:.2f}'
		             + ('' if mth == 'ward_area' else f' · ARI {r.ari_ref:.2f}'), fontsize=8.5, loc='left')
		ax.set_axis_off()
	fig.suptitle(f'{k} HFTs con cada método (umbral 0.{t})', x=0.01, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.01, 0.935, 'Cada HFT tiene el color del HFT del método actual con el que comparte más área. '
	         'ARI: acuerdo con el método actual (1 = el mismo mapa)', color=INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.93), h_pad=3)
	fig.savefig(os.path.join(HERE, f'fig4_mapas_0{t}_k{k}.png'), dpi=130)
	plt.close(fig)


if __name__ == '__main__':
	tabs, keep = [], {}
	for t in THR:
		s, tab, parts, cuts = run(t)
		tabs.append(tab)
		keep[t] = (s, tab, parts, cuts)
		fig_fixed(tab, t, len(s.n), cuts)
		fig_auto(tab, t, len(s.n), cuts)
	tab = pd.concat(tabs)
	dif = compare(tab)
	summ = summary(dif)
	summ.to_csv(os.path.join(HERE, 'resumen.csv'), index=False)
	pd.set_option('display.width', 200)
	print('\nMedian difference with ward_area at the same number of HFTs and share of cases better than it,')
	print(f'thresholds {THR}, 2 to {KMAX} HFTs:')
	print(summ.round(3).to_string(index=False))
	fig_summary(dif)
	if TMAP in keep:
		s, tab, parts, cuts = keep[TMAP]
		k = cuts['p25']
		fig_maps(s, tab, parts, TMAP, k, ['ward_area', 'ward_area_kmeans', 'kmeans', 'hsom'])
		sel = tab[tab.k.isin([cuts['p50'], k]) & (tab.params == '')]
		print(f'\nThreshold 0.{TMAP}, cuts p50 and p25 of ward_area ({cuts["p50"]} and {k} HFTs):')
		print(sel[['method', 'k'] + METRICS + ['ari_ref']].round(3).to_string(index=False))
	print('\nFiguras y tablas en', HERE)
