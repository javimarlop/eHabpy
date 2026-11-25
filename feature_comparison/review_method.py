"""Revisión del método de compare_features.py: ¿cambian sus conclusiones con medidas que no favorecen a la media
por construcción, con otros pesos del bloque de dispersión o con otra forma de la SD?

    /opt/anaconda3/envs/ehabpy/bin/python -s review_method.py

Reutiliza las funciones de compare_features.py (mismos datos, Ward ponderado, bootstrap con las mismas semillas).
Añade:
  gmi      información mutua gaussiana entre el HFT y los 7 valores del píxel (nats):
           0.5 (logdet Sigma_total - sum_h p_h logdet Sigma_h). Es el criterio de verosimilitud de una mezcla de
           normales con covarianza propia de cada HFT: premia separar medias Y agrupar segmentos de heterogeneidad
           parecida (un HFT homogéneo tiene un determinante pequeño). No favorece a la media por construcción.
  mi_bin   información mutua no paramétrica HFT-variable (64 clases por cuantiles de los píxeles), dividida por la
           entropía de la variable, media de las 7: premia cualquier diferencia de distribución (posición,
           dispersión, forma) entre HFTs.
Conjuntos nuevos (bloque de dispersión centrado y escalado para que su peso, suma de varianzas ponderadas, sea r
veces el de las medias):
  mean+sd_rX       SD re-tipificada (como mean+sd) con r = X
  mean+sdnat_rX    SD en unidades de píxel sin re-tipificar (pesos relativos naturales entre variables), r = X
  mean+logsd_rX    log(SD) re-tipificado, r = X
  mean+logcv4_r05  log(CV) solo de las variables de escala de razón positivas (NDVImax, pendiente, precip., TWI)
  median_pipe      mediana tipificada como std_stats() con STAT_AGG = 'median' (SD de mediana + varianza)
Escribe review_medidas.csv, review_resumen.csv, review_bloque_sd.csv y review_boot_semilla.csv en esta carpeta.
"""

import os
import sys
import time
from itertools import combinations

sys.argv = sys.argv[:1]  # compare_features lee los umbrales de sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_features as cf  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402

THR = [3, 4, 5, 6, 7]
KS = cf.KS
NB = 64
NBOOT = 30
POS4 = ['NDVImax', 'Pendiente', 'PrecipAnual', 'TWI']


def wvar(X, w):
	return cf.wstd(X, w)[1] ** 2


def features(s, name, names):
	"""Descriptores de los conjuntos (los de compare_features.py y los nuevos)."""
	if name in cf.SETS:
		return s.features(name)[0]
	if name == 'median_pipe':
		mu = np.average(s.median, axis=0, weights=s.n)
		sp = np.sqrt((np.sum(s.n[:, None] * s.sd ** 2, 0) + np.sum(s.n[:, None] * (s.median - mu) ** 2, 0)) / s.N)
		return (s.median - mu) / sp
	loc = s.mean
	wl = wvar(loc, s.n).sum()
	kind, r = name.split('+')[1].split('_r')
	r = float(r.replace('0', '0.', 1)) if r.startswith('0') else float(r)
	if kind in ('sd', 'sdnat'):
		d = s.sd.copy()
	elif kind == 'logsd':
		d = np.log(np.maximum(s.sd, np.quantile(s.sd[s.sd > 0], 0.01)))
	elif kind == 'logcv4':
		j = [names.index(v) for v in POS4]
		d = np.log(s.cv[:, j])
	d = d - np.average(d, axis=0, weights=s.n)
	if kind != 'sdnat':
		d = d / np.sqrt(wvar(d, s.n))
	d = d * np.sqrt(r * wl / wvar(d, s.n).sum())
	return np.hstack([loc, d])


class Extra:
	"""Histogramas de cada segmento por clases de cuantiles de cada variable (para mi_bin) y logdet total."""
	def __init__(self, s):
		self.H = np.zeros((s.ns, NB, s.nv))
		self.hx = []
		for j in range(s.nv):
			e = np.quantile(s.Z[:, j], np.linspace(0, 1, NB + 1)[1:-1])
			b = np.searchsorted(e, s.Z[:, j], side='right')
			self.H[:, :, j] = np.bincount(s.g * NB + b, minlength=s.ns * NB).reshape(s.ns, NB)
			p = self.H[:, :, j].sum(0) / s.N
			p = p[p > 0]
			self.hx.append(-np.sum(p * np.log(p)))
		self.hx = np.array(self.hx)
		C = np.zeros((s.nv, s.nv))
		Qt = s.Q.sum(0)
		C[s.pairs[:, 0], s.pairs[:, 1]] = Qt
		C[s.pairs[:, 1], s.pairs[:, 0]] = Qt
		mu = s.S.sum(0) / s.N
		self.ld_tot = np.linalg.slogdet(C / s.N - np.outer(mu, mu) + 1e-6 * np.eye(s.nv))[1]


def measures(cl, s, ex):
	_, cl = np.unique(cl, return_inverse=True)
	k = cl.max() + 1
	Nh, Sh, Qh = cf.rowsum(s.n, cl, k), cf.rowsum(s.S, cl, k), cf.rowsum(s.Q, cl, k)
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
	sil = cf.silhouette(s.D, cl)
	gmi = 0.5 * (ex.ld_tot - np.sum(p * ld))
	Hh = cf.rowsum(ex.H, cl, k) / s.N  # k x NB x nv
	px = Hh.sum(0, keepdims=True)
	ph = Hh.sum(1, keepdims=True)
	with np.errstate(divide='ignore', invalid='ignore'):
		t = np.where(Hh > 0, Hh * np.log(Hh / (ph * px)), 0)
	mi = t.sum((0, 1)) / ex.hx
	msd = cf.rowsum(s.n[:, None] * s.sd, cl, k) / Nh[:, None]
	gsd = np.average(s.sd, axis=0, weights=s.n)
	tsd = (s.n[:, None] * (s.sd - gsd) ** 2).sum(0)
	sd_eta2 = np.mean((Nh[:, None] * (msd - gsd) ** 2).sum(0) / tsd)
	return {'r2': r2, 'jm_mean': jm.mean(), 'jm_min': jm.min(), 'sil': np.sum(sil * s.n) / s.N,
	        'n_eff': np.exp(-np.sum(p * np.log(p))), 'gmi': gmi, 'mi_bin': mi.mean(), 'sd_eta2': sd_eta2}


SETS = ['mean', 'median', 'median_pipe', 'mean+sd', 'mean+sd_half', 'mean+sd_pix',
        'mean+sd_r01', 'mean+sd_r025', 'mean+sd_r1',
        'mean+sdnat_r025', 'mean+sdnat_r05', 'mean+sdnat_r1',
        'mean+logsd_r025', 'mean+logsd_r05', 'mean+logsd_r2', 'mean+logcv4_r05']
METRICS = ['r2', 'jm_mean', 'jm_min', 'sil', 'n_eff', 'gmi', 'mi_bin', 'sd_eta2', 'boot']


def wcorr(x, y, w):
	c = np.cov(x, y, aweights=w)
	return c[0, 1] / np.sqrt(c[0, 0] * c[1, 1])


if __name__ == '__main__':
	t0 = time.time()
	P, names, shape = cf.read_pixels()
	rows, keep, blk, bs = [], {}, [], []
	for t in THR:
		s = cf.Segments(t, P, shape)
		ex = Extra(s)
		parts = {}
		# estructura del bloque de SD: correlación con el tamaño del segmento y entre las SD
		Zsd = (s.sd - np.average(s.sd, 0, weights=s.n)) / cf.wstd(s.sd, s.n)[1]
		R = np.cov(Zsd.T, aweights=s.n)
		R = R / np.sqrt(np.outer(np.diag(R), np.diag(R)))
		ev = np.sort(np.linalg.eigvalsh(R))[::-1]
		lognz = np.log(s.n)
		for j, v in enumerate(names):
			blk.append({'threshold': t, 'var': v, 'r_sd_logn': wcorr(s.sd[:, j], lognz, s.n),
			            'r_sd_logn_unw': np.corrcoef(s.sd[:, j], lognz)[0, 1],
			            'r_sd_absmean': wcorr(s.sd[:, j], np.abs(s.mean[:, j]), s.n),
			            'pc1_share_sdblock': ev[0] / ev.sum(), 'mean_abs_r_between_sds':
			            np.abs(R[np.triu_indices(s.nv, 1)]).mean()})
		for name in SETS:
			F = features(s, name, names)
			tree = cf.ward(F, s.n)
			pk = cf.cut_tree(tree, s.ns, KS)
			boot, _ = cf.bootstrap(F, s, pk, NBOOT, seed=1000 + t)
			p25, p50 = cf.cut_height(tree, s.ns, .25), cf.cut_height(tree, s.ns, .5)
			for k in KS:
				parts[(name, k)] = pk[k]
				rows.append({'threshold': t, 'set': name, 'k': k, **measures(pk[k], s, ex), 'boot': boot[k],
				             'ari_ref': cf.ari(pk[k], parts[('mean', k)], s.n), 'p25': p25, 'p50': p50})
			# bootstrap con otra semilla y más réplicas, solo media y mediana
			if name in ('mean', 'median', 'median_pipe'):
				b2, _ = cf.bootstrap(F, s, pk, 100, seed=7000 + t)
				bs += [{'threshold': t, 'set': name, 'k': k, 'boot100': b2[k]} for k in KS]
		keep[t] = (s, parts)
		print(f'threshold 0.{t}: {time.time() - t0:.0f} s', flush=True)
	tab = pd.DataFrame(rows)
	# estabilidad entre umbrales
	st = []
	for name in SETS:
		for k in KS:
			v = []
			for a, b in combinations(THR, 2):
				sa, sb = keep[a][0], keep[b][0]
				_, ia, ib = np.intersect1d(sa.pix, sb.pix, assume_unique=True, return_indices=True)
				v.append(cf.ari_pix(keep[a][1][(name, k)][sa.g[ia]], keep[b][1][(name, k)][sb.g[ib]]))
			st.append({'set': name, 'k': k, 'ari_thr': np.mean(v)})
	st = pd.DataFrame(st)
	tab.to_csv(os.path.join(HERE, 'review_medidas.csv'), index=False)
	pd.DataFrame(blk).to_csv(os.path.join(HERE, 'review_bloque_sd.csv'), index=False)
	bs = pd.DataFrame(bs)
	bs.to_csv(os.path.join(HERE, 'review_boot_semilla.csv'), index=False)
	# resumen: diferencia mediana con 'mean' al mismo umbral y k, todos los k y k 3-10 (cortes p25-p50 de 'mean')
	ref = tab[tab.set == 'mean'].set_index(['threshold', 'k'])[METRICS]
	out = []
	for name, g in tab[tab.set != 'mean'].groupby('set', sort=False):
		r = ref.loc[list(zip(g.threshold, g.k))].to_numpy()
		d = pd.DataFrame(g[METRICS].to_numpy() - r, columns=METRICS)
		d['threshold'], d['k'] = g.threshold.to_numpy(), g.k.to_numpy()
		sref = st[st.set == 'mean'].set_index('k').ari_thr
		dthr = st[st.set == name].set_index('k').ari_thr - sref
		for lab, sel in [('k2-20', d.k >= 2), ('k3-10', (d.k >= 3) & (d.k <= 10))]:
			row = {'set': name, 'range': lab, 'ari_ref': g.ari_ref[sel.to_numpy()].median(),
			       'ratio_p25': g.groupby('threshold').p25.first().median()}
			for m in METRICS:
				row['d_' + m] = d.loc[sel, m].median()
				row['pct_' + m] = np.mean(d.loc[sel, m] > 1e-9)
				# nº de umbrales (de 5) con mediana por umbral > 0
				row['thr_' + m] = int((d[sel].groupby('threshold')[m].median() > 0).sum())
			kk = dthr.index[(dthr.index >= (3 if lab == 'k3-10' else 2)) & (dthr.index <= (10 if lab == 'k3-10' else 20))]
			row['d_ari_thr'] = dthr.loc[kk].median()
			out.append(row)
	summ = pd.DataFrame(out)
	summ.to_csv(os.path.join(HERE, 'review_resumen.csv'), index=False)
	pd.set_option('display.width', 250)
	pd.set_option('display.max_columns', 60)
	cols = ['set', 'range', 'ari_ref', 'ratio_p25'] + ['d_' + m for m in METRICS] + ['d_ari_thr']
	print(summ[cols].round(3).to_string(index=False))
	print('\nshare of cases better / thresholds (of 5) with median diff > 0:')
	print(summ[['set', 'range'] + ['pct_' + m for m in METRICS]].round(2).to_string(index=False))
	print(summ[['set', 'range'] + ['thr_' + m for m in METRICS]].to_string(index=False))
	print('\nSD block structure (median over thresholds):')
	print(pd.DataFrame(blk).groupby('var').median(numeric_only=True).drop(columns='threshold').round(3).to_string())
	b = bs.pivot_table(index=['threshold', 'k'], columns='set', values='boot100')
	b30 = tab[tab.set.isin(['mean', 'median', 'median_pipe'])].pivot_table(index=['threshold', 'k'], columns='set',
	                                                                       values='boot')
	print('\nBootstrap median - mean, seed 1000+t / 30 reps vs seed 7000+t / 100 reps:')
	for nm in ['median', 'median_pipe']:
		d30, d100 = b30[nm] - b30['mean'], b[nm] - b['mean']
		print(nm, 'all k: 30 reps %.3f (%.0f%%), 100 reps %.3f (%.0f%%)' % (d30.median(), 100 * (d30 > 0).mean(),
		      d100.median(), 100 * (d100 > 0).mean()),
		      '| per threshold (100 reps):', d100.groupby(level=0).median().round(3).to_dict(),
		      '| k3-10 (100):', round(d100[(d100.index.get_level_values(1) >= 3) &
		                                    (d100.index.get_level_values(1) <= 10)].median(), 3))
	print(f'\nDone in {time.time() - t0:.0f} s')
