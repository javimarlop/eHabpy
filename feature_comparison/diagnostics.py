"""Diagnóstico de los datos: por qué media o mediana, y la desviación típica (SD) o el CV de cada variable, cambiarían
(o no) la clasificación de los segmentos en HFTs.

    /opt/anaconda3/envs/ehabpy/bin/python -s diagnostics.py [UMBRALES]

UMBRALES: índices separados por comas (por defecto 3,5,7). Lee ../algorithm_comparison/data (seg_0k.bin, grid.csv,
vars.csv, de segment_sums.R) y los 7 rásters ../inVars/Guajares_<i>_*.tif (NoData -9999). Píxeles usados: seg > 0 y
datos en las 7 variables. Escribe en diag/ (por umbral k):
  diag_seg_0k.csv     por segmento y variable: n, borde (fracción de píxeles con un vecino-4 de otro segmento o sin
                      datos), media, mediana, sd (unidades originales), asimetría, curtosis (exceso), dif = (media -
                      mediana) / sd, z de media y mediana (tipificadas con la media y SD de los píxeles, como std_stats)
  diag_var_0k.csv     resumen por variable (las medidas del informe)
y diag/log.txt con el informe.
"""

import glob
import os
import sys

import numpy as np
import pandas as pd
from osgeo import gdal
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'algorithm_comparison', 'data')
INV = os.path.join(HERE, '..', 'inVars')
OUT = os.path.join(HERE, 'diag')
os.makedirs(OUT, exist_ok=True)
THR = [int(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else [3, 5, 7]
LOG = open(os.path.join(OUT, 'log.txt'), 'w')


def say(*a):
	s = ' '.join(str(x) for x in a)
	print(s)
	LOG.write(s + '\n')


def wcorr(x, y, w=None):
	w = np.ones_like(x) if w is None else w
	mx, my = np.average(x, weights=w), np.average(y, weights=w)
	return np.sum(w * (x - mx) * (y - my)) / np.sqrt(np.sum(w * (x - mx) ** 2) * np.sum(w * (y - my) ** 2))


def spear(x, y):
	return stats.spearmanr(x, y).statistic


def r2_ols(y, X):
	X = np.column_stack([np.ones(len(y)), X])
	b, *_ = np.linalg.lstsq(X, y, rcond=None)
	return 1 - np.sum((y - X @ b) ** 2) / np.sum((y - y.mean()) ** 2)


# --- Rásters ---
grid = pd.read_csv(os.path.join(DATA, 'grid.csv'))
nr, nc = int(grid['nrow'][0]), int(grid['ncol'][0])
vars_ = pd.read_csv(os.path.join(DATA, 'vars.csv'))
short = list(vars_['short'])
files = sorted(glob.glob(os.path.join(INV, 'Guajares_[0-9]_*.tif')))
assert len(files) == 7, files
P = []
for i, (f, v) in enumerate(zip(files, vars_['var'])):
	# vars.csv cambia '-' por '_' (nombres de ENV_VARS); comprobar el orden por el nombre
	assert os.path.basename(f).replace('-', '_').startswith(v), (f, v)
	ds = gdal.Open(f)
	assert (ds.RasterYSize, ds.RasterXSize) == (nr, nc), (f, ds.RasterYSize, ds.RasterXSize)
	b = ds.GetRasterBand(1)
	nd = b.GetNoDataValue()
	x = b.ReadAsArray().astype(float).ravel()
	x[(x == nd) | ~np.isfinite(x)] = np.nan
	P.append(x)
	say(f'{short[i]:15s} {os.path.basename(f)[:60]}  NoData={nd}  min={np.nanmin(x):.4g} max={np.nanmax(x):.4g} '
	    f'neg={np.mean(x[np.isfinite(x)] < 0) * 100:.2f}%')
P = np.column_stack(P)
valid = np.all(np.isfinite(P), axis=1)
NV = P.shape[1]


def edge_mask(seg2):
	"""Píxeles con algún vecino-4 de otro segmento (o sin datos / fuera de la malla)."""
	e = np.zeros_like(seg2, dtype=bool)
	p = np.pad(seg2, 1, constant_values=0)
	c = p[1:-1, 1:-1]
	for dy, dx in ((0, 1), (2, 1), (1, 0), (1, 2)):
		e |= p[dy:dy + nr, dx:dx + nc] != c
	return e.ravel()


def per_segment(t):
	seg = np.fromfile(os.path.join(DATA, f'seg_0{t}.bin'), dtype=np.int32)
	assert seg.size == nr * nc
	bad = (seg > 0) & ~valid
	say(f'\n=== Threshold 0.{t} ===')
	say(f'pixels seg>0 without data in some variable: {bad.sum()} (excluded)')
	ok = (seg > 0) & valid
	edge = edge_mask(np.where(ok, seg, 0).reshape(nr, nc))[ok]
	g = seg[ok]
	X = P[ok]
	# tipificación global de los píxeles (como segment_sums.R / std_stats: media y SD de los píxeles)
	mu, sdp = X.mean(0), X.std(0, ddof=1)
	segs, g = np.unique(g, return_inverse=True)
	ns = len(segs)
	n = np.bincount(g).astype(float)
	ef = np.bincount(g, weights=edge) / n
	rows = []
	order0 = np.argsort(g, kind='stable')
	starts = np.r_[0, np.cumsum(n).astype(int)[:-1]]
	cnt = n.astype(int)
	res = {}
	for j in range(NV):
		x = X[:, j]
		o = np.lexsort((x, g))
		xs = x[o]
		med = (xs[starts + (cnt - 1) // 2] + xs[starts + cnt // 2]) / 2
		m = np.bincount(g, weights=x) / n
		d = x - m[g]
		m2 = np.bincount(g, weights=d ** 2) / n
		m3 = np.bincount(g, weights=d ** 3) / n
		m4 = np.bincount(g, weights=d ** 4) / n
		sd = np.sqrt(m2 * n / np.maximum(n - 1, 1))
		with np.errstate(invalid='ignore', divide='ignore'):
			sk = m3 / m2 ** 1.5
			ku = m4 / m2 ** 2 - 3
		npos = np.bincount(g, weights=x > 0)
		nneg = np.bincount(g, weights=x < 0)
		# media de los píxeles interiores (sin borde)
		ni = np.bincount(g, weights=~edge)
		with np.errstate(invalid='ignore', divide='ignore'):
			mint = np.bincount(g, weights=x * ~edge) / ni
		res[short[j]] = dict(mean=m, median=med, sd=sd, skew=sk, kurt=ku, npos=npos, nneg=nneg, mint=mint,
		                     zmean=(m - mu[j]) / sdp[j], zmed=(med - mu[j]) / sdp[j], sdp=sdp[j], mup=mu[j],
		                     wss=np.sum(d ** 2), tss=np.sum((x - mu[j]) ** 2),
		                     pix_skew=stats.skew(x), pix_kurt=stats.kurtosis(x),
		                     # asimetría de los residuos dentro del segmento (todos los píxeles juntos)
		                     res_skew=stats.skew(d), res_kurt=stats.kurtosis(d))
	say(f'{ns} segments, {int(n.sum())} pixels; segment size: min {n.min():.0f}, median {np.median(n):.0f}, '
	    f'max {n.max():.0f}; edge fraction median {np.median(ef):.3f}')
	# tabla por segmento
	tab = pd.DataFrame({'seg_row': segs, 'n': n, 'edge_frac': ef})
	for v in short:
		r = res[v]
		with np.errstate(invalid='ignore', divide='ignore'):
			dif = (r['mean'] - r['median']) / r['sd']
		for k in ('mean', 'median', 'sd', 'skew', 'kurt', 'zmean', 'zmed'):
			tab[f'{v}_{k}'] = r[k]
		tab[f'{v}_dif'] = dif
	tab.to_csv(os.path.join(OUT, f'diag_seg_0{t}.csv'), sep=' ', index=False)
	return n, ef, res, tab


def report(t):
	n, ef, res, tab = per_segment(t)
	w = n
	big = n >= 30  # momentos estables
	logn = np.log10(n)
	out = []
	say('\n-- 1. Mean vs median --')
	say(f'{"var":15s} {"pixSkew":>7s} {"resSkew":>7s} {"resKurt":>7s} {"segSkew50":>9s} {"segSkewIQR":>13s} '
	    f'{"|skew|>1%":>9s} {"dif>0.1%seg":>11s} {"dif>0.1%pix":>11s} {"r(z)":>6s} {"rw(z)":>6s} '
	    f'{"medAbsDz":>8s} {"sdZmean":>7s} {"rhoDifLogn":>10s} {"rhoDifEdge":>10s} {"rhoDifSkew":>10s}')
	Zm = np.column_stack([res[v]['zmean'] for v in short])
	Zd = np.column_stack([res[v]['zmed'] for v in short])
	for v in short:
		r = res[v]
		with np.errstate(invalid='ignore', divide='ignore'):
			dif = (r['mean'] - r['median']) / r['sd']
		okd = np.isfinite(dif)
		adif = np.abs(dif)
		sk = r['skew'][big & np.isfinite(r['skew'])]
		dz = np.abs(r['zmean'] - r['zmed'])
		row = dict(var=v, pix_skew=r['pix_skew'], pix_kurt=r['pix_kurt'], res_skew=r['res_skew'], res_kurt=r['res_kurt'],
		           seg_skew_med=np.median(sk), seg_skew_q25=np.quantile(sk, .25), seg_skew_q75=np.quantile(sk, .75),
		           seg_kurt_med=np.median(r['kurt'][big & np.isfinite(r['kurt'])]),
		           seg_abs_skew_gt1=np.mean(np.abs(sk) > 1),
		           dif_gt01_seg=np.mean(adif[okd] > 0.1), dif_gt01_pix=np.sum(w[okd] * (adif[okd] > 0.1)) / w[okd].sum(),
		           dif_med=np.median(adif[okd]), dif_q90=np.quantile(adif[okd], .9),
		           r_mean_median=np.corrcoef(r['zmean'], r['zmed'])[0, 1],
		           rw_mean_median=wcorr(r['zmean'], r['zmed'], w),
		           med_abs_dz=np.median(dz), q90_abs_dz=np.quantile(dz, .9), sd_zmean=np.sqrt(np.cov(r['zmean'], aweights=w)),
		           rho_dif_logn=spear(adif[okd], logn[okd]), rho_dif_edge=spear(adif[okd], ef[okd]),
		           rho_dif_abs_skew=spear(adif[okd & np.isfinite(r['skew'])], np.abs(r['skew'][okd & np.isfinite(r['skew'])])),
		           # media de los píxeles interiores frente a la media de todos y a la mediana (en z global)
		           med_abs_dz_interior=np.nanmedian(np.abs(r['mint'] - r['mean']) / r['sdp']))
		out.append(row)
		say(f'{v:15s} {row["pix_skew"]:7.2f} {row["res_skew"]:7.2f} {row["res_kurt"]:7.2f} {row["seg_skew_med"]:9.2f} '
		    f'[{row["seg_skew_q25"]:5.2f},{row["seg_skew_q75"]:5.2f}] {row["seg_abs_skew_gt1"] * 100:9.1f} '
		    f'{row["dif_gt01_seg"] * 100:11.1f} {row["dif_gt01_pix"] * 100:11.1f} {row["r_mean_median"]:6.3f} '
		    f'{row["rw_mean_median"]:6.3f} {row["med_abs_dz"]:8.3f} {row["sd_zmean"]:7.3f} {row["rho_dif_logn"]:10.2f} '
		    f'{row["rho_dif_edge"]:10.2f} {row["rho_dif_abs_skew"]:10.2f}')
	# |dif| > 0.1 por clase de tamaño
	cls = pd.qcut(n, 4, labels=['Q1 small', 'Q2', 'Q3', 'Q4 large'])
	say('share of segments with |mean-median| > 0.1 SD, by size quartile (n range):')
	for c in cls.categories:
		s = np.asarray(cls == c)
		shares = [np.nanmean(np.abs(tab.loc[s, f'{v}_dif']) > 0.1) * 100 for v in short]
		say(f'  {c:9s} n {n[s].min():6.0f}-{n[s].max():7.0f}: ' + ' '.join(f'{v[:7]}={x:4.0f}%' for v, x in zip(short, shares)))
	# efecto en las distancias entre segmentos (las de ward_area)
	from scipy.spatial.distance import pdist
	dm, dd = pdist(Zm), pdist(Zd)
	say(f'distances between segments (7 z-variables): corr(mean-based, median-based) = {np.corrcoef(dm, dd)[0, 1]:.4f}; '
	    f'median |d_mean - d_median| / median d_mean = {np.median(np.abs(dm - dd)) / np.median(dm):.3f}')
	say(f'z(mean) - z(median) as share of between-segment spread of z(mean): ' +
	    ' '.join(f'{r["var"][:7]}={r["med_abs_dz"] / r["sd_zmean"]:.3f}' for r in out))

	say('\n-- 2. Within-segment SD --')
	say(f'{"var":15s} {"within%":>7s} {"rhoSD,mean":>10s} {"rhoSD,logn":>10s} {"rhoSD,edge":>10s} {"CV(SD)":>7s} '
	    f'{"CVw(SD)":>7s} {"R2 logSD~logn":>13s} {"R2 ~7means":>10s} {"R2 ~7means+logn":>15s} {"R2 ~means+means2+logn":>21s} '
	    f'{"sdSD/sdMean":>11s}')
	big2 = n >= 10
	LS = np.column_stack([np.log(res[v]['sd'][big2]) for v in short])
	M7 = np.column_stack([res[v]['zmean'][big2] for v in short])
	for k, v in enumerate(short):
		r = res[v]
		sd = r['sd']
		within = r['wss'] / r['tss']
		cvsd = np.std(sd[big2]) / np.mean(sd[big2])
		mw = np.average(sd, weights=w)
		cvw = np.sqrt(np.average((sd - mw) ** 2, weights=w)) / mw
		y = LS[:, k]
		r2n = r2_ols(y, logn[big2])
		r2m = r2_ols(y, M7)
		r2mn = r2_ols(y, np.column_stack([M7, logn[big2]]))
		r2q = r2_ols(y, np.column_stack([M7, M7 ** 2, logn[big2]]))
		# variación entre segmentos de la SD tipificada frente a la de las medias tipificadas
		ratio = np.std(sd[big2] / r['sdp']) / np.std(r['zmean'][big2])
		row = out[k]
		row.update(within_share=within, rho_sd_mean=spear(sd, r['mean']), rho_sd_logn=spear(sd, logn),
		           rho_sd_edge=spear(sd, ef), cv_sd=cvsd, cvw_sd=cvw, r2_logsd_logn=r2n, r2_logsd_means=r2m,
		           r2_logsd_means_logn=r2mn, r2_logsd_quad=r2q, sdSD_over_sdMean=ratio)
		say(f'{v:15s} {within * 100:7.1f} {row["rho_sd_mean"]:10.2f} {row["rho_sd_logn"]:10.2f} {row["rho_sd_edge"]:10.2f} '
		    f'{cvsd:7.2f} {cvw:7.2f} {r2n:13.2f} {r2m:10.2f} {r2mn:15.2f} {r2q:21.2f} {ratio:11.2f}')
	C = np.corrcoef(LS.T)
	ev = np.linalg.eigvalsh(C)[::-1]
	say(f'log SDs of the 7 variables (n>=10, {big2.sum()} segs): mean |r| between variables = '
	    f'{np.mean(np.abs(C[np.triu_indices(NV, 1)])):.2f}; PC1 explains {ev[0] / NV * 100:.0f}%, PC1-2 {ev[:2].sum() / NV * 100:.0f}%')
	# SD de cada segmento tipificada con la SD global (pooled within-sd a la potencia de la media)
	say('pixel-weighted mean within-segment SD / pixel SD: ' +
	    ' '.join(f'{v[:7]}={np.average(res[v]["sd"], weights=w) / res[v]["sdp"]:.2f}' for v in short))

	say('\n-- 3. CV = SD / |mean| --')
	say(f'{"var":15s} {"pixNeg%":>7s} {"segSign±%":>9s} {"pixInSign±%":>11s} {"|m|<.5sd seg%":>13s} {"pix%":>6s} '
	    f'{"CV q05":>7s} {"q25":>7s} {"q50":>7s} {"q75":>7s} {"q95":>7s} {"max":>9s} {"CVmax/CVmed":>11s}')
	for k, v in enumerate(short):
		r = res[v]
		m, sd = r['mean'], r['sd']
		near0 = np.abs(m) < 0.5 * sd
		sign = (r['npos'] > 0) & (r['nneg'] > 0)
		with np.errstate(divide='ignore', invalid='ignore'):
			cv = sd / np.abs(m)
		cvf = cv[np.isfinite(cv)]
		q = np.quantile(cvf, [.05, .25, .5, .75, .95])
		pixneg = np.sum(r['nneg']) / n.sum()
		row = out[k]
		row.update(pix_neg=pixneg, seg_sign_change=sign.mean(), pix_in_sign_change_segs=np.sum(w * sign) / w.sum(),
		           near0_seg=near0.mean(), near0_pix=np.sum(w * near0) / w.sum(), cv_q05=q[0], cv_q25=q[1], cv_q50=q[2],
		           cv_q75=q[3], cv_q95=q[4], cv_max=cvf.max(), cv_inf=np.sum(~np.isfinite(cv)),
		           rho_cv_sd=spear(cvf, sd[np.isfinite(cv)]), rho_cv_mean=spear(cvf, np.abs(m)[np.isfinite(cv)]))
		say(f'{v:15s} {pixneg * 100:7.2f} {sign.mean() * 100:9.1f} {row["pix_in_sign_change_segs"] * 100:11.1f} '
		    f'{near0.mean() * 100:13.1f} {row["near0_pix"] * 100:6.1f} {q[0]:7.3f} {q[1]:7.3f} {q[2]:7.3f} {q[3]:7.3f} '
		    f'{q[4]:7.3f} {cvf.max():9.3g} {cvf.max() / q[2]:11.1f}')
	say('Spearman of CV with SD and with |mean| (which one drives the CV):')
	say('  ' + ' '.join(f'{r["var"][:7]}: sd {r["rho_cv_sd"]:.2f} |m| {r["rho_cv_mean"]:.2f};' for r in out))
	# temperatura: °C frente a K (escala de intervalo)
	tv = [v for v in short if v.startswith('Temp')][0]
	r = res[tv]
	cvc, cvk = r['sd'] / np.abs(r['mean']), r['sd'] / (r['mean'] + 273.15)
	say(f'{tv}: CV in degC vs K: Spearman r = {spear(cvc, cvk):.3f}, Pearson r(sd, CV_K) = {np.corrcoef(r["sd"], cvk)[0, 1]:.4f}; '
	    f'CV_degC range {cvc.min():.4f}-{cvc.max():.4f}, CV_K range {cvk.min():.5f}-{cvk.max():.5f} '
	    f'(mean range {r["mean"].min():.2f}-{r["mean"].max():.2f} degC)')
	pd.DataFrame(out).to_csv(os.path.join(OUT, f'diag_var_0{t}.csv'), sep=' ', index=False)


for t in THR:
	report(t)
LOG.close()
