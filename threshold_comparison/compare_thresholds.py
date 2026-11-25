"""Compara criterios para elegir el umbral de segmentación con los datos de la última ejecución.

Lee (sin modificarlos) csv/park_<id>_hri_results<k>.csv y shp/park_segm_<id>_<k>.shp de segm_optim,
recalcula para cada umbral el I de Moran de las medias de los segmentos y la varianza dentro de
los segmentos de cada variable (como moranvar.py) y aplica:

  - actual (optim = 1 y 2): primer cruce de las curvas normalizadas de moranvar_plots.R
  - GS de Espindola et al. (2006): mínimo de V_norm + MI_norm (normalización min-max)
  - F de Johnson et al. (2015), como en i.segment.uspo, con alpha 0.5, 1 y 2 (máximo)
  - GS con límites fijos: mínimo de V / V_total + |MI| (no depende de los umbrales probados)

para el rango completo de umbrales, para rangos más cortos y quitando cada variable.

    PYTHONNOUSERSITE=1 python compare_thresholds.py
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import libpysal
import esda
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from osgeo import gdal

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
SEGM = os.path.normpath(os.path.join(HERE, '..', 'segm_optim'))
sys.path.insert(0, SEGM)
from config import ENV_VARS  # noqa: E402

os.chdir(SEGM)
PARK = int(np.atleast_1d(np.genfromtxt('csv/segm_done.csv', dtype=int))[0])
VARS = list(ENV_VARS)
SHORT = [v.split('_')[2] for v in VARS]
NV = len(VARS)
KS = np.arange(1, 10)

# Paleta (dataviz: orden categórico fijo, texto en tinta neutra)
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6']
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.titleweight': 'bold', 'axes.titlesize': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})


def clean_id(v):
	s = str(v).strip()
	return s[:-2] if s.endswith('.0') else s


# --- Estadísticos por umbral y variable ---
MI = np.full((9, NV), np.nan)    # I de Moran de las medias de los segmentos (vecinos tipo torre, como moranvar.py)
MIQ = np.full((9, NV), np.nan)   # con vecinos tipo reina (también los que se tocan solo en una esquina)
MIK = np.full((9, NV), np.nan)   # con los 4 segmentos más cercanos (centroides): sin segmentos aislados
ISL = np.zeros(9, dtype=int)     # segmentos sin ningún vecino tipo torre
SE = np.full((9, NV), np.nan)    # su error típico (hipótesis de normalidad)
VW = np.full((9, NV), np.nan)    # varianza dentro de los segmentos ponderada por área (moranvar.py)
VT = np.full((9, NV), np.nan)    # varianza total de los píxeles = dentro + entre segmentos
NSEG = np.zeros(9, dtype=int)
for k in KS:
	hri = pd.read_csv(f'csv/park_{PARK}_hri_results{k}.csv', sep=' ', dtype={'segm_id': str})
	hri['segm_id'] = hri['segm_id'].map(clean_id)
	gdf = gpd.read_file(f'shp/park_segm_{PARK}_{k}.shp')
	gdf['segm_id'] = gdf['segm_id'].map(clean_id)
	diss = gdf.dissolve(by='segm_id', as_index=False, aggfunc='first')
	hri = hri[hri['segm_id'].isin(set(diss['segm_id']))]
	diss = diss.set_index('segm_id').loc[hri['segm_id']].reset_index()
	w = libpysal.weights.Rook.from_dataframe(diss, use_index=False, silence_warnings=True)
	wq = libpysal.weights.Queen.from_dataframe(diss, use_index=False, silence_warnings=True)
	wk = libpysal.weights.KNN.from_dataframe(diss, k=4)
	n = hri['sumpamask'].to_numpy(float)
	NSEG[k - 1], ISL[k - 1] = len(hri), len(w.islands)
	for j, v in enumerate(VARS):
		m = pd.to_numeric(hri[v + 'pamean'], errors='coerce').to_numpy()
		var = pd.to_numeric(hri[v + 'pavar'], errors='coerce').to_numpy()
		mo = esda.Moran(m, w, permutations=0)
		MI[k - 1, j], SE[k - 1, j] = abs(mo.I), mo.seI_norm
		MIQ[k - 1, j] = abs(esda.Moran(m, wq, permutations=0).I)
		MIK[k - 1, j] = abs(esda.Moran(m, wk, permutations=0).I)
		VW[k - 1, j] = np.sum(var * n) / n.sum()
		grand = np.sum(m * n) / n.sum()
		VT[k - 1, j] = VW[k - 1, j] + np.sum(n * (m - grand) ** 2) / n.sum()


# --- Criterios: ks son los umbrales probados (1-9) y vs las variables usadas ---
def minmax(x, high_good=False):
	r = x.max(0) - x.min(0)
	r[r == 0] = np.inf
	return (x.max(0) - x) / r if high_good else (x - x.min(0)) / r


def crossing(thr, x1, x2):
	"""Cruces de dos curvas interpolados sobre los umbrales, como en moranvar_plots.R."""
	above = (x1 > x2).astype(int)
	ip = np.where(np.diff(above) != 0)[0]
	s1, s2 = x1[ip + 1] - x1[ip], x2[ip + 1] - x2[ip]
	t = (x2[ip] - x1[ip]) / (s1 - s2)
	return thr[ip] + t * (thr[ip + 1] - thr[ip])


def current(ks, vs, optim, mi=MI):
	M, V = mi[np.ix_(ks - 1, vs)], VW[np.ix_(ks - 1, vs)]
	x1, x2 = minmax(M).mean(1), minmax(V, high_good=True).mean(1)
	xp = crossing(ks.astype(float), x1, x2)
	if len(xp) == 0:
		return 5
	res = np.round(xp[0]) - (1 if optim == 2 else 0)
	return int(ks[np.argmin(np.abs(ks - res))])


def gs_scores(ks, vs, mi=MI):
	M, V = mi[np.ix_(ks - 1, vs)], VW[np.ix_(ks - 1, vs)]
	return (minmax(M) + minmax(V)).mean(1)  # menor = mejor


def f_scores(ks, vs, alpha, mi=MI):
	M, V = mi[np.ix_(ks - 1, vs)], VW[np.ix_(ks - 1, vs)]
	nv, nm = minmax(V, high_good=True).mean(1), minmax(M, high_good=True).mean(1)
	den = alpha ** 2 * nm + nv
	return np.where(den > 0, (1 + alpha ** 2) * nv * nm / np.where(den > 0, den, 1), 0)  # mayor = mejor


def gsfix_scores(ks, vs, mi=MI):
	return (VW[np.ix_(ks - 1, vs)] / VT[np.ix_(ks - 1, vs)] + mi[np.ix_(ks - 1, vs)]).mean(1)  # menor = mejor


CRITERIA = [
	('Actual (optim = 2): 1.er cruce - 1', lambda ks, vs, mi=MI: current(ks, vs, 2, mi)),
	('Actual (optim = 1): 1.er cruce', lambda ks, vs, mi=MI: current(ks, vs, 1, mi)),
	('GS Espindola (mín. V + MI)', lambda ks, vs, mi=MI: int(ks[np.argmin(gs_scores(ks, vs, mi))])),
	('F Johnson, alfa = 0.5', lambda ks, vs, mi=MI: int(ks[np.argmax(f_scores(ks, vs, 0.5, mi))])),
	('F Johnson, alfa = 1', lambda ks, vs, mi=MI: int(ks[np.argmax(f_scores(ks, vs, 1, mi))])),
	('F Johnson, alfa = 2', lambda ks, vs, mi=MI: int(ks[np.argmax(f_scores(ks, vs, 2, mi))])),
	('GS límites fijos (mín. V/Vtot + |MI|)', lambda ks, vs, mi=MI: int(ks[np.argmin(gsfix_scores(ks, vs, mi))])),
]
ALLV = np.arange(NV)
RANGES = [('0.1-0.9', KS), ('0.1-0.8', KS[:8]), ('0.2-0.9', KS[1:]), ('0.1-0.7', KS[:7]), ('0.3-0.9', KS[2:])]

sel = pd.DataFrame({name: [f(ks, ALLV) for _, ks in RANGES] for name, f in CRITERIA}, index=[r for r, _ in RANGES]).T
lovo = pd.DataFrame({name: [f(KS, np.delete(ALLV, j)) for j in range(NV)] for name, f in CRITERIA},
                    index=['sin ' + s for s in SHORT]).T
wsel = pd.DataFrame({name: [f(KS, ALLV, m) for m in (MIQ, MIK)] for name, f in CRITERIA},
                    index=['vecinos reina', '4 vecinos más cercanos']).T

out = HERE
tab = pd.DataFrame({'umbral': KS / 10, 'segmentos': NSEG, 'aislados_torre': ISL, 'MI_medio': MI.mean(1),
                    'MI_reina': MIQ.mean(1), 'MI_4vecinos': MIK.mean(1), 'V_sobre_Vtot_medio': (VW / VT).mean(1),
                    'GS_Espindola': gs_scores(KS, ALLV), 'F_alfa0.5': f_scores(KS, ALLV, 0.5), 'F_alfa1': f_scores(KS, ALLV, 1),
                    'F_alfa2': f_scores(KS, ALLV, 2), 'GS_limites_fijos': gsfix_scores(KS, ALLV)})
tab.to_csv(os.path.join(out, 'criterios_por_umbral.csv'), index=False, float_format='%.4f')
pd.concat([sel, lovo, wsel], axis=1).to_csv(os.path.join(out, 'umbral_elegido.csv'))
pd.DataFrame(MI, index=KS / 10, columns=SHORT).to_csv(os.path.join(out, 'moran_por_variable.csv'), float_format='%.4f')
pd.DataFrame(VW / VT, index=KS / 10, columns=SHORT).to_csv(os.path.join(out, 'V_sobre_Vtot_por_variable.csv'), float_format='%.4f')
with pd.option_context('display.width', 200, 'display.max_columns', 20):
	print(tab.round(3).to_string(index=False))
	print('\nUmbral elegido (x 0.1) según el rango de umbrales probados:\n', sel.to_string())
	print('\nQuitando una variable (rango 0.1-0.9):\n', lovo.to_string())
	print('\nCon otros vecinos para el I de Moran (rango 0.1-0.9):\n', wsel.to_string())

# --- Figura 1: curvas y criterios ---
def label_right(a, x, ys, names, gap):
	"""Etiquetas a la derecha de las líneas, separadas al menos gap en y para que no se solapen."""
	order = np.argsort(ys)
	pos = np.asarray(ys, float)[order]
	for i in range(1, len(pos)):
		pos[i] = max(pos[i], pos[i - 1] + gap)
	for p, i in zip(pos, order):
		a.text(x, p, names[i], fontsize=7, color=INK2, va='center')


fig, ax = plt.subplots(2, 2, figsize=(12, 8.8))
thr = KS / 10
for a in ax.flat:
	a.grid(True, color=GRID, lw=0.8)
	a.set_axisbelow(True)
	a.set_xticks(thr)
	a.set_xlabel('Umbral de similitud')

a = ax[0, 0]
x1, x2 = minmax(MI).mean(1), minmax(VW, high_good=True).mean(1)
a.plot(thr, x1, '-o', color=C[0], lw=2, ms=6, label='I de Moran normalizado (bajo = mejor)')
a.plot(thr, x2, '-o', color=C[1], lw=2, ms=6, label='1 - varianza normalizada (alto = mejor)')
xp = crossing(KS.astype(float), x1, x2)
a.plot(xp / 10, np.interp(xp, KS, x1), 'D', color=INK, ms=7, label=f'1.er cruce ({xp[0] / 10:.2f})')
for k, lab, ha in [(current(KS, ALLV, 1), ' optim = 1', 'left'), (current(KS, ALLV, 2), 'optim = 2 ', 'right')]:
	a.axvline(k / 10, color=INK2, ls='--', lw=1)
	a.text(k / 10, 1.03, lab, ha=ha, va='bottom', color=INK2, fontsize=8)
a.set_ylim(-0.05, 1.1)
a.set_title('a) Método actual: cruce de las curvas normalizadas', loc='left')
a.legend(frameon=False, loc='lower left', fontsize=8)

a = ax[0, 1]
for j in range(NV):
	a.plot(thr, MI[:, j], color='#b9b8b2', lw=1)
label_right(a, 0.915, MI[-1], SHORT, 0.045)
mse = np.sqrt((SE ** 2).mean(1))
a.fill_between(thr, MI.mean(1) - 1.96 * mse, MI.mean(1) + 1.96 * mse, color=C[0], alpha=0.15, lw=0)
a.plot(thr, MI.mean(1), '-o', color=C[0], lw=2, ms=6, label='Media, vecinos torre (actual)')
a.plot(thr, MIQ.mean(1), '-s', color=C[2], lw=2, ms=5, label='Media, vecinos reina')
a.plot(thr, MIK.mean(1), '-^', color=C[6], lw=2, ms=5, label='Media, 4 vecinos más cercanos')
a.legend(frameon=False, loc='upper center', fontsize=8)
a.set_xlim(0.05, 1.0)
a.set_ylabel('|I de Moran| entre segmentos vecinos')
a.set_title('b) I de Moran: cada variable (gris) y media (banda ± 1.96 e.t.)', loc='left')
a.set_xticks(thr, [f'{t:.1f}\n{n}\n{i}' for t, n, i in zip(thr, NSEG, ISL)])
a.annotate('umbral\nsegmentos\nsin vecino torre', xy=(0, 0), xycoords='axes fraction', xytext=(-8, -7),
           textcoords='offset points', ha='right', va='top', fontsize=plt.rcParams['font.size'], color=INK2)
a.set_xlabel('')

a = ax[1, 0]
for j in range(NV):
	a.plot(thr, VW[:, j] / VT[:, j], color='#b9b8b2', lw=1)
label_right(a, 0.915, VW[-1] / VT[-1], SHORT, 0.035)
a.plot(thr, (VW / VT).mean(1), '-o', color=C[1], lw=2, ms=6)
a.set_xlim(0.05, 1.0)
a.set_ylabel('Varianza dentro de los segmentos / total')
a.set_title('c) Varianza no explicada por los segmentos (1 - R²)', loc='left')

a = ax[1, 1]
series = [('GS Espindola', 1 - minmax(gs_scores(KS, ALLV)[:, None])[:, 0], C[1]),
          ('F Johnson (alfa = 1)', minmax(f_scores(KS, ALLV, 1)[:, None], high_good=False)[:, 0], C[2]),
          ('GS límites fijos', 1 - minmax(gsfix_scores(KS, ALLV)[:, None])[:, 0], C[3])]
for name, s, col in series:
	a.plot(thr, s, '-o', color=col, lw=2, ms=6, label=name)
	b = np.argmax(s)
	a.plot(thr[b], s[b], 'o', ms=12, mfc='none', mec=INK, mew=1.5)
a.axvline(current(KS, ALLV, 2) / 10, color=C[0], ls='--', lw=1.5, label='Actual (optim = 2)')
a.set_ylim(-0.05, 1.12)
a.set_ylabel('Puntuación reescalada (1 = mejor)')
a.set_title('d) Criterios alternativos (círculo = umbral elegido)', loc='left')
a.legend(frameon=False, loc='lower left', fontsize=8)
fig.suptitle(f'Elección del umbral de segmentación, área {PARK} (Guajares, PA_BUFFER = 500)', x=0.01, ha='left',
             fontweight='bold', color=INK)
fig.tight_layout()
fig.savefig(os.path.join(out, 'fig1_criterios.png'), dpi=150)

# --- Figura 2: umbral elegido según el rango probado y quitando cada variable ---
both = pd.concat([sel, lovo, wsel], axis=1)
cmap = LinearSegmentedColormap.from_list('b', ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
fig, a = plt.subplots(figsize=(12, 3.6))
a.imshow(both.to_numpy(), cmap=cmap, vmin=1, vmax=9, aspect='auto')
for i in range(both.shape[0]):
	for j in range(both.shape[1]):
		v = both.iat[i, j]
		a.text(j, i, f'0.{v}', ha='center', va='center', fontsize=8, color='white' if v >= 5 else INK)
a.set_xticks(range(both.shape[1]), both.columns, rotation=35, ha='right')
a.set_yticks(range(both.shape[0]), both.index)
a.axvline(len(RANGES) - 0.5, color='white', lw=4)
a.axvline(len(RANGES) + NV - 0.5, color='white', lw=4)
a.tick_params(length=0)
for s in a.spines.values():
	s.set_visible(False)
a.set_title('Umbral elegido por cada criterio\nrangos de umbrales probados  |  quitando una variable  |  '
            'vecinos del I de Moran', loc='left')
fig.tight_layout()
fig.savefig(os.path.join(out, 'fig2_estabilidad.png'), dpi=150)

# --- Figura 3: segmentaciones de los umbrales elegidos sobre NDVImax ---
chosen = [4, 5, 7]  # GS / F, F alfa 0.5 / límites fijos y el método actual con el rango completo
ds = gdal.Open(os.path.join('..', 'inVars', ENV_VARS[VARS[0]]['file']))
gt, band = ds.GetGeoTransform(), ds.GetRasterBand(1)
x0, y0, x1b, y1b = gpd.read_file(f'shp/park_segm_{PARK}_1.shp').total_bounds
c0, r0 = int((x0 - gt[0]) / gt[1]), int((gt[3] - y1b) / -gt[5])
nc, nr = int(round((x1b - x0) / gt[1])), int(round((y1b - y0) / -gt[5]))
img = band.ReadAsArray(c0, r0, nc, nr).astype(float)
img[img == ENV_VARS[VARS[0]]['nodata']] = np.nan
ext = (gt[0] + c0 * gt[1], gt[0] + (c0 + nc) * gt[1], gt[3] + (r0 + nr) * gt[5], gt[3] + r0 * gt[5])
pa = gpd.read_file(os.path.join('pa', [f for f in os.listdir('pa') if f.endswith('.gpkg')][0]))
who = {k: [n for n in sel.index if sel.loc[n, '0.1-0.9'] == k] for k in chosen}
fig, axs = plt.subplots(1, len(chosen), figsize=(6 * len(chosen), 4.6))
for a, k in zip(np.atleast_1d(axs), chosen):
	a.imshow(img, extent=ext, cmap='Greys', vmin=np.nanpercentile(img, 2), vmax=np.nanpercentile(img, 98))
	g = gpd.read_file(f'shp/park_segm_{PARK}_{k}.shp').dissolve(by='segm_id')
	g.boundary.plot(ax=a, color=C[1], lw=0.6)
	pa.to_crs(g.crs).boundary.plot(ax=a, color=C[0], lw=1.2)
	a.set_title(f'Umbral 0.{k}: {NSEG[k - 1]} segmentos', loc='left')
	a.text(0.01, -0.02, '\n'.join(who[k]), transform=a.transAxes, va='top', fontsize=7.5, color=INK2)
	a.set_xticks([])
	a.set_yticks([])
fig.suptitle('Segmentos (naranja) de los umbrales elegidos sobre NDVImax; en azul, el área quemada', x=0.01, ha='left',
             fontweight='bold')
fig.tight_layout()
fig.savefig(os.path.join(out, 'fig3_segmentos.png'), dpi=150)
print('\nFiguras y tablas en', out)
