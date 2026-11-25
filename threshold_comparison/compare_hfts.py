"""Compara los HFT obtenidos con distintos umbrales de segmentación.

Los HFT de cada umbral salen de moranvar_plots.R con el umbral fijado (EHAB_FORCE_RES), ejecutado en
una copia de segm_optim, con SCENARIO = umbral0<k>. Para cada umbral y corte (p25, p50, p75) calcula:

  - el número de HFTs;
  - la proporción de la varianza de los píxeles que explican los HFTs (R², media de las variables),
    a partir de las medias, varianzas y tamaños de los segmentos (csv/park_<id>_hri_results<k>.csv);
  - el acuerdo entre los mapas de HFTs de dos umbrales (índice de Rand ajustado, por píxeles).

    PYTHONNOUSERSITE=1 python compare_hfts.py RESULTADOS   # results/<STUDY_AREA> de la copia
"""

import glob
import os
import shutil
import sys
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.ticker

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
SEGM = os.path.normpath(os.path.join(HERE, '..', 'segm_optim'))
sys.path.insert(0, SEGM)
from config import ENV_VARS, RESOLUTION  # noqa: E402

ROOT = os.path.abspath(sys.argv[1])
os.chdir(SEGM)
PARK = int(np.atleast_1d(np.genfromtxt('csv/segm_done.csv', dtype=int))[0])
VARS = list(ENV_VARS)
CUTS = ['p25', 'p50', 'p75']
RUNS = sorted(int(os.path.basename(d)[-1]) for d in glob.glob(os.path.join(ROOT, 'umbral0[1-9]')))

INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
C = {4: '#eb6834', 5: '#1baf7a', 7: '#2a78d6'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.titleweight': 'bold', 'axes.titlesize': 10, 'axes.spines.top': False,
                     'axes.spines.right': False})


def clean_id(v):
	s = str(v).strip()
	return s[:-2] if s.endswith('.0') else s


def class_shp(k, cut):
	return os.path.join(ROOT, f'umbral0{k}', cut, f'park_segm_{PARK}_{k}_{cut}_class.shp')


def r2(k, cut):
	"""R² de los HFTs y de los segmentos (media de las variables) y número de HFTs."""
	hri = pd.read_csv(f'csv/park_{PARK}_hri_results{k}.csv', sep=' ', dtype={'segm_id': str})
	hri['segm_id'] = hri['segm_id'].map(clean_id)
	cls = gpd.read_file(class_shp(k, cut))
	cls['segm_id'] = cls['segm_id'].map(clean_id)
	hft = cls.drop_duplicates('segm_id').set_index('segm_id')['hclst_m']
	hri = hri[hri['segm_id'].isin(hft.index)]
	h = hft.loc[hri['segm_id']].to_numpy()
	n = hri['sumpamask'].to_numpy(float)
	out_h, out_s = [], []
	for v in VARS:
		m = pd.to_numeric(hri[v + 'pamean'], errors='coerce').to_numpy()
		var = pd.to_numeric(hri[v + 'pavar'], errors='coerce').to_numpy()
		grand = np.sum(n * m) / n.sum()
		within = np.sum(n * var) / n.sum()
		between_s = np.sum(n * (m - grand) ** 2) / n.sum()
		tot = within + between_s
		between_h = sum(n[h == c].sum() * (np.sum(n[h == c] * m[h == c]) / n[h == c].sum() - grand) ** 2
		                for c in np.unique(h)) / n.sum()
		out_h.append(between_h / tot)
		out_s.append(between_s / tot)
	return np.mean(out_h), np.mean(out_s), len(np.unique(h))


def ari(k1, k2, cut):
	"""Índice de Rand ajustado entre los mapas de HFTs de dos umbrales, contando píxeles."""
	a = gpd.read_file(class_shp(k1, cut)).dissolve(by='hclst_m').reset_index()[['hclst_m', 'geometry']]
	b = gpd.read_file(class_shp(k2, cut)).dissolve(by='hclst_m').reset_index()[['hclst_m', 'geometry']]
	ov = gpd.overlay(a, b, how='intersection', keep_geom_type=True)
	ov['px'] = ov.area / RESOLUTION ** 2
	tab = ov.pivot_table(index='hclst_m_1', columns='hclst_m_2', values='px', aggfunc='sum', fill_value=0).to_numpy()
	comb = lambda x: x * (x - 1) / 2
	n = tab.sum()
	sij, sa, sb = comb(tab).sum(), comb(tab.sum(1)).sum(), comb(tab.sum(0)).sum()
	exp = sa * sb / comb(n)
	return (sij - exp) / ((sa + sb) / 2 - exp)


rows = []
for k in RUNS:
	for cut in CUTS:
		rh, rs, nh = r2(k, cut)
		rows.append({'umbral': k / 10, 'corte': cut, 'HFTs': nh, 'R2_HFTs': rh, 'R2_segmentos': rs})
tab = pd.DataFrame(rows)
pairs = [(a, b) for i, a in enumerate(RUNS) for b in RUNS[i + 1:]]
agree = pd.DataFrame({f'0.{a} vs 0.{b}': [ari(a, b, cut) for cut in CUTS] for a, b in pairs}, index=CUTS)
tab.to_csv(os.path.join(HERE, 'hfts_por_umbral.csv'), index=False, float_format='%.4f')
agree.to_csv(os.path.join(HERE, 'hfts_acuerdo_ari.csv'), float_format='%.3f')
print(tab.round(3).to_string(index=False))
print('\nÍndice de Rand ajustado entre mapas de HFTs (1 = idénticos, 0 = como al azar):\n', agree.round(3).to_string())

# --- Figura 4: varianza explicada por los HFTs según el número de HFTs ---
fig, a = plt.subplots(figsize=(6.4, 4.4))
a.grid(True, color=GRID, lw=0.8)
a.set_axisbelow(True)
for k in RUNS:
	t = tab[tab['umbral'] == k / 10].sort_values('HFTs')
	a.plot(t['HFTs'], t['R2_HFTs'], '-o', color=C.get(k, INK2), lw=2, ms=7,
	       label=f"umbral 0.{k} ({t['R2_segmentos'].iat[0]:.2f} con los segmentos)")
	for _, r in t.iterrows():
		a.annotate(r['corte'], (r['HFTs'], r['R2_HFTs']), xytext=(5, -10), textcoords='offset points', fontsize=7,
		           color=INK2)
a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
a.set_xlabel('Número de HFTs')
a.set_ylabel('Varianza de los píxeles explicada por los HFTs (R²)')
a.set_title('Varianza explicada por los HFTs (media de las 7 variables)', loc='left')
a.legend(frameon=False, fontsize=8, loc='lower right')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'fig4_hfts_r2.png'), dpi=150)

# --- Figura 5: mapas de HFTs (p50) de moranvar_plots.R, uno al lado del otro ---
maps = [(k, os.path.join(ROOT, f'umbral0{k}', 'p50', f'map_{PARK}_{k}_p50_segms_hclust.png')) for k in RUNS]
fig, axs = plt.subplots(1, len(maps), figsize=(6.5 * len(maps), 5.6))
for a, (k, f) in zip(np.atleast_1d(axs), maps):
	a.imshow(mpimg.imread(f))
	t = tab[(tab['umbral'] == k / 10) & (tab['corte'] == 'p50')].iloc[0]
	a.set_title(f"Umbral 0.{k}: {int(t['HFTs'])} HFTs, R² = {t['R2_HFTs']:.2f}", loc='left')
	a.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'fig5_mapas_hfts_p50.png'), dpi=110, bbox_inches='tight')

# Resultados completos de moranvar_plots.R de cada umbral, para verlos con detalle
for k in RUNS:
	dst = os.path.join(HERE, 'hfts', f'umbral0{k}')
	shutil.rmtree(dst, ignore_errors=True)
	shutil.copytree(os.path.join(ROOT, f'umbral0{k}'), dst)
print('\nFiguras y tablas en', HERE)
