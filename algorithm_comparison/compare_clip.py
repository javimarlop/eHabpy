"""Compara los resultados del modelo sin y con CLIP_PERCENTILE (ejecuciones completas en copias de segm_optim).

    Rscript segment_sums.R <segm_optim> 1,2,3,4,5,6,7,8,9 data_clip/<run> <results/AREA/SCENARIO>   # por ejecución
    python compare_clip.py

Ejecuciones (carpetas de data_clip, con sus sumas de píxeles, umbrales y HFTs de cada corte):
  A  sin recorte (CLIP_PERCENTILE = 0): debe dar los mismos resultados que antes de añadir la opción
  B  CLIP_PERCENTILE = 0.5 en la segmentación y la clasificación
  C  la segmentación de A con CLIP_PERCENTILE = 0.5 solo en la clasificación (moranvar_plots.R), para separar
     el efecto del recorte en cada paso
Todas las medidas usan los valores originales de los píxeles (sin recortar), como las gráficas del modelo:
r2, jm_mean, jm_min, sil, n_eff y eta2 de cada variable (compare_algorithms.py y compare_pca.py), r2_segm (parte
de la varianza de los píxeles entre los segmentos del umbral elegido: lo más que pueden explicar sus HFTs) y
ari_A (acuerdo con los HFTs de A del mismo corte, por píxeles, aunque los segmentos sean distintos).
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import compare_algorithms as ca
import compare_pca as cp

HERE = ca.HERE
DATA = os.path.join(HERE, 'data_clip')
RUNS = {'A': 'Sin recorte', 'B': 'Recorte 0.5-99.5 % (segmentación y clasificación)',
        'C': 'Recorte solo en la clasificación'}
COL = {'A': ca.INK, 'B': ca.SLOTS[0], 'C': ca.SLOTS[1]}
MARK = {'A': 'o', 'B': 's', 'C': '^'}
CUTS = ['p10', 'p25', 'p50', 'p75', 'p90']
VARS = cp.VARS


def load(run):
	"""Segmentos del umbral elegido, tabla de umbrales y HFTs de cada corte de una ejecución."""
	d = os.path.join(DATA, run)
	thr = pd.read_csv(os.path.join(d, 'thresholds.csv'), sep=' ')
	t = int(thr.loc[thr.selected == 1, 'threshold'].iloc[0])
	ca.DATA = d  # Segments lee de ca.DATA
	s = ca.Segments(t)
	ids = pd.read_csv(os.path.join(d, f'sums_0{t}.csv'), dtype={'segm_id': str}).segm_id.to_numpy()[s.rows]
	hfts = {}
	for c in CUTS:
		f = os.path.join(d, f'hfts_{c}.csv')
		if os.path.exists(f):
			h = pd.read_csv(f, dtype={'segm_id': str}).set_index('segm_id').hft
			hfts[c] = h.reindex(ids).to_numpy()
	nr, nc = pd.read_csv(os.path.join(d, 'grid.csv')).iloc[0]
	seg = np.fromfile(os.path.join(d, f'seg_0{t}.bin'), dtype=np.int32).reshape(nr, nc)
	return {'t': t, 's': s, 'thr': thr, 'hfts': hfts, 'seg': seg}


def pixel_labels(r, cut):
	"""HFT de cada píxel de la malla (-1 fuera de los segmentos)."""
	rowmap = np.full(r['seg'].max() + 1, -1)
	rowmap[r['s'].rows + 1] = np.arange(len(r['s'].rows))
	idx = rowmap[r['seg']]
	return np.where(idx >= 0, r['hfts'][cut][np.maximum(idx, 0)], -1)


def fig_thresholds(runs):
	fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
	for x in ['A', 'B']:
		th = runs[x]['thr']
		axs[0].plot(th.threshold / 10, th.segments, color=COL[x], marker=MARK[x], lw=1.8, label=RUNS[x])
		axs[1].plot(th.threshold / 10, th.score, color=COL[x], marker=MARK[x], lw=1.8, label=RUNS[x])
		sel = th[th.selected == 1]
		for ax, c in zip(axs, ['segments', 'score']):
			ax.scatter(sel.threshold / 10, sel[c], s=150, facecolors='none', edgecolors=COL[x], linewidths=1.8, zorder=3)
	axs[0].set_title('Número de segmentos', loc='left')
	axs[1].set_title('Puntuación del umbral (THRESHOLD_METHOD, mayor = mejor)', loc='left')
	for ax in axs:
		ax.set_xlabel('Umbral de segmentación')
	axs[0].legend(fontsize=8.5)
	fig.suptitle('Segmentación sin y con recorte de las capas (círculo: umbral elegido)', x=0.01, ha='left',
	             fontweight='bold', fontsize=11)
	fig.tight_layout(rect=(0, 0, 1, 0.93))
	fig.savefig(os.path.join(HERE, 'figC1_umbrales.png'), dpi=150)
	plt.close(fig)


def fig_r2(runs):
	"""R² de los píxeles (valores originales) de los HFTs del umbral elegido de cada ejecución según su número."""
	fig, ax = plt.subplots(figsize=(8.5, 5))
	for x, r in runs.items():
		c = pd.read_csv(os.path.join(DATA, x, 'r2_hfts.csv'), sep=' ')
		c = c[np.isclose(c.threshold, r['t'] / 10) & (c.hfts <= 30)]
		ax.plot(c.hfts, c.r2, color=COL[x], marker=MARK[x], ms=3.5, lw=1.8, label=f'{RUNS[x]} (umbral 0.{r["t"]})')
	ax.set_xlabel('Número de HFTs')
	ax.set_title('R² de los píxeles (media de las variables, valores originales)', loc='left')
	ax.legend(fontsize=8.5, loc='lower right')
	fig.suptitle('HFTs del umbral elegido de cada ejecución', x=0.01, ha='left', fontweight='bold', fontsize=11)
	fig.tight_layout(rect=(0, 0, 1, 0.94))
	fig.savefig(os.path.join(HERE, 'figC4_r2.png'), dpi=150)
	plt.close(fig)


def fig_eta2(tab, cut):
	g = tab[tab.cut == cut].set_index('run')
	rows = ['r2_segm'] + [f'eta2_{v}' for v in VARS]
	labels = ['Segmentos (techo)'] + [f'{v} (HFTs)' for v in VARS]
	fig, ax = plt.subplots(figsize=(9.5, 5.2))
	y = np.arange(len(rows), dtype=float)
	y[1:] += 0.6
	for i, r in enumerate(rows):
		vals = [g.loc[x, r] for x in RUNS if x in g.index]
		ax.plot([min(vals), max(vals)], [y[i]] * 2, color=ca.GRID, lw=3, solid_capstyle='round', zorder=1)
	for x in RUNS:
		if x in g.index:
			ax.scatter([g.loc[x, r] for r in rows], y, s=55, color=COL[x], marker=MARK[x], edgecolors='#fcfcfb',
			           linewidths=1.2, zorder=3, label=f'{RUNS[x]} ({int(g.loc[x, "hfts"])} HFTs)')
	ax.set_yticks(y, labels)
	ax.invert_yaxis()
	ax.set_xlim(0, 1)
	ax.grid(axis='y', visible=False)
	ax.set_xlabel('Parte de la varianza de los píxeles (valores originales) entre segmentos o entre HFTs')
	ax.legend(loc='upper right', fontsize=8.5)
	fig.suptitle(f'Qué explican los segmentos del umbral elegido y los HFTs del corte {cut}', x=0.01, ha='left',
	             fontweight='bold', fontsize=11)
	fig.text(0.01, 0.915, 'Arriba, la media de las variables con los segmentos (lo más que pueden explicar sus HFTs); '
	         'abajo, cada variable con los HFTs', color=ca.INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.91))
	fig.savefig(os.path.join(HERE, f'figC2_eta2_{cut}.png'), dpi=150)
	plt.close(fig)


def fig_maps(runs, tab, cut):
	ref = pixel_labels(runs['A'], cut)
	cols = ca.SLOTS + ['#86b6ef', '#f4a582', '#7fd8b4', '#f6cd73', '#9085e9', '#c3c2b7']
	ok = np.argwhere(ref >= 0)
	(r0, c0), (r1, c1) = ok.min(0), ok.max(0) + 1
	fig, axs = plt.subplots(1, 3, figsize=(17, 4.6))
	for ax, x in zip(axs, RUNS):
		lab = pixel_labels(runs[x], cut)
		# renumeración: cada HFT con el número del HFT de A con el que comparte más píxeles
		both = (lab >= 0) & (ref >= 0)
		u, inv = np.unique(lab[both], return_inverse=True)
		v, invr = np.unique(ref[both], return_inverse=True)
		T = np.zeros((len(u), len(v)))
		np.add.at(T, (inv, invr), 1)
		new, used = {}, set()
		for i, j in sorted(np.ndindex(T.shape), key=lambda ij: -T[ij]):
			if u[i] not in new and j not in used:
				new[u[i]] = j
				used.add(j)
		for k in u:
			if k not in new:
				new[k] = len(v) + len([q for q in new.values() if q >= len(v)])
		img = np.full(lab.shape, np.nan)
		img[lab >= 0] = np.vectorize(new.get)(lab[lab >= 0])
		n = int(np.nanmax(img)) + 1
		cmap = ListedColormap(cols[:max(n, 2)] if n <= len(cols) else plt.get_cmap('tab20').colors[:n])
		ax.imshow(img[r0:r1, c0:c1], cmap=cmap, vmin=-0.5, vmax=cmap.N - 0.5, interpolation='nearest')
		r = tab[(tab.run == x) & (tab.cut == cut)].iloc[0]
		ax.set_title(f'{RUNS[x]}\numbral 0.{runs[x]["t"]} · {int(r.hfts)} HFTs · R² {r.r2:.2f} · JM medio {r.jm_mean:.2f}, '
		             f'mín. {r.jm_min:.2f}' + ('' if x == 'A' else f' · ARI {r.ari_A:.2f}'), fontsize=8.5, loc='left')
		ax.set_axis_off()
	fig.suptitle(f'HFTs del corte {cut} de cada ejecución', x=0.01, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.01, 0.9, 'Cada HFT tiene el color del HFT sin recorte con el que comparte más píxeles. '
	         'ARI: acuerdo con el mapa sin recorte (1 = el mismo mapa)', color=ca.INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.9))
	fig.savefig(os.path.join(HERE, f'figC3_mapas_{cut}.png'), dpi=130)
	plt.close(fig)


if __name__ == '__main__':
	pd.set_option('display.width', 220)
	runs = {x: load(x) for x in RUNS if os.path.isdir(os.path.join(DATA, x))}
	rows = []
	for x, r in runs.items():
		s = r['s']
		r2_segm = cp.eta2(np.arange(len(s.n)), s, np.eye(s.nv)).mean()
		for c, cl in r['hfts'].items():
			ari = np.nan
			if x != 'A' and c in runs['A']['hfts']:
				a, b = pixel_labels(runs['A'], c).ravel(), pixel_labels(r, c).ravel()
				both = (a >= 0) & (b >= 0)
				ari = ca.ari(a[both], b[both], np.ones(both.sum()))
			ev = cp.eta2(cl, s, np.eye(s.nv))
			rows.append({'run': x, 'threshold': r['t'] / 10, 'segments': len(s.n), 'r2_segm': r2_segm, 'cut': c,
			             **ca.measures(cl, s), 'ari_A': ari, **{f'eta2_{v}': e for v, e in zip(VARS, ev)}})
	tab = pd.DataFrame(rows)
	tab.to_csv(os.path.join(HERE, 'clip_medidas.csv'), index=False)
	print(tab[['run', 'threshold', 'segments', 'r2_segm', 'cut', 'hfts', 'r2', 'jm_mean', 'jm_min', 'sil', 'n_eff', 'ari_A']]
	      .round(3).to_string(index=False))
	print('\neta2 of every variable:')
	print(tab[['run', 'cut', 'hfts'] + [f'eta2_{v}' for v in VARS]].round(2).to_string(index=False))
	fig_thresholds(runs)
	fig_r2(runs)
	for c in ['p25', 'p50']:
		if all(c in r['hfts'] for r in runs.values()):
			fig_eta2(tab, c)
			fig_maps(runs, tab, c)
	print('\nFiguras y tablas en', HERE)
