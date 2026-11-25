"""Compara la clasificación de los segmentos con Ward ponderado sobre los tres primeros ejes de una PCA de las
variables (como la entrada de la segmentación con SEGM_PCA = True) con la actual (CLASS_METHOD = 'ward_area',
las siete variables tipificadas).

    Rscript segment_sums.R            # una vez: sumas de los píxeles de cada segmento y PCA en data/
    python compare_pca.py [UMBRALES] [UMBRAL_FIGURAS]

La PCA es la de los píxeles tipificados (data/pca_0k.csv, como i.pca -n); cada segmento se sitúa en ella por
su media (la media de las puntuaciones de sus píxeles). Métodos, todos con Ward ponderado por el nº de píxeles:
  ward_area   referencia: las siete variables tipificadas (cada variable el mismo peso)
  pca3        PC1-PC3 con su varianza: quita los ejes 4 a 7, pero PC1 sigue pesando más
  pca3_equal  PC1-PC3 con la misma desviación típica: cada uno de los tres gradientes el mismo peso
  pca3_range  PC1-PC3 con el mismo rango de los píxeles, como en la segmentación (i.pca los reescala a 0-255
              e i.segment escala cada capa a su rango)
Medidas (con las de compare_algorithms.py): r2 (las siete variables), r2_pc (los tres ejes, con el mismo peso),
jm_mean y jm_min (las siete variables), n_eff, ari_ref (con ward_area con el mismo nº de HFTs) y eta2 de cada
variable y de cada eje (parte de su varianza de los píxeles que explican los HFTs).
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import compare_algorithms as ca

HERE, DATA, KMAX = ca.HERE, ca.DATA, ca.KMAX
METHODS = ['ward_area', 'pca3', 'pca3_equal', 'pca3_range']
LABEL = {'ward_area': 'Actual (ward_area): 7 variables tipificadas', 'pca3': 'PC1-PC3 con su varianza',
         'pca3_equal': 'PC1-PC3 con el mismo peso (misma desv. típica)',
         'pca3_range': 'PC1-PC3 con el mismo rango (como la segmentación)'}
COL = {'ward_area': ca.INK, 'pca3': ca.SLOTS[0], 'pca3_equal': ca.SLOTS[1], 'pca3_range': ca.SLOTS[2]}
MARK = {'ward_area': 'o', 'pca3': 's', 'pca3_equal': '^', 'pca3_range': 'D'}
CUTS = [('p10', .1), ('p25', .25), ('p50', .5), ('p75', .75), ('p90', .9)]


def cross(s):
	"""Matriz de sumas de productos cruzados de los píxeles tipificados de toda el área."""
	C = np.zeros((s.nv, s.nv))
	q = s.Q.sum(0)
	C[s.pairs[:, 0], s.pairs[:, 1]] = q
	C[s.pairs[:, 1], s.pairs[:, 0]] = q
	return C


def eta2(cl, s, W):
	"""Parte de la varianza de los píxeles de cada combinación lineal (columnas de W) de las variables
	tipificadas que explican los grupos cl de los segmentos: W = identidad, cada variable; W = cargas, cada eje."""
	_, cl = np.unique(cl, return_inverse=True)
	k = cl.max() + 1
	Sh, Nh = ca.rowsum(s.S, cl, k) @ W, np.bincount(cl, s.n)
	St = s.S.sum(0) @ W
	tot = np.einsum('ij,ik,kj->j', W, cross(s), W) - St ** 2 / s.N
	return ((Sh ** 2 / Nh[:, None]).sum(0) - St ** 2 / s.N) / tot


def run(t):
	s = ca.Segments(t)
	p = pd.read_csv(os.path.join(DATA, f'pca_0{t}.csv'))
	V = p[[f'load{i}' for i in range(1, s.nv + 1)]].to_numpy().T[:, :3]  # cargas de PC1-PC3 (columnas)
	Z = s.M @ V  # puntuación media de los píxeles de cada segmento en PC1-PC3
	feats = {'ward_area': s.M, 'pca3': Z, 'pca3_equal': Z / p['sd'].to_numpy()[:3],
	         'pca3_range': Z / p['range'].to_numpy()[:3]}
	ks, ns = list(range(2, KMAX + 1)), len(s.n)
	parts, cuts, rows = {}, {}, []
	for mth, X in feats.items():
		tree = ca.agnes(X, s.n, 'ward')
		for k, cl in ca.cut_tree(tree, ns, ks).items():
			parts[(mth, k)] = cl
		cuts[mth] = {c: int(ca.cut_height(tree, ns, f)) for c, f in CUTS}
	for (mth, k), cl in parts.items():
		ev, ep = eta2(cl, s, np.eye(s.nv)), eta2(cl, s, V)
		rows.append({'threshold': t, 'method': mth, 'k': k, **ca.measures(cl, s), 'r2_pc': ep.mean(),
		             'ari_ref': ca.ari(cl, parts[('ward_area', k)], s.n),
		             **{f'eta2_{v}': e for v, e in zip(VARS, ev)}, **{f'eta2_PC{i}': e for i, e in enumerate(ep, 1)}})
	tab = pd.DataFrame(rows)
	tab.to_csv(os.path.join(HERE, f'pca_medidas_0{t}.csv'), index=False)
	return s, p, tab, parts, cuts


# --- Figuras ---

def mark_cuts(ax, cuts):
	for c in ['p25', 'p50']:
		ax.axvline(cuts['ward_area'][c], color=ca.MUTED, lw=0.8, ls='--', zorder=0)
		ax.text(cuts['ward_area'][c], 1, f' {c}', transform=ax.get_xaxis_transform(), color=ca.MUTED, fontsize=7.5, va='top')


def fig_k(tab, t, ns, cuts):
	title = {'r2': 'R² de los píxeles: las 7 variables', 'r2_pc': 'R² de los píxeles: PC1-PC3 con el mismo peso',
	         'jm_mean': 'JM medio entre pares de HFTs (máx. 2)', 'jm_min': 'JM mínimo: el par de HFTs menos separable',
	         'n_eff': 'Nº efectivo de HFTs (por su área)', 'ari_ref': 'Acuerdo con el método actual (ARI por área)'}
	fig, axs = plt.subplots(3, 2, figsize=(11, 11.5))
	for ax, m in zip(axs.flat, title):
		for mth in METHODS:
			g = tab[tab.method == mth].sort_values('k')
			ax.plot(g.k, g[m], color=COL[mth], lw=2.4 if mth == 'ward_area' else 1.5, marker=MARK[mth], ms=4,
			        label=LABEL[mth], zorder=3 if mth == 'ward_area' else 2)
		mark_cuts(ax, cuts)
		ax.set_title(title[m], loc='left')
		ax.set_xticks(range(2, KMAX + 1, 2))
		ax.set_xlabel('Número de HFTs')
	h, l = axs[0, 0].get_legend_handles_labels()
	fig.legend(h, l, loc='lower center', ncol=2, fontsize=9)
	fig.suptitle(f'Clasificación con las 7 variables o con PC1-PC3 (umbral 0.{t}, {ns} segmentos)',
	             x=0.06, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.06, 0.952, 'Ward ponderado por área en todos los casos. Líneas discontinuas: cortes p25 y p50 del método actual',
	         color=ca.INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0.05, 1, 0.95))
	fig.savefig(os.path.join(HERE, f'figP1_k_0{t}.png'), dpi=150)
	plt.close(fig)


def fig_eta2(tab, p, t, k):
	"""eta2 de cada variable y de cada eje con k HFTs, por método."""
	load = p[[f'load{i}' for i in range(1, len(VARS) + 1)]].to_numpy()[:3] * p['sd'].to_numpy()[:3, None]
	axes_lab = []
	for i in range(3):  # cada eje con sus variables de mayor carga (correlación)
		top = [f'{VARS[j]} {load[i, j]:+.2f}' for j in np.argsort(-abs(load[i]))[:2]]
		axes_lab.append(f'PC{i + 1} ({p.eigenvalue[i] / p.eigenvalue.sum():.0%}): ' + ', '.join(top))
	rows = [f'eta2_{v}' for v in VARS] + [f'eta2_PC{i}' for i in (1, 2, 3)]
	labels = list(VARS) + axes_lab
	fig, ax = plt.subplots(figsize=(10, 6.4))
	g = tab[tab.k == k].set_index('method')
	y = np.arange(len(rows), dtype=float)
	y[len(VARS):] += 0.8  # separa las variables de los ejes
	for i, r in enumerate(rows):
		vals = [g.loc[m, r] for m in METHODS]
		ax.plot([min(vals), max(vals)], [y[i]] * 2, color=ca.GRID, lw=3, solid_capstyle='round', zorder=1)
	for j, mth in enumerate(METHODS):
		ax.scatter([g.loc[mth, r] for r in rows], y, s=60 if mth == 'ward_area' else 46, color=COL[mth], marker=MARK[mth],
		           edgecolors='#fcfcfb', linewidths=1.2, zorder=3 if mth == 'ward_area' else 2, label=LABEL[mth])
	ax.set_yticks(y, labels)
	ax.invert_yaxis()
	ax.set_xlim(0, 1)
	ax.grid(axis='y', visible=False)
	ax.set_xlabel('eta² (parte de la varianza de los píxeles que explican los HFTs)')
	ax.legend(loc='upper right', fontsize=8.5)
	fig.suptitle(f'Qué explican {k} HFTs de cada método (umbral 0.{t})', x=0.01, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.01, 0.925, 'Arriba, cada variable; abajo, cada eje de la PCA de los píxeles (entre paréntesis, su varianza; '
	         'después, sus dos variables más correlacionadas)', color=ca.INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.92))
	fig.savefig(os.path.join(HERE, f'figP2_eta2_0{t}_k{k}.png'), dpi=150)
	plt.close(fig)


def fig_maps(s, tab, parts, t, k):
	"""Mapas de los HFTs de cada método con k HFTs, cada HFT con el color del HFT de ward_area con el que comparte
	más área (como fig_maps de compare_algorithms.py)."""
	nr, nc = pd.read_csv(os.path.join(DATA, 'grid.csv')).iloc[0]
	seg = np.fromfile(os.path.join(DATA, f'seg_0{t}.bin'), dtype=np.int32).reshape(nr, nc)
	rowmap = np.full(seg.max() + 1, -1)
	rowmap[s.rows + 1] = np.arange(len(s.rows))
	ok = np.argwhere(seg > 0)
	(r0, c0), (r1, c1) = ok.min(0), ok.max(0) + 1
	seg = seg[r0:r1, c0:c1]
	ref = parts[('ward_area', k)]
	cols = ca.SLOTS + ['#86b6ef', '#f4a582', '#7fd8b4', '#f6cd73', '#9085e9', '#c3c2b7']
	cmap = ListedColormap(cols[:max(k, 2)] if k <= len(cols) else plt.get_cmap('tab20').colors[:k])
	fig, axs = plt.subplots(2, 2, figsize=(11.5, 9.2))
	for ax, mth in zip(axs.flat, METHODS):
		p = parts[(mth, k)]
		T = np.zeros((p.max() + 1, ref.max() + 1))
		np.add.at(T, (p, ref), s.n)
		lab, used = np.full(p.max() + 1, -1), set()
		for i, j in sorted(np.ndindex(T.shape), key=lambda ij: -T[ij]):
			if lab[i] < 0 and j not in used:
				lab[i] = j
				used.add(j)
		for i in np.flatnonzero(lab < 0):
			lab[i] = min(set(range(k)) - used)
			used.add(lab[i])
		img = lab[p[np.maximum(rowmap[seg], 0)]].astype(float)
		img[(seg == 0) | (rowmap[seg] < 0)] = np.nan
		ax.imshow(img, cmap=cmap, vmin=-0.5, vmax=cmap.N - 0.5, interpolation='nearest')
		r = tab[(tab.method == mth) & (tab.k == k)].iloc[0]
		ax.set_title(f'{LABEL[mth]}\nR² {r.r2:.2f} (PC1-PC3 {r.r2_pc:.2f}) · JM medio {r.jm_mean:.2f}, mín. {r.jm_min:.2f}'
		             + ('' if mth == 'ward_area' else f' · ARI {r.ari_ref:.2f}'), fontsize=8.5, loc='left')
		ax.set_axis_off()
	fig.suptitle(f'{k} HFTs con cada método (umbral 0.{t})', x=0.01, ha='left', fontweight='bold', fontsize=11)
	fig.text(0.01, 0.935, 'Cada HFT tiene el color del HFT del método actual con el que comparte más área. '
	         'ARI: acuerdo con el método actual (1 = el mismo mapa)', color=ca.INK2, fontsize=8.5)
	fig.tight_layout(rect=(0, 0, 1, 0.93), h_pad=3)
	fig.savefig(os.path.join(HERE, f'figP3_mapas_0{t}_k{k}.png'), dpi=130)
	plt.close(fig)


VARS = list(pd.read_csv(os.path.join(DATA, 'vars.csv')).short)

if __name__ == '__main__':
	pd.set_option('display.width', 220)
	tabs, keep = [], {}
	for t in ca.THR:
		s, p, tab, parts, cuts = run(t)
		tabs.append(tab)
		keep[t] = (s, p, tab, parts, cuts)
		print(f'\n=== Threshold 0.{t}: {len(s.n)} segments | HFTs of the cuts p10, p25, p50, p75, p90 of each tree:')
		for mth in METHODS:
			print(f'  {mth:<11} {list(cuts[mth].values())}')
		fig_k(tab, t, len(s.n), cuts)
	tab = pd.concat(tabs)

	# Diferencias con ward_area con el mismo nº de HFTs (umbrales y 2 a KMAX HFTs)
	cols = ['r2', 'r2_pc', 'jm_mean', 'jm_min', 'n_eff'] + [f'eta2_{v}' for v in VARS]
	ref = tab[tab.method == 'ward_area'].set_index(['threshold', 'k'])[cols]
	out = []
	for mth in METHODS[1:]:
		g = tab[tab.method == mth].set_index(['threshold', 'k'])[cols]
		d = g - ref.loc[g.index]
		row = {'method': mth, 'ari_ref': tab[tab.method == mth].ari_ref.median()}
		for c in cols:
			row['d_' + c] = d[c].median()
		for c in ['r2', 'r2_pc', 'jm_mean', 'jm_min']:
			row['better_' + c] = np.mean(d[c] > 1e-9)
		out.append(row)
	summ = pd.DataFrame(out)
	summ.to_csv(os.path.join(HERE, 'pca_resumen.csv'), index=False)
	print(f'\nMedian difference with ward_area at the same number of HFTs (thresholds {ca.THR}, 2 to {KMAX} HFTs)')
	print(summ[['method', 'ari_ref'] + ['d_' + c for c in ['r2', 'r2_pc', 'jm_mean', 'jm_min', 'n_eff']] +
	           ['better_' + c for c in ['r2', 'r2_pc', 'jm_mean', 'jm_min']]].round(3).to_string(index=False))
	print('\nMedian difference of the eta2 of each variable:')
	print(summ[['method'] + [f'd_eta2_{v}' for v in VARS]].round(3).to_string(index=False))

	if ca.TMAP in keep:
		s, p, tab, parts, cuts = keep[ca.TMAP]
		for k in sorted({cuts['ward_area']['p25'], cuts['ward_area']['p50']}):
			fig_eta2(tab, p, ca.TMAP, k)
			print(f'\nThreshold 0.{ca.TMAP}, {k} HFTs:')
			print(tab[tab.k == k][['method', 'r2', 'r2_pc', 'jm_mean', 'jm_min', 'n_eff', 'ari_ref'] +
			                      [f'eta2_{v}' for v in VARS]].round(2).to_string(index=False))
		fig_maps(s, tab, parts, ca.TMAP, cuts['ward_area']['p25'])
	print('\nFiguras y tablas en', HERE)
