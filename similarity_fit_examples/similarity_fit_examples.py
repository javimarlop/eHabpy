"""
Examples of the fit of the similarity in the study area itself (similarity.py, ECOSIM_COVARIANCE 'mcd' and 'classic'):
Guajares_zona_quemada_dNBR2022_10m / gee_guajares / p25, 100 m, compared with its ecoregion. The two runs were made in
an isolated copy of the repository (scratchpad simfit/), not in segm_optim/results.

  1. fit_maps.png: the HFTs of the study area on the 100 m grid of the similarity (the HFT that covers at least half of
     each pixel, as similarity.py does) and, for three HFTs, their similarity map (Mahalanobis, 0-1) in the study area
     and what the fit counts: pixels of the HFT similar to it (>= its threshold T), pixels of the HFT not similar, and
     pixels of other HFTs similar to it.
  2. fit_indicator.png: SelfSim_pct, StudySim_pct and SelfFit_pct of every HFT, robust (MCD) and classic covariance.
Each PNG goes with a CSV of the data drawn. The counts of the maps are checked against fit_matrix.csv of similarity.py.

    /opt/anaconda3/envs/ehabpy/bin/python -s similarity_fit_examples.py <scratchpad>/simfit
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from osgeo import gdal

gdal.UseExceptions()
SIMFIT = sys.argv[1]
REPO = os.path.join(SIMFIT, 'repo')
sys.path.insert(0, os.path.join(REPO, 'segm_optim'))
import similarity as S  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(REPO, 'segm_optim/results/Guajares_zona_quemada_dNBR2022_10m/gee_guajares/p25')
SHP = os.path.join(P, 'park_segm_1_7_p25_optim_class.shp')
COLOURS = pd.read_csv(os.path.join(P, 'hft_colours_1_7_p25_optim.csv'), sep=' ')
SHOW = [5, 2, 1]  # the best fit (but HFT 6, 76 pixels), the worst and one in between
INK, INK2 = '#0b0b0b', '#52514e'
CAT = ['#008300', '#2a78d6', '#eb6834']  # HFT similar, HFT not similar, other HFT similar
plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2})


def read(path, win=None):
    ds = gdal.Open(path)
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray(*win) if win else b.ReadAsArray()
    return a.astype(float), b.GetNoDataValue(), ds


def load(cov):
    d = os.path.join(SIMFIT, cov)
    hri = pd.read_csv(os.path.join(d, 'hri_results.csv'), sep=' ', na_values=['NA', 'None'])
    fitm = pd.read_csv(os.path.join(d, 'fit_matrix.csv'), sep=' ')
    return d, hri, fitm


def main():
    d, hri, fitm = load('mcd')
    _, _, ds = read(os.path.join(d, 'sim_1_hft1.tif'))
    grid = S.Grid(ds.GetGeoTransform(), ds.RasterXSize, ds.RasterYSize, ds.GetProjection())
    cov = S.hft_cover(SHP, grid, 10.0, 0.5)
    r0, c0, nr, nc = cov['window']
    win = (c0, r0, nc, nr)
    assign = cov['assign']
    hfts = [h for h in cov['hfts']]
    gt = grid.gt
    ext = (gt[0] + c0 * gt[1], gt[0] + (c0 + nc) * gt[1], gt[3] - (r0 + nr) * gt[1], gt[3] - r0 * gt[1])
    T = dict(zip(hri.hft, hri.AveHFTSim))

    # check: the pixels of each HFT k similar to h (sim >= T_h) are those of fit_matrix.csv
    sims, worst = {}, 0
    for h in hfts:
        s, nd, _ = read(os.path.join(d, f'sim_1_hft{h}.tif'), win)
        s[s == nd] = np.nan
        sims[h] = s
        for k in hfts:
            m = assign == k
            ok = m & np.isfinite(s)
            n_sim = int((s[ok] >= T[h]).sum())
            ref = fitm[(fitm.hft_ref == h) & (fitm.hft == k)]
            if len(ref):
                worst = max(worst, abs(n_sim - int(ref.NrSimPix.iloc[0])), abs(int(ok.sum()) - int(ref.NrPixHFT.iloc[0])))
    print(f'largest difference with fit_matrix.csv (pixels): {worst}')

    # 1. maps
    col = dict(zip(COLOURS.HFT, COLOURS.colour))
    order = list(COLOURS.HFT)
    fig, axs = plt.subplots(len(SHOW) + 1, 2, figsize=(10, 3.25 * (len(SHOW) + 1)), constrained_layout=True)
    inside = assign > 0
    lut = ListedColormap([col[h] for h in order])
    idx = np.full(assign.shape, np.nan)
    for i, h in enumerate(order):
        idx[assign == h] = i
    ax = axs[0, 0]
    ax.imshow(idx, cmap=lut, vmin=-0.5, vmax=len(order) - 0.5, extent=ext, interpolation='nearest')
    ax.set_title('HFTs of the study area (100 m pixels covered at least half by an HFT)', loc='left', fontsize=9.5)
    hand = [plt.Rectangle((0, 0), 1, 1, color=col[h]) for h in order]
    ax.legend(hand, [f'HFT {h}' for h in order], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False,
              fontsize=8.5)
    axs[0, 1].axis('off')
    rows = []
    for i, h in enumerate(SHOW, start=1):
        r = hri[hri.hft == h].iloc[0]
        s = sims[h]
        ax = axs[i, 0]
        im = ax.imshow(np.where(inside, s, np.nan), cmap='Blues', vmin=0, vmax=1, extent=ext, interpolation='nearest')
        ax.contour(np.flipud(assign == h).astype(float), levels=[0.5], colors=INK, linewidths=0.8,
                   extent=ext, origin='lower')
        ax.set_title(f'HFT {h}: similarity (Mahalanobis, 0-1) in the study area\nblack line: the HFT; threshold '
                     f'T = {T[h]:.2f}', loc='left', fontsize=9.5)
        fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, label='similarity')
        cat = np.full(assign.shape, np.nan)
        sim_ok = np.isfinite(s) & (s >= T[h])
        cat[(assign == h) & sim_ok] = 0
        cat[(assign == h) & ~sim_ok] = 1
        cat[inside & (assign != h) & sim_ok] = 2
        ax = axs[i, 1]
        ax.imshow(np.where(inside, 0, np.nan), cmap=ListedColormap(['#e4e3df']), extent=ext, interpolation='nearest')
        ax.imshow(cat, cmap=ListedColormap(CAT), vmin=-0.5, vmax=2.5, extent=ext, interpolation='nearest')
        ax.set_title(f'HFT {h}: SelfSim {r.SelfSim_pct:.0f} %, StudySim {r.StudySim_pct:.0f} %, '
                     f'SelfFit {r.SelfFit_pct:.0f} %', loc='left', fontsize=9.5)
        n = [int(np.nansum(cat == v)) for v in range(3)]
        rows += [{'hft': h, 'threshold_T': T[h], 'category': c, 'pixels': v, 'ha': v} for c, v in
                 zip(['HFT similar to itself', 'HFT not similar', 'other HFT similar'], n)]
    hand = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CAT] + [plt.Rectangle((0, 0), 1, 1, color='#e4e3df')]
    axs[0, 1].legend(hand, ['pixels of the HFT similar to it', 'pixels of the HFT not similar',
                            'pixels of other HFTs similar to it', 'other pixels of the study area'],
                     loc='upper left', frameon=False, fontsize=9, title='Right-hand maps (similarity >= T)',
                     title_fontsize=9.5, alignment='left')
    axs[0, 1].text(0.0, 0.42, 'SelfSim = green / (green + blue)\nStudySim = (green + orange) / study area\n'
                   'SelfFit = green / (green + orange)', transform=axs[0, 1].transAxes, va='top', fontsize=9,
                   color=INK, linespacing=1.6)
    for ax in axs.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
    x0, y0 = ext[0] + 300, ext[2] + 300
    axs[0, 0].plot([x0, x0 + 2000], [y0, y0], color=INK, lw=3, solid_capstyle='butt')
    axs[0, 0].text(x0 + 1000, y0 + 200, '2 km', ha='center', fontsize=8.5, color=INK)
    fig.suptitle('Similarity of each HFT in its own study area (Guajares, p25, 100 m, robust covariance MCD)',
                 x=0.01, ha='left', fontsize=10.5, color=INK)
    fig.savefig(os.path.join(OUT, 'fit_maps.png'), dpi=170)
    plt.close(fig)
    pd.DataFrame(rows).assign(ha=lambda t: t.pixels * (gt[1] ** 2) / 1e4).to_csv(
        os.path.join(OUT, 'fit_maps.csv'), sep=' ', index=False, float_format='%.4f')

    # 2. indicator per HFT, MCD and classic
    meas = [('SelfSim_pct', 'SelfSim: % of the HFT\nsimilar to itself'),
            ('StudySim_pct', 'StudySim: % of the study area\nsimilar to the HFT'),
            ('SelfFit_pct', 'SelfFit: % of that similar area\nthat is the HFT itself')]
    tab = []
    for cv, c in (('mcd', '#2a78d6'), ('classic', '#eb6834')):
        _, h2, _ = load(cv)
        h2 = h2.assign(covariance=cv)[['hft', 'covariance', 'NrPixHFT', 'AveHFTSim', 'Covariance'] + [m for m, _ in meas]]
        tab.append(h2)
    tab = pd.concat(tab)
    fig, axs = plt.subplots(1, 3, figsize=(11, 4.2), sharey=True, constrained_layout=True)
    ypos = {h: i for i, h in enumerate(order[::-1])}
    for ax, (m, lab) in zip(axs, meas):
        for cv, c, dy in (('mcd', '#2a78d6', 0.13), ('classic', '#eb6834', -0.13)):
            t = tab[tab.covariance == cv]
            ax.scatter(t[m], [ypos[h] + dy for h in t.hft], s=42, color=c, edgecolor='white', linewidth=1.5,
                       zorder=3, label='robust (MCD)' if cv == 'mcd' else 'classic')
        ax.set_title(lab, loc='left', fontsize=9.5)
        ax.set_xlim(0, 100)
        ax.grid(axis='x', color='#e4e3df', lw=0.8)
        ax.set_axisbelow(True)
        for s_ in ('top', 'right'):
            ax.spines[s_].set_visible(False)
        ax.set_xlabel('%')
    axs[0].set_yticks(list(ypos.values()))
    axs[0].set_yticklabels([f'HFT {h} ({int(tab[(tab.hft == h)].NrPixHFT.iloc[0])} px)' for h in ypos])
    axs[2].legend(loc='lower left', frameon=False, fontsize=9)
    fig.suptitle('Fit of the similarity in the study area, per HFT (Guajares, p25, 100 m)', x=0.01, ha='left',
                 fontsize=10.5, color=INK)
    fig.savefig(os.path.join(OUT, 'fit_indicator.png'), dpi=170)
    plt.close(fig)
    tab.to_csv(os.path.join(OUT, 'fit_indicator.csv'), sep=' ', index=False, float_format='%.4f')
    print('done')


if __name__ == '__main__':
    main()
