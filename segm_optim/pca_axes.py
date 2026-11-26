#### eHabitat+: the axes of the PCA that the segmentation uses (SEGM_PCA, SEGM_PCA_AXES, SEGM_PCA_VARIANCE)
#### License: CC BY-SA 3.0

"""
With SEGM_PCA = True, segmentation_pca_par.py segments axes of a standardized PCA of the variables of each area
(i.pca -n): the first SEGM_PCA_AXES when it is a number of axes (1 or more; all of them if the area has fewer), else
the first DEFAULT_AXES (3) when SEGM_PCA_VARIANCE is 0 (the default), or else the first ones that together explain at
least SEGM_PCA_VARIANCE % of the variance (from MIN_VARIANCE to 100; 100 = all of them). SEGM_PCA_AXES 0 (the default)
leaves the choice to SEGM_PCA_VARIANCE, which a number of axes overrides. It writes csv/pca_axes_<area>.csv with every
axis (eigenvalue, % of the variance, cumulative %, whether it is used), and says in its log how much of the variance
the axes used explain; moranvar_plots.R copies it to the results and writes it in processing_metadata.txt.

Before the segmentation, pca_preview.py computes the same PCA from the GeoTIFFs, without GRASS, with how much each
variable loads on each axis, to choose the axes first (the app, in HFT mapping, and its command line).
"""

import csv
import math
import re

DEFAULT_AXES = 3
MIN_VARIANCE = 50
TOL = 1e-6  # % printed by i.pca with two decimals: a cumulative % within this of the target reaches it
LINE = re.compile(r'^\s*PC(\d+)\s+(\S+)\s+\(.*\)\s+\[\s*(-?[0-9.]+)%\]\s*$')


def variance_option(cfg):
	"""SEGM_PCA_VARIANCE of config.py (CFG, a dict) as a number: 0 (the first DEFAULT_AXES axes, also without the option
	or with None) or from MIN_VARIANCE to 100; ValueError saying what is wrong otherwise."""
	v = cfg.get('SEGM_PCA_VARIANCE', 0)
	if v is None:
		return 0
	if isinstance(v, bool) or not isinstance(v, (int, float)) or v != v:
		raise ValueError(f'SEGM_PCA_VARIANCE in config.py must be 0 or a % from {MIN_VARIANCE} to 100, not {v!r}')
	if v != 0 and not MIN_VARIANCE <= v <= 100:
		raise ValueError(f'SEGM_PCA_VARIANCE in config.py must be 0 (the first {DEFAULT_AXES} axes) or a % from '
		                 f'{MIN_VARIANCE} to 100, not {v!r}')
	return v


def axes_option(cfg):
	"""SEGM_PCA_AXES of config.py (CFG, a dict) as a whole number: 0 (the axes of SEGM_PCA_VARIANCE, also without the
	option or with None) or a number of axes, 1 or more; ValueError saying what is wrong otherwise."""
	v = cfg.get('SEGM_PCA_AXES', 0)
	if v is None:
		return 0
	if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 or v != int(v):
		raise ValueError(f'SEGM_PCA_AXES in config.py must be 0 (the axes of SEGM_PCA_VARIANCE) or a number of axes, 1 or '
		                 f'more, not {v!r}')
	return int(v)


def parse(text):
	"""The axes of the output of i.pca: [(axis 1.., eigenvalue, % of the variance)], in order."""
	out = []
	for line in str(text).splitlines():
		m = LINE.match(line)
		if m:
			try:
				eig = float(m.group(2))
			except ValueError:
				eig = float('nan')
			out.append((int(m.group(1)), eig, float(m.group(3))))
	return sorted(out)


def choose(percents, variance=0, axes=0):
	"""How many axes to use: the first AXES (or all, if there are fewer) when AXES is a number of axes, else the first
	DEFAULT_AXES (or all, if there are fewer) with VARIANCE 0, else the fewest first axes whose cumulative % of the
	variance reaches VARIANCE (all of them if it is never reached)."""
	n = len(percents)
	if axes:
		return min(int(axes), n)
	if not variance:
		return min(DEFAULT_AXES, n)
	if variance >= 100:  # all of them (the % of i.pca, rounded to two decimals, may add up to 99.99 or 100.01)
		return n
	cum = 0.0
	for k, p in enumerate(percents, 1):
		cum += p
		if cum >= variance - TOL:
			return k
	return n


def summary(axes, k, variance=0, n_axes=0):
	"""The sentence of the log and of processing_metadata.txt about the axes used (N_AXES: SEGM_PCA_AXES)."""
	cum = sum(p for _, _, p in axes[:k])
	used = 'PC1' if k == 1 else f'PC1-PC{k}'
	why = (f'SEGM_PCA_AXES = {n_axes}' + (f', but the area has {len(axes)} axes' if k < n_axes else '') if n_axes
	       else f'the first {DEFAULT_AXES} axes, SEGM_PCA_VARIANCE = 0' if not variance and k == DEFAULT_AXES
	       else f'all the axes, fewer than {DEFAULT_AXES}' if not variance
	       else f'the first ones that explain at least {variance:g} %')
	return f'{used} of {len(axes)} axes, {cum:.1f} % of the variance ({why})'


def write(path, axes, k):
	"""csv/pca_axes_<area>.csv: axis, eigenvalue, % of the variance, cumulative % and used (1/0)."""
	cum = 0.0
	with open(path, 'w', newline='') as f:
		w = csv.writer(f, delimiter=' ')
		w.writerow(['axis', 'eigenvalue', 'variance_pct', 'cumulative_pct', 'used'])
		for i, (n, eig, p) in enumerate(axes, 1):
			cum += p
			w.writerow([f'PC{n}', eig, p, round(cum, 2), int(i <= k)])

