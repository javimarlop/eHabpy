#### eHabitat+: landscape metrics of a map of classes (patch, class and landscape levels), as the R package
#### landscapemetrics computes them (Hesselbarth et al. 2019, Ecography 42: 1648-1657), which follows FRAGSTATS
#### (McGarigal et al. 2023, https://www.fragstats.org)
#### License: CC BY-SA 3.0

"""
The 133 metrics of landscapemetrics 2.2 (list_lsm()), in numpy and scipy so that they need no R package (landscapemetrics
is not in conda-forge) and run on large grids: the similar areas of the HFTs in their ecoregion (similarity.py) and the
patches of the user (landscape_metrics.py: rasters or polygons, e.g. a habitat at several dates).

    import lsm
    out = lsm.compute(classes, valid, res=100)              # every metric; classes: integers, valid: cells with data
    out = lsm.compute(classes, valid, 100, ['c_pland', 'c_enn_mn', 'l_shdi'], directions=4, edge_depth=2)

The defaults are those of landscapemetrics (calculate_lsm): patches of 8 neighbours (directions), an edge depth of 1
cell for the core areas, cells at the edge of the landscape (or next to cells without data) taken as edge for the core
areas (consider_boundary = FALSE), the edges with cells without data or outside the map not counted in the total edge
(count_boundary = FALSE), the 4 neighbours of each cell for the adjacencies (neighbourhood = 4, ordered pairs) and
entropies in bits (base = 'log2'). The landscape is every cell with data: PLAND, PD, ED, LPI, DIVISION, MESH... are
relative to it. Tests: segm_optim/tests/test_lsm.py (values of landscapemetrics 2.2.1 for small maps).

Each metric is level_metric: p_ (one value per patch), c_ (one per class) or l_ (one for the landscape), as the
functions lsm_p_area, lsm_c_ai, lsm_l_shdi... of landscapemetrics. CATALOGUE describes them (name, type, unit, what it
measures and how to read it), for the tables and the guides of the app.
"""

import math

import numpy as np

# Catalogue ------------------------------------------------------------------------------------------------------------

# The metrics of every patch, with their statistics of the classes and the landscape (_mn, _sd, _cv). dcore: the
# statistics of the disjunct core areas, which landscapemetrics computes from ncore (the number of core areas of each
# patch), as done here.
PATCH = ('area', 'cai', 'circle', 'contig', 'core', 'enn', 'frac', 'gyrate', 'ncore', 'para', 'perim', 'shape')
STATS = ('area', 'cai', 'circle', 'contig', 'core', 'dcore', 'enn', 'frac', 'gyrate', 'para', 'shape')
CLASS = ('ai', 'ca', 'clumpy', 'cohesion', 'cpland', 'dcad', 'division', 'ed', 'iji', 'lpi', 'lsi', 'mesh', 'ndca',
         'nlsi', 'np', 'pafrac', 'pd', 'pladj', 'pland', 'split', 'tca', 'te')
LANDSCAPE = ('ai', 'cohesion', 'condent', 'contag', 'dcad', 'division', 'ed', 'ent', 'iji', 'joinent', 'lpi', 'lsi',
             'mesh', 'msidi', 'msiei', 'mutinf', 'ndca', 'np', 'pafrac', 'pd', 'pladj', 'pr', 'prd', 'relmutinf', 'rpr',
             'shdi', 'shei', 'sidi', 'siei', 'split', 'ta', 'tca', 'te')
SUFFIX = {'mn': 'mean', 'sd': 'standard deviation', 'cv': 'coefficient of variation'}

# metric: (name, type, unit, what it measures, how to read it, fragmentation): fragmentation '+' when larger values
# mean a more fragmented (or more scattered, less connected) habitat, '-' when smaller values do, '' when neither
_BASE = {
	'area': ('Patch area', 'area and edge', 'ha', 'The area of a patch.',
	         'Smaller patches hold fewer individuals and more edge: a falling mean patch area is a sign of fragmentation.',
	         '-'),
	'cai': ('Core area index', 'core area', '%', 'The share of a patch that is core area: farther from its edge than the '
	        'edge depth.', 'Small or narrow patches have little core (0 %); large compact ones approach 100 %.', '-'),
	'circle': ('Related circumscribing circle', 'shape', '0–1', '1 − patch area / area of the smallest circle that '
	           'contains the patch (all the corners of its cells).',
	           '0 for a round patch, about 0.36 for a square one (or one cell), near 1 for long, narrow ones (corridors, '
	           'strips).', ''),
	'contig': ('Contiguity index', 'shape', '0–1', 'How connected the cells of a patch are, with a 3 × 3 template that '
	           'counts its side neighbours twice and its corner neighbours once.',
	           '0 for one cell, near 1 for compact patches; low values for thin, branching or perforated patches.', '-'),
	'core': ('Core area', 'core area', 'ha', 'The area of a patch farther from its edge than the edge depth.',
	         'The habitat free of edge effects: it shrinks faster than the area when patches break up.', '-'),
	'dcore': ('Disjunct core areas', 'core area', 'number', 'The number of separate core areas of each patch (as '
	          'landscapemetrics computes the dcore statistics, from ncore).',
	          'Patches whose core is split in several parts by narrow necks or holes have more than one.', ''),
	'enn': ('Euclidean nearest-neighbour distance', 'aggregation', 'm', 'The distance from a patch to the nearest patch '
	        'of the same class, edge to edge (between the centres of their cells). NA for a class with one patch.',
	        'Larger distances mean more isolated patches, harder to reach for organisms that disperse between them.',
	        '+'),
	'frac': ('Fractal dimension index', 'shape', '1–2', '2 ln(0.25 perimeter) / ln(area): the complexity of the edge of a '
	         'patch at every scale.', '1 for simple shapes (squares), towards 2 for very convoluted edges.', ''),
	'gyrate': ('Radius of gyration', 'area and edge', 'm', 'The mean distance between the cells of a patch and its '
	           'centroid: how far the patch extends.',
	           'Large for large or elongated patches; the mean of the class falls when the patches break up.', '-'),
	'ncore': ('Number of core areas', 'core area', 'number', 'The number of separate core areas of a patch.',
	          'More than 1 when narrow necks or holes split its core.', ''),
	'para': ('Perimeter-area ratio', 'shape', 'm/m²', 'Perimeter / area of a patch.',
	         'High for small or elongated patches (it depends on their size).', '+'),
	'perim': ('Patch perimeter', 'area and edge', 'm', 'The length of the edge of a patch, its cells next to any other '
	          'class, to cells without data or to the edge of the map included.', '', ''),
	'shape': ('Shape index', 'shape', '≥ 1', '0.25 perimeter / √area: the perimeter compared with that of a square of '
	          'the same area.', '1 for a square, larger for more complex or elongated patches.', ''),
}
_CLASS = {
	'ai': ('Aggregation index', 'aggregation', '%', 'Like adjacencies (shared sides between cells of the class) / their '
	       'maximum possible number for its area (He et al. 2000).',
	       '100 when the class is one compact patch; low when its cells are scattered.', '-'),
	'ca': ('Total (class) area', 'area and edge', 'ha', 'The area of the class.', 'Habitat loss: the amount, not its '
	       'configuration.', '-'),
	'clumpy': ('Clumpiness index', 'aggregation', '−1 to 1', 'How much more (or less) clumped the cells of the class '
	           'are than at random, given its share of the landscape.',
	           '−1 maximally dispersed, 0 random, 1 maximally clumped.', '-'),
	'cohesion': ('Patch cohesion index', 'aggregation', '0–100', 'The physical connectedness of the patches of the class, '
	             'from their perimeters and areas (Schumaker 1996).',
	             'Rises as the class becomes clumped in large, compact patches; falls as it breaks up.', '-'),
	'cpland': ('Core area percentage of landscape', 'core area', '%', 'The core area of the class / the landscape.',
	           'The share of the landscape that is interior habitat, away from edges.', '-'),
	'dcad': ('Disjunct core area density', 'core area', 'per 100 ha', 'The number of separate core areas of the class '
	         'per 100 ha of landscape.', 'More core areas per area: a more subdivided interior habitat.', '+'),
	'division': ('Division index', 'aggregation', '0–1', 'The probability that two cells chosen at random in the landscape '
	             'are not in the same patch of the class (Jaeger 2000).',
	             '0 when the class is one patch covering the landscape; near 1 when it is scattered in small patches.',
	             '+'),
	'ed': ('Edge density', 'area and edge', 'm/ha', 'The edge of the class with other classes / the landscape.',
	       'More edge per area: more patches or more irregular ones.', '+'),
	'iji': ('Interspersion and juxtaposition index', 'aggregation', '0–100', 'How evenly the class is next to every '
	        'other class (its edges with each one). NA with fewer than 3 classes.',
	        '100 when it touches every other class equally, low when it touches only a few.', ''),
	'lpi': ('Largest patch index', 'area and edge', '%', 'The largest patch of the class / the landscape.',
	        'Falls when the largest habitat block breaks up.', '-'),
	'lsi': ('Landscape shape index', 'aggregation', '≥ 1', 'The edge of the class / the minimum edge for its area (a '
	        'square).', '1 for one square patch; grows as the class splits into more patches or more irregular ones.',
	        '+'),
	'mesh': ('Effective mesh size', 'aggregation', 'ha', 'Σ patch area² / landscape area: the size of the patch in which '
	         'a cell chosen at random in the landscape lies, on average (Jaeger 2000).',
	         'The best-known measure of fragmentation: it falls when the patches break up or shrink.', '-'),
	'ndca': ('Number of disjunct core areas', 'core area', 'number', 'The number of separate core areas of the class.',
	         '', ''),
	'nlsi': ('Normalized landscape shape index', 'aggregation', '0–1', 'LSI rescaled between its minimum and maximum '
	         'for the area of the class.', '0 for one compact patch, 1 for a maximally disaggregated class.', '+'),
	'np': ('Number of patches', 'aggregation', 'number', 'The number of patches of the class.',
	       'More patches of the same habitat: subdivision (but also a new patch where there was none).', '+'),
	'pafrac': ('Perimeter-area fractal dimension', 'shape', '1–2', '2 / the slope of ln(area) on ln(perimeter) of the '
	           'patches of the class. NA with fewer than 10 patches.',
	           '1 for simple shapes, towards 2 for convoluted ones, at every patch size.', ''),
	'pd': ('Patch density', 'aggregation', 'per 100 ha', 'The number of patches of the class per 100 ha of landscape.',
	       'More patches per area: more subdivided.', '+'),
	'pladj': ('Percentage of like adjacencies', 'aggregation', '%', 'The share of the cell sides of the class that touch '
	          'the same class.', 'High when the class is clumped; low when it is scattered.', '-'),
	'pland': ('Percentage of landscape', 'area and edge', '%', 'The area of the class / the landscape.',
	          'Habitat amount: below 10–30 % the configuration matters most for many species.', '-'),
	'split': ('Splitting index', 'aggregation', '≥ 1', 'Landscape area² / Σ patch area²: the number of patches of equal '
	          'size that would give the same division (Jaeger 2000).',
	          '1 for one patch covering the landscape; grows as the class breaks up.', '+'),
	'tca': ('Total core area', 'core area', 'ha', 'The core area of the class.', '', '-'),
	'te': ('Total edge', 'area and edge', 'm', 'The length of the edges of the class with other classes (not with cells '
	       'without data nor with the edge of the map).', '', '+'),
}
_LANDSCAPE = {
	'ai': ('Aggregation index', 'aggregation', '%', 'The aggregation index of each class weighted by its share of the '
	       'landscape.', '', '-'),
	'cohesion': ('Patch cohesion index', 'aggregation', '0–100', 'The connectedness of all the patches.', '', '-'),
	'condent': ('Conditional entropy', 'complexity', 'bits', 'The entropy of the class of a cell when the class of its '
	            'neighbour is known (Nowosad & Stepinski 2019).',
	            'The configurational complexity: high when the classes are finely intermixed.', '+'),
	'contag': ('Contagion', 'aggregation', '0–100', 'How much the classes are aggregated, from the adjacencies between '
	           'them (Li & Reynolds 1993). NA with one class.',
	           'High with few large, compact patches; low with many small, intermixed ones.', '-'),
	'dcad': ('Disjunct core area density', 'core area', 'per 100 ha', 'The separate core areas per 100 ha.', '', ''),
	'division': ('Division index', 'aggregation', '0–1', 'The probability that two cells chosen at random are not in '
	             'the same patch.', '', '+'),
	'ed': ('Edge density', 'area and edge', 'm/ha', 'The edges between classes / the landscape.', '', '+'),
	'ent': ('Shannon entropy', 'complexity', 'bits', 'The entropy of the classes (their composition, in bits).',
	        'The diversity of classes, as SHDI in bits.', ''),
	'iji': ('Interspersion and juxtaposition index', 'aggregation', '0–100', 'How evenly every pair of classes is '
	        'adjacent. NA with fewer than 3 classes.', '', ''),
	'joinent': ('Joint entropy', 'complexity', 'bits', 'The entropy of the pairs of classes of adjacent cells: the '
	            'overall complexity of the landscape (composition and configuration).', '', ''),
	'lpi': ('Largest patch index', 'area and edge', '%', 'The largest patch / the landscape.', '', '-'),
	'lsi': ('Landscape shape index', 'aggregation', '≥ 1', 'All the edges / the minimum edge for the landscape area.',
	        '', '+'),
	'mesh': ('Effective mesh size', 'aggregation', 'ha', 'Σ patch area² / landscape area, over all the patches.', '',
	         '-'),
	'msidi': ('Modified Simpson\'s diversity index', 'diversity', '≥ 0', '−ln Σ P², P = share of each class.', '', ''),
	'msiei': ('Modified Simpson\'s evenness index', 'diversity', '0–1', 'MSIDI / ln(number of classes).', '', ''),
	'mutinf': ('Mutual information', 'complexity', 'bits', 'How much the class of a cell tells about the class of its '
	           'neighbour: ENT − CONDENT (Nowosad & Stepinski 2019).',
	           'High when the classes form large aggregates (spatial autocorrelation).', '-'),
	'ndca': ('Number of disjunct core areas', 'core area', 'number', 'The separate core areas of all the classes.', '',
	         ''),
	'np': ('Number of patches', 'aggregation', 'number', 'The patches of all the classes.', '', '+'),
	'pafrac': ('Perimeter-area fractal dimension', 'shape', '1–2', 'As the class PAFRAC, with every patch. NA with '
	           'fewer than 10 patches.', '', ''),
	'pd': ('Patch density', 'aggregation', 'per 100 ha', 'The patches per 100 ha.', '', '+'),
	'pladj': ('Percentage of like adjacencies', 'aggregation', '%', 'The share of the cell sides that touch the same '
	          'class.', '', '-'),
	'pr': ('Patch richness', 'diversity', 'number', 'The number of classes.', '', ''),
	'prd': ('Patch richness density', 'diversity', 'per 100 ha', 'The number of classes per 100 ha.', '', ''),
	'relmutinf': ('Relative mutual information', 'complexity', '0–1', 'MUTINF / ENT: the mutual information '
	              'independent of the number and share of the classes.', '', '-'),
	'rpr': ('Relative patch richness', 'diversity', '%', 'The number of classes / the maximum possible: in '
	        'landscape_metrics.py and the app, the classes of all the maps of the run.',
	        'Below 100 % in a map that lacks some of the classes of the others.', ''),
	'shdi': ('Shannon\'s diversity index', 'diversity', '≥ 0', '−Σ P ln P, P = share of each class.',
	         '0 with one class; grows with more classes and more even shares.', ''),
	'shei': ('Shannon\'s evenness index', 'diversity', '0–1', 'SHDI / ln(number of classes).',
	         '1 when every class covers the same area.', ''),
	'sidi': ('Simpson\'s diversity index', 'diversity', '0–1', '1 − Σ P²: the probability that two cells chosen at '
	         'random are of different classes.', '', ''),
	'siei': ('Simpson\'s evenness index', 'diversity', '0–1', 'SIDI / (1 − 1 / number of classes).', '', ''),
	'split': ('Splitting index', 'aggregation', '≥ 1', 'Landscape area² / Σ patch area², over all the patches.', '', '+'),
	'ta': ('Total area', 'area and edge', 'ha', 'The area of the landscape: its cells with data.', '', ''),
	'tca': ('Total core area', 'core area', 'ha', 'The core area of all the classes.', '', ''),
	'te': ('Total edge', 'area and edge', 'm', 'The length of the edges between classes.', '', '+'),
}


# What the statistics of the patches (_mn, _sd, _cv) summarize: the value of every patch, with what it is when the
# name does not say it
_OF_PATCHES = {
	'area': ('the patch area', ''),
	'cai': ('the core area index', 'the share of each patch farther from its edge than the edge depth'),
	'circle': ('the related circumscribing circle', '0 for round patches, about 0.36 for squares, near 1 for long, '
	           'narrow ones'),
	'contig': ('the contiguity index', '0 for one cell, near 1 for compact patches'),
	'core': ('the core area', 'the area of each patch farther from its edge than the edge depth'),
	'dcore': ('the number of core areas', 'as landscapemetrics computes the dcore statistics'),
	'enn': ('the distance to the nearest patch of the same class', 'edge to edge, between the centres of their cells; NA '
	        'when a class has one patch'),
	'frac': ('the fractal dimension index', '1 for squares, towards 2 for convoluted edges'),
	'gyrate': ('the radius of gyration', 'the mean distance between the cells of each patch and its centroid'),
	'para': ('the perimeter-area ratio', ''),
	'shape': ('the shape index', '0.25 perimeter / √area: 1 for a square'),
}
# Metrics of the amount of a class, not of its configuration: their change is loss or gain, not fragmentation
AMOUNT = ('ca', 'pland', 'tca', 'cpland', 'area')
# The class metrics that are 0 for a class absent from a map (no cells, no patches); the rest are NA
ZERO_WHEN_ABSENT = ('ca', 'pland', 'np', 'pd', 'lpi', 'te', 'ed', 'tca', 'cpland', 'ndca', 'dcad', 'mesh')


def absent_class(value, metrics):
	"""The row of the class VALUE in a map that does not have it, with the class metrics METRICS (keys c_...)."""
	out = {'class': int(value)}
	for k in metrics:
		if k.startswith('c_'):
			m = k[2:]
			out[m] = 0.0 if m in ZERO_WHEN_ABSENT else None
	return out


def _entry(level, metric):
	stat = metric.rsplit('_', 1)[1] if '_' in metric and metric.rsplit('_', 1)[1] in SUFFIX else None
	base = metric.rsplit('_', 1)[0] if stat else metric
	if stat:
		name, typ, unit, what, read, frag = _BASE[base]
		noun, detail = _OF_PATCHES[base]
		of = 'over the patches of the class' if level == 'class' else 'over all the patches of the landscape'
		what = f'The {SUFFIX[stat]} of {noun}, {of}' + (f' ({detail}).' if detail else '.')
		if stat == 'cv':
			unit, frag, read = '%', '', 'The variability of the patches: 0 when they are all alike.'
		elif stat == 'sd':
			frag, read = '', ''
	else:
		name, typ, unit, what, read, frag = (_BASE if level == 'patch' else _CLASS if level == 'class' else _LANDSCAPE)[
			metric]
	label = metric.upper()
	return {'key': f'{level[0]}_{metric}', 'level': level, 'metric': metric, 'label': label,
	        'name': name + (f' ({SUFFIX[stat]})' if stat else ''), 'type': typ + ' metric', 'unit': unit, 'what': what,
	        'read': read, 'fragmentation': frag, 'amount': metric in AMOUNT, 'function': f'lsm_{level[0]}_{metric}'}


def _catalogue():
	out = [_entry('patch', m) for m in PATCH]
	for level, special in (('class', CLASS), ('landscape', LANDSCAPE)):
		names = list(special) + [f'{b}_{s}' for b in STATS for s in SUFFIX]
		out += [_entry(level, m) for m in sorted(names)]
	return out


CATALOGUE = _catalogue()
KEYS = {m['key']: m for m in CATALOGUE}
DOCS = 'https://r-spatialecology.github.io/landscapemetrics/reference/{function}.html'
REFERENCE = ('Hesselbarth, M.H.K., Sciaini, M., With, K.A., Wiegand, K., Nowosad, J. (2019). landscapemetrics: an '
             'open-source R tool to calculate landscape metrics. Ecography 42: 1648-1657. '
             'https://doi.org/10.1111/ecog.04617')
FRAGSTATS = ('McGarigal, K., Cushman, S.A., Ene, E. (2023). FRAGSTATS v4.3: Spatial Pattern Analysis Program for '
             'Categorical Maps. University of Massachusetts, Amherst. https://www.fragstats.org')


def keys(level=None):
	"""The keys of the metrics (p_area, c_ai, l_shdi...), of one LEVEL ('patch', 'class', 'landscape') or all."""
	return [m['key'] for m in CATALOGUE if level is None or m['level'] == level]


# Helpers --------------------------------------------------------------------------------------------------------------

def _min_edge(n):
	"""Minimum perimeter (cell sides) of N cells (a square, or as close to one as possible)."""
	n = np.asarray(n, dtype=np.float64)
	s = np.floor(np.sqrt(n))
	m = n - s * s
	return np.where(m == 0, 4 * s, np.where(n <= s * (1 + s), 4 * s + 2, 4 * s + 4))


def _max_adjacencies(n):
	"""Maximum like adjacencies (shared sides) of N cells."""
	n = np.asarray(n, dtype=np.float64)
	s = np.floor(np.sqrt(n))
	m = n - s * s
	return np.where(m == 0, 2 * s * (s - 1), np.where(m <= s, 2 * s * (s - 1) + 2 * m - 1, 2 * s * (s - 1) + 2 * m - 2))


def _entropy(counts):
	c = np.asarray(counts, dtype=np.float64).ravel()
	c = c[c > 0]
	if not c.size:
		return 0.0
	p = c / c.sum()
	return float(-(p * np.log2(p)).sum())


def _value(v):
	"""A number of the results: None for NA (NaN or infinite)."""
	if v is None:
		return None
	v = float(v)
	return v if math.isfinite(v) else None


def _stats(values, groups, n_groups):
	"""Mean, standard deviation (n - 1) and coefficient of variation (%) of VALUES (NaN = NA) by GROUPS (0 ..
	n_groups - 1), as R: a group with an NA value gives NA, one value gives an NA standard deviation."""
	values = np.asarray(values, dtype=np.float64)
	cnt = np.bincount(groups, minlength=n_groups).astype(np.float64)
	bad = np.bincount(groups, weights=~np.isfinite(values), minlength=n_groups) > 0
	v = np.where(np.isfinite(values), values, 0.0)
	with np.errstate(invalid='ignore', divide='ignore'):
		mn = np.bincount(groups, weights=v, minlength=n_groups) / cnt
		dev = (v - mn[groups]) ** 2
		sd = np.sqrt(np.bincount(groups, weights=dev, minlength=n_groups) / (cnt - 1))
		cv = sd / mn * 100
	sd[cnt < 2] = np.nan
	for a in (mn, sd, cv):
		a[bad | (cnt == 0)] = np.nan
	return mn, sd, cv


def _shift_pairs(a, diagonal=False):
	"""The pairs of neighbouring cells of A (views): right, down and, with DIAGONAL, down-right and down-left."""
	out = [(a[:, :-1], a[:, 1:]), (a[:-1, :], a[1:, :])]
	if diagonal:
		out += [(a[:-1, :-1], a[1:, 1:]), (a[:-1, 1:], a[1:, :-1])]
	return out


# The landscape --------------------------------------------------------------------------------------------------------

class Landscape:
	"""A map of classes (integers) with its cells without data, and the parts of its metrics, computed when first
	needed. PATCH_CLASSES: the classes whose patches are needed (None = all; fewer is faster for a class-level metric
	of one class, as the similar areas of similarity.py)."""

	def __init__(self, classes, valid, res, directions=8, edge_depth=1, patch_classes=None):
		classes = np.asarray(classes)
		valid = np.asarray(valid, dtype=bool) & (np.isfinite(classes) if classes.dtype.kind == 'f' else True)
		if classes.ndim != 2:
			raise ValueError('the map of classes must have rows and columns')
		if directions not in (4, 8):
			raise ValueError(f'directions must be 4 or 8, not {directions}')
		edge_depth = int(edge_depth)
		if edge_depth < 1:
			raise ValueError('the edge depth must be at least 1 cell')
		self.res, self.directions, self.edge_depth = float(res), int(directions), edge_depth
		self.shape = classes.shape
		values = np.unique(classes[valid])
		if classes.dtype.kind == 'f':
			if values.size and np.any(values != np.round(values)):
				raise ValueError('the classes must be whole numbers')
			values = values.astype(np.int64)
		self.values = values.astype(np.int64)
		self.k = len(values)
		h, w = classes.shape
		# padded with a ring of cells without data, so that the neighbours of every cell exist
		idx = np.full((h + 2, w + 2), -1, dtype=np.int32)
		if self.k:
			idx[1:-1, 1:-1][valid] = np.searchsorted(values, classes[valid].astype(values.dtype))
		self.idx = idx
		self.inside = np.zeros(idx.shape, dtype=bool)
		self.inside[1:-1, 1:-1] = True
		self.n_cells = np.bincount(idx[idx >= 0], minlength=self.k).astype(np.float64)
		self.N = float(self.n_cells.sum())
		self.cell_ha = self.res * self.res / 1e4
		self.area_ha = self.N * self.cell_ha
		if patch_classes is None:
			self.want = np.arange(self.k)
		else:
			self.want = np.flatnonzero(np.isin(self.values, np.asarray(list(patch_classes), dtype=np.int64)))
		self._c = {}

	def _cached(self, name, fn):
		if name not in self._c:
			self._c[name] = fn()
		return self._c[name]

	# adjacencies (4 neighbours) ----------------------------------------------------------------------------------------
	@property
	def adjacency(self):
		"""(C, E_nodata, E_outside): C[i, j] the cell sides of class i next to class j, both ways (the co-occurrence
		matrix of landscapemetrics, its diagonal twice the like adjacencies); E_nodata[i] and E_outside[i] the cell sides
		of class i next to cells without data and to the edge of the map."""
		def make():
			k = self.k
			C = np.zeros((k, k), dtype=np.float64)
			e_nd = np.zeros(k)
			e_out = np.zeros(k)
			for (a, b), (ia, ib) in zip(_shift_pairs(self.idx), _shift_pairs(self.inside)):
				m = (a >= 0) & (b >= 0)
				if k:
					cnt = np.bincount(a[m].astype(np.int64) * k + b[m], minlength=k * k).reshape(k, k)
					C += cnt + cnt.T
				for x, y, iy in ((a, b, ib), (b, a, ia)):  # x has data, y not
					m = (x >= 0) & (y < 0)
					e_nd += np.bincount(x[m & iy], minlength=k)
					e_out += np.bincount(x[m & ~iy], minlength=k)
			return C, e_nd, e_out
		return self._cached('adjacency', make)

	# patches ----------------------------------------------------------------------------------------------------------
	@property
	def patches(self):
		"""(labels, class index of every patch, cells of every patch): the patches of the classes wanted, numbered from
		1 class after class (0 = none)."""
		def make():
			from scipy import ndimage
			st = ndimage.generate_binary_structure(2, 2 if self.directions == 8 else 1)
			lab = np.zeros(self.idx.shape, dtype=np.int32)
			pclass = []
			boxes = ndimage.find_objects(self.idx + 1)  # the box of every class (value + 1)
			n0 = 0
			for ci in self.want:
				box = boxes[ci] if ci < len(boxes) else None
				if box is None:
					continue
				m = self.idx[box] == ci
				lb, n = ndimage.label(m, structure=st, output=np.int32)
				sub = lab[box]
				sub[m] = lb[m] + n0
				pclass += [ci] * n
				n0 += n
			pclass = np.asarray(pclass, dtype=np.int64)
			cells = np.bincount(lab.ravel(), minlength=n0 + 1)[1:].astype(np.float64)
			return lab, pclass, cells
		return self._cached('patches', make)

	@property
	def n_patches(self):
		return len(self.patches[1])

	@property
	def perimeter(self):
		"""Cell sides of every patch next to anything that is not the patch (other classes, cells without data or the
		edge of the map)."""
		def make():
			lab = self.patches[0]
			n = self.n_patches
			p = np.zeros(n + 1)
			for a, b in _shift_pairs(lab):
				m = a != b
				p += np.bincount(a[m], minlength=n + 1) + np.bincount(b[m], minlength=n + 1)
			return p[1:]
		return self._cached('perimeter', make)

	@property
	def edge_cells(self):
		"""The cells with data with a neighbour (4) of another class, without data or outside the map."""
		def make():
			idx = self.idx
			e = np.zeros(idx.shape, dtype=bool)
			core = idx[1:-1, 1:-1]
			e[1:-1, 1:-1] = (core >= 0) & ((idx[:-2, 1:-1] != core) | (idx[2:, 1:-1] != core) |
			                               (idx[1:-1, :-2] != core) | (idx[1:-1, 2:] != core))
			return e
		return self._cached('edge_cells', make)

	@property
	def core_cells(self):
		"""The core cells: farther from the edge of their class (other classes, cells without data, the edge of the
		map) than the edge depth (in cells, 4 neighbours), as landscapemetrics with consider_boundary = FALSE."""
		def make():
			from scipy import ndimage
			if self.edge_depth == 1:
				return (self.idx >= 0) & ~self.edge_cells
			d = ndimage.distance_transform_cdt(np.where(self.edge_cells | (self.idx < 0), 0, 1).astype(np.int8),
			                                   metric='taxicab')
			return (self.idx >= 0) & (d >= self.edge_depth)
		return self._cached('core_cells', make)

	@property
	def core(self):
		def make():
			lab = self.patches[0]
			return np.bincount(lab[self.core_cells], minlength=self.n_patches + 1)[1:].astype(np.float64)
		return self._cached('core', make)

	@property
	def ncore(self):
		"""The separate core areas of every patch (their cells joined with the neighbours of the patches)."""
		def make():
			from scipy import ndimage
			lab = self.patches[0]
			m = self.core_cells & (lab > 0)
			st = ndimage.generate_binary_structure(2, 2 if self.directions == 8 else 1)
			cl, nc = ndimage.label(m, structure=st, output=np.int32)
			owner = np.zeros(nc + 1, dtype=np.int64)
			owner[cl[m]] = lab[m]
			return np.bincount(owner[1:], minlength=self.n_patches + 1)[1:].astype(np.float64)
		return self._cached('ncore', make)

	def _cells_by_rows(self, mask, fn):
		"""FN(rows, columns, labels) on the cells of MASK with a patch, in blocks of rows (to bound the memory)."""
		lab = self.patches[0]
		h, w = lab.shape
		step = max(1, int(4e6 // max(w, 1)))
		for r0 in range(0, h, step):
			sub = lab[r0:r0 + step]
			m = sub > 0 if mask is None else (sub > 0) & mask[r0:r0 + step]
			rr, cc = np.nonzero(m)
			if rr.size:
				fn(rr + r0, cc, sub[rr, cc])

	@property
	def gyrate(self):
		"""The mean distance (m) between the centres of the cells of every patch and its centroid."""
		def make():
			n = self.n_patches
			sx, sy, sd = np.zeros(n + 1), np.zeros(n + 1), np.zeros(n + 1)
			cells = np.concatenate([[1.0], self.patches[2]])

			def sums(rr, cc, ll):
				sx[:] += np.bincount(ll, weights=cc, minlength=n + 1)
				sy[:] += np.bincount(ll, weights=rr, minlength=n + 1)
			self._cells_by_rows(None, sums)
			cx, cy = sx / cells, sy / cells

			def dist(rr, cc, ll):
				sd[:] += np.bincount(ll, weights=np.hypot(cc - cx[ll], rr - cy[ll]), minlength=n + 1)
			self._cells_by_rows(None, dist)
			return sd[1:] / cells[1:] * self.res
		return self._cached('gyrate', make)

	@property
	def contig(self):
		"""The contiguity index of every patch: ((cells + 2 × side neighbours + corner neighbours in the patch, each
		cell with each of its neighbours) / cells - 1) / 12."""
		def make():
			lab = self.patches[0]
			n = self.n_patches
			side, corner = np.zeros(n + 1), np.zeros(n + 1)
			for i, (a, b) in enumerate(_shift_pairs(lab, diagonal=True)):
				m = (a > 0) & (a == b)
				c = 2 * np.bincount(a[m], minlength=n + 1)  # both cells of every pair
				if i < 2:
					side += c
				else:
					corner += c
			cells = self.patches[2]
			return ((corner[1:] + 2 * side[1:] + cells) / cells - 1) / 12
		return self._cached('contig', make)

	@property
	def enn(self):
		"""The distance (m) from every patch to the nearest patch of its class, between the centres of their edge cells
		(NaN for the patches of a class with one patch)."""
		def make():
			lab, pclass, _ = self.patches
			n = self.n_patches
			best = np.full(n + 1, np.inf)
			per_class = np.bincount(pclass, minlength=self.k)
			rows, cols, labs = [], [], []

			def take(rr, cc, ll):
				keep = per_class[pclass[ll - 1]] > 1
				rows.append(rr[keep].astype(np.int32))
				cols.append(cc[keep].astype(np.int32))
				labs.append(ll[keep].astype(np.int32))
			self._cells_by_rows(self.edge_cells, take)
			if rows:
				rr, cc, ll = np.concatenate(rows), np.concatenate(cols), np.concatenate(labs)
				cl = pclass[ll - 1]
				for ci in np.unique(cl):
					sel = cl == ci
					_nearest_other(np.column_stack([rr[sel], cc[sel]]).astype(np.float64), ll[sel], best)
			out = best[1:] * self.res
			out[per_class[pclass] < 2] = np.nan
			return out
		return self._cached('enn', make)

	@property
	def circle(self):
		"""1 - area / area of the smallest circle that contains every corner of the cells of every patch."""
		def make():
			lab, pclass, cells = self.patches
			rows, cols, labs = [], [], []

			def take(rr, cc, ll):
				rows.append(rr.astype(np.int32))
				cols.append(cc.astype(np.int32))
				labs.append(ll.astype(np.int32))
			self._cells_by_rows(self.edge_cells, take)  # the cells at the ends of every row of a patch are edge cells
			diam = np.zeros(self.n_patches)
			if rows:
				diam = _diameters(np.concatenate(rows), np.concatenate(cols), np.concatenate(labs), cells)
			with np.errstate(invalid='ignore', divide='ignore'):
				return 1 - cells / (math.pi * (diam / 2) ** 2)
		return self._cached('circle', make)

	# per patch, in the units of landscapemetrics ----------------------------------------------------------------------
	def patch_values(self, metric):
		cells = self.patches[2]
		a_m2 = cells * self.res * self.res
		p_m = self.perimeter * self.res
		with np.errstate(invalid='ignore', divide='ignore'):
			if metric == 'area':
				return cells * self.cell_ha
			if metric == 'perim':
				return p_m
			if metric == 'core':
				return self.core * self.cell_ha
			if metric == 'cai':
				return self.core / cells * 100
			if metric in ('ncore', 'dcore'):
				return self.ncore
			if metric == 'para':
				return p_m / a_m2
			if metric == 'shape':
				return 0.25 * p_m / np.sqrt(a_m2)
			if metric == 'frac':
				f = 2 * np.log(0.25 * p_m) / np.log(a_m2)
				return np.where(np.isnan(f), 1.0, f)
			if metric == 'gyrate':
				return self.gyrate
			if metric == 'contig':
				return self.contig
			if metric == 'enn':
				return self.enn
			if metric == 'circle':
				return self.circle
		raise ValueError(f'unknown metric of the patches: {metric}')


def _nearest_other(pts, labels, best):
	"""For the points PTS (rows, columns) of the patches LABELS: the distance from every patch to the nearest point of
	another patch, into BEST (indexed by label; only lowered). The k nearest points of each point, k growing until one
	is of another patch or cannot be nearer than the best distance found for its patch; the few points left (of large
	isolated patches): for each of their patches, the exact distance of one of its points (its k nearest points, k
	the points of the patch + 1, include one of another patch) bounds the search to the points of the other patches
	within that distance of the patch."""
	from scipy.spatial import cKDTree
	n = len(pts)
	if n < 2 or len(np.unique(labels)) < 2:
		return
	tree = cKDTree(pts)
	todo = np.arange(n)
	k = 8
	while todo.size:
		kk = min(k, n)
		left = []
		step = max(1, int(4e6 // kk))
		for s in range(0, todo.size, step):
			ch = todo[s:s + step]
			d, j = tree.query(pts[ch], k=kk)
			if kk == 1:
				d, j = d[:, None], j[:, None]
			other = labels[j] != labels[ch][:, None]
			has = other.any(axis=1)
			first = other.argmax(axis=1)
			dist = d[np.arange(len(ch)), first]
			np.minimum.at(best, labels[ch][has], dist[has])
			un = ~has
			left.append((ch[un], d[un, -1]))
		ch = np.concatenate([x[0] for x in left]) if left else np.zeros(0, dtype=np.int64)
		lb = np.concatenate([x[1] for x in left]) if left else np.zeros(0)
		keep = lb < best[labels[ch]]
		todo = ch[keep]
		if kk >= n:
			break
		k *= 4
		if k > 512 and todo.size:
			_nearest_left(pts, labels, best, np.unique(labels[todo]), tree)
			break


def _nearest_left(pts, labels, best, left, tree):
	"""The exact distance from each patch LEFT to the nearest point of another patch, into BEST."""
	from scipy.spatial import cKDTree
	n = len(pts)
	by_label = np.argsort(labels, kind='stable')
	sl = labels[by_label]
	by_row = np.argsort(pts[:, 0], kind='stable')
	rows = pts[by_row, 0]
	for lab_ in left:
		mine = by_label[np.searchsorted(sl, lab_, 'left'):np.searchsorted(sl, lab_, 'right')]
		ub = best[lab_]
		if not np.isfinite(ub):
			d, j = tree.query(pts[mine[0]], k=min(n, len(mine) + 1))
			d, j = np.atleast_1d(d), np.atleast_1d(j)
			other = labels[j] != lab_
			if not other.any():
				continue
			ub = float(d[other.argmax()])
		p = pts[mine]
		r0, r1 = p[:, 0].min() - ub, p[:, 0].max() + ub
		c0, c1 = p[:, 1].min() - ub, p[:, 1].max() + ub
		cand = by_row[np.searchsorted(rows, r0, 'left'):np.searchsorted(rows, r1, 'right')]
		cand = cand[(pts[cand, 1] >= c0) & (pts[cand, 1] <= c1) & (labels[cand] != lab_)]
		if cand.size:
			d, _ = cKDTree(p).query(pts[cand], k=1, distance_upper_bound=ub * (1 + 1e-9) + 1e-9)
			d = d[np.isfinite(d)]
			if d.size:
				ub = min(ub, float(d.min()))
		best[lab_] = min(best[lab_], ub)


def _circumcircle(a, b, c):
	"""Centre and radius of the circle through A, B and C (None if they are on a line)."""
	ax, ay = a
	bx, by = b
	cx, cy = c
	d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
	if abs(d) < 1e-12:
		return None
	ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
	uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
	return (ux, uy), math.hypot(ax - ux, ay - uy)


def _enclosing_diameter(px, py):
	"""The diameter of the smallest circle that contains the points (PX, PY) (Welzl's algorithm, iterative, the points
	in a fixed pseudo-random order so that the result does not change between runs)."""
	pts = list(zip(px.tolist(), py.tolist()))
	if len(pts) == 1:
		return 0.0
	order = np.random.default_rng(0).permutation(len(pts))
	pts = [pts[i] for i in order]
	eps = 1e-9

	def out(p, c, r):
		return math.hypot(p[0] - c[0], p[1] - c[1]) > r + eps

	c, r = pts[0], 0.0
	for i in range(1, len(pts)):
		if not out(pts[i], c, r):
			continue
		c, r = pts[i], 0.0
		for j in range(i):
			if not out(pts[j], c, r):
				continue
			c = ((pts[i][0] + pts[j][0]) / 2, (pts[i][1] + pts[j][1]) / 2)
			r = math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]) / 2
			for k in range(j):
				if not out(pts[k], c, r):
					continue
				cc = _circumcircle(pts[i], pts[j], pts[k])
				if cc is None:  # on a line: the circle of the two farthest of them
					a, b = max(((pts[i], pts[j]), (pts[i], pts[k]), (pts[j], pts[k])),
					           key=lambda q: math.hypot(q[0][0] - q[1][0], q[0][1] - q[1][1]))
					cc = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2), math.hypot(a[0] - b[0], a[1] - b[1]) / 2
				c, r = cc
	return 2 * r


def _diameters(rr, cc, ll, cells):
	"""The diameter (cells) of the smallest circle around the corners of the cells of every patch (CELLS: the cells of
	each patch): around the outer corners of the first and last cell of every row of the patch (its convex hull is made
	of them); for a patch that fills its bounding rectangle, the diagonal of the rectangle."""
	n = len(cells)
	order = np.lexsort((cc, rr, ll))
	rr, cc, ll = rr[order], cc[order], ll[order]
	key_change = np.flatnonzero((np.diff(ll) != 0) | (np.diff(rr) != 0)) + 1
	starts = np.concatenate([[0], key_change])
	ends = np.concatenate([key_change, [len(ll)]]) - 1
	r, lft, rgt, lb = rr[starts], cc[starts], cc[ends], ll[starts]
	# corners (x = column, y = row) of the first and last cell of every row
	x = np.concatenate([lft, lft, rgt + 1, rgt + 1]).astype(np.float64)
	y = np.concatenate([r, r + 1, r, r + 1]).astype(np.float64)
	lbl = np.concatenate([lb, lb, lb, lb])
	o = np.argsort(lbl, kind='stable')
	x, y, lbl = x[o], y[o], lbl[o]
	bounds = np.flatnonzero(np.diff(lbl)) + 1
	starts = np.concatenate([[0], bounds])
	ends = np.concatenate([bounds, [len(lbl)]])
	w = np.maximum.reduceat(x, starts) - np.minimum.reduceat(x, starts)
	h = np.maximum.reduceat(y, starts) - np.minimum.reduceat(y, starts)
	ids = lbl[starts] - 1
	out = np.zeros(n)
	full = cells[ids] == w * h  # rectangles (one cell included): their diagonal
	out[ids[full]] = np.hypot(w[full], h[full])
	for s, e in zip(starts[~full], ends[~full]):
		pts = np.unique(np.column_stack([x[s:e], y[s:e]]), axis=0)
		if len(pts) > 8:
			try:
				from scipy.spatial import ConvexHull
				pts = pts[ConvexHull(pts).vertices]
			except Exception:  # noqa: BLE001 - points on a line: all of them
				pass
		out[lbl[s] - 1] = _enclosing_diameter(pts[:, 0], pts[:, 1])
	return out


# Metrics --------------------------------------------------------------------------------------------------------------

def _class_metric(L, metric, needs):
	"""The values of a class-level metric for every class (array of L.k; NaN = NA)."""
	k = L.k
	nc = L.n_cells
	A = L.area_ha
	ncl = np.bincount(L.patches[1], minlength=k).astype(np.float64) if needs else None
	with np.errstate(invalid='ignore', divide='ignore'):
		if metric == 'ca':
			return nc * L.cell_ha
		if metric == 'pland':
			return nc / L.N * 100
		if metric == 'np':
			return ncl
		if metric == 'pd':
			return ncl / A * 100
		if metric in ('te', 'ed'):
			C = L.adjacency[0]
			te = (C.sum(axis=1) - np.diag(C)) * L.res
			return te if metric == 'te' else te / A
		if metric in ('lsi', 'nlsi'):
			C, e_nd, e_out = L.adjacency
			perim = C.sum(axis=1) - np.diag(C) + e_nd + e_out
			pmin = _min_edge(nc)
			if metric == 'lsi':
				return perim / pmin
			h, w = L.shape
			cb = 2 * h + 2 * w
			pi = nc / L.N
			lim = (0.5 * L.N + 0.5 * cb) / L.N
			even = L.N % 2 == 0
			pmax = np.where(pi <= 0.5, 4 * nc,
			                np.where(pi <= lim, 3 * L.N - 2 * nc + (0 if even else 3), perim + 4 * (L.N - nc)))
			out = (perim - pmin) / (pmax - pmin)
			return np.where(np.isfinite(out), out, np.nan)
		if metric == 'ai':
			C = L.adjacency[0]
			return (np.diag(C) / 2) / _max_adjacencies(nc) * 100
		if metric == 'pladj':
			C = L.adjacency[0]
			return np.diag(C) / (4 * nc) * 100
		if metric == 'clumpy':
			C = L.adjacency[0]
			g = np.diag(C) / (4 * nc - _min_edge(nc))
			p = nc / L.N
			out = np.where(g >= p, (g - p) / (1 - p), np.where(p >= 0.5, (g - p) / (1 - p), (g - p) / -p))
			return np.where(~np.isfinite(g) | (p == 1), np.nan, out)
		if metric == 'iji':
			C = L.adjacency[0].copy()
			if k < 3:
				return np.full(k, np.nan)
			np.fill_diagonal(C, 0)
			r = C / C.sum(axis=1, keepdims=True)
			t = np.where(r > 0, r * np.log(np.where(r > 0, r, 1)), 0.0)
			return -t.sum(axis=1) / math.log(k - 1) * 100
		if metric in ('lpi', 'division', 'mesh', 'split', 'cohesion'):
			pc, cells = L.patches[1], L.patches[2]
			a = cells / L.N  # share of the landscape of every patch
			if metric == 'lpi':
				out = np.full(k, np.nan)
				np.fmax.at(out, pc, a * 100)
				return out
			s2 = np.bincount(pc, weights=a * a, minlength=k)
			if metric == 'division':
				return np.where(ncl > 0, 1 - s2, np.nan)
			if metric == 'mesh':
				return np.where(ncl > 0, s2 * A, np.nan)
			if metric == 'split':
				return np.where(ncl > 0, 1 / s2, np.nan)
			p = L.perimeter
			num = np.bincount(pc, weights=p, minlength=k)
			den = np.bincount(pc, weights=p * np.sqrt(cells), minlength=k)
			return (1 - num / den) / (1 - 1 / math.sqrt(L.N)) * 100
		if metric == 'pafrac':
			pc = L.patches[1]
			a = np.log(L.patch_values('area') * 1e4)
			p = np.log(L.patch_values('perim'))
			out = np.full(k, np.nan)
			for ci in np.flatnonzero(ncl >= 10):
				s = pc == ci
				out[ci] = _pafrac(a[s], p[s])
			return out
		if metric in ('tca', 'cpland'):
			t = np.bincount(L.patches[1], weights=L.core, minlength=k) * L.cell_ha
			return t if metric == 'tca' else t / A * 100
		if metric in ('ndca', 'dcad'):
			t = np.bincount(L.patches[1], weights=L.ncore, minlength=k)
			return t if metric == 'ndca' else t / A * 100
	raise ValueError(f'unknown class metric: {metric}')


def _pafrac(log_area, log_perim):
	"""2 / the slope of ln(area) on ln(perimeter) (least squares), NaN if the perimeters are all equal."""
	x = log_perim - log_perim.mean()
	sxx = float((x * x).sum())
	if sxx <= 0:
		return np.nan
	return 2 / (float((x * (log_area - log_area.mean())).sum()) / sxx)


def _landscape_metric(L, metric, classes_max=None):
	k = L.k
	A = L.area_ha
	P = L.n_cells / L.N if L.N else np.zeros(0)
	with np.errstate(invalid='ignore', divide='ignore'):
		if metric == 'ta':
			return A
		if metric == 'pr':
			return float(k)
		if metric == 'prd':
			return k / A * 100
		if metric == 'rpr':
			return k / classes_max * 100 if classes_max else np.nan
		if metric == 'shdi':
			return float(-(P * np.log(P)).sum())
		if metric == 'shei':
			return 0.0 if k == 1 else float(-(P * np.log(P)).sum() / math.log(k))
		if metric == 'sidi':
			return float(1 - (P * P).sum())
		if metric == 'siei':
			return float(1 - (P * P).sum()) / (1 - 1 / k) if k > 1 else np.nan
		if metric == 'msidi':
			return float(-math.log((P * P).sum()))
		if metric == 'msiei':
			return float(-math.log((P * P).sum())) / math.log(k) if k > 1 else np.nan
		if metric in ('np', 'pd'):
			n = float(L.n_patches)
			return n if metric == 'np' else n / A * 100
		if metric in ('te', 'ed'):
			C = L.adjacency[0]
			te = float(np.tril(C, -1).sum()) * L.res
			return te if metric == 'te' else te / A
		if metric == 'lsi':
			C, e_nd, e_out = L.adjacency
			return float(np.tril(C, -1).sum() + e_nd.sum() + e_out.sum()) / float(_min_edge(L.N))
		if metric == 'pladj':
			C, e_nd, _ = L.adjacency
			return float(np.trace(C)) / float((4 * L.n_cells - e_nd).sum()) * 100
		if metric == 'ai':
			ai = _class_metric(L, 'ai', False)
			return float(np.nansum(ai * P))
		if metric == 'contag':
			if k < 2:
				return np.nan
			C = L.adjacency[0]
			q = C / C.sum()
			q = q[q > 0]
			return (1 + float((q * np.log(q)).sum()) / (2 * math.log(k))) * 100
		if metric == 'iji':
			if k < 3:
				return np.nan
			C = L.adjacency[0]
			e = np.tril(C, -1)
			e = e[np.tril_indices(k, -1)]
			q = e / e.sum()
			q = q[q > 0]
			return -float((q * np.log(q)).sum()) / math.log(0.5 * k * (k - 1)) * 100
		if metric in ('ent', 'condent', 'joinent', 'mutinf', 'relmutinf'):
			C = L.adjacency[0]
			comp = _entropy(C.sum(axis=0))
			cplx = _entropy(C)
			if metric == 'ent':
				return comp
			if metric == 'joinent':
				return cplx
			if metric == 'condent':
				return cplx - comp
			mut = comp - (cplx - comp)
			if metric == 'mutinf':
				return mut
			return 1.0 if mut == 0 else mut / comp
		if metric in ('lpi', 'division', 'mesh', 'split', 'cohesion'):
			cells = L.patches[2]
			a = cells / L.N
			if metric == 'lpi':
				return float(a.max() * 100)
			s2 = float((a * a).sum())
			if metric == 'division':
				return 1 - s2
			if metric == 'mesh':
				return s2 * A
			if metric == 'split':
				return 1 / s2
			p = L.perimeter
			return (1 - p.sum() / (p * np.sqrt(cells)).sum()) / (1 - 1 / math.sqrt(L.N)) * 100
		if metric == 'pafrac':
			if L.n_patches < 10:
				return np.nan
			return _pafrac(np.log(L.patch_values('area')), np.log(L.patch_values('perim')))
		if metric in ('tca', 'ndca', 'dcad'):
			if metric == 'tca':
				return float(L.core.sum()) * L.cell_ha
			n = float(L.ncore.sum())
			return n if metric == 'ndca' else n / A * 100
	raise ValueError(f'unknown landscape metric: {metric}')


def compute(classes, valid, res, metrics=None, directions=8, edge_depth=1, classes_max=None, only_classes=None,
            keep=False):
	"""The landscape metrics of the map CLASSES (rows x columns of whole numbers) on its cells VALID (the rest: no data,
	outside the landscape), with cells of RES metres. METRICS: keys of CATALOGUE (p_area, c_pland, l_shdi...; None =
	all). DIRECTIONS: neighbours of the patches (8 or 4); EDGE_DEPTH: cells; CLASSES_MAX: for RPR (None = NA). ONLY_CLASSES:
	the class-level (and patch-level) metrics of these classes only (faster; no landscape-level metrics then).
	Returns {'patch': [{'class', 'id', 'cells', metric: value...}], 'class': [{'class', metric: value...}],
	'landscape': {metric: value}, 'classes': [...], 'cells': cells with data, 'res': res}, the values None for NA; with
	KEEP, also 'landscape_object': the Landscape (its patches, numbered as the ids, in .patches[0] with a ring of 0)."""
	wanted = keys() if metrics is None else list(dict.fromkeys(metrics))
	unknown = [m for m in wanted if m not in KEYS]
	if unknown:
		raise ValueError('unknown metrics: ' + ', '.join(unknown))
	if only_classes is not None:
		wanted = [m for m in wanted if not m.startswith('l_')]
	L = Landscape(classes, valid, res, directions, edge_depth, patch_classes=only_classes)
	out = {'classes': [int(v) for v in L.values], 'cells': int(L.N), 'res': L.res, 'directions': L.directions,
	       'edge_depth': L.edge_depth, 'patch': [], 'class': [], 'landscape': {}}
	if keep:
		out['landscape_object'] = L
	if not L.k:
		return out
	want_cls = L.want
	p_metrics = [KEYS[m]['metric'] for m in wanted if m.startswith('p_')]
	c_metrics = [KEYS[m]['metric'] for m in wanted if m.startswith('c_')]
	l_metrics = [KEYS[m]['metric'] for m in wanted if m.startswith('l_')]
	needs_patches = bool(p_metrics) or any(m not in ('ca', 'pland', 'te', 'ed', 'lsi', 'nlsi', 'ai', 'pladj', 'clumpy',
	                                                  'iji') for m in c_metrics) or bool(l_metrics)
	pvals = {}
	if needs_patches:
		lab, pclass, cells = L.patches
		for m in set(p_metrics) | {x.rsplit('_', 1)[0] for x in c_metrics + l_metrics if x.rsplit('_', 1)[-1] in SUFFIX}:
			pvals[m] = L.patch_values(m)
		for i in range(L.n_patches if p_metrics else 0):
			row = {'class': int(L.values[pclass[i]]), 'id': i + 1, 'cells': int(cells[i])}
			row.update({m: _value(pvals[m][i]) for m in p_metrics})
			out['patch'].append(row)
	cls_rows = {ci: {'class': int(L.values[ci])} for ci in want_cls}
	for m in c_metrics:
		base, _, stat = m.rpartition('_')
		if stat in SUFFIX and base in STATS:
			mn, sd, cv = _stats(pvals[base], L.patches[1], L.k)
			v = {'mn': mn, 'sd': sd, 'cv': cv}[stat]
		else:
			v = _class_metric(L, m, needs_patches)
		for ci in want_cls:
			cls_rows[ci][m] = _value(v[ci])
	out['class'] = [cls_rows[ci] for ci in want_cls]
	for m in l_metrics:
		base, _, stat = m.rpartition('_')
		if stat in SUFFIX and base in STATS:
			v = pvals[base]
			if v.size == 0 or not np.all(np.isfinite(v)):
				out['landscape'][m] = None
				continue
			mn = float(v.mean())
			sd = float(v.std(ddof=1)) if v.size > 1 else np.nan
			out['landscape'][m] = _value({'mn': mn, 'sd': sd, 'cv': sd / mn * 100 if mn else np.nan}[stat])
		else:
			out['landscape'][m] = _value(_landscape_metric(L, m, classes_max))
	return out
