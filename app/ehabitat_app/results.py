"""Results of the HFT model in segm_optim/results/<study area>/<scenario>: the plots and tables of the scenario
(threshold, R2, tree of the cuts, correlation of the variables) and, for each cut of the tree (p10 ... p90 with
SIMIL = 3, or the scenario folder itself), the plots of the HFTs and their shapefile, shown as an interactive map."""

import csv
import os
import re

from . import inputs, paths

CUT_FILES = [
	# kind, label, regular expression (after <id>_<k><tag>_)
	('map', 'HFT map', r'map_{id}_(?P<k>\d+){tag}_segms_hclust\.png'),
	('radar', 'HFT profiles (radar plot)', r'radarplot_{id}_(?P<k>\d+){tag}_segms_mean\.png'),
	('boxplot', 'Pixel values by HFT', r'boxplot_{id}_(?P<k>\d+){tag}_hft_pixels\.png'),
	('heatmap', 'Heatmap of the HFT profiles', r'heatmap_{id}_(?P<k>\d+){tag}_hft_profiles\.png'),
	('biplot', 'PCA biplot of the HFTs', r'pca_biplot_{id}_(?P<k>\d+){tag}_hft\.png'),
	('hclust', 'Classification tree of the segments', r'hclust_{id}_(?P<k>\d+){tag}_segms_mean\.png'),
]
SCENARIO_FILES = [
	('threshold', "Moran's I and variance by threshold", r'{id}_moran_var_mean\.png'),
	('r2', 'R² of the HFTs and cuts of the tree', r'r2_{id}_hfts\.png'),
	('tree', 'Tree of the HFTs of the cuts', r'hft_tree_{id}_\d+\.png'),
	('cor', 'Correlation between the variables', r'cor_{id}_vars\.png'),
]
CUT_RE = re.compile(r'^p\d+$')
SHP_RE = re.compile(r'^park_segm_(?P<id>-?\d+)_(?P<k>\d+)(?:_(?P<cut>p\d+))?(?P<optim>_optim)?_class\.shp$')
# Set1 of RColorBrewer, the colours of moranvar_plots.R when it does not use the colours of the tree
SET1 = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF']


def list_results():
	"""{study area: [scenarios]} with results (newest first)."""
	out = []
	if not os.path.isdir(paths.RESULTS):
		return out
	for sa in sorted(os.listdir(paths.RESULTS)):
		d = os.path.join(paths.RESULTS, sa)
		if not os.path.isdir(d) or sa.startswith('.'):
			continue
		scen = []
		for sc in os.listdir(d):
			p = os.path.join(d, sc)
			if os.path.isdir(p) and not sc.startswith('.') and _has_results(p):
				scen.append((os.path.getmtime(p), sc))
		if scen:
			out.append({'study_area': sa, 'scenarios': [s for _, s in sorted(scen, reverse=True)]})
	return out


def _has_results(folder):
	for root, _, files in os.walk(folder):
		if any(f.endswith(('.png', '.shp')) for f in files):
			return True
	return False


def _url(path):
	"""Address of a file of the results in the server of the app (with its time, so that a new run is not cached)."""
	import urllib.parse
	return 'files/' + urllib.parse.quote(paths.rel(path)) + f'?v={int(os.path.getmtime(path))}'


def _table(path, limit=200):
	"""A table written by moranvar_plots.R (separated by spaces, with a header)."""
	with open(path, newline='', encoding='utf-8', errors='replace') as f:
		rows = [r for r in csv.reader(f, delimiter=' ') if r]
	return {'file': os.path.basename(path), 'header': rows[0] if rows else [], 'rows': rows[1:limit + 1],
	        'url': _url(path)}


def scenario(study_area, scen):
	root = os.path.join(paths.RESULTS, study_area, scen)
	if not paths.inside(root, paths.RESULTS) or not os.path.isdir(root):
		raise ValueError(f'no results in {study_area}/{scen}')
	files = os.listdir(root)
	cut_dirs = sorted((d for d in files if CUT_RE.match(d) and os.path.isdir(os.path.join(root, d))),
	                  key=lambda d: int(d[1:]))
	ids = set()
	for f in files:
		for pat in (r'^(-?\d+)_moran_var_mean\.png$', r'^(-?\d+)_optim_thresholds\.csv$', r'^r2_(-?\d+)_hfts\.(png|csv)$'):
			m = re.match(pat, f)
			if m:
				ids.add(m.group(1))
	shps = []
	for d in ([os.path.join(root, c) for c in cut_dirs] or [root]):
		for f in sorted(os.listdir(d)):
			m = SHP_RE.match(f)
			if m:
				ids.add(m.group('id'))
				shps.append((d, f, m))
	out = {'study_area': study_area, 'scenario': scen, 'path': root, 'rel': paths.rel(root), 'ids': sorted(ids, key=int),
	       'per_id': {}, 'labels': area_labels(study_area)}
	meta = os.path.join(root, 'processing_metadata.txt')
	out['metadata'] = open(meta, encoding='utf-8', errors='replace').read() if os.path.isfile(meta) else None
	overall = os.path.join(root, 'overall_optim_thresholds.csv')
	thresholds = {}
	if os.path.isfile(overall):
		for line in open(overall):
			p = line.split()
			if len(p) >= 2:
				thresholds[p[0]] = p[1]
	for i in out['ids']:
		rec = {'id': i, 'threshold': thresholds.get(i), 'plots': [], 'tables': [], 'cuts': []}
		for kind, label, pat in SCENARIO_FILES:
			rx = re.compile('^' + pat.format(id=re.escape(i)) + '$')
			for f in sorted(files):
				if rx.match(f):
					rec['plots'].append({'kind': kind, 'label': label, 'file': f, 'url': _url(os.path.join(root, f))})
		for name, label in ((f'r2_{i}_cuts.csv', 'R² of each cut'), (f'{i}_optim_thresholds.csv', 'Segmentation thresholds'),
		                    (f'r2_{i}_hfts.csv', 'R² by number of HFTs')):
			p = os.path.join(root, name)
			if os.path.isfile(p):
				rec['tables'].append(dict(_table(p), label=label))
		# the HFTs of every cut with the mean of their variables (hft_stats_<id>_<k>.csv, with more statistics,
		# is only in the folder)
		for f in sorted(files):
			if re.match(rf'^hft_means_{re.escape(i)}_\d+\.csv$', f):
				rec['tables'].append(dict(_table(os.path.join(root, f)),
				                          label='HFTs of each cut: area and mean of the variables (to compare study areas)'))
		for d, shp, m in shps:
			if m.group('id') != i:
				continue
			cut = m.group('cut') or 'results'
			tag = (f"_{m.group('cut')}" if m.group('cut') else '') + (m.group('optim') or '')
			plots = []
			dfiles = sorted(os.listdir(d))
			for kind, label, pat in CUT_FILES:
				rx = re.compile('^' + pat.format(id=re.escape(i), tag=re.escape(tag)) + '$')
				for f in dfiles:
					if rx.match(f):
						plots.append({'kind': kind, 'label': label, 'file': f, 'url': _url(os.path.join(d, f))})
			stem = shp[:-len('_class.shp')].replace('park_segm_', '', 1)
			rec['cuts'].append({'cut': cut, 'suggested': bool(m.group('optim')), 'k': m.group('k'),
			                    'shp': os.path.join(d, shp), 'plots': plots,
			                    'colours': _first([os.path.join(d, f'hft_colours_{stem}.csv')])})
		# the cuts from more HFTs to fewer, and the counts of HFTs of each cut from r2_<id>_cuts.csv
		counts = {}
		for t in rec['tables']:
			if t['file'] == f'r2_{i}_cuts.csv' and 'cut' in t['header'] and 'hfts' in t['header']:
				ci, hi = t['header'].index('cut'), t['header'].index('hfts')
				counts = {r[ci]: r[hi] for r in t['rows'] if len(r) > max(ci, hi)}
		for c in rec['cuts']:
			c['hfts'] = counts.get(c['cut'])
		rec['cuts'].sort(key=lambda c: int(c['cut'][1:]) if CUT_RE.match(c['cut']) else 0)
		out['per_id'][i] = rec
	return out


def _first(cands):
	for c in cands:
		if os.path.isfile(c):
			return c
	return None


def _hex(c):
	c = str(c).strip()
	return c if c.startswith('#') else '#' + c


def _interp(cols, n):
	"""N colours interpolated along COLS (as colorRampPalette of R)."""
	if n <= len(cols):
		return cols[:n]
	rgb = [tuple(int(c[i:i + 2], 16) for i in (1, 3, 5)) for c in cols]
	out = []
	for j in range(n):
		t = j / (n - 1) * (len(rgb) - 1)
		a = int(t)
		b = min(a + 1, len(rgb) - 1)
		f = t - a
		out.append('#%02X%02X%02X' % tuple(round(rgb[a][k] + (rgb[b][k] - rgb[a][k]) * f) for k in range(3)))
	return out


def _variable_names(meta, n):
	"""Names of the variables (ENV_VARS of processing_metadata.txt), if there are N."""
	if not meta:
		return None
	m = re.search(r'^ENV_VARS = \{(.*?)^\}', meta, re.M | re.S)
	if not m:
		return None
	names = re.findall(r"^\s*,?'([^']+)':", m.group(1), re.M)
	return names if len(names) == n else None


def hft_layer(shp, colours=None, metadata=None):
	"""The HFTs of a shapefile of the results for the map: GeoJSON in WGS 84, colours (those of the plots, from
	hft_colours_*.csv of moranvar_plots.R, or Set1 as the plots without the colours of the tree), legend with the area
	of each HFT, and the profile of each HFT (its mean of each variable, rescaled 0-1 among the HFTs)."""
	import geopandas as gpd
	if not paths.inside(shp, paths.RESULTS):
		raise ValueError('file outside the results')
	g = gpd.read_file(shp)
	if 'hclst_m' not in g.columns:
		raise ValueError(f'{os.path.basename(shp)} has no column of the HFTs (hclst_m)')
	g['hft'] = g['hclst_m'].astype(int)
	area = g.geometry.area if g.crs is not None and g.crs.is_projected else None
	if area is not None:
		g['ha'] = (area * (g.crs.axis_info[0].unit_conversion_factor ** 2) / 1e4).round(2)
	hfts = sorted(g['hft'].unique().tolist())
	order, cmap = None, {}
	if colours and os.path.isfile(colours):
		with open(colours, newline='') as f:
			rows = [r for r in csv.reader(f, delimiter=' ') if r]
		if rows and rows[0][:2] == ['HFT', 'colour']:
			order = [int(r[0]) for r in rows[1:]]
			cmap = {int(r[0]): _hex(r[1]) for r in rows[1:]}
	if not cmap:
		pal = _interp(SET1, max(hfts))
		cmap = {h: pal[h - 1] for h in hfts}
	order = [h for h in (order or []) if h in hfts]
	order += [h for h in hfts if h not in order]
	# the variables: the columns between cat_pa and wdpaid (as moranvar_plots.R writes them)
	cols = list(g.columns)
	try:
		var_cols = cols[cols.index('cat_pa') + 1:cols.index('wdpaid')]
	except ValueError:
		var_cols = [c for c in cols if c not in ('hclst_m', 'segm_id', 'cat', 'label', 'cat_pa', 'wdpaid', 'geometry',
		                                         'hft', 'ha')]
	names = _variable_names(metadata, len(var_cols)) or var_cols
	total = float(g['ha'].sum()) if 'ha' in g else None
	legend = []
	profiles = {}
	for h in order:
		sub = g[g['hft'] == h]
		ha = float(sub['ha'].sum()) if 'ha' in g else None
		legend.append({'hft': h, 'colour': cmap.get(h, '#888888'), 'segments': int(len(sub)),
		               'ha': round(ha, 1) if ha is not None else None,
		               'pct': round(100 * ha / total, 1) if ha is not None and total else None})
		if var_cols:
			profiles[h] = [round(float(sub[c].iloc[0]), 3) if len(sub) else None for c in var_cols]
	w = g.to_crs(4326)
	keep = ['hft', 'segm_id'] + (['ha'] if 'ha' in w else [])
	w = w[keep + ['geometry']]
	gj = w.__geo_interface__
	for f in gj['features']:
		f.pop('bbox', None)
		f['properties'] = {k: (int(v) if k in ('hft', 'segm_id') else v) for k, v in f['properties'].items()}
	_round(gj)
	minx, miny, maxx, maxy = w.total_bounds
	return {'geojson': gj, 'legend': legend, 'bounds': [[float(miny), float(minx)], [float(maxy), float(maxx)]],
	        'variables': names, 'profiles': {str(k): v for k, v in profiles.items()},
	        'colours_from': 'tree' if colours and os.path.isfile(colours) else 'set1', 'file': paths.rel(shp)}


def _round(gj, nd=6):
	def r(c):
		if isinstance(c, (list, tuple)):
			if c and isinstance(c[0], (int, float)):
				return [round(x, nd) for x in c]
			return [r(x) for x in c]
		return c
	for f in gj['features']:
		geom = f.get('geometry')
		if geom and 'coordinates' in geom:
			geom['coordinates'] = r(geom['coordinates'])


def _study_file(study_area):
	"""(file, layer, ID column, name column) of the study areas of the results of STUDY_AREA: those chosen in the app
	if they are that file, else the file of segm_optim/pa with its name (and the ID column of its last setup.py)."""
	settings = paths.load_settings().get('study_area') or {}
	if settings.get('file') and inputs.grass_name(settings['file']) == study_area and os.path.isfile(settings['file']):
		f, layer, col = settings['file'], settings.get('layer'), settings.get('col_id') or 'cat'
		label = settings['label_col'] if 'label_col' in settings else None
	else:
		if not os.path.isdir(paths.PA_DIR):
			return None
		key = lambda s: re.sub(r'[^A-Za-z0-9]', '', s)
		target = {key(study_area), key(re.sub(r'^v_', '', study_area))}
		found = [os.path.join(paths.PA_DIR, x) for x in sorted(os.listdir(paths.PA_DIR))
		         if x.lower().endswith(inputs.VECTOR_EXT) and key(os.path.splitext(x)[0]) in target]
		if not found:
			return None
		import hft_state
		f, layer = found[0], None
		col = (hft_state.get_setup(study_area) or {}).get('col_id') or 'cat'
		label = 'guess'
	if label == 'guess':
		try:
			label = inputs.guess_label_col(inputs.vector_info(f, layer)['fields'])
		except Exception:
			label = None
	return f, layer, col, label


def area_labels(study_area):
	"""{ID (text): name} of the study areas of STUDY_AREA (empty without a column of names)."""
	sf = _study_file(study_area)
	if not sf or not sf[3]:
		return {}
	try:
		return {str(a['id']): a['label'] for a in inputs.vector_areas(*sf) if a['label']}
	except Exception:
		return {}


def study_outline(study_area):
	"""The polygons of the study areas (the file of segm_optim/pa with its name) for the map, with the ID of each."""
	sf = _study_file(study_area)
	if not sf:
		return None
	f, layer, col, label = sf
	try:
		return inputs.vector_geojson(f, layer, col_id=col, label_col=label)
	except Exception:
		return None
