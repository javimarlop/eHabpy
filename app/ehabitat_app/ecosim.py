"""Similarity of the HFTs in their ecoregion (segm_optim/similarity.py): its options in config.py (ECOSIM_*), the
areas the HFTs can be compared with (the ecoregion of their set of variables, or another folder of ecoregion variables
such as a comparison area downloaded from MonitorEO), the plan of a run (the pairs of variables of both scales, the
grid of the similarity and the pixels of each HFT on it), the run itself (a job) and the runs already made, with their
tables, plots and maps. A run is the folder similarity_<pixel size>m (similarity_<pixel size>m_<tag of the area> with
another area than the default one) of a cut of the results (or of the scenario folder with SIMIL 1 or 2), finished
when it has similarity_run.json."""

import csv
import inspect
import importlib
import json
import os
import re
import sys
import threading
import time

from . import inputs, jobs, model, paths, results

import hft_state  # segm_optim (paths puts it in sys.path)

RUN_JSON = 'similarity_run.json'
# _ for the decimal point (12_5m); _<tag>: compared with another area than the ecoregion of the set (similarity.py)
RUN_RE = re.compile(r'^similarity_(?P<res>\d+(?:[._]\d+)?)m(?:_(?P<tag>[A-Za-z0-9_-]+))?$')
CUT_VALUE_RE = re.compile(r'^(p\d+|results)?$')
# the warning of similarity.py when it does not use the study_area folder of the area
STUDY_AT_AREA_RE = re.compile(r'(?:^|[/\\])study_area is not used: (.+)$', re.S)
ENGINE = os.path.join(paths.SEGM, 'similarity.py')
PLOTS = [
	# regular expression, label
	(r'^sim_maps_(?P<id>-?\d+)\.png$', 'Similarity of each HFT in the ecoregion'),
	(r'^sim_metrics_(?P<id>-?\d+)\.png$', 'Metrics of the similar areas of each HFT'),
]
# Colours of the map of the similar areas (similar_<id>_hft<h>.tif): 0 not similar and 255 NoData are transparent
SIMILAR = [(1, '#f4a261', 'Similar, in a patch smaller than the HFT'),
           (2, '#d62828', 'Similar, in a patch at least as large as the HFT'),
           (3, '#1d3557', 'The HFT')]
ZERO = 'd7d7d7'  # similarity 0 (p <= 0.001) in the map of the similarity

# The options of config.py of similarity.py, with their default values (those of config.py in the repository). They
# are not in model.PARAMS, so that the page of the HFT model does not show them. type: as in model.PARAMS, and cut (a
# cut of the tree) and file (a vector file, '' = none or the default one).
PARAMS = [
	{'name': 'ECOSIM_CUT', 'type': 'cut', 'default': '', 'label': 'Cut of the tree',
	 'help': 'The cut whose HFTs are compared with their ecoregion: empty = the suggested cut of each study area (_optim '
	         'in its files), or p10 … p90. With SIMIL 1 or 2 the only cut is the folder of the scenario ("results").'},
	{'name': 'ECOSIM_AREA', 'type': 'area', 'default': '', 'label': 'Ecoregion or comparison area',
	 'help': 'The area whose pixels are compared with the HFTs: empty = the ecoregion grid of the set of variables of the '
	         'scenario (inVars/<set>/ecoregion), or another folder of ecoregion variables of inVars with the same variables '
	         '(e.g. a comparison area downloaded from MonitorEO, inside or outside the study area). Its pixel size is that '
	         'of its variables. The statistics of the HFTs come from its study_area grid (the study area at the pixel of '
	         'the area) when it has one, else from the area grid itself.'},
	{'name': 'ECOSIM_RESOLUTION', 'type': 'number', 'default': 0, 'min': 0, 'step': 1, 'unit': 'm', 'label': 'Pixel size',
	 'help': 'Pixel size of the similarity. 0 = the pixel of the ecoregion variables; otherwise a whole multiple of it '
	         '(e.g. 200, 500 or 1000 m with 100 m variables): each coarse pixel is the mean of its pixels with data, and '
	         'needs data in at least half of its area. The pixels of the HFTs are taken on this same grid.'},
	{'name': 'ECOSIM_STAT', 'type': 'select', 'default': 'mean', 'label': 'Similarity threshold',
	 'choices': [['mean', 'Mean similarity of the HFT pixels (AveHFTSim)'],
	             ['median', 'Median similarity of the HFT pixels (MedianHFTSim)']],
	 'help': 'The pixels of the ecoregion with a similarity at least this statistic of the pixels of the HFT (those with '
	         'a similarity above 0) are the similar areas of the HFT.'},
	{'name': 'ECOSIM_NEIGHBOURS', 'type': 'select', 'cast': 'int', 'default': 8, 'label': 'Neighbourhood of the patches',
	 'choices': [[8, '8 neighbours (queen, as eHabitat+)'], [4, '4 neighbours (rook)']],
	 'help': 'Similar pixels in the same patch: those that touch at a side or a corner (8), or only at a side (4). The '
	         'patches give HSR, the number of patches and the aggregation.'},
	{'name': 'ECOSIM_MIN_COVER', 'type': 'number', 'default': 0.5, 'min': 0.01, 'max': 1, 'step': 0.05,
	 'label': 'Minimum cover of an HFT pixel',
	 'help': 'Fraction of a pixel of the similarity that an HFT has to cover for the pixel to belong to the HFT (each '
	         'pixel goes to the HFT that covers most of it).'},
	{'name': 'ECOSIM_ECOREGIONS', 'type': 'file', 'default': '', 'label': 'Ecoregions',
	 'help': 'Polygons of the ecoregions: empty = ecoregion/ecoregions.gpkg of the set of variables (MonitorEO), or a '
	         'file of your own (ID column ECO_ID, or else its first integer column). Each HFT is compared with the rest of '
	         'its own ecoregion.'},
	{'name': 'ECOSIM_MARINE', 'type': 'file', 'default': '', 'label': 'Marine ecoregions',
	 'help': 'Optional polygons of marine ecoregions (e.g. MEOW), only to count them in MarEco of the THDI; empty = none '
	         '(MarEco = 0).'},
]
BY_NAME = {p['name']: p for p in PARAMS}


# Options --------------------------------------------------------------------------------------------------------------

def defaults():
	return {p['name']: p['default'] for p in PARAMS}


def options():
	"""The options with their values in config.py and their defaults."""
	cfg = model.read_config()
	out = []
	for p in PARAMS:
		q = dict(p)
		q['value'] = cfg.get(p['name'], p['default'])
		out.append(q)
	return {'params': out, 'defaults': defaults()}


def area_dir(value):
	"""Absolute folder of a value of ECOSIM_AREA (relative to inVars, or absolute); None for ''."""
	value = os.path.expanduser(str(value or '').strip())
	if not value:
		return None
	return os.path.normpath(value if os.path.isabs(value) else os.path.join(paths.INVARS, value))


def _literal(p, value):
	"""Python literal of VALUE for config.py, checked."""
	if p['type'] == 'area':
		value = str(value or '').strip().replace('\\', '/').rstrip('/')
		if value and not os.path.isdir(area_dir(value)):
			raise ValueError(f"{p['label']}: the folder {value} does not exist (a folder of inVars with ecoregion "
			                 'variables, e.g. gee_guajares_snevada/ecoregion)')
		return repr(value)
	if p['type'] == 'cut':
		value = str(value or '').strip()
		if not CUT_VALUE_RE.match(value):
			raise ValueError(f"{p['label']}: {value!r} is not a cut of the tree (p10 … p90, results or empty)")
		return repr(value)
	if p['type'] == 'file':
		value = str(value or '').strip()
		if value:
			value = os.path.abspath(os.path.expanduser(value))
			if not os.path.isfile(value):
				raise ValueError(f"{p['label']}: {value} not found")
			if not value.lower().endswith(inputs.VECTOR_EXT):
				raise ValueError(f"{p['label']}: {os.path.basename(value)} is not a vector file (.gpkg, .shp, .geojson)")
		return repr(value)
	return model._literal(p, value)


def save(values):
	"""Write VALUES ({name: value} of the ECOSIM_* options) to config.py, keeping its comments; returns the names
	changed."""
	lits = {name: _literal(BY_NAME[name], value) for name, value in (values or {}).items() if name in BY_NAME}
	cfg = model.read_config()
	changed = []
	for name, lit in lits.items():
		if name not in cfg or repr(cfg.get(name)) != lit:
			hft_state.write_option(name, lit)
			changed.append(name)
	return changed


# Plan -----------------------------------------------------------------------------------------------------------------

_engine_lock = threading.Lock()
_engine_mtime = {'t': None}


def _engine():
	"""The module similarity of segm_optim (imported again when the file changed, e.g. after an update of eHabpy), or a
	ValueError that says why it cannot be used."""
	with _engine_lock:
		try:
			mtime = os.path.getmtime(ENGINE)
		except OSError:
			mtime = None
		try:
			import similarity
			if mtime is not None and _engine_mtime['t'] is not None and mtime != _engine_mtime['t']:
				similarity = importlib.reload(similarity)
			_engine_mtime['t'] = mtime
		except ImportError as e:
			if getattr(e, 'name', None) == 'similarity':
				raise ValueError('this copy of eHabpy has no segm_optim/similarity.py: update it (Check for updates)')
			raise ValueError(f'segm_optim/similarity.py cannot be loaded: {e}')
		except SyntaxError as e:
			raise ValueError(f'segm_optim/similarity.py cannot be loaded: {e}')
	return similarity


def _engine_has_area():
	"""Whether similarity.py can compare the HFTs with another area than the ecoregion of their set (--area)."""
	try:
		with open(ENGINE, encoding='utf-8') as f:
			return "'--area'" in f.read()
	except OSError:
		return False


def _plain(x):
	"""X with the numbers of numpy as numbers of Python (for JSON)."""
	if isinstance(x, dict):
		return {str(k): _plain(v) for k, v in x.items()}
	if isinstance(x, (list, tuple)):
		return [_plain(v) for v in x]
	if hasattr(x, 'item') and callable(x.item) and not isinstance(x, (str, bytes)):
		try:
			return x.item()
		except (TypeError, ValueError):
			return str(x)
	return x


def _abs(path):
	"""Absolute path of a file named as similarity.py names them: relative to the copy of eHabpy (or to segm_optim),
	or absolute."""
	if not path:
		return None
	path = os.path.expanduser(str(path))
	if os.path.isabs(path):
		return os.path.normpath(path)
	for base in (paths.ROOT, paths.SEGM):
		p = os.path.normpath(os.path.join(base, path))
		if os.path.exists(p):
			return p
	return os.path.normpath(os.path.join(paths.ROOT, path))


def _num(v):
	if v is None or v == '':
		return None
	try:
		return float(v)
	except (TypeError, ValueError):
		raise ValueError(f'{v!r} is not a number')


def _res_name(res):
	"""Name of the folder of a run with pixels of RES metres (as similarity.py)."""
	res = float(res)
	return f'similarity_{res:.0f}m' if abs(res - round(res)) < 1e-6 else f"similarity_{res:g}m".replace('.', '_')


def colours(shp, hfts, colours_csv=None):
	"""{HFT: colour} of the HFTs of a shapefile of the results: those of the plots (hft_colours_<stem>.csv of
	moranvar_plots.R) or, without it, Set1 as the plots without the colours of the tree (as results.hft_layer)."""
	if not colours_csv and shp:
		stem = os.path.basename(shp)[:-len('_class.shp')].replace('park_segm_', '', 1)
		colours_csv = os.path.join(os.path.dirname(shp), f'hft_colours_{stem}.csv')
	cmap = {}
	if colours_csv and os.path.isfile(colours_csv):
		with open(colours_csv, newline='') as f:
			rows = [r for r in csv.reader(f, delimiter=' ') if r]
		if rows and rows[0][:2] == ['HFT', 'colour']:
			cmap = {int(r[0]): results._hex(r[1]) for r in rows[1:] if len(r) > 1}
	hfts = [int(h) for h in hfts if h is not None]
	if hfts and any(h not in cmap for h in hfts):
		pal = results._interp(results.SET1, max(hfts))
		for h in hfts:
			cmap.setdefault(h, pal[h - 1] if 0 < h <= len(pal) else '#888888')
	return cmap


def areas(study_area, scenario):
	"""The areas the HFTs of a scenario can be compared with (similarity.comparison_areas): the ecoregion grid of
	their set of variables first (default), then the other folders of ecoregion variables of inVars, each with its
	kind, name, pixel, the variables it lacks and its study_area grid. [] when similarity.py cannot compare with other
	areas."""
	sim = _engine()
	if not hasattr(sim, 'comparison_areas'):
		return []
	out = _plain(sim.comparison_areas(study_area=study_area, scenario=scenario, cfg=model.read_config()) or [])
	out.sort(key=lambda a: not a.get('default'))
	for a in out:
		p = _abs(a.get('polygons'))
		a['polygons_path'] = p if p and os.path.isfile(p) else None
	return out


def _area_record(area):
	"""The area of a plan or a run with the absolute path of its polygons (for the map)."""
	if not isinstance(area, dict):
		return area
	p = _abs(area.get('polygons'))
	return dict(area, path=p if p and os.path.isfile(p) else None)


def _reference_rejected(p):
	"""Why similarity.py does not use the study_area folder of the area of a plan (other CRS, pixel, lattice, variables
	or specs), or None when it uses it or the area has none: area.reference_note of the plan, or, in a similarity.py
	without it, the warning '<area>/study_area is not used: <why>'."""
	area = p.get('area')
	if isinstance(area, dict) and 'reference_note' in area:
		return area['reference_note'] or None
	for w in p.get('warnings') or []:
		m = STUDY_AT_AREA_RE.search(str(w))
		if m:
			return m.group(1)
	return None


def plan(study_area, scenario, cut=None, res=None, min_cover=None, ecoregions=None, area=None):
	"""What similarity.py would do (similarity.plan): the variables of both scales, the grid of the similarity and the
	pixels of each HFT on it, with the colours of the HFTs and the folder of the run of each study area. None = the
	option of config.py; AREA '' = the ecoregion of the set of variables, or a folder of ecoregion variables."""
	sim = _engine()
	if cut is not None and not CUT_VALUE_RE.match(str(cut)):
		raise ValueError(f'{cut!r} is not a cut of the tree')
	kw = {}
	if area is not None:
		area = str(area).strip().replace('\\', '/').rstrip('/')
		if 'area' in inspect.signature(sim.plan).parameters:
			kw['area'] = area
		elif area:
			raise ValueError('this copy of similarity.py cannot compare the HFTs with another area: update eHabpy '
			                 '(Check for updates)')
	p = _plain(sim.plan(study_area=study_area, scenario=scenario, cut=cut, res=_num(res), min_cover=_num(min_cover),
	                    ecoregions=ecoregions, **kw))
	p['area'] = _area_record(p.get('area'))
	targets = p.get('targets') or []
	pixels = sum(int(h.get('ref_pixels') or 0) for t in targets for h in t.get('hfts') or [])
	# no HFT pixel on the grid of their statistics: the area does not cover the study area and has no study_area grid
	# (similarity.py says so in its errors and leaves those targets out; the page explains what to do)
	p['no_hft_pixels'] = bool(p.get('grid')) and pixels == 0 and (any(t.get('hfts') for t in targets)
	                                                               or bool(p.get('errors')))
	p['reference_rejected'] = _reference_rejected(p)
	for t in targets:
		shp = _abs(t.get('shp'))
		hfts = t.get('hfts') or []
		cmap = colours(shp, [h.get('hft') for h in hfts])
		for h in hfts:
			h['colour'] = cmap.get(h.get('hft'), '#888888')
		run_dir = _abs(t['folder']) if t.get('folder') else \
			os.path.join(os.path.dirname(shp), _res_name(p['res'])) if shp and p.get('res') else None
		t['run_folder'] = paths.rel(run_dir) if run_dir else None
		t['run_exists'] = bool(run_dir and os.path.isdir(run_dir) and os.listdir(run_dir))
	eco = p.get('ecoregions')
	if eco and eco.get('file'):
		eco['path'] = _abs(eco['file'])
	return p


def sources(study_area, scenario):
	"""What the similarity of a scenario starts from: its study areas (IDs) and cuts with their HFTs, the cut chosen in
	config.py (if the scenario has it), the options and the plan with them (or the error that stops it)."""
	r = results.scenario(study_area, scenario)
	cuts, suggested = {}, {}
	for i, rec in r['per_id'].items():
		for c in rec['cuts']:
			d = cuts.setdefault(c['cut'], {'cut': c['cut'], 'hfts': {}, 'suggested': []})
			if c.get('hfts'):
				d['hfts'][i] = c['hfts']
			if c['suggested']:
				d['suggested'].append(i)
				suggested[i] = c['cut']
	ordered = sorted(cuts.values(), key=lambda c: int(c['cut'][1:]) if results.CUT_RE.match(c['cut']) else 0)
	names = [c['cut'] for c in ordered]
	cfg = model.read_config()
	cut = str(cfg.get('ECOSIM_CUT', '') or '')
	if cut not in names and not (cut == '' and suggested):
		cut = '' if suggested else (names[0] if names else '')
	out = {'study_area': study_area, 'scenario': scenario, 'rel': r['rel'], 'ids': r['ids'], 'labels': r['labels'],
	       'cuts': ordered, 'suggested': suggested, 'cut': cut, 'options': options(), 'plan': None, 'error': None,
	       'engine': os.path.isfile(ENGINE), 'areas': [], 'areas_error': None,
	       'area': str(cfg.get('ECOSIM_AREA', '') or '')}
	try:
		out['areas'] = areas(study_area, scenario)
	except (ValueError, OSError, RuntimeError) as e:
		out['areas_error'] = str(e)
	if not r['ids'] or not names:
		out['error'] = f'{study_area}/{scenario} has no HFTs (no shapefiles of the HFTs)'
		return out
	try:
		out['plan'] = plan(study_area, scenario, cut)
	except ValueError as e:
		out['error'] = str(e)
	return out


# Run ------------------------------------------------------------------------------------------------------------------

def run(values, study_area, scenario, cut=''):
	"""Save the options and run similarity.py for the HFTs of a scenario and cut (a job); its result is the folder of
	the run made."""
	if jobs.running('model'):
		raise ValueError('the HFT model is running: compute the similarity when it finishes')
	if jobs.running('ecosim'):
		raise ValueError('the similarity is already being computed')
	root = os.path.join(paths.RESULTS, str(study_area or ''), str(scenario or ''))
	if not study_area or not scenario or not paths.inside(root, paths.RESULTS) or not os.path.isdir(root):
		raise ValueError(f'no results in {study_area}/{scenario}')
	if not os.path.isfile(os.path.join(paths.SEGM, 'similarity.py')):
		raise ValueError('this copy of eHabpy has no segm_optim/similarity.py: update it (Check for updates)')
	cut = str(cut or '').strip()
	if not CUT_VALUE_RE.match(cut):
		raise ValueError(f'{cut!r} is not a cut of the tree')
	area = str((values or {}).get('ECOSIM_AREA') or '').strip()
	if area and not _engine_has_area():
		raise ValueError('this copy of similarity.py cannot compare the HFTs with another area: update eHabpy (Check '
		                 'for updates)')
	save(dict(values or {}, ECOSIM_CUT=cut))
	cmd = [sys.executable, 'similarity.py', '--study-area', study_area, '--scenario', scenario, '--cut', cut]
	if _engine_has_area():
		area = str(model.read_config().get('ECOSIM_AREA', '') or '')
		if area:
			cmd += ['--area', area]
	started = time.time()

	def after(job):
		made = [x for x in runs(study_area, scenario) if x['mtime'] >= started - 1]
		if made:
			job.say('\nResults in ' + ', '.join(x['folder'] for x in made))
		else:
			job.say('\nWARNING: no finished run (similarity_run.json) was found')
		return {'study_area': study_area, 'scenario': scenario, 'cut': made[0]['cut'] if made else cut,
		        'folder': made[0]['folder'] if made else None, 'folders': [x['folder'] for x in made]}

	job = jobs.start_commands('ecosim', 'Similarity in the ecoregion', [cmd], after=after)
	return job.id


# Results --------------------------------------------------------------------------------------------------------------

def _read_run(folder):
	with open(os.path.join(folder, RUN_JSON), encoding='utf-8') as f:
		rec = json.load(f)
	if not isinstance(rec, dict):
		raise ValueError(f'{RUN_JSON} is not a record of a run')
	return rec


def runs(study_area, scenario):
	"""The finished runs of a scenario (folders similarity_<res>m with similarity_run.json in its cuts, or in the
	scenario folder with SIMIL 1 or 2), newest first."""
	root = os.path.join(paths.RESULTS, str(study_area or ''), str(scenario or ''))
	if not paths.inside(root, paths.RESULTS) or not os.path.isdir(root):
		return []
	dirs = [(c, os.path.join(root, c)) for c in sorted(os.listdir(root))
	        if results.CUT_RE.match(c) and os.path.isdir(os.path.join(root, c))] + [('results', root)]
	out = []
	for cut, d in dirs:
		for name in sorted(os.listdir(d)):
			m = RUN_RE.match(name)
			folder = os.path.join(d, name)
			if not m or not os.path.isfile(os.path.join(folder, RUN_JSON)):
				continue
			try:
				rec = _read_run(folder)
			except (OSError, ValueError):
				continue
			hfts = rec.get('hfts') or []
			counts = {}
			for h in hfts:
				counts[h.get('status') or 'unknown'] = counts.get(h.get('status') or 'unknown', 0) + 1
			res = rec.get('res') or float(m.group('res').replace('_', '.'))
			area = rec.get('area') if isinstance(rec.get('area'), dict) else {}
			out.append({'cut': cut, 'res': res, 'date': rec.get('date'), 'ids': [str(i) for i in rec.get('ids') or []],
			            'folder': paths.rel(folder), 'hfts': len(hfts),
			            'counts': counts, 'stat': (rec.get('options') or {}).get('stat'), 'seconds': rec.get('seconds'),
			            'tag': area.get('tag') or m.group('tag'), 'area': {k: area.get(k) for k in
			                                                                ('folder', 'tag', 'kind', 'name')} if area else None,
			            'mtime': os.path.getmtime(os.path.join(folder, RUN_JSON))})
	out.sort(key=lambda x: x['mtime'], reverse=True)
	return _plain(out)


def _run_dir(folder):
	d = _abs(folder)
	if not d or not paths.inside(d, paths.RESULTS) or not os.path.isdir(d):
		raise ValueError(f'{folder} is not a folder of the results')
	return d


def _rows(table):
	"""The rows of a table of results._table as {column: value} ('NA' as None)."""
	if not table:
		return []
	h = table['header']
	return [{k: (v if v != 'NA' else None) for k, v in zip(h, r)} for r in table['rows']]


def run_results(folder):
	"""A finished run: its record (similarity_run.json), the tables hri_results.csv and thdi_res.csv, the plots, the
	maps of each HFT with its colour, the ecoregions and similarity_metadata.txt."""
	d = _run_dir(folder)
	if not os.path.isfile(os.path.join(d, RUN_JSON)):
		raise ValueError(f'{paths.rel(d)} has no finished run ({RUN_JSON})')
	rec = _plain(_read_run(d))
	names = rec.get('tables') or {}
	tables = {}
	for key, default, label in (('hri', 'hri_results.csv', 'Similarity and landscape metrics of each HFT'),
	                            ('thdi', 'thdi_res.csv', 'Total Habitat Diversity and Replaceability (THDI)')):
		p = os.path.join(d, names.get(key) or default)
		tables[key] = dict(results._table(p, limit=2000), label=label) if os.path.isfile(p) else None
	# plots: those of the record first, then any other PNG of the folder
	pngs = [f for f in (rec.get('plots') or []) if os.path.isfile(os.path.join(d, f))]
	pngs += sorted(f for f in os.listdir(d) if f.lower().endswith('.png') and f not in pngs)
	plots = []
	for f in pngs:
		label = f
		for pat, lab in PLOTS:
			m = re.match(pat, f)
			if m:
				label = lab + (f" · ID {m.group('id')}" if len(rec.get('ids') or []) > 1 else '')
		plots.append({'label': label, 'file': f, 'url': results._url(os.path.join(d, f))})
	# maps of each HFT, with the colours of the HFT map: the files of the HFTs next to the run first (a scenario moved
	# to <scenario>_previous keeps its own files, while the recorded paths are those of the new scenario)
	def local(path):
		if not path:
			return None
		cand = os.path.join(os.path.dirname(d), os.path.basename(str(path)))
		return cand if os.path.isfile(cand) else _abs(path)

	by_id = {}
	for h in rec.get('hfts') or []:
		by_id.setdefault(str(h.get('id')), []).append(h.get('hft'))
	cmaps = {}
	for i, hfts in by_id.items():
		shp = local((rec.get('shapefiles') or {}).get(i))
		cmaps[i] = colours(shp, hfts, local((rec.get('colours') or {}).get(i)))
	maps = []
	for h in rec.get('hfts') or []:
		i = str(h.get('id'))
		m = {k: h.get(k) for k in ('hft', 'status', 'threshold', 'area_ha', 'ref_pixels', 'vars_used', 'ecoregion')}
		m['id'] = i
		m['colour'] = cmaps.get(i, {}).get(h.get('hft'), '#888888')
		for kind in ('sim', 'similar'):
			p = os.path.join(d, h[kind]) if h.get(kind) else None
			m[kind] = paths.rel(p) if p and os.path.isfile(p) else None
		maps.append(m)
	eco = rec.get('ecoregions') or None
	if eco and eco.get('file'):
		p = _abs(eco['file'])
		eco = dict(eco, path=p if os.path.isfile(p) else None)
	meta = os.path.join(d, 'similarity_metadata.txt')
	return {'folder': paths.rel(d), 'path': d, 'run': rec, 'tables': tables, 'hri': _rows(tables['hri']),
	        'thdi': _rows(tables['thdi']), 'plots': plots, 'maps': maps, 'ecoregions': eco,
	        'area': _area_record(rec.get('area')),
	        'metadata': open(meta, encoding='utf-8', errors='replace').read() if os.path.isfile(meta) else None}


def preview(tif, kind='sim', bbox=None, max_size=1200):
	"""A map of a run for the map of the app (PNG in Web Mercator): KIND 'sim', the similarity from 0 to 1 (viridis, 0
	in light grey), or 'similar', the similar areas and the HFT (categories). BBOX ([west, south, east, north] in WGS
	84): only that part (the detail of the view when the map is zoomed in), else the whole map."""
	import numpy as np
	path = _abs(tif)
	if not path or not paths.inside(path, paths.RESULTS) or not os.path.isfile(path):
		raise ValueError('only the maps of the results')
	if kind not in ('sim', 'similar'):
		raise ValueError(f'unknown kind of map {kind!r}')
	from osgeo import gdal
	if bbox is not None and (len(bbox) != 4 or not all(isinstance(v, (int, float)) for v in bbox)):
		raise ValueError('the rectangle of the view must be [west, south, east, north]')
	coarse = not bbox and gdal.Open(path).RasterXSize > max_size
	# coarser than the raster: the mean similarity, and the highest class so that small patches and HFTs are seen
	a, valid, bbox = inputs.warp_preview(path, max_size, 'average' if kind == 'sim' else 'max', overviews=kind == 'sim',
	                                     bbox=bbox)
	note = None
	if coarse:
		note = 'each pixel of the view is the mean of the pixels it covers' if kind == 'sim' else \
			'each pixel of the view has the highest class of the pixels it covers'
	rgba = np.zeros((4,) + a.shape, dtype='uint8')
	if kind == 'sim':
		v = np.clip(np.where(valid, a, 0), 0, 1)
		lut = inputs._ramp(inputs.VIRIDIS)
		idx = np.rint(v * 255).astype('uint8')
		for i in range(3):
			rgba[i] = lut[idx, i]
		rgba[3] = np.where(valid, 230, 0)
		zero = valid & (v < 0.001)
		for i in range(3):
			rgba[i][zero] = int(ZERO[2 * i:2 * i + 2], 16)
		rgba[3][zero] = 120
		return {'image': inputs.png_uri(rgba), 'bbox': bbox, 'width': int(a.shape[1]), 'kind': kind, 'min': 0, 'max': 1,
		        'ramp': ['#' + c for c in inputs.VIRIDIS], 'zero': '#' + ZERO, 'empty': not valid.any(), 'note': note}
	v = np.rint(np.where(valid, a, 0)).astype('int64')
	for value, colour, _ in SIMILAR:
		m = valid & (v == value)
		for i in range(3):
			rgba[i][m] = int(colour[1 + 2 * i:3 + 2 * i], 16)
		rgba[3][m] = 235
	return {'image': inputs.png_uri(rgba), 'bbox': bbox, 'width': int(a.shape[1]), 'kind': kind, 'empty': not valid.any(),
	        'note': note,
	        'legend': [{'value': value, 'colour': colour, 'label': label} for value, colour, label in SIMILAR]}


def _same(a, b):
	try:
		return os.path.samefile(a, b)
	except OSError:
		return os.path.normpath(a) == os.path.normpath(b)


def outline(path, field=None, ids=None, folder=None):
	"""Polygons of the ecoregions or of a comparison area for the map (GeoJSON in WGS 84), only those with FIELD in IDS
	when given: the files of inVars, segm_optim/pa and the results, those of the options ECOSIM_ECOREGIONS and
	ECOSIM_MARINE and, with FOLDER (a finished run), the files recorded by that run."""
	from osgeo import gdal
	p = _abs(path)
	cfg = model.read_config()
	own = {os.path.abspath(os.path.expanduser(str(cfg.get(n) or ''))) for n in ('ECOSIM_ECOREGIONS', 'ECOSIM_MARINE')
	       if cfg.get(n)}
	ok = bool(p) and os.path.isfile(p) and (any(paths.inside(p, f) for f in (paths.INVARS, paths.PA_DIR, paths.RESULTS))
	                                         or p in own)
	if p and os.path.isfile(p) and not ok and folder:
		try:
			rec = _read_run(_run_dir(folder))
		except (OSError, ValueError):
			rec = {}
		recorded = [(rec.get('ecoregions') or {}).get('file'), rec.get('marine'), (rec.get('area') or {}).get('polygons')]
		ok = any(_same(p, _abs(f)) for f in recorded if f)
	if not ok:
		raise ValueError('only the ecoregions of inVars, segm_optim/pa, the results, the options or the run')
	info = inputs.vector_info(p)
	label = inputs.guess_label_col(info['fields'])
	names = {f['name'] for f in info['fields']}
	ids = [int(i) for i in (ids or []) if str(i).lstrip('-').isdigit()]
	if not (field and field in names and ids):
		return inputs.vector_geojson(p, col_id=field if field in names else None, label_col=label)
	tmp = f'/vsimem/ecoregions_{os.getpid()}_{time.time_ns()}.geojson'
	try:
		gdal.VectorTranslate(tmp, p, options=gdal.VectorTranslateOptions(
			format='GeoJSON', where=f'"{field}" IN ({", ".join(map(str, ids))})'))
		return inputs.vector_geojson(tmp, col_id=field, label_col=label)
	finally:
		gdal.Unlink(tmp)
