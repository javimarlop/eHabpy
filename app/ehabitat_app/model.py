"""Parameters of the HFT model (config.py), the IDs to process (palist.csv) and the steps of the model: setup.py
imports the variables and the study area into GRASS, run_hft.py runs the steps that are out of date (hft_state.py)."""

import math
import os
import re
import sys
import multiprocessing

from . import inputs, jobs, paths

import hft_state  # segm_optim (paths puts it in sys.path)
import var_order  # segm_optim

THRESHOLD_METHODS = ('f', 'gs', 'gs_fixed', 'crossing', 'crossing_first')

# The options of config.py that the app changes, with their default values (those of config.py in the repository)
# and the help of config.py. type: bool, int, number, select, text, threshold (a method or a fixed threshold), folder,
# cuts (a list of percentages) and thresholds (a list of thresholds of the segmentation).
PARAMS = [
	# Scenario and run
	{'name': 'SCENARIO', 'group': 'run', 'type': 'text', 'default': '', 'label': 'Scenario',
	 'help': 'Subfolder of the results, results/<study area>/<scenario> (letters, digits, - and _), to keep apart the '
	         'results of different options. Empty = the name of the folder of variables (INVARS_FOLDER), or "base" for '
	         'inVars itself. If the scenario already has results made with other options, they are moved to '
	         '<scenario>_previous.'},
	{'name': 'FORCE_RESTART', 'group': 'run', 'type': 'bool', 'default': False, 'label': 'Restart the segmentation from scratch',
	 'help': 'Ticked: deletes the control files and csv/shp/tiffs and runs the segmentation again from the start. '
	         'Unticked: resumes an interrupted segmentation, or does nothing if it finished with the same options (if the '
	         'options, palist.csv or the data changed, it starts from scratch).'},
	{'name': 'GRASSDB', 'group': 'run', 'type': 'folder', 'default': '~/grassdata', 'label': 'GRASS database',
	 'help': 'Folder of the GRASS GIS database, where setup.py creates the project of the study area (not the folder '
	         'where GRASS is installed).'},
	# Segmentation
	{'name': 'PA_BUFFER', 'group': 'segm', 'type': 'number', 'default': 500, 'min': 0, 'step': 10, 'unit': 'm',
	 'label': 'Buffer around the study area',
	 'help': 'Distance in map units (metres) added to the rectangle of the study area for the working region (rounded up '
	         'to a whole number of pixels). 0 = the rectangle of the study area.'},
	{'name': 'CLIP_TO_PA', 'group': 'segm', 'type': 'bool', 'default': False, 'label': 'Clip to the study area',
	 'help': 'Ticked: only the pixels inside the polygons of the study area are analysed. Unticked: the whole rectangle '
	         'of the working region (with the buffer).'},
	{'name': 'SEGM_PCA', 'group': 'segm', 'type': 'bool', 'default': True, 'label': 'Segment the first 3 axes of a PCA',
	 'help': 'Ticked: segments the first three axes of a standardized PCA of the variables (needs at least 3). '
	         'Unticked: segments the variables themselves, all with the same weight (highly correlated variables weigh '
	         'more). In both cases i.segment scales each layer to its range, so the units do not matter.'},
	{'name': 'CLIP_PERCENTILE', 'group': 'segm', 'type': 'number', 'default': 0, 'min': 0, 'max': 49.9, 'step': 0.1,
	 'unit': 'percentile', 'label': 'Clipping of extreme values',
	 'help': '0 = no clipping. > 0: before the segmentation and the classification, each layer is clipped to its '
	         'percentiles P and 100 - P (e.g. 0.5), so that a variable with a long tail (e.g. TWI) does not end up with '
	         'almost no weight. The plots and statistics of the HFTs always use the original values.'},
	{'name': 'SEGM_THRESHOLDS', 'group': 'segm', 'type': 'thresholds', 'default': list(hft_state.DEFAULT_THRESHOLDS),
	 'label': 'Segmentation thresholds',
	 'help': 'Thresholds of i.segment (greater than 0 and smaller than 1, separated by commas, at most 99; they are sorted '
	         'from lowest to highest): lower thresholds give smaller segments. Each threshold starts from the segments of the '
	         'previous one, so the list itself changes the segments. E.g. 0.02, 0.04, 0.06, 0.08, 0.1 for very small habitats, '
	         'with a smaller minimum size of a segment: no segment is smaller than it, whatever the thresholds. Every '
	         'threshold adds time, and the classification into HFTs compares every pair of segments of every threshold '
	         'tested: a threshold with up to 7,500 segments needs about 70 bytes × (number of segments)² (3.9 GB with '
	         '7,500), and one with more is classified with little memory, in a time that grows with (number of '
	         'segments)²: about a minute with 55,000 segments and four with 166,000. Changing them runs '
	         'the segmentation again. With other thresholds than the nine of 0.1 to 0.9, the segm_id of a segment has 2 '
	         'digits for the threshold and 7 for the segment (1 and 5 with those nine).'},
	{'name': 'SMALL_SEGM_RADIUS', 'group': 'segm', 'type': 'int', 'default': 10, 'min': 1, 'step': 1, 'unit': 'pixels',
	 'label': 'Radius to merge small segments',
	 'help': 'i.segment already merges every segment smaller than the minimum size into its most similar neighbour; this '
	         'only acts on the pieces left smaller (cut by pixels without data, or by the outline of the study area): they '
	         'are merged with a larger segment within this distance (see the rule below), in repeated passes; those '
	         'without any so close, with the nearest one. In pixels: 10 = 100 m with 10 m pixels.'},
	{'name': 'SMALL_SEGM_MERGE', 'group': 'segm', 'type': 'select', 'default': 'frequent',
	 'label': 'Merge small segments into',
	 'choices': [['frequent', 'The most frequent larger segment within the radius (as before)'],
	             ['similar', 'The most similar larger segment within the radius']],
	 'help': '"frequent": the larger segment with the most pixels within the radius, whether it is similar or not (the '
	         'rule up to now). "similar": the larger segment whose means in the layers of the segmentation (PC1-PC3, or the '
	         'variables), each scaled to its range as i.segment does, are the closest.'},
	{'name': 'MIN_SEGM_SIZE', 'group': 'segm', 'type': 'int', 'default': 0, 'min': 0, 'step': 1, 'unit': 'pixels',
	 'label': 'Minimum size of a segment',
	 'help': 'The minimum size of i.segment: segments with fewer pixels are merged into their most similar neighbour '
	         '(and the pieces left, see the radius above). 0 = automatic: the square root of the number of pixels with data of each area, so it grows '
	         'with its size. A number sets the same minimum in every area (100 = 1 ha with 10 m pixels).'},
	# Threshold
	{'name': 'THRESHOLD_METHOD', 'group': 'thr', 'type': 'threshold', 'default': 'f', 'label': 'Choice of the segmentation threshold',
	 'choices': [['f', 'F-function (Johnson et al. 2015, as i.segment.uspo)'],
	             ['gs', 'Global Score (Espindola et al. 2006)'],
	             ['gs_fixed', 'Global Score without rescaling between thresholds'],
	             ['crossing', 'First crossing of MI and V minus one threshold, the previous one (original eHabitat+)'],
	             ['crossing_first', 'First crossing of MI and V'],
	             ['fixed', 'Fixed threshold']],
	 'help': "Among the segmentation thresholds (by default the nine from 0.1 to 0.9), with the Moran's I of the means of "
	         'the segments (MI, low = segments different from their neighbours) and the variance within them (V, low = '
	         'homogeneous segments). In the test with Guajares, "f" with alpha 0.5 was the most stable. A fixed threshold '
	         'must be one of the segmentation thresholds: the segmentation stops at it and only that threshold is summarised '
	         'and shown, which saves time (the higher the threshold, the less).'},
	{'name': 'THRESHOLD_ALPHA', 'group': 'thr', 'type': 'number', 'default': 0.5, 'min': 0.05, 'step': 0.05,
	 'label': 'Weight of V against MI (alpha of F)', 'only_if': {'THRESHOLD_METHOD': 'f'},
	 'help': '1 = the same weight; 0.5 = V weighs half (higher thresholds, fewer segments); 2 = V weighs double (lower '
	         'thresholds, more segments). Only with the F-function.'},
	{'name': 'MORAN_NEIGHBOURS', 'group': 'thr', 'type': 'select', 'default': 'rook', 'label': "Neighbours for Moran's I",
	 'choices': [['rook', 'Rook: segments that share an edge'],
	             ['queen', 'Queen: also those that only touch at a corner']],
	 'help': '"queen" leaves fewer segments without neighbours when there are many pixels without data between them.'},
	# Classification
	{'name': 'CLASS_METHOD', 'group': 'class', 'type': 'select', 'default': 'ward_area', 'label': 'Classification of the segments',
	 'choices': [['ward_area', 'Ward weighted by area (default)'],
	             ['ward_area_kmeans', 'Ward weighted by area, then weighted k-means'],
	             ['ward', 'Unweighted Ward (original eHabitat+)']],
	 'help': '"ward_area": Ward weighted by the number of pixels of each segment, on its mean of each variable '
	         'standardized with the standard deviation of the pixels (the HFTs minimize the variance of the pixels within '
	         'them). "ward_area_kmeans" improves the separability a little, but the HFTs of the cuts are no longer exactly '
	         'nested. "ward": on the means and variances rescaled to 0-1.'},
	{'name': 'SIMIL', 'group': 'class', 'type': 'select', 'default': 3, 'cast': 'int', 'label': 'Number of HFTs (cut of the tree)',
	 'choices': [[3, 'Cuts at 10, 25, 50, 75 and 90 % of the height of the tree (p10 ... p90)'],
	             [1, 'Cut at the 25th percentile of the merge heights, with an upper limit'],
	             [2, 'No grouping of segments: set by the upper limits']],
	 'help': 'With 3 the results of each cut are saved in p10, p25, p50, p75 and p90 by default (from more HFTs to fewer; the cuts can be changed below) and the '
	         'cut closest to the elbow of the R² curve is suggested (_optim in its file names).'},
	{'name': 'SIMIL_CUTS', 'group': 'class', 'type': 'cuts', 'default': [10, 25, 50, 75, 90], 'only_if': {'SIMIL': 3},
	 'label': 'Cuts of the tree (% of its height)',
	 'help': 'Percentages of the height of the tree (1 to 99, separated by commas) where it is cut; each gives a classification '
	         'in its folder (p15, p40...), from more HFTs (the lowest) to fewer. After a run, change them and run again from '
	         '"plots": only the classification is repeated, and the results of the scenario are replaced (the previous ones go '
	         'to <scenario>_previous).'},
	{'name': 'HFT_SUBCLASSES', 'group': 'class', 'type': 'bool', 'default': False, 'label': 'Subclasses of the pixels of the HFTs',
	 'help': 'The HFTs group whole segments by the mean of their pixels, so they do not separate what changes inside a segment '
	         '(e.g. a wet meadow that is only part of one). Ticked, the pixels of each HFT are also divided into subclasses by '
	         'k-means: a subclass is divided in two, again and again, while the division separates its pixels at least as much '
	         'as the two most similar HFTs of the cut are separated, so the number of subclasses follows from the variability '
	         'of the HFTs of each cut, without choosing it. Each cut gets a map of the subclasses (GeoTIFF '
	         'hft_subclasses_<id>_<k>_<cut>.tif, code = 100 x HFT + subclass) and their table, and its shapefile of the HFTs '
	         'has every segment split into the subclasses of its pixels (fields subclass, hft_<n>_sub_<y>, and sub_code), '
	         'which the similarity can compare instead of the HFTs. The subclasses follow the '
	         'pixels, not the segments, so their map has more isolated pixels than that of the HFTs.'},
	{'name': 'STAT_AGG', 'group': 'class', 'type': 'select', 'default': 'mean', 'label': 'Statistic of the segments',
	 'choices': [['mean', 'Mean'], ['median', 'Median']],
	 'help': 'Statistic of the pixels of each segment and variable used to classify the segments.'},
]
BY_NAME = {p['name']: p for p in PARAMS}
GROUPS = [['run', 'Scenario and run'], ['segm', 'Segmentation'], ['thr', 'Segmentation threshold'],
          ['class', 'Classification into HFTs']]
STEPS = [{'name': n, 'label': lab} for n, lab in (
	('segm', 'Segmentation (GRASS i.segment, {n} thresholds)'), ('stats', 'Statistics of the segments'),
	('moranvar', "Moran's I, variance and threshold"), ('plots', 'HFTs and plots (R)'))]


def thresholds_of(cfg):
	"""The thresholds of the segmentation in CFG (SEGM_THRESHOLDS), or the default ones if it is not valid."""
	try:
		return hft_state.segm_thresholds(cfg)
	except ValueError:
		return list(hft_state.DEFAULT_THRESHOLDS)


def step_label(name, cfg):
	label = next(s['label'] for s in STEPS if s['name'] == name)
	return label.replace('{n}', str(len(thresholds_of(cfg))))


def read_config():
	return hft_state.read_config()


def pa_buffer():
	"""PA_BUFFER of config.py as the model uses it (0 when it is not there), in m."""
	v = read_config().get('PA_BUFFER', 0)
	return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else 0.0


def default_nproc():
	return max(1, multiprocessing.cpu_count() - 1)


def params():
	"""The parameters with their values in config.py and their defaults."""
	cfg = read_config()
	ths = thresholds_of(cfg)
	out = []
	for p in PARAMS:
		q = dict(p)
		v = cfg.get(p['name'], p['default'])
		if p['type'] == 'threshold':
			q['value'] = 'fixed' if isinstance(v, (int, float)) and not isinstance(v, bool) else v
			# a fixed threshold: by default the one in the middle of the thresholds (0.5 of 0.1-0.9)
			q['fixed'] = float(v) if q['value'] == 'fixed' else ths[(len(ths) - 1) // 2]
			q['thresholds'] = ths
		elif p['type'] == 'thresholds':
			try:
				q['value'] = hft_state.parse_thresholds(v)
			except ValueError:
				q['value'] = v  # shown as it is in config.py; saving it sorts it (_thresholds)
		elif p['name'] == 'SMALL_SEGM_MERGE':
			try:
				q['value'] = hft_state.merge_rule(cfg)  # 'Similar' or None in config.py: an option of the list
			except ValueError:
				q['value'] = v
		else:
			q['value'] = v
		out.append(q)
	s = paths.load_settings()
	return {'params': out, 'groups': GROUPS, 'nproc': s.get('nproc') or default_nproc(), 'nproc_max': multiprocessing.cpu_count(),
	        'nproc_default': default_nproc()}


def _thresholds(p, value):
	"""The thresholds of VALUE (a list, or a text with numbers separated by commas, or by semicolons with decimal
	commas: '0,1; 0,2') sorted from the lowest to the highest and without repeated ones, as config.py needs them, checked."""
	text = '' if isinstance(value, (list, tuple)) else str(value or '').strip()
	if isinstance(value, (list, tuple)):
		items = list(value)
	elif ';' in text:  # '0,1; 0,2': decimal commas
		items = [x.replace(',', '.') for x in re.split(r'[;\s]+', text) if x]
	else:
		items = [x for x in re.split(r'[,\s]+', text) if x]
	hint = ' (write the decimals with a point: 0.1, 0.2, or separate the thresholds with semicolons: 0,1; 0,2)' \
		if ';' not in text and re.search(r'\d,\d', text) else ''
	try:
		nums = sorted(float(x) for x in items)
	except (TypeError, ValueError):
		raise ValueError(f"{p['label']}: write numbers between 0 and 1 separated by commas (e.g. 0.1, 0.2, 0.3)")
	nums = [v for i, v in enumerate(nums) if i == 0 or abs(v - nums[i - 1]) > hft_state.THRESHOLD_TOL]
	try:
		return hft_state.parse_thresholds(nums)
	except ValueError as e:
		raise ValueError(f"{p['label']}: {str(e).replace('SEGM_THRESHOLDS: ', '')}{hint}")


def _literal(p, value, fixed=None, thresholds=None):
	"""Python literal of VALUE for config.py, checked. THRESHOLDS: those of the segmentation, for a fixed threshold."""
	name, t = p['name'], p['type']
	if t == 'bool':
		return repr(bool(value))
	if t in ('int', 'number'):
		try:
			v = float(value)
		except (TypeError, ValueError):
			raise ValueError(f"{p['label']}: {value!r} is not a number")
		if not math.isfinite(v):
			raise ValueError(f"{p['label']}: {value!r} is not a number")
		if 'min' in p and v < p['min']:
			raise ValueError(f"{p['label']}: must be at least {p['min']}")
		if 'max' in p and v > p['max']:
			raise ValueError(f"{p['label']}: must be at most {p['max']}")
		if t == 'int' or v.is_integer():
			return str(int(round(v)))
		return repr(v)
	if t == 'select':
		allowed = [c[0] for c in p['choices']]
		if p.get('cast') == 'int':
			value = int(value)
		if value not in allowed:
			raise ValueError(f"{p['label']}: invalid value {value!r}")
		return repr(value)
	if t == 'cuts':
		items = value if isinstance(value, (list, tuple)) else re.split(r'[,;\s]+', str(value or '').strip())
		try:
			cuts = sorted({int(float(x)) for x in items if str(x).strip() != ''})
		except ValueError:
			raise ValueError(f"{p['label']}: write whole percentages separated by commas (e.g. 10, 25, 50, 75, 90)")
		if not cuts or cuts[0] < 1 or cuts[-1] > 99 or any(float(x) != int(float(x)) for x in items if str(x).strip() != ''):
			raise ValueError(f"{p['label']}: write whole percentages between 1 and 99 (e.g. 10, 25, 50, 75, 90)")
		return repr(cuts)
	if t == 'thresholds':
		return repr(_thresholds(p, value))
	if t == 'threshold':
		if value == 'fixed':
			try:
				v = float(fixed)
			except (TypeError, ValueError):
				raise ValueError('Fixed threshold: not a number')
			ths = thresholds or list(hft_state.DEFAULT_THRESHOLDS)
			same = [x for x in ths if abs(x - v) <= hft_state.THRESHOLD_TOL]
			if not same:
				raise ValueError('Fixed threshold: must be one of the segmentation thresholds: '
				                 + ', '.join(map(hft_state.threshold_label, ths)))
			return repr(same[0])
		if value not in THRESHOLD_METHODS:
			raise ValueError(f'Invalid threshold method: {value!r}')
		return repr(value)
	if t in ('text', 'folder'):
		value = str(value or '').strip()
		if name == 'SCENARIO' and value:
			err = hft_state.scenario_error(value)
			if err:
				raise ValueError(f'Scenario {value!r}: {err}')
		if name == 'GRASSDB' and not value:
			raise ValueError('The GRASS database cannot be empty')
		return repr(value)
	raise ValueError(f'unknown type {t}')


def save_params(values, nproc=None):
	"""Write VALUES ({name: value}, and THRESHOLD_FIXED for a fixed threshold) to config.py, keeping its comments."""
	cfg = read_config()
	# the thresholds of the segmentation after saving, which a fixed threshold must be one of
	ths = _thresholds(BY_NAME['SEGM_THRESHOLDS'], values['SEGM_THRESHOLDS']) if 'SEGM_THRESHOLDS' in values else thresholds_of(cfg)
	lits = {}
	for name, value in values.items():
		if name not in BY_NAME:
			continue
		lits[name] = _literal(BY_NAME[name], value, values.get('THRESHOLD_FIXED'), ths)
	if 'THRESHOLD_METHOD' not in lits:  # a fixed threshold of config.py must be one of the new thresholds
		err = hft_state.thresholds_error(dict(cfg, SEGM_THRESHOLDS=ths))
		if err:
			raise ValueError(err)
	changed = []
	for name, lit in lits.items():
		if repr(cfg.get(name)) != lit or name not in cfg:
			hft_state.write_option(name, lit)
			changed.append(name)
	if nproc is not None:
		paths.save_settings(nproc=max(1, int(nproc)))
	return changed


def defaults():
	return {p['name']: p['default'] for p in PARAMS}


# Study area and IDs --------------------------------------------------------------------------------------------------

def study_area():
	"""The study area chosen in the app or, the first time, that of the last setup.py (run_state.json) or the only
	vector file of segm_optim/pa: {file, layer, col_id, ids (to process)}."""
	s = paths.load_settings().get('study_area')
	if s and s.get('file') and os.path.isfile(s['file']):
		return s
	cfg = read_config()
	rec = hft_state.get_setup(cfg.get('STUDY_AREA')) or {}
	f = rec.get('study_area_file')
	if not (f and os.path.isfile(f)):
		cands = sorted(os.path.join(paths.PA_DIR, x) for x in os.listdir(paths.PA_DIR)
		               if x.lower().endswith(inputs.VECTOR_EXT)) if os.path.isdir(paths.PA_DIR) else []
		same = [c for c in cands if inputs.grass_name(c) == cfg.get('STUDY_AREA')]
		f = (same or cands or [None])[0] if (same or len(cands) == 1) else None
	if not f:
		return None
	ids = None
	palist = os.path.join(paths.SEGM, 'palist.csv')
	if os.path.isfile(palist):
		ids = [int(x) for x in open(palist).read().split() if x.lstrip('-').isdigit()] or None
	return {'file': f, 'layer': rec.get('layer'), 'col_id': rec.get('col_id') or cfg.get('COL_ID') or 'cat', 'ids': ids}


def write_palist(ids):
	ids = sorted({int(i) for i in ids})
	if not ids:
		raise ValueError('choose at least one ID of the study area')
	palist = os.path.join(paths.SEGM, 'palist.csv')
	text = ''.join(f'{i}\n' for i in ids)
	old = open(palist).read() if os.path.isfile(palist) else None
	if old is None or sorted(old.split()) != sorted(text.split()):
		with open(palist, 'w') as f:
			f.write(text)
	return ids


# Status --------------------------------------------------------------------------------------------------------------

def setup_needed(cfg=None, state=None):
	"""Why setup.py has to be run with the choices of the app (empty if not): no GRASS project, another study area,
	another set of variables or input files changed."""
	cfg = read_config() if cfg is None else cfg
	state = hft_state.load() if state is None else state
	sa = study_area()
	out = []
	if not sa:
		return ['choose the study area (Input variables)']
	name = inputs.grass_name(sa['file'])
	grassdb = os.path.expanduser(str(cfg.get('GRASSDB') or '~/grassdata'))
	if not os.path.isdir(os.path.join(grassdb, name, 'PERMANENT')):
		out.append(f'no GRASS project {name} in {grassdb}')
	if cfg.get('STUDY_AREA') != name:
		out.append(f"the study area changed ({cfg.get('STUDY_AREA')} -> {name})")
	rec = hft_state.get_setup(name, state) or {}
	if rec:
		if rec.get('study_area_file') and os.path.abspath(rec['study_area_file']) != os.path.abspath(sa['file']):
			out.append('another file of the study area')
		if (rec.get('layer') or None) != (sa.get('layer') or None):
			out.append('another layer of the study area')
		if rec.get('col_id', 'cat') != sa.get('col_id', 'cat'):
			out.append(f"another column of IDs ({rec.get('col_id')} -> {sa.get('col_id')})")
		if os.path.expanduser(rec.get('grassdb', '')) != grassdb:
			out.append('another GRASS database')
	if cfg.get('STUDY_AREA') == name:
		out += hft_state.setup_changes(cfg, state)
	return out


def status():
	"""Study area, variables, scenario and the steps that are up to date or have to run."""
	try:
		cfg, state = read_config(), hft_state.load()
	except Exception as e:
		return {'error': f'config.py cannot be read: {e}'}
	sa = study_area()
	out = {
		'study_area': cfg.get('STUDY_AREA'), 'invars_folder': hft_state.invars_folder(cfg),
		'scenario': hft_state.scenario_name(cfg), 'results_dir': hft_state.results_dir(cfg).replace(os.sep, '/'),
		'variables': list(cfg.get('ENV_VARS', {})), 'resolution': cfg.get('RESOLUTION'),
		'study_area_file': paths.rel(sa['file']) if sa else None,
		'selected_study_area': inputs.grass_name(sa['file']) if sa else None,
		'ids': sa.get('ids') if sa else None,
	}
	try:
		out['all_ids'] = inputs.vector_ids(sa['file'], sa.get('layer'), sa.get('col_id') or 'cat') if sa else []
	except Exception:
		out['all_ids'] = []
	if sa and not out['ids']:
		out['ids'] = out['all_ids']
	try:
		out['setup'] = setup_needed(cfg, state)
	except Exception as e:
		out['setup'] = [f'cannot be checked: {e}']
	steps = []
	try:
		setup = out['setup'] if isinstance(out['setup'], list) else []
		for name, reasons in hft_state.plan(cfg, state, setup):
			steps.append({'name': name, 'label': step_label(name, cfg), 'reasons': reasons})
	except Exception as e:
		out['error'] = f'the state of the steps cannot be checked: {e}'
	out['steps'] = steps
	out['running'] = bool(jobs.running('model'))
	return out


def scenario_check(name):
	"""Whether the scenario NAME already has results made with other options (they would be moved to _previous)."""
	cfg, state = read_config(), hft_state.load()
	name = (name or '').strip() or hft_state.scenario_name(dict(cfg, SCENARIO=''))
	err = hft_state.scenario_error(name)
	if err:
		return {'name': name, 'error': err}
	this = dict(cfg, SCENARIO=name)
	folder = os.path.join(paths.SEGM, hft_state.results_dir(this))
	has = os.path.isdir(folder) and bool(os.listdir(folder))
	try:
		changes = hft_state.scenario_changes(this, state)
	except Exception:
		changes = None
	return {'name': name, 'exists': has, 'changes': changes, 'folder': paths.rel(folder)}


# Run -----------------------------------------------------------------------------------------------------------------

def run(values=None, nproc=None, first=None, scenario=None):
	"""Save the parameters, write palist.csv, run setup.py if needed and run_hft.py (from the step FIRST)."""
	if jobs.running('model'):
		raise ValueError('the model is already running')
	if jobs.running('ecosim'):
		raise ValueError('the similarity is being computed with the results: run the model when it finishes')
	if values:
		save_params(values, nproc)
	sa = study_area()
	if not sa:
		raise ValueError('choose the study area first in "Input variables"')
	cfg = read_config()
	ids = sa.get('ids') or inputs.vector_ids(sa['file'], sa.get('layer'), sa.get('col_id') or 'cat')
	write_palist(ids)
	py = sys.executable
	commands = []
	need = setup_needed(cfg)
	if need:
		cmd = [py, 'setup.py', '--grassdb', str(cfg.get('GRASSDB') or '~/grassdata'), '--study-area', sa['file'],
		       '--col-id', sa.get('col_id') or 'cat']
		if sa.get('layer'):
			cmd += ['--layer', sa['layer']]
		commands.append(cmd)
	scen = (scenario if scenario is not None else cfg.get('SCENARIO')) or ''
	cmd = [py, 'run_hft.py']
	if scen:
		cmd += ['--scenario', scen]
	if first:
		if first not in hft_state.NAMES:
			raise ValueError(f'unknown step {first}')
		cmd += ['--from', first]
	commands.append(cmd)
	env = dict(os.environ, EHAB_PROGRESS='1', EHAB_NPROC=str(int(nproc or paths.load_settings().get('nproc') or default_nproc())))
	title = 'HFT mapping' + (f' (from {first})' if first else '')

	def after(job):
		cfg2 = read_config()
		job.say(f"\nResults in segm_optim/{hft_state.results_dir(cfg2)}".replace(os.sep, '/'))
		return {'study_area': cfg2.get('STUDY_AREA'), 'scenario': hft_state.scenario_name(cfg2)}

	intro = (('setup.py writes the new order of the variables in config.py: ' if var_order.only_order(need) else
	          'The data have to be imported into GRASS (setup.py): ') + '; '.join(need)) if need else None
	job = jobs.start_commands('model', title, commands, env=env, after=after, intro=intro)
	return job.id
