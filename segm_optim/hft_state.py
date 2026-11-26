#### eHabitat+ HFT module: record of the steps run, to know which ones must be run again
#### License: CC BY-SA 3.0

"""
Every step of the HFT module records in run_state.json the inputs it was run
with: the options of config.py it uses, palist.csv, the input files and its
own code. Before running, a step checks that the steps before it are up to
date with the current inputs and, if they are not, stops saying which script
to run. run_hft.py uses the same record to run only the steps that are needed.
It also records the inputs that produced the results of every scenario
(SCENARIO in config.py, or the name of INVARS_FOLDER; results in
results/<STUDY_AREA>/<SCENARIO>), to warn when a scenario name is used again
with other options.

    python3 hft_state.py              # status of every step

and, used by moranvar_plots.R:

    python3 hft_state.py check STEP   # exit code 3 if a step before STEP is out of date
    python3 hft_state.py start STEP   # STEP starts with the current inputs
    python3 hft_state.py done STEP    # STEP finished
"""

import ast
import glob
import hashlib
import json
import os
import posixpath
import re
import runpy
import sys

import pca_axes
import var_order

HERE = os.path.dirname(os.path.abspath(__file__))
STATE = os.path.join(HERE, 'run_state.json')
CONFIG = os.path.join(HERE, 'config.py')
PALIST = os.path.join(HERE, 'palist.csv')
INVARS = os.path.normpath(os.path.join(HERE, '..', 'inVars'))

RASTER_EXT = ('.tif', '.tiff')
SHP_PARTS = ('.shp', '.dbf', '.shx', '.prj', '.cpg')
NOT_READY = 3  # exit code of 'check' when a step before is out of date
DEFAULT_SCENARIO = 'base'
SCENARIO_RE = re.compile(r'[A-Za-z0-9_-]+')  # also checked in moranvar_plots.R

# Steps in order: name, what it does, command that runs it, files with its code
STEPS = [
	('segm', 'segmentation', 'python3 segmentation_pca_par.py', ['segmentation_pca_par.py', 'pca_axes.py']),
	('stats', 'statistics of the segments', 'python3 get_mean_median_var.py',
	 ['get_mean_median_var.py', 'ehab_segm_core.py']),
	('moranvar', "Moran's I, variance and threshold", 'python3 moranvar.py', ['moranvar.py']),
	('plots', 'HFTs and plots', 'Rscript moranvar_plots.R', ['moranvar_plots.R', 'plot_utils.R', 'CreateRadialPlot.R',
	                                                         'subclasses.R', 'ward_lowmem.py']),
]
NAMES = [s[0] for s in STEPS]
LABEL = {s[0]: s[1] for s in STEPS}
COMMAND = {s[0]: s[2] for s in STEPS}
CODE = {s[0]: s[3] for s in STEPS}


def load():
	try:
		with open(STATE) as f:
			return json.load(f)
	except (OSError, ValueError):
		return {}


def save(state):
	tmp = STATE + '.tmp'
	with open(tmp, 'w') as f:
		json.dump(state, f, indent=1, sort_keys=True)
	os.replace(tmp, STATE)


def read_config():
	"""Options of config.py as they are now (not as imported when the script started)."""
	return runpy.run_path(CONFIG)


def _assignment_span(text, name):
	"""(start, end) in TEXT of the last statement NAME = value at the top level of config.py, from the start of its line
	to the end of the value (so a value written over several lines is replaced whole, and the comment after it is
	kept), or None if there is none or config.py cannot be parsed."""
	try:
		tree = ast.parse(text)
	except (SyntaxError, ValueError):
		return None
	node = None
	for n in tree.body:
		if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name) and n.targets[0].id == name:
			node = n  # the last one, which is the value that Python keeps
	if node is None or node.end_lineno is None:
		return None
	lines = text.splitlines(keepends=True)
	offset = lambda i: sum(len(x) for x in lines[:i - 1])  # first character of line i (from 1)
	end_line = lines[node.end_lineno - 1]
	end_col = len(end_line.encode('utf-8')[:node.end_col_offset].decode('utf-8', errors='ignore'))  # ast counts bytes
	return offset(node.lineno), offset(node.end_lineno) + end_col


def set_option(text, name, value):
	"""Set NAME = value in the text of config.py keeping its comment, or append it. The whole value is replaced, also
	when it is written over several lines (e.g. a long list)."""
	line = f'{name} = {value}'
	span = _assignment_span(text, name)
	if span:
		return text[:span[0]] + line + text[span[1]:]
	pattern = re.compile(rf'^{name}\s*=.*?(?=\s+#|$)', re.M)  # config.py does not parse: a value on one line
	if pattern.search(text):
		return pattern.sub(lambda m: line, text, count=1)
	return text.rstrip('\n') + '\n' + line + '\n'


def write_option(name, value):
	"""Set NAME = value (a Python literal) in config.py."""
	with open(CONFIG) as f:
		text = f.read()
	with open(CONFIG, 'w') as f:
		f.write(set_option(text, name, value))


def invars_folder(cfg):
	"""Subfolder of ../inVars with the variables (INVARS_FOLDER in config.py); '' for ../inVars itself."""
	return str(cfg.get('INVARS_FOLDER') or '').strip().strip('/\\').replace('\\', '/')


def scenario_name(cfg):
	"""SCENARIO of config.py; if empty, the name of INVARS_FOLDER (as in moranvar_plots.R), or 'base'."""
	return str(cfg.get('SCENARIO') or re.sub(r'[^A-Za-z0-9_-]', '_', invars_folder(cfg)) or DEFAULT_SCENARIO)


def scenario_error(name):
	"""Why NAME cannot be a scenario name (None if it can)."""
	if not SCENARIO_RE.fullmatch(name):
		return 'use only letters, digits, - and _'
	if name.endswith('_previous'):
		return "it cannot end in '_previous' (used for the previous results of a scenario)"
	return None


def results_dir(cfg):
	return os.path.join('results', str(cfg.get('STUDY_AREA')), scenario_name(cfg))


def file_sig(path):
	"""Size and modification time of a file: they change when the file is replaced or edited."""
	try:
		st = os.stat(path)
	except OSError:
		return 'missing'
	return f'{st.st_size} {st.st_mtime_ns}'


def vector_sigs(path):
	"""Signatures of a vector file and, for a shapefile, of its .dbf, .shx, .prj and .cpg."""
	stem, ext = os.path.splitext(path)
	parts = SHP_PARTS if ext.lower() == '.shp' else (ext,)
	return {os.path.basename(stem + e): file_sig(stem + e) for e in parts if os.path.exists(stem + e)}


def raster_files(folder=''):
	"""GeoTIFFs of ../inVars/FOLDER, as paths relative to ../inVars ('/' as separator)."""
	return sorted(posixpath.join(folder, os.path.basename(f)) for f in glob.glob(os.path.join(INVARS, folder, '*'))
	              if f.lower().endswith(RASTER_EXT) and os.path.isfile(f))


def in_folder(files, folder):
	"""Entries of FILES (paths relative to ../inVars) that are directly in ../inVars/FOLDER."""
	return {f: v for f, v in files.items() if posixpath.dirname(f) == folder}


def where(folder):
	return f'inVars/{folder}' if folder else 'inVars'


def _hash(obj):
	return hashlib.sha1(json.dumps(obj, sort_keys=True).encode()).hexdigest()[:16]


def _code_sig(files):
	h = hashlib.sha1()
	for f in files:
		with open(os.path.join(HERE, f), 'rb') as fh:
			h.update(fh.read())
	return h.hexdigest()[:12]


def get_setup(area, state=None):
	return (load() if state is None else state).get('setup', {}).get(area)


def set_setup(area, rec):
	state = load()
	state.setdefault('setup', {})[area] = rec
	save(state)


def setup_digest(cfg=None, state=None):
	"""Hash of what setup.py imported for the variables of INVARS_FOLDER and the study area."""
	cfg = read_config() if cfg is None else cfg
	rec = get_setup(cfg.get('STUDY_AREA'), state)
	if not rec:
		return 'unknown'
	# the record keeps the files of every folder imported; only those of INVARS_FOLDER count
	rec = {k: v for k, v in rec.items() if k != 'invars_folder'}
	rec['variables'] = in_folder(rec.get('variables', {}), invars_folder(cfg))
	return _hash(rec)


def setup_changes(cfg=None, state=None):
	"""Why setup.py has to be run again: INVARS_FOLDER changed, or input files changed since setup.py
	imported them (only the first when setup.py left no record)."""
	cfg = read_config() if cfg is None else cfg
	folder = invars_folder(cfg)
	# the variables of config.py (ENV_VARS, filled by setup.py) must be those of INVARS_FOLDER
	used = {posixpath.dirname(p['file']) for p in cfg.get('ENV_VARS', {}).values()}
	if used and used != {folder}:
		return [f"INVARS_FOLDER changed: the variables of config.py are those of {', '.join(map(where, sorted(used)))}, "
		        f"not of {where(folder)}"]
	rec = get_setup(cfg.get('STUDY_AREA'), state)
	if not rec:
		return []
	if rec.get('invars_folder', '') != folder:
		return [f"INVARS_FOLDER changed ({where(rec.get('invars_folder', ''))} -> {where(folder)})"]
	if not os.path.isdir(os.path.join(INVARS, folder)):
		return [f'{where(folder)} not found']
	out = []
	now = {f: file_sig(os.path.join(INVARS, f)) for f in raster_files(folder)}
	old = in_folder(rec.get('variables', {}), folder)
	for f in sorted(set(now) | set(old)):
		if f not in old:
			out.append(f'new file inVars/{f}')
		elif f not in now:
			out.append(f'inVars/{f} removed')
		elif now[f] != old[f]:
			out.append(f'inVars/{f} changed')
	if rec.get('study_area_file') and vector_sigs(rec['study_area_file']) != rec.get('study_area'):
		out.append(f"{os.path.relpath(rec['study_area_file'], HERE)} changed")
	# the order of the variables chosen in the app (variables_order.txt of the folder), not yet in config.py
	indir = os.path.join(INVARS, folder)
	prev = var_order.previous_order(cfg.get('ENV_VARS', {}), folder)
	if prev and var_order.read_order(indir):
		want = [os.path.basename(p) for p in var_order.ordered([os.path.join(indir, posixpath.basename(f))
		                                                        for f in now], indir, prev)]
		if [n for n in want if n in prev] != prev:
			out.append(f'{var_order.ORDER_CHANGED} ({where(folder)}/{var_order.ORDER_FILE})')
	return out


def progress(kind, value=''):
	"""Tell the desktop app how far a step is (it sets EHAB_PROGRESS; elsewhere nothing is printed). One line each:
	STEPS a,b,c (the steps to run), STEP name (the one that starts), UNITS n (adds n units of work to the step) and
	DONE n (n units are done). Scripts of the same step in parallel add up."""
	if os.environ.get('EHAB_PROGRESS'):
		os.write(1, f'##{kind} {value}\n'.encode())


# Thresholds of the segmentation (SEGM_THRESHOLDS in config.py). Every threshold is known by its position k = 1..n in
# the sorted list, which names its files (shp/park_segm_<id>_<k>.shp, csv/park_<id>_hri_results<k>.csv...). The
# segmentation writes the value of every k to csv/segm_thresholds.csv, which the later steps read for the labels.
DEFAULT_THRESHOLDS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
MAX_THRESHOLDS = 99  # k has at most two digits (it is part of segm_id)
THRESHOLD_TOL = 1e-9  # a fixed THRESHOLD_METHOD is the threshold of SEGM_THRESHOLDS within this distance
THRESHOLDS_CSV = os.path.join(HERE, 'csv', 'segm_thresholds.csv')
THRESHOLD_METHODS = ('f', 'gs', 'gs_fixed', 'crossing', 'crossing_first')  # THRESHOLD_METHOD, or a fixed threshold
SMALL_SEGM_MERGES = ('frequent', 'similar')  # SMALL_SEGM_MERGE: 'frequent', the rule up to eHabpy 1.5, is the default
SEGM_ID_MAX = 2 ** 63 - 1  # segm_id is a 64-bit integer in SQLite (GRASS) and NumPy: larger ones lose digits
# Memory of the classification into HFTs of moranvar_plots.R per pair of segments of the threshold with the most
# segments: dense matrices of the distances between them (class_tree, the silhouette) of two thresholds at once (the R2
# curves classify every threshold tested), measured: 2.5 GB with 6,000 segments (about 70 bytes x 6,000²)
CLASS_BYTES = 70
# A threshold with more segments is classified by ward_lowmem.py, without those matrices (the same results, with memory
# proportional to the segments; its time grows with segments²: about a minute with 55,000, four with 166,000): with
# this limit the matrices need at most 70 bytes x 7,500² = 3.9 GB. The same number is LOWMEM_SEGMENTS in moranvar_plots.R
LOWMEM_SEGMENTS = 7500


def parse_thresholds(value):
	"""A value of SEGM_THRESHOLDS as the list of its thresholds (floats); ValueError saying what is wrong if it is not a
	list of 1 to MAX_THRESHOLDS numbers greater than 0 and smaller than 1, from the lowest to the highest without
	repeating any (the order is that of the chain of the segmentation and gives the index k of each threshold, so it is
	not changed silently: the app sorts the thresholds typed before writing them)."""
	if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
		raise ValueError(f'SEGM_THRESHOLDS: it must be a list of thresholds between 0 and 1, e.g. {DEFAULT_THRESHOLDS}, not {value!r}')
	out = []
	for v in value:
		if isinstance(v, bool) or not isinstance(v, (int, float)) or v != v:
			raise ValueError(f'SEGM_THRESHOLDS: {v!r} is not a number')
		if not 0 < v < 1:
			raise ValueError(f'SEGM_THRESHOLDS: every threshold must be greater than 0 and smaller than 1, not {v!r}')
		v = float(v)
		if out and abs(v - out[-1]) <= THRESHOLD_TOL:
			raise ValueError(f'SEGM_THRESHOLDS: {threshold_label(out[-1])} is repeated')
		if out and v < out[-1]:
			raise ValueError(f'SEGM_THRESHOLDS: the thresholds must go from the lowest to the highest, '
			                 f'and {threshold_label(v)} comes after {threshold_label(out[-1])}')
		out.append(v)
	if not out:
		raise ValueError('SEGM_THRESHOLDS: the list cannot be empty')
	if len(out) > MAX_THRESHOLDS:
		raise ValueError(f'SEGM_THRESHOLDS: at most {MAX_THRESHOLDS} thresholds, not {len(out)}')
	return out


def segm_thresholds(cfg=None):
	"""The thresholds of the segmentation (SEGM_THRESHOLDS in config.py, DEFAULT_THRESHOLDS without it), from the lowest
	to the highest: k is the position of each one, from 1. ValueError if the option is not valid."""
	cfg = read_config() if cfg is None else cfg
	return parse_thresholds(cfg.get('SEGM_THRESHOLDS', DEFAULT_THRESHOLDS))


def is_default_thresholds(ths):
	"""Whether the thresholds THS (a sorted list) are the nine of 0.1-0.9, DEFAULT_THRESHOLDS: the layout of segm_id
	and the inputs recorded are then those of eHabpy up to 1.5."""
	return (isinstance(ths, (list, tuple)) and len(ths) == len(DEFAULT_THRESHOLDS)
	        and all(isinstance(a, (int, float)) and abs(a - b) <= THRESHOLD_TOL for a, b in zip(ths, DEFAULT_THRESHOLDS)))


def threshold_label(value):
	"""A threshold as text for messages, plots and tables: 0.1, 0.05, 0.125 (no exponent, no trailing zeros)."""
	return f'{float(value):.15f}'.rstrip('0').rstrip('.')


def write_run_thresholds(ths, path=None):
	"""csv/segm_thresholds.csv: 'k value', the value of every threshold k of the segmentation (written by it)."""
	with open(path or THRESHOLDS_CSV, 'w') as f:
		f.write('k value\n' + ''.join(f'{k} {threshold_label(v)}\n' for k, v in enumerate(ths, 1)))


def run_thresholds(path=None, cfg=None):
	"""{k: value} of the thresholds of the last segmentation (csv/segm_thresholds.csv). Without that file (a
	segmentation of eHabpy 1.5 or before), or if it cannot be read, those of SEGM_THRESHOLDS in config.py (CFG), which
	the steps check to be those of the segmentation (check_ready): the nine of 0.1-0.9 of those segmentations, value =
	k / 10, or the default ones if the option is not valid."""
	path = path or THRESHOLDS_CSV
	try:
		with open(path) as f:
			rows = [line.split() for line in f.read().splitlines()[1:] if line.strip()]
		out = {}
		for row in rows:
			k, v = row  # ValueError if a row does not have two fields
			out[int(k)] = float(v)
		if out:
			return dict(sorted(out.items()))
		raise ValueError('no thresholds')
	except OSError:
		pass
	except ValueError as e:
		print(f'WARNING: {os.path.basename(path)} cannot be read ({e}): the thresholds of config.py are used', file=sys.stderr)
	try:
		ths = segm_thresholds(cfg)
	except ValueError:
		ths = DEFAULT_THRESHOLDS
	return dict(enumerate(ths, 1))


def fixed_level(cfg=None):
	"""Level k (1..n: the position of the threshold in SEGM_THRESHOLDS) when THRESHOLD_METHOD in config.py is a fixed
	threshold, else None. The segmentation then goes up the chain of thresholds only as far as that one, and the later
	steps calculate and show only that one. ValueError when the number is not one of the thresholds of SEGM_THRESHOLDS
	(or that option is not valid)."""
	cfg = read_config() if cfg is None else cfg
	t = cfg.get('THRESHOLD_METHOD')
	if isinstance(t, (int, float)) and not isinstance(t, bool):
		ths = segm_thresholds(cfg)
		for k, v in enumerate(ths, 1):
			if abs(t - v) <= THRESHOLD_TOL:
				return k
		raise ValueError(f'THRESHOLD_METHOD = {t!r} in config.py is not one of the thresholds of the segmentation '
		                 f"(SEGM_THRESHOLDS): {', '.join(map(threshold_label, ths))}")
	return None


def segm_id_layout(ths):
	"""(digits of k, digits of the segment number) in segm_id = area ID || k || segment number, both zero-padded:
	(1, 5) with the nine thresholds of 0.1-0.9, as up to eHabpy 1.5 (at most 99,999 segments per threshold), else
	(2, 7) (at most 9,999,999). The same width for all the segments of a run, so a segm_id is unique among the areas
	and the thresholds."""
	return (1, 5) if is_default_thresholds(ths) else (2, 7)


def thresholds_error(cfg=None):
	"""Why SEGM_THRESHOLDS or a fixed THRESHOLD_METHOD of config.py are not valid (None if they are)."""
	cfg = read_config() if cfg is None else cfg
	try:
		segm_thresholds(cfg)
		fixed_level(cfg)
	except ValueError as e:
		return str(e)
	return None


def merge_rule(cfg=None):
	"""SMALL_SEGM_MERGE of config.py as the segmentation uses it, 'frequent' or 'similar' (without the option, or with
	None or '', 'frequent'; in capitals too); ValueError for any other value."""
	cfg = read_config() if cfg is None else cfg
	v = cfg.get('SMALL_SEGM_MERGE', 'frequent')
	rule = 'frequent' if v is None or v == '' else v.strip().lower() if isinstance(v, str) else v
	if rule not in SMALL_SEGM_MERGES:
		raise ValueError(f"SMALL_SEGM_MERGE in config.py must be {' or '.join(map(repr, SMALL_SEGM_MERGES))}, not {v!r}")
	return rule


def segm_id_error(area, ths, kmax=None):
	"""Why the segm_id of the segments of the area with ID AREA could be larger than SEGM_ID_MAX with the thresholds THS
	(segmented up to the threshold KMAX, by default all), else None. The largest one is the area ID, the largest k and
	the largest segment number allowed (segm_id_layout), so it is known before the segmentation. IDs that are not whole
	numbers are not checked."""
	s = str(area).strip()
	digits = s[1:] if s.startswith('-') else s
	if not digits.isdigit():
		return None
	kd, cd = segm_id_layout(ths)
	largest = int(f'{digits}{(kmax or len(ths)):0{kd}d}{10 ** cd - 1}')
	if largest <= SEGM_ID_MAX:
		return None
	return (f'area {s}: the segm_id of its segments could reach {largest} (the area ID followed by {kd} digit(s) for '
	        f'the threshold and {cd} for the segment), larger than {SEGM_ID_MAX}, the largest whole number that GRASS '
	        'and NumPy keep exactly: use a column of shorter IDs for the study areas (COL_ID)'
	        + (', or the nine thresholds of 0.1-0.9 (1 + 5 digits)' if (kd, cd) != (1, 5) else ''))


def options_error(cfg=None):
	"""Why SEGM_THRESHOLDS, THRESHOLD_METHOD, SMALL_SEGM_MERGE, HFT_SUBCLASSES or SEGM_PCA_VARIANCE of config.py are
	not valid, or why
	the segm_id of an area of palist.csv would not fit (segm_id_error), None if all is right: checked before any step
	runs (run_hft.py, check_ready), so that a wrong value does not stop the run after the steps before the one that
	uses it."""
	cfg = read_config() if cfg is None else cfg
	err = thresholds_error(cfg)
	if err:
		return err
	t = cfg.get('THRESHOLD_METHOD', 'f')
	if not (isinstance(t, str) and t in THRESHOLD_METHODS) and (isinstance(t, bool) or not isinstance(t, (int, float))):
		return (f"THRESHOLD_METHOD in config.py must be {', '.join(map(repr, THRESHOLD_METHODS))} or one of the "
		        f'thresholds of SEGM_THRESHOLDS (a number), not {t!r}')
	try:
		merge_rule(cfg)
	except ValueError as e:
		return str(e)
	if not isinstance(cfg.get('HFT_SUBCLASSES', False), bool):
		return f"HFT_SUBCLASSES in config.py must be True or False, not {cfg.get('HFT_SUBCLASSES')!r}"
	try:
		pca_axes.variance_option(cfg)
	except ValueError as e:
		return str(e)
	if os.path.isfile(PALIST):
		ths = segm_thresholds(cfg)
		errs = [e for e in (segm_id_error(a, ths, fixed_level(cfg)) for a in sorted(set(open(PALIST).read().split()))) if e]
		if errs:
			return 'the IDs of palist.csv are too long for segm_id: ' + '; '.join(errs[:3]) + (f' (and {len(errs) - 3} more)' if len(errs) > 3 else '')
	return None


def class_memory_gb(nseg):
	"""Memory, in GB, that the classification into HFTs (moranvar_plots.R) needs when the threshold with the most
	segments that it classifies has NSEG segments (CLASS_BYTES per pair): with more than LOWMEM_SEGMENTS, that of the
	thresholds classified with the matrices of the distances, LOWMEM_SEGMENTS at most."""
	return CLASS_BYTES * min(nseg, LOWMEM_SEGMENTS) ** 2 / 1e9


def step_inputs(step, cfg=None, state=None):
	"""What a step depends on, besides the steps before it."""
	cfg = read_config() if cfg is None else cfg

	def opts(*names):
		return {n: repr(cfg.get(n, d)) for n, d in names}

	# options added later, only recorded when they are not at their default (so that adding them did not
	# make the steps already run out of date)
	def later(*names):
		return {n: repr(cfg.get(n)) for n, d in names if cfg.get(n, d) != d}

	if step == 'segm':
		inputs = opts(('GRASSDB', None), ('GRASSLOC', None), ('STUDY_AREA', None), ('COL_ID', None),
		              ('RESOLUTION', None), ('ENV_VARS0', None), ('CLIP_TO_PA', False), ('PA_BUFFER', 0),
		              ('SEGM_PCA', True), ('SMALL_SEGM_RADIUS', 5))
		inputs.update(later(('CLIP_PERCENTILE', 0), ('MIN_SEGM_SIZE', 0)))
		# the axes of the PCA (SEGM_PCA_VARIANCE; 0 or None, the first 3 axes as before the option, not recorded): only
		# with SEGM_PCA, which ignores it otherwise
		try:
			pv = pca_axes.variance_option(cfg)
		except ValueError:
			pv = cfg.get('SEGM_PCA_VARIANCE')
		if cfg.get('SEGM_PCA', True) and pv:
			inputs['SEGM_PCA_VARIANCE'] = repr(pv)
		# SMALL_SEGM_MERGE as the segmentation reads it (merge_rule: 'Frequent' or None are 'frequent', the default, not
		# recorded); an invalid value is recorded as it is (the scripts stop saying why: check_ready)
		try:
			rule = merge_rule(cfg)
		except ValueError:
			rule = cfg.get('SMALL_SEGM_MERGE')
		if rule != 'frequent':
			inputs['SMALL_SEGM_MERGE'] = repr(rule)
		# the thresholds, as an option of later(): not recorded when they are the nine of 0.1-0.9, so the segmentations made before the option existed stay up to date; an invalid value is recorded
		# as it is (the scripts stop saying why: check_ready)
		try:
			ths = segm_thresholds(cfg)
		except ValueError:
			ths = cfg.get('SEGM_THRESHOLDS')
		if not is_default_thresholds(ths):
			inputs['SEGM_THRESHOLDS'] = repr(ths)
		try:
			fixed = fixed_level(cfg)
		except ValueError:
			fixed = None
		if fixed:
			inputs['thresholds segmented'] = f'1-{fixed}'  # a fixed threshold stops the chain there
		ids = open(PALIST).read().split() if os.path.isfile(PALIST) else []
		inputs['palist.csv'] = ' '.join(sorted(set(ids)))
		inputs['imported data'] = setup_digest(cfg, state)
	elif step == 'stats':
		inputs = opts(('STAT_AGG', 'mean'), ('RESOLUTION', None), ('ENV_VARS', {}))
		for p in cfg.get('ENV_VARS', {}).values():
			inputs['inVars/' + p['file']] = file_sig(os.path.join(INVARS, p['file']))
	elif step == 'moranvar':
		inputs = opts(('THRESHOLD_METHOD', 'f'), ('THRESHOLD_ALPHA', 0.5), ('MORAN_NEIGHBOURS', 'rook'))
	else:
		inputs = opts(('SIMIL', 2), ('CLASS_METHOD', 'ward_area'))
		inputs.update(later(('CLIP_PERCENTILE', 0), ('SIMIL_CUTS', [10, 25, 50, 75, 90]), ('HFT_SUBCLASSES', False)))
		inputs['SCENARIO'] = repr(scenario_name(cfg))  # where the results are written
	inputs['code'] = _code_sig(CODE[step])
	return inputs


def chain(step, cfg=None, state=None):
	"""Hash of the inputs of a step and of all the steps before it."""
	cfg = read_config() if cfg is None else cfg
	h = ''
	for name in NAMES[:NAMES.index(step) + 1]:
		h = _hash([h, step_inputs(name, cfg, state)])
	return h


def record(step, state=None):
	return (load() if state is None else state).get('steps', {}).get(step)


def _short(v, n=40):
	v = str(v)
	return v if len(v) <= n else v[:n - 3] + '...'


def _diff(step, old, new):
	"""Differences between two sets of inputs of a step, as messages."""
	out = []
	for k in sorted(set(old) | set(new)):
		if old.get(k) == new.get(k):
			continue
		if k == 'code':
			out.append(f"{' or '.join(CODE[step])} changed")
		elif k == 'imported data':
			out.append('input data imported again by setup.py')
		elif k.startswith('inVars/'):
			if k in old and k in new:  # files added or left out: already in 'ENV_VARS changed'
				out.append(f'{k} changed')
		else:
			show = lambda v: 'not set' if v is None else _short(v)  # options of later() at their default
			out.append(f'{k} changed ({show(old.get(k))} -> {show(new.get(k))})')
	return out


def reasons(step, cfg=None, state=None):
	"""Why a step has to be run (empty list if it is up to date)."""
	cfg = read_config() if cfg is None else cfg
	state = load() if state is None else state
	rec = record(step, state)
	if not rec:
		return ['not run yet']
	out = []
	if rec['chain'] != chain(step, cfg, state):
		out = _diff(step, rec['inputs'], step_inputs(step, cfg, state))
		if not out:
			out.append('the inputs of a step before it changed')
	if not rec.get('done'):
		out.append('did not finish (interrupted, failed or still running)')
	return out


def plan(cfg=None, state=None, setup=None):
	"""(step, reasons) for every step: once a step has to run, all the steps after it too. SETUP: why setup.py has to
	import the data again (by default setup_changes): then the segmentation has to run too, since the record of the
	steps still has the data imported before."""
	cfg = read_config() if cfg is None else cfg
	state = load() if state is None else state
	setup = setup_changes(cfg, state) if setup is None else setup
	out, needed = [], False
	for name in NAMES:
		r = reasons(name, cfg, state)
		if setup and name == NAMES[0] and r == []:
			r = ['setup.py will write the new order of the variables in config.py' if var_order.only_order(setup)
			     else 'setup.py will import the input data again']
		if needed and not r:
			r = ['a step before it has to run']
		needed = needed or bool(r)
		out.append((name, r))
	return out


def not_ready(step, cfg=None, state=None):
	"""Message saying what to run when the input data or a step before STEP are out of date, else None."""
	cfg = read_config() if cfg is None else cfg
	state = load() if state is None else state
	changes = setup_changes(cfg, state)
	if changes:
		return (f"ERROR: the input data changed since setup.py imported them: {'; '.join(changes)}.\n"
		        "Run 'python3 run_hft.py' (it runs setup.py and every step that is needed), "
		        "or 'python3 setup.py' first.")
	for name, r in plan(cfg, state):
		if name == step:
			return None
		if r:
			return (f"ERROR: run the {LABEL[name]} first ({'; '.join(r)}).\n"
			        f"Run 'python3 run_hft.py' to run every step that is needed, or '{COMMAND[name]}'.")


def check_ready(step):
	"""Stop with a message saying what to run when a step before STEP is out of date, or what is wrong in the options
	of config.py checked before every step (options_error: SEGM_THRESHOLDS, THRESHOLD_METHOD, SMALL_SEGM_MERGE,
	HFT_SUBCLASSES)."""
	err = options_error()
	msg = f'ERROR: {err}' if err else not_ready(step)
	if msg:
		print(msg, file=sys.stderr)
		sys.exit(NOT_READY)


def start(step):
	"""STEP starts with the current inputs (it is out of date until done() is called)."""
	cfg, state = read_config(), load()
	state.setdefault('steps', {})[step] = {'inputs': step_inputs(step, cfg, state),
	                                       'chain': chain(step, cfg, state), 'done': False}
	save(state)


def done(step):
	state = load()
	rec = record(step, state)
	if rec:
		rec['done'] = True
		if step == 'plots':
			# Inputs of all the steps that produced the results of this scenario (without its name)
			cfg = read_config()
			inputs = {s: dict(record(s, state)['inputs']) for s in NAMES if record(s, state)}
			inputs.get('plots', {}).pop('SCENARIO', None)
			state.setdefault('scenarios', {}).setdefault(str(cfg.get('STUDY_AREA')), {})[scenario_name(cfg)] = inputs
		save(state)


def scenario_changes(cfg=None, state=None):
	"""How the current inputs differ from those that produced the results of the scenario:
	None if there is no record of them, an empty list if they are the same."""
	cfg = read_config() if cfg is None else cfg
	state = load() if state is None else state
	rec = state.get('scenarios', {}).get(str(cfg.get('STUDY_AREA')), {}).get(scenario_name(cfg))
	if rec is None:
		return None
	out = []
	for name in NAMES:
		now = step_inputs(name, cfg, state)
		now.pop('SCENARIO', None)
		out += _diff(name, rec.get(name, {}), now)
	return out


def reset(first='segm'):
	"""Forget the steps from FIRST on (their files have been deleted)."""
	state = load()
	for name in NAMES[NAMES.index(first):]:
		state.get('steps', {}).pop(name, None)
	save(state)


def print_status(cfg=None, state=None):
	cfg = read_config() if cfg is None else cfg
	state = load() if state is None else state
	print(f"Study area: {cfg.get('STUDY_AREA')} | variables: {where(invars_folder(cfg))} | "
	      f"scenario: {scenario_name(cfg)} (results in {results_dir(cfg)})")
	changes = setup_changes(cfg, state)
	if changes:
		print(f"   input data: to import again with setup.py ({'; '.join(changes)})")
	for i, (name, r) in enumerate(plan(cfg, state), 1):
		print(f"{i:>3}. {LABEL[name]:<42} {'to run: ' + '; '.join(r) if r else 'up to date'}")


if __name__ == '__main__':
	os.chdir(HERE)
	if len(sys.argv) == 1:
		print_status()
	elif len(sys.argv) == 3 and sys.argv[1] in ('check', 'start', 'done') and sys.argv[2] in NAMES:
		{'check': check_ready, 'start': start, 'done': done}[sys.argv[1]](sys.argv[2])
	else:
		sys.exit(f"usage: python3 hft_state.py [check|start|done {'|'.join(NAMES)}]")
