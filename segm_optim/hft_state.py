#### eHabitat+ HFT module: record of the steps run, to know which ones must be run again
#### License: CC BY-SA 3.0

"""
Every step of the HFT module records in run_state.json the inputs it was run
with: the options of config.py it uses, palist.csv, the input files and its
own code. Before running, a step checks that the steps before it are up to
date with the current inputs and, if they are not, stops saying which script
to run. run_hft.py uses the same record to run only the steps that are needed.
It also records the inputs that produced the results of every scenario
(SCENARIO in config.py, results in results/<STUDY_AREA>/<SCENARIO>), to warn
when a scenario name is used again with other options.

    python3 hft_state.py              # status of every step

and, used by moranvar_plots.R:

    python3 hft_state.py check STEP   # exit code 3 if a step before STEP is out of date
    python3 hft_state.py start STEP   # STEP starts with the current inputs
    python3 hft_state.py done STEP    # STEP finished
"""

import glob
import hashlib
import json
import os
import re
import runpy
import sys

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
	('segm', 'segmentation', 'python3 segmentation_pca_par.py', ['segmentation_pca_par.py']),
	('stats', 'statistics of the segments', 'python3 get_mean_median_var.py',
	 ['get_mean_median_var.py', 'ehab_segm_core.py']),
	('moranvar', "Moran's I, variance and threshold", 'python3 moranvar.py', ['moranvar.py']),
	('plots', 'HFTs and plots', 'Rscript moranvar_plots.R', ['moranvar_plots.R', 'CreateRadialPlot.R']),
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


def set_option(text, name, value):
	"""Set NAME = value in the text of config.py keeping its comment, or append it."""
	pattern = re.compile(rf'^{name}\s*=.*?(?=\s+#|$)', re.M)
	line = f'{name} = {value}'
	if pattern.search(text):
		return pattern.sub(lambda m: line, text, count=1)
	return text.rstrip('\n') + '\n' + line + '\n'


def write_option(name, value):
	"""Set NAME = value (a Python literal) in config.py."""
	with open(CONFIG) as f:
		text = f.read()
	with open(CONFIG, 'w') as f:
		f.write(set_option(text, name, value))


def scenario_name(cfg):
	return str(cfg.get('SCENARIO') or DEFAULT_SCENARIO)


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


def raster_files():
	return sorted(os.path.basename(f) for f in glob.glob(os.path.join(INVARS, '*')) if f.lower().endswith(RASTER_EXT))


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


def setup_changes(cfg=None, state=None):
	"""Input files changed since setup.py imported them (empty when setup.py left no record)."""
	cfg = read_config() if cfg is None else cfg
	rec = get_setup(cfg.get('STUDY_AREA'), state)
	if not rec:
		return []
	out = []
	now = {f: file_sig(os.path.join(INVARS, f)) for f in raster_files()}
	old = rec.get('variables', {})
	for f in sorted(set(now) | set(old)):
		if f not in old:
			out.append(f'new file inVars/{f}')
		elif f not in now:
			out.append(f'inVars/{f} removed')
		elif now[f] != old[f]:
			out.append(f'inVars/{f} changed')
	if rec.get('study_area_file') and vector_sigs(rec['study_area_file']) != rec.get('study_area'):
		out.append(f"{os.path.relpath(rec['study_area_file'], HERE)} changed")
	return out


def step_inputs(step, cfg=None, state=None):
	"""What a step depends on, besides the steps before it."""
	cfg = read_config() if cfg is None else cfg

	def opts(*names):
		return {n: repr(cfg.get(n, d)) for n, d in names}

	if step == 'segm':
		inputs = opts(('GRASSDB', None), ('GRASSLOC', None), ('STUDY_AREA', None), ('COL_ID', None),
		              ('RESOLUTION', None), ('ENV_VARS0', None), ('CLIP_TO_PA', False), ('PA_BUFFER', 0),
		              ('SEGM_PCA', True), ('SMALL_SEGM_RADIUS', 5))
		ids = open(PALIST).read().split() if os.path.isfile(PALIST) else []
		inputs['palist.csv'] = ' '.join(sorted(set(ids)))
		setup = get_setup(cfg.get('STUDY_AREA'), state)
		inputs['imported data'] = _hash(setup) if setup else 'unknown'
	elif step == 'stats':
		inputs = opts(('STAT_AGG', 'mean'), ('RESOLUTION', None), ('ENV_VARS', {}))
		for p in cfg.get('ENV_VARS', {}).values():
			inputs['inVars/' + p['file']] = file_sig(os.path.join(INVARS, p['file']))
	elif step == 'moranvar':
		inputs = opts(('THRESHOLD_METHOD', 'f'), ('THRESHOLD_ALPHA', 0.5), ('MORAN_NEIGHBOURS', 'rook'))
	else:
		inputs = opts(('SIMIL', 2), ('CLASS_METHOD', 'ward_area'))
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
			out.append(f'{k} changed')
		else:
			out.append(f'{k} changed ({_short(old.get(k))} -> {_short(new.get(k))})')
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


def plan(cfg=None, state=None):
	"""(step, reasons) for every step: once a step has to run, all the steps after it too."""
	cfg = read_config() if cfg is None else cfg
	state = load() if state is None else state
	out, needed = [], False
	for name in NAMES:
		r = reasons(name, cfg, state)
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
	"""Stop with a message saying what to run when a step before STEP is out of date."""
	msg = not_ready(step)
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
	print(f"Study area: {cfg.get('STUDY_AREA')} | scenario: {scenario_name(cfg)} (results in {results_dir(cfg)})")
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
