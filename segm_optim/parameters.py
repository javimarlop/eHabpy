#### eHabitat+: the parameters of a run in one file, to use them in another study area or share them
#### License: CC BY-SA 3.0

"""
A parameters file (YAML) keeps everything that defines how eHabitat+ treats a study area, but not the study area
itself, so that another user (or you, elsewhere) can load it and run the same method on other data:

- model: the options of the HFT model in config.py (segmentation, threshold, classification, subclasses);
- similarity: the options of the similarity (ECOSIM_*) that do not point to files;
- recipe: the recipe of the variables of MonitorEO (period, season, pixel filters, variables and metrics), that of the
  set of variables of the model (inVars/<set>/gee_vars.yaml) or else the last one of the app, without the comparison
  area;
- variables_order: the order of the variables of the model (the first one defines the pixels analysed).

What belongs to a computer or a study area is left out (LOCAL): the GRASS database and project, the study area, its
IDs and pixel size, the variables imported, the set of variables, the scenario, FORCE_RESTART, and the folders or files
of the comparison area and the ecoregions. Loading the file changes the other options of config.py, writes the recipe
as the last recipe of the app (app_recipe.yaml; the app shows it in Variables from MonitorEO, and get_vars_gee.py
downloads it with --recipe app_recipe.yaml) and the order of the variables of the set of the model, if it has
variables of the same names.

run_hft.py writes one for every run, parameters.yaml in the folder of the results of the scenario; the app saves and
loads them in HFT mapping (Save to a file… and Load from a file…).

    python3 parameters.py save FILE    # the parameters of config.py (and the recipe of the set) to FILE
    python3 parameters.py show FILE    # what loading FILE would change
    python3 parameters.py load FILE    # FILE into config.py, app_recipe.yaml and the order of the variables
"""

import argparse
import copy
import datetime
import math
import os
import re
import sys

import yaml

import hft_state
import pca_axes
import var_order

FORMAT = 'eHabitat+ parameters'
VERSION = 1
FILE = 'parameters.yaml'  # the name of the one of each scenario (run_hft.py)
APP_RECIPE = os.path.join(hft_state.HERE, 'app_recipe.yaml')  # the last recipe of the app (app/ehabitat_app/paths.py)

# The options of config.py of each section, as their values (type) are checked when loaded: bool, int, number (both
# finite, within BOUNDS), pca_variance (0 or 50-100), pca_axes (0 or a number of axes), choice (one of the list), cuts
# (whole percentages 1-99), thresholds, threshold (a method or a number of the thresholds) and cut (of the tree)
MODEL = {
	'PA_BUFFER': 'number', 'CLIP_TO_PA': 'bool', 'SEGM_PCA': 'bool', 'SEGM_PCA_VARIANCE': 'pca_variance',
	'SEGM_PCA_AXES': 'pca_axes', 'CLIP_PERCENTILE': 'number',
	'SEGM_THRESHOLDS': 'thresholds', 'SMALL_SEGM_RADIUS': 'int', 'SMALL_SEGM_MERGE': ('frequent', 'similar'),
	'MIN_SEGM_SIZE': 'int', 'THRESHOLD_METHOD': 'threshold', 'THRESHOLD_ALPHA': 'number',
	'MORAN_NEIGHBOURS': ('rook', 'queen'), 'CLASS_METHOD': ('ward_area', 'ward_area_kmeans', 'ward'),
	'SIMIL': (1, 2, 3), 'SIMIL_CUTS': 'cuts', 'HFT_SUBCLASSES': 'bool', 'STAT_AGG': ('mean', 'median'),
}
SIMILARITY = {
	'ECOSIM_CUT': 'cut', 'ECOSIM_RESOLUTION': 'number', 'ECOSIM_CLASSES': ('hft', 'subclass'),
	'ECOSIM_STAT': ('mean', 'median'), 'ECOSIM_COVARIANCE': ('mcd', 'classic'), 'ECOSIM_NEIGHBOURS': (8, 4),
	'ECOSIM_MIN_COVER': 'number',
}
SECTIONS = {'model': MODEL, 'similarity': SIMILARITY}
# (minimum, maximum) of the numbers, those of the forms of the app (app/ehabitat_app/model.py and ecosim.py)
BOUNDS = {'PA_BUFFER': (0, None), 'CLIP_PERCENTILE': (0, 49.9), 'SMALL_SEGM_RADIUS': (1, None), 'MIN_SEGM_SIZE': (0, None),
          'THRESHOLD_ALPHA': (0.05, None), 'ECOSIM_RESOLUTION': (0, None), 'ECOSIM_MIN_COVER': (0.01, 1)}
# The value that the scripts use when config.py does not have the option (an older config.py, kept by the updates):
# saved as such, so that the file gives the method that was run (getattr of segmentation_pca_par.py, moranvar.py and
# get_mean_median_var.py, read_config of moranvar_plots.R, DEFAULTS of similarity.py)
USED = {'PA_BUFFER': 0, 'CLIP_TO_PA': False, 'SEGM_PCA': True, 'SEGM_PCA_VARIANCE': 0, 'SEGM_PCA_AXES': 0,
        'CLIP_PERCENTILE': 0,
        'SEGM_THRESHOLDS': list(hft_state.DEFAULT_THRESHOLDS), 'SMALL_SEGM_RADIUS': 5, 'SMALL_SEGM_MERGE': 'frequent',
        'MIN_SEGM_SIZE': 0, 'THRESHOLD_METHOD': 'f', 'THRESHOLD_ALPHA': 0.5, 'MORAN_NEIGHBOURS': 'rook',
        'CLASS_METHOD': 'ward_area', 'SIMIL': 2, 'SIMIL_CUTS': [10, 25, 50, 75, 90], 'HFT_SUBCLASSES': False,
        'STAT_AGG': 'mean', 'ECOSIM_CUT': '', 'ECOSIM_RESOLUTION': 0, 'ECOSIM_CLASSES': 'hft', 'ECOSIM_STAT': 'mean',
        'ECOSIM_COVARIANCE': 'mcd', 'ECOSIM_NEIGHBOURS': 8, 'ECOSIM_MIN_COVER': 0.5}
# Options of the model added after the first parameters files: a file without one of them (saved before it existed)
# was made with its value of USED, which a load sets (so that, e.g., the % of the variance of an older file is not
# overridden by the number of axes of the PCA of config.py)
LATER = ('SEGM_PCA_VARIANCE', 'SEGM_PCA_AXES')
LOCAL = ('GRASSDB', 'GRASSLOC', 'STUDY_AREA', 'COL_ID', 'RESOLUTION', 'ENV_VARS0', 'ENV_VARS', 'INVARS_FOLDER',
         'SCENARIO', 'FORCE_RESTART', 'ECOSIM_AREA', 'ECOSIM_ECOREGIONS', 'ECOSIM_MARINE')
THRESHOLD_METHODS = ('f', 'gs', 'gs_fixed', 'crossing', 'crossing_first')
# of a study area or a computer; the folder of the set too, so that a download in another study area does not replace
# the variables of the first one (it is gee_<study area> when loaded)
RECIPE_LOCAL = ('comparison_area', 'comparison_area_layer', 'drive_dir', 'project', 'folder')
CUT_RE = re.compile(r'^(p\d+|results)?$')  # as CUT_VALUE_RE of app/ehabitat_app/ecosim.py

HEADER = ('# Parameters of eHabitat+: the options of the HFT model and of the similarity, the recipe of the variables of\n'
          '# MonitorEO and the order of the variables, without the study area. Load them in the app (HFT mapping, Load\n'
          '# from a file...) or with: cd segm_optim && python3 parameters.py load <this file>\n')


def _key(name):
	"""The name of a variable without its order prefix (1_, 2_...), extension and capitals, as similarity.py pairs them."""
	return re.sub(r'^\d+_', '', os.path.splitext(os.path.basename(str(name)))[0]).lower()


def model_value(name, cfg):
	"""The value of the option NAME of CFG as the scripts use it (USED when it is not there or None)."""
	if name == 'SMALL_SEGM_MERGE':
		try:
			return hft_state.merge_rule(cfg)
		except ValueError:
			return cfg.get(name)
	if name == 'SEGM_THRESHOLDS':
		try:
			return hft_state.segm_thresholds(cfg)
		except ValueError:
			return cfg.get(name)
	v = cfg.get(name)
	return USED.get(name) if v is None else v


def set_recipe(cfg, app_recipe=True):
	"""(recipe, its file) of the set of variables of the model, inVars/<set>/gee_vars.yaml, or else (APP_RECIPE) the last
	recipe of the app; (None, None) without any."""
	folder = os.path.join(hft_state.INVARS, hft_state.invars_folder(cfg))
	for path in (os.path.join(folder, 'gee_vars.yaml'),) + ((APP_RECIPE,) if app_recipe else ()):
		if os.path.isfile(path):
			try:
				with open(path, encoding='utf-8') as f:
					rec = yaml.safe_load(f)
			except (OSError, yaml.YAMLError):
				continue
			if isinstance(rec, dict) and rec.get('variables'):
				return rec, path
	return None, None


def portable_recipe(rec):
	"""The recipe REC without what belongs to a study area or a computer."""
	if not isinstance(rec, dict):
		return None
	rec = copy.deepcopy(rec)
	for k in RECIPE_LOCAL:
		rec.pop(k, None)
	return rec


def variables_order(cfg):
	"""The files of the variables of the model in their order: those of ENV_VARS of config.py (the last setup.py)."""
	ev = cfg.get('ENV_VARS') or {}
	out = []
	for v in ev.values() if isinstance(ev, dict) else []:
		f = v.get('file') if isinstance(v, dict) else None
		if f:
			out.append(os.path.basename(str(f)))
	return out


def collect(cfg=None, recipe=None, recipe_source=None, order=None, app_recipe=True):
	"""The parameters of CFG (config.py by default) as a dict: the sections model and similarity (every option, with
	the value the scripts use), the recipe (that of the set of the model by default, or else with APP_RECIPE the last one
	of the app; False for none) and the order of the variables (that of config.py by default)."""
	cfg = hft_state.read_config() if cfg is None else cfg
	if recipe is None:
		recipe, recipe_source = set_recipe(cfg, app_recipe)
	out = {'format': FORMAT, 'version': VERSION,
	       'saved': datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
	       'from': {'study_area': cfg.get('STUDY_AREA'), 'scenario': hft_state.scenario_name(cfg),
	                'variables': hft_state.invars_folder(cfg) or 'inVars',
	                'recipe': os.path.relpath(recipe_source, hft_state.HERE) if recipe_source else None}}
	for sec, names in SECTIONS.items():
		out[sec] = {n: model_value(n, cfg) for n in names}
	out['recipe'] = portable_recipe(recipe) if recipe else None
	out['variables_order'] = list(variables_order(cfg) if order is None else order)
	return out


def problems(data):
	"""The options of DATA (collect()) that a load would refuse (e.g. a value of a config.py edited by hand outside the
	limits of the forms): [why], [] if none."""
	out = []
	model = data.get('model') or {}
	try:
		ths = hft_state.parse_thresholds(model.get('SEGM_THRESHOLDS', hft_state.DEFAULT_THRESHOLDS))
	except ValueError:
		ths = None
	for sec, names in SECTIONS.items():
		for name, value in (data.get(sec) or {}).items():
			if name in names:
				try:
					check(name, value, ths)
				except (ValueError, OverflowError) as e:
					out.append(str(e))
	return out


def save(path, data=None):
	"""Write DATA (collect() by default) to the YAML file PATH; returns PATH."""
	data = collect() if data is None else data
	tmp = path + '.tmp'
	with open(tmp, 'w', encoding='utf-8') as f:
		f.write(HEADER)
		yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
	os.replace(tmp, path)
	return path


def read(path):
	"""The parameters of the file PATH, checked: ValueError saying what is wrong with it."""
	try:
		with open(path, encoding='utf-8-sig') as f:
			data = yaml.safe_load(f)
	except yaml.YAMLError as e:
		raise ValueError(f'{os.path.basename(path)} is not valid YAML: {e}') from e
	if not isinstance(data, dict) or data.get('format') != FORMAT:
		raise ValueError(f'{os.path.basename(path)} is not a parameters file of eHabitat+ (it has no format: {FORMAT})')
	v = data.get('version')
	if not isinstance(v, int) or v < 1:
		raise ValueError(f'{os.path.basename(path)}: invalid version {v!r}')
	if v > VERSION:
		raise ValueError(f'{os.path.basename(path)} was saved by a newer eHabitat+ (version {v} of the file): update '
		                 'eHabitat+ to load it')
	for sec in SECTIONS:
		if data.get(sec) is not None and not isinstance(data[sec], dict):
			raise ValueError(f'{os.path.basename(path)}: {sec} must be a list of options')
	err = recipe_error(data.get('recipe'))
	if err:
		raise ValueError(f'{os.path.basename(path)}: the recipe of MonitorEO is not valid: {err}')
	order = data.get('variables_order')
	if order is not None and not (isinstance(order, list) and all(isinstance(x, str) for x in order)):
		raise ValueError(f'{os.path.basename(path)}: variables_order must be a list of file names')
	return data


def recipe_error(rec):
	"""Why the recipe REC of a parameters file cannot be the recipe of the app (None if it can, or if there is none)."""
	if rec is None:
		return None
	if not isinstance(rec, dict):
		return 'it is not a recipe'
	if not all(isinstance(k, str) for k in rec):
		return 'its keys must be names'
	v = rec.get('variables')
	if not isinstance(v, list) or not v or not all(isinstance(x, dict) and all(isinstance(k, str) for k in x) for x in v):
		return 'it has no list of variables'
	for k, kind in (('where', dict), ('rgb', dict), ('exclude_landcover', list)):
		if rec.get(k) is not None and not isinstance(rec[k], kind):
			return f'{k} is not valid'
	for k in ('study', 'quality'):  # the app reads them as lists of options (gee._normal); ecoregion can be empty
		if k in rec and not isinstance(rec[k], dict):
			return f'{k} must be a list of options'
	if rec.get('ecoregion') is not None and not isinstance(rec['ecoregion'], dict):
		return 'ecoregion must be a list of options'
	for k in ('period', 'season'):
		x = rec.get(k)
		if x is not None and not (isinstance(x, list) and len(x) == 2 and all(isinstance(d, str) for d in x)):
			return f'{k} must be two dates'
	return None


def check(name, value, thresholds=None):
	"""The Python literal of VALUE for the option NAME of config.py; ValueError if it is not valid. THRESHOLDS: those of
	the segmentation, which a fixed THRESHOLD_METHOD must be one of."""
	kind = MODEL.get(name) or SIMILARITY.get(name)
	bad = ValueError(f'{name}: invalid value {value!r}')
	if kind is None:
		raise ValueError(f'{name} is not an option that a parameters file can change')
	if isinstance(kind, tuple):
		if isinstance(value, float) and math.isfinite(value) and value.is_integer():
			value = int(value)  # 1.0 for SIMIL 1
		if isinstance(value, bool) or value not in kind:
			raise ValueError(f"{name}: {value!r} is not one of {', '.join(map(repr, kind))}")
		return repr(value)
	if kind == 'bool':
		if not isinstance(value, bool):
			raise bad
		return repr(value)
	if kind in ('int', 'number'):
		if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
			raise ValueError(f'{name}: {value!r} is not a number')
		lo, hi = BOUNDS.get(name, (0, None))
		if value < lo or (hi is not None and value > hi):
			raise ValueError(f'{name}: {value!r} must be ' + (f'from {lo} to {hi}' if hi is not None else f'at least {lo}'))
		if kind == 'int' and value != int(value):
			raise ValueError(f'{name}: {value!r} is not a whole number')
		return repr(int(value)) if float(value).is_integer() else repr(float(value))
	if kind == 'pca_variance':
		try:
			v = pca_axes.variance_option({name: value})
			return repr(int(v) if float(v).is_integer() else v)  # 90.0 as 90, as config.py has it
		except ValueError as e:
			raise ValueError(str(e).replace(' in config.py', '')) from e
	if kind == 'pca_axes':
		try:
			return repr(pca_axes.axes_option({name: value}))  # 4.0 as 4
		except ValueError as e:
			raise ValueError(str(e).replace(' in config.py', '')) from e
	if kind == 'cuts':
		if not (isinstance(value, list) and value and all(isinstance(x, int) and not isinstance(x, bool) and 1 <= x <= 99
		                                                  for x in value)):
			raise ValueError(f'{name}: {value!r} is not a list of whole percentages from 1 to 99')
		return repr(sorted(set(value)))
	if kind == 'thresholds':
		try:
			return repr(hft_state.parse_thresholds(value))
		except ValueError as e:
			raise ValueError(str(e)) from e
	if kind == 'threshold':
		if isinstance(value, str):
			if value not in THRESHOLD_METHODS:
				raise ValueError(f"{name}: {value!r} is not a method ({', '.join(THRESHOLD_METHODS)}) nor a threshold")
			return repr(value)
		if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
			raise bad
		ths = thresholds or list(hft_state.DEFAULT_THRESHOLDS)
		same = [x for x in ths if abs(x - value) <= hft_state.THRESHOLD_TOL]
		if not same:
			raise ValueError(f'{name}: the fixed threshold {value} is not one of the segmentation thresholds')
		return repr(same[0])
	if kind == 'cut':
		if not isinstance(value, str) or not CUT_RE.match(value):
			raise ValueError(f'{name}: {value!r} is not a cut of the tree (p10 … p90, results or empty)')
		return repr(value)
	raise bad


def order_for_set(order, folder):
	"""The files of the set FOLDER (a folder of inVars) in the order ORDER (names of files of another area: they are
	paired by name, without their order prefix and capitals), then the rest as they are; (files, the names of ORDER
	without a file in the set). [] if no file of the set is in ORDER."""
	try:
		files = sorted(f for f in os.listdir(folder) if f.lower().endswith(('.tif', '.tiff')))
	except OSError:
		files = []
	by_key = {}
	for f in files:
		by_key.setdefault(_key(f), f)
	ordered = [by_key[_key(n)] for n in order if _key(n) in by_key]
	ordered = list(dict.fromkeys(ordered))
	missing = [n for n in order if _key(n) not in by_key]
	if not ordered:
		return [], missing
	return ordered + [f for f in files if f not in ordered], missing


def plan(data, cfg=None, study_area=None):
	"""What loading DATA (read()) would do with config.py CFG: {'options': [{name, section, old, new, literal}] (only
	the changes), 'errors': [why an option cannot be loaded], 'ignored': [options of other areas or unknown], 'recipe':
	{variables, folder (that of the download in this study area: recipe_folder), crs (a fixed CRS of the recipe, or
	None)} or None, 'order': {folder, files, missing} or None}. STUDY_AREA: the name of the study area chosen (else that
	of config.py)."""
	cfg = hft_state.read_config() if cfg is None else cfg
	out = {'options': [], 'errors': [], 'ignored': [], 'recipe': None, 'order': None}
	model = data.get('model') or {}
	try:
		ths = hft_state.parse_thresholds(model['SEGM_THRESHOLDS']) if 'SEGM_THRESHOLDS' in model else \
			hft_state.segm_thresholds(cfg)
	except ValueError:
		ths = None
	for sec, names in SECTIONS.items():
		values = dict(data.get(sec) or {})
		filled = {n for n in LATER if n not in values} if sec == 'model' and values else set()
		values.update({n: USED[n] for n in filled})
		for name, value in values.items():
			if name not in names:
				out['ignored'].append(f'{name} ({"of a study area or a computer" if name in LOCAL else "unknown"})')
				continue
			try:
				lit = check(name, value, ths)
			except ValueError as e:
				out['errors'].append(str(e))
				continue
			if (name not in cfg and name not in filled) or repr(model_value(name, cfg)) != lit:
				out['options'].append({'name': name, 'section': sec, 'old': model_value(name, cfg), 'new': value,
				                       'literal': lit})
	rec = portable_recipe(data.get('recipe'))
	err = recipe_error(rec)
	if err:
		out['errors'].append(f'the recipe of MonitorEO is not valid: {err}')
	elif rec:
		crs = [str(c) for c in (rec.get('crs'), (rec.get('study') or {}).get('crs'), (rec.get('ecoregion') or {}).get('crs'))
		       if c not in (None, '', 'auto')]
		out['recipe'] = {'variables': len(rec['variables']), 'folder': recipe_folder(cfg, study_area),
		                 'crs': ', '.join(dict.fromkeys(crs)) or None}
	order = data.get('variables_order') or []
	if order:
		folder = os.path.join(hft_state.INVARS, hft_state.invars_folder(cfg))
		files, missing = order_for_set(order, folder)
		out['order'] = {'folder': hft_state.invars_folder(cfg) or 'inVars', 'files': files, 'missing': missing}
	return out


def recipe_folder(cfg, study_area=None):
	"""The folder of inVars where the loaded recipe downloads its variables: gee_<study area> (STUDY_AREA, the study
	area chosen; else that of config.py, which is that of the last setup.py), with _2, _3... if it exists already, so
	that the download neither replaces the variables of another study area nor mixes them with those of another recipe."""
	sa = re.sub(r'[^A-Za-z0-9_-]', '_', str(study_area or cfg.get('STUDY_AREA') or '')).strip('_')
	base = f'gee_{sa}' if sa else 'gee'
	name, i = base, 1
	while os.path.exists(os.path.join(hft_state.INVARS, name)):
		i += 1
		name = f'{base}_{i}'
	return name


def apply(data, cfg=None, write_recipe=True, write_order=True, study_area=None):
	"""Load DATA (read()) into config.py, app_recipe.yaml and the order of the variables of the set of the model;
	ValueError (nothing written) if an option is not valid. STUDY_AREA: as in plan(). Returns plan()."""
	cfg = hft_state.read_config() if cfg is None else cfg
	p = plan(data, cfg, study_area)
	if p['errors']:
		raise ValueError('; '.join(p['errors']))
	new = dict(cfg)
	new.update({o['name']: o['new'] for o in p['options']})
	err = hft_state.thresholds_error(new)
	if err:
		raise ValueError(err)
	for o in p['options']:
		hft_state.write_option(o['name'], o['literal'])
	if write_recipe and p['recipe']:
		rec = {'folder': p['recipe']['folder'], **portable_recipe(data['recipe'])}
		tmp = APP_RECIPE + '.tmp'
		with open(tmp, 'w', encoding='utf-8') as f:
			f.write('# Recipe of variables of MonitorEO for eHabitat+ (get_vars_gee.py), loaded from a parameters file\n')
			yaml.safe_dump(rec, f, sort_keys=False, allow_unicode=True, default_flow_style=None)
		os.replace(tmp, APP_RECIPE)
	if write_order and p['order'] and p['order']['files']:
		folder = os.path.join(hft_state.INVARS, hft_state.invars_folder(cfg))
		if var_order.read_order(folder) != p['order']['files']:
			var_order.write_order(folder, p['order']['files'])
		else:
			p['order']['unchanged'] = True
	return p


def describe(p):
	"""The lines that tell what plan() P does."""
	lines = [f"  {o['name']}: {o['old']!r} -> {o['new']!r}" for o in p['options']] or ['  no option of config.py changes']
	lines += [f'  not valid: {e}' for e in p['errors']]
	lines += [f'  left out: {x}' for x in p['ignored']]
	if p['recipe']:
		lines.append(f"  recipe of MonitorEO with {p['recipe']['variables']} variables -> app_recipe.yaml, to download into "
		             f"inVars/{p['recipe']['folder']} (python3 get_vars_gee.py --recipe app_recipe.yaml)")
		if p['recipe']['crs']:
			lines.append(f"  the recipe fixes the CRS {p['recipe']['crs']}: check that it suits this study area")
	if p['order']:
		o = p['order']
		if o['files']:
			lines.append(f"  order of the variables of {o['folder']}: {', '.join(o['files'])}")
		else:
			lines.append(f"  order of the variables: none of them is in {o['folder']} (it does not change; load the file "
			             'again after downloading its variables to apply it)')
		if o['missing'] and o['files']:
			lines.append(f"  variables of the file that {o['folder']} does not have: {', '.join(o['missing'])}")
	return lines


def main(argv=None):
	ap = argparse.ArgumentParser(description='Save the parameters of eHabitat+ to a file, or load them from one.')
	ap.add_argument('action', choices=('save', 'show', 'load'))
	ap.add_argument('file')
	a = ap.parse_args(argv)
	a.file = os.path.abspath(os.path.join(START, a.file))
	os.chdir(hft_state.HERE)
	if a.action == 'save':
		data = collect()
		save(a.file, data)
		for w in problems(data):
			print(f'WARNING: a load of this file would refuse it: {w}')
		print(f"Parameters saved in {a.file}: {len(data['model'])} options of the model, {len(data['similarity'])} of the "
		      f"similarity, {'the recipe ' + data['from']['recipe'] if data['recipe'] else 'no recipe'} and the order of "
		      f"{len(data['variables_order'])} variables")
		return
	try:
		data = read(a.file)
		p = plan(data) if a.action == 'show' else apply(data)
	except (ValueError, OSError, OverflowError) as e:
		sys.exit(f'ERROR: {e}')
	print(('Loading ' if a.action == 'show' else 'Loaded ') + a.file + ':')
	print('\n'.join(describe(p)))
	if a.action == 'load' and p['options']:
		print('Run setup.py if the variables changed, then run_hft.py')


START = os.getcwd()  # the folder of a relative FILE of the command line

if __name__ == '__main__':
	main()
