"""Parameters files (segm_optim/parameters.py): the options of the HFT model and of the similarity, the recipe of
MonitorEO and the order of the variables in one file, saved to use them in another study area or share them, and
loaded with the checks of the forms of the app."""

import os

from . import ecosim, inputs, jobs, model

import parameters  # segm_optim (paths puts it in sys.path)


def _label(name):
	p = model.BY_NAME.get(name) or ecosim.BY_NAME.get(name)
	return p['label'] if p else name


def _text(v):
	if isinstance(v, bool):
		return 'yes' if v else 'no'
	if isinstance(v, list):
		return ', '.join(str(x) for x in v)
	return '(empty)' if v in ('', None) else str(v)


def _study_area_name():
	"""The name of the study area chosen in the app (that of its file, as setup.py names it), or None."""
	try:
		sa = model.study_area()
		return inputs.grass_name(sa['file']) if sa and sa.get('file') else None
	except (OSError, ValueError, KeyError):
		return None


def export(path):
	"""Write the parameters of config.py (and the recipe of the set of the model) to PATH (.yaml added if it has no
	extension of YAML); what the file has."""
	path = str(path or '').strip()
	if not path:
		raise ValueError('choose the file')
	if not path.lower().endswith(('.yaml', '.yml')):
		path += '.yaml'
	path = os.path.abspath(os.path.expanduser(path))
	data = parameters.collect()
	parameters.save(path, data)
	return {'path': path, 'model': len(data['model']), 'similarity': len(data['similarity']),
	        'recipe': data['from']['recipe'], 'variables': len(data['variables_order']), 'problems': parameters.problems(data)}


def preview(path):
	"""What loading the parameters file PATH would change, for the confirmation of the app."""
	data = parameters.read(path)
	p = parameters.plan(data, study_area=_study_area_name())
	return {'file': os.path.basename(path), 'saved': data.get('saved'), 'from': data.get('from') or {},
	        'options': [{'name': o['name'], 'section': o['section'], 'label': _label(o['name']),
	                     'old': _text(o['old']), 'new': _text(o['new'])} for o in p['options']],
	        'errors': p['errors'], 'ignored': p['ignored'], 'recipe': p['recipe'], 'order': p['order']}


def load(path):
	"""Load the parameters file PATH: its options into config.py (checked as in the forms, nothing is written if one is
	not valid), its recipe as the last recipe of the app and the order of the variables of the set of the model."""
	if jobs.running('model') or jobs.running('ecosim') or jobs.running('gee'):
		raise ValueError('the model, the similarity or a download of MonitorEO is running: load the parameters when it '
		                 'finishes')
	data = parameters.read(path)
	sa = _study_area_name()
	p = parameters.plan(data, study_area=sa)
	if p['errors']:
		raise ValueError('; '.join(p['errors']))
	mvals = {o['name']: o['new'] for o in p['options'] if o['section'] == 'model'}
	evals = {o['name']: o['new'] for o in p['options'] if o['section'] == 'similarity'}
	if 'THRESHOLD_METHOD' in mvals and not isinstance(mvals['THRESHOLD_METHOD'], str):
		mvals['THRESHOLD_FIXED'] = mvals['THRESHOLD_METHOD']
		mvals['THRESHOLD_METHOD'] = 'fixed'
	# every option checked before writing any
	cfg = model.read_config()
	ths = model._thresholds(model.BY_NAME['SEGM_THRESHOLDS'], mvals['SEGM_THRESHOLDS']) if 'SEGM_THRESHOLDS' in mvals \
		else model.thresholds_of(cfg)
	for name, value in mvals.items():
		if name in model.BY_NAME:
			model._literal(model.BY_NAME[name], value, mvals.get('THRESHOLD_FIXED'), ths)
	for name, value in evals.items():
		ecosim._literal(ecosim.BY_NAME[name], value)
	changed = (model.save_params(mvals) if mvals else []) + (ecosim.save(evals) if evals else [])
	done = parameters.apply(data, study_area=sa)  # the recipe and the order of the variables (the options are written)
	return {'changed': changed, 'recipe': done['recipe'], 'order': done['order'], 'ignored': p['ignored']}
