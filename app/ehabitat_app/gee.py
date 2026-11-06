"""MonitorEO: the catalog of variables, the connection with Google Earth Engine (the shared account of MonitorEO or the
user's own) and the download of a recipe of variables over the study area and its ecoregions with get_vars_gee.py,
which keeps the package monitoreo at its last version of GitHub and computes the variables with it."""

import copy
import importlib.metadata
import os
import shutil
import sys

import yaml

from . import jobs, model, paths

RECIPE = os.path.join(paths.SEGM, 'gee_vars.yaml')  # the recipe of the repository (with its comments)
CREDENTIALS = os.path.join('~', '.config', 'earthengine', 'credentials')

SENSORS = ['Sentinel-2', 'Landsat', 'Modis']
METRICS = ['mean', 'median', 'min', 'max', 'std', 'sum', 'count', 'mode', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95']
TEMPORAL = [['period', 'Whole period'], ['annual', 'Each year, then between years'],
            ['monthly', 'Each month, between years and mean of the months'], ['16-day', 'Each 16 days, between years and mean']]
WORLDCOVER = [[10, 'Tree cover'], [20, 'Shrubland'], [30, 'Grassland'], [40, 'Cropland'], [50, 'Built-up'],
              [60, 'Bare / sparse vegetation'], [70, 'Snow and ice'], [80, 'Permanent water bodies'],
              [90, 'Herbaceous wetland'], [95, 'Mangroves'], [100, 'Moss and lichen']]


def _monitoreo():
	try:
		import monitoreo  # noqa: F401
		from monitoreo import catalog, services
		return catalog, services
	except ImportError as e:
		raise ValueError(f'the package monitoreo is not installed ({e}): press "Update MonitorEO"') from e


def version():
	try:
		return importlib.metadata.version('monitoreo')
	except importlib.metadata.PackageNotFoundError:
		return None


def status():
	out = {'version': version(), 'credentials': os.path.isfile(os.path.expanduser(CREDENTIALS)),
	       'connected': False, 'project': None, 'shared': None, 'can_share': False}
	s = paths.load_settings()
	out['account'] = s.get('ee_account') or ('own' if s.get('ee_project') else 'shared')
	out['ee_project'] = s.get('ee_project') or os.environ.get('EE_PROJECT') or ''
	try:
		_, services = _monitoreo()
	except ValueError:
		return out
	out['can_share'] = hasattr(services, 'shared')
	out['connected'] = bool(services.project())
	out['project'] = services.project()
	out['shared'] = services.shared() if hasattr(services, 'shared') else False
	return out


def catalog():
	"""Variables of the catalog of MonitorEO and the options of a recipe."""
	cat, services = _monitoreo()
	try:
		rows = services.list_variables(lang='en')
	except TypeError:
		rows = services.list_variables()
	src = cat.load().get('source', {})
	return {'variables': rows, 'source': src, 'sensors': SENSORS, 'metrics': METRICS, 'temporal': TEMPORAL,
	        'worldcover': WORLDCOVER}


def describe(ref):
	_, services = _monitoreo()
	d = services.describe_variable(None, ref=int(ref), lang='en')
	return {k: (v if isinstance(v, (str, int, float, bool, list, dict, type(None))) else str(v)) for k, v in d.items()}


def connect(account, project=None):
	"""Connect to Earth Engine (in a job: it may take a while): ACCOUNT 'shared' (the shared account of MonitorEO) or
	'own' (PROJECT, with the credentials of the Google account)."""
	project = (project or '').strip()
	if account == 'own' and not project:
		raise ValueError('write the Google Cloud project registered in Earth Engine')
	paths.save_settings(ee_account=account, ee_project=project or None)

	def work(job):
		_, services = _monitoreo()
		if account == 'shared':
			if not hasattr(services, 'shared'):
				raise ValueError('this version of monitoreo has no shared account: press "Update MonitorEO"')
			job.say('Connecting with the shared account of MonitorEO...')
			os.environ.pop('EE_PROJECT', None)
			got = services.init(None, shared=True)
		else:
			job.say(f'Connecting to Earth Engine with the project {project}...')
			got = services.init(project)
		out = {'project': got, 'shared': services.shared() if hasattr(services, 'shared') else False}
		if hasattr(services, 'eecu_usage'):
			try:
				u = services.eecu_usage(lang='en')
				out['eecu'] = {k: u.get(k) for k in ('text', 'left_hours', 'limit_hours', 'used_percent', 'left_percent',
				                                     'warning', 'account', 'resets')}
			except Exception as e:
				out['eecu'] = {'text': f'The EECUs of this month cannot be read ({e}).'}
		job.say('Connected.' + (f" {out['eecu']['text']}" if out.get('eecu', {}).get('text') else ''))
		return out

	return jobs.start_function('gee-connect', 'Connection to Earth Engine', work).id


def authenticate():
	"""Sign in with a Google account in the browser (once on each computer), as 'earthengine authenticate'."""
	def work(job):
		import ee
		job.say('The web browser opens to sign in with your Google account and give permission to Earth Engine...')
		ee.Authenticate(auth_mode='localhost', force=True)
		job.say('Authenticated: the credentials are kept in ~/.config/earthengine.')
		return True
	return jobs.start_function('gee-auth', 'Authentication with Earth Engine', work).id


def update():
	"""Install the last version of monitoreo of GitHub (as get_vars_gee.py does on every run)."""
	code = 'import get_vars_gee as g; g.update_monitoreo(True)'
	before = version()

	def after(job):
		import importlib
		importlib.invalidate_caches()
		now = None
		try:
			import subprocess
			r = subprocess.run([sys.executable, '-c', 'import importlib.metadata as m; print(m.version("monitoreo"))'],
			                   capture_output=True, text=True, creationflags=jobs.NO_WINDOW)
			now = r.stdout.strip() or None
		except OSError:
			pass
		if now and now != before and 'monitoreo' in sys.modules:
			job.say(f'\nmonitoreo {before} -> {now}: close and open the app again to use the new version in the '
			        'catalog (the downloads already use it).')
		return {'before': before, 'now': now}

	return jobs.start_commands('gee-update', 'Update MonitorEO', [[sys.executable, '-c', code]], after=after).id


# Recipes -------------------------------------------------------------------------------------------------------------

def _read_yaml(path):
	with open(path, encoding='utf-8') as f:
		rec = yaml.safe_load(f) or {}
	if not isinstance(rec, dict):
		raise ValueError(f'{path} is not a recipe')
	return rec


def default_recipe():
	"""The recipe of the repository (gee_vars.yaml)."""
	return _normal(_read_yaml(RECIPE))


def last_recipe():
	"""The last recipe used in the app, or that of the repository."""
	if os.path.isfile(paths.APP_RECIPE):
		try:
			return {'recipe': _normal(_read_yaml(paths.APP_RECIPE)), 'source': paths.rel(paths.APP_RECIPE)}
		except (OSError, ValueError, yaml.YAMLError):
			pass
	return {'recipe': default_recipe(), 'source': paths.rel(RECIPE)}


def load_recipe(path):
	try:
		return {'recipe': _normal(_read_yaml(path)), 'source': path}
	except yaml.YAMLError as e:
		raise ValueError(f'{os.path.basename(path)} is not valid YAML: {e}') from e


def _normal(rec):
	"""The recipe with the keys that the form shows."""
	rec = copy.deepcopy(rec)
	rec.setdefault('folder', 'gee')
	rec.setdefault('crs', 'auto')
	rec.setdefault('study', {})
	rec['study'].setdefault('resolution', 10)
	rec['study'].setdefault('buffer', None)
	rec['study'].setdefault('drive', False)
	rec.setdefault('ecoregion', {'resolution': None})
	if rec['ecoregion'] is None:
		rec['ecoregion'] = {'resolution': None}
	rec['ecoregion'].setdefault('drive', False)
	rec['ecoregion'].setdefault('crs', None)
	rec.setdefault('period', ['2018-01-01', '2021-12-31'])
	rec.setdefault('season', None)
	rec.setdefault('quality', {})
	q = rec['quality']
	for k in ('clouds', 'water', 'snow', 'aerosols'):
		q.setdefault(k, True)
	q.setdefault('cloud_score', 0.6)
	q.setdefault('scene_clouds', None)
	rec.setdefault('exclude_landcover', [])
	rec.setdefault('drive_folder', 'eHabitat_GEE')
	rec.setdefault('drive_dir', None)
	rec.setdefault('variables', [])
	for k in ('period',):
		rec[k] = [str(x) for x in rec[k]] if rec[k] else rec[k]
	return rec


def clean_recipe(rec):
	"""The recipe as get_vars_gee.py reads it: without empty optional keys."""
	rec = copy.deepcopy(rec)
	for k in ('season', 'drive_dir', 'exclude_landcover'):
		if rec.get(k) in (None, '', []):
			rec.pop(k, None)
	if rec.get('ecoregion') and not rec['ecoregion'].get('resolution'):
		rec['ecoregion'] = {'resolution': None}
	elif rec.get('ecoregion') and not rec['ecoregion'].get('crs'):
		rec['ecoregion'].pop('crs', None)
	q = rec.get('quality') or {}
	for k in ('cloud_score', 'scene_clouds'):
		if q.get(k) in ('', None) and k == 'scene_clouds':
			q.pop(k, None)
	vars_ = []
	for v in rec.get('variables') or []:
		v = {k: x for k, x in v.items() if x not in (None, '', [])}
		if isinstance(v.get('metric'), list) and len(v['metric']) == 1:
			v['metric'] = v['metric'][0]
		vars_.append(v)
	rec['variables'] = vars_
	if not vars_:
		raise ValueError('the recipe has no variables')
	if not rec.get('folder') or not str(rec['folder']).strip():
		raise ValueError('the destination folder is missing')
	return rec


def save_recipe(rec, path):
	rec = clean_recipe(rec)
	with open(path, 'w', encoding='utf-8') as f:
		f.write('# Recipe of variables of MonitorEO for eHabitat+ (get_vars_gee.py), saved by the app\n')
		yaml.safe_dump(rec, f, sort_keys=False, allow_unicode=True, default_flow_style=None)
	return path


def run(rec, opts):
	"""Plan (dry_run) or download the recipe REC with get_vars_gee.py. OPTS: dry_run, only, overwrite, use, no_wait,
	drive_dir, no_update."""
	if jobs.running('gee'):
		raise ValueError('a download is already running')
	sa = model.study_area()
	if not sa:
		raise ValueError('choose the study area first in "Input variables"')
	s = paths.load_settings()
	account = s.get('ee_account') or 'shared'
	project = (s.get('ee_project') or '').strip() if account == 'own' else ''
	rec = clean_recipe(rec)
	rec['project'] = project
	if account == 'shared':
		for g in ('study', 'ecoregion'):
			if (rec.get(g) or {}).get('drive') and not opts.get('dry_run'):
				raise ValueError('the shared account of MonitorEO cannot export to Google Drive: untick "Export through '
				                 'Google Drive" or use your own Earth Engine account')
	save_recipe(rec, paths.APP_RECIPE)
	folder = str(rec['folder']).strip('/\\')
	cmd = [sys.executable, 'get_vars_gee.py', '--recipe', paths.APP_RECIPE, '--folder', folder,
	       '--study-area', sa['file']]
	if sa.get('layer'):
		cmd += ['--layer', sa['layer']]
	if project:
		cmd += ['--project', project]
	if opts.get('dry_run'):
		cmd.append('--dry-run')
	if opts.get('only') in ('study', 'ecoregion'):
		cmd += ['--only', opts['only']]
	if opts.get('overwrite'):
		cmd.append('--overwrite')
	if opts.get('use') and not opts.get('dry_run'):
		cmd.append('--use')
	if opts.get('no_wait'):
		cmd.append('--no-wait')
	if opts.get('no_update'):
		cmd.append('--no-update')
	env = dict(os.environ)
	if account == 'shared':
		env.pop('EE_PROJECT', None)
		env.pop('GOOGLE_CLOUD_PROJECT', None)
	title = 'Plan of the download from MonitorEO' if opts.get('dry_run') else 'Download from MonitorEO'

	def after(job):
		if not opts.get('dry_run'):
			dest = os.path.join(paths.INVARS, folder)
			if os.path.isdir(dest):
				shutil.copy2(paths.APP_RECIPE, os.path.join(dest, 'gee_vars.yaml'))
		return {'folder': folder}

	intro = ('Shared account of MonitorEO' if account == 'shared' else f'Own Earth Engine account, project {project}')
	return jobs.start_commands('gee', title, [cmd], env=env, after=after, intro=intro).id
