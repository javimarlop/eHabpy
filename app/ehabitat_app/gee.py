"""MonitorEO: the catalog of variables, the connection with Google Earth Engine (the shared account of MonitorEO or the
user's own) and the download of a recipe of variables over the study area and its ecoregions (or a comparison area of
the user) with get_vars_gee.py, which keeps the package monitoreo at its last version of GitHub and computes the
variables with it. The app updates monitoreo in the background when it opens (as it updates itself) and its catalog
is read again when the installed package changes, without opening the app again."""

import copy
import importlib
import importlib.machinery
import importlib.metadata
import json
import os
import re
import shutil
import sys
import threading
import urllib.parse
import urllib.request

import yaml

from . import inputs, jobs, model, paths

RECIPE = os.path.join(paths.SEGM, 'gee_vars.yaml')  # the recipe of the repository (with its comments)
CREDENTIALS = os.path.join('~', '.config', 'earthengine', 'credentials')
SLUG_RE = re.compile(r'[^A-Za-z0-9_-]+')

SENSORS = ['Sentinel-2', 'Landsat', 'Modis']
METRICS = ['mean', 'median', 'min', 'max', 'std', 'sum', 'count', 'mode', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95']
TEMPORAL = [['period', 'Whole period'], ['annual', 'Each year, then between years'],
            ['monthly', 'Each month, between years and mean of the months'], ['16-day', 'Each 16 days, between years and mean']]
WORLDCOVER = [[10, 'Tree cover'], [20, 'Shrubland'], [30, 'Grassland'], [40, 'Cropland'], [50, 'Built-up'],
              [60, 'Bare / sparse vegetation'], [70, 'Snow and ice'], [80, 'Permanent water bodies'],
              [90, 'Herbaceous wetland'], [95, 'Mangroves'], [100, 'Moss and lichen']]


# The package monitoreo, read again when it changes ------------------------------------------------------------------

_lock = threading.RLock()
_loaded = {'sig': None}  # what identified the package when it was imported (see _signature)
_auto = {'started': False}  # the update of monitoreo when the app opens: only once


def installed():
	"""(version, commit of GitHub or None, folder of an editable copy or None) of the installed package monitoreo, or
	None if it is not installed (as get_vars_gee.installed_monitoreo)."""
	importlib.invalidate_caches()
	try:
		dist = importlib.metadata.distribution('monitoreo')
	except importlib.metadata.PackageNotFoundError:
		return None
	try:
		url = json.loads(dist.read_text('direct_url.json') or '{}')
	except ValueError:
		url = {}
	commit = (url.get('vcs_info') or {}).get('commit_id')
	m = re.search(r'/archive/([0-9a-f]{40})\.zip$', url.get('url', ''))
	if m:
		commit = m.group(1)
	folder = None
	if (url.get('dir_info') or {}).get('editable') and url.get('url', '').startswith('file:'):
		folder = urllib.request.url2pathname(urllib.parse.urlparse(url['url']).path)
	return dist.version, commit, folder


def version():
	got = installed()
	return got[0] if got else None


def _package_dir():
	"""The folder of the package monitoreo that an import would load now (not the one already imported)."""
	spec = importlib.machinery.PathFinder.find_spec('monitoreo')
	if spec and spec.submodule_search_locations:
		return list(spec.submodule_search_locations)[0]
	m = sys.modules.get('monitoreo')
	return os.path.dirname(m.__file__) if m is not None and getattr(m, '__file__', None) else None


def _signature():
	"""What identifies the installed package monitoreo: its folder, version and commit, and the size and time of
	its catalog (catalog.json and catalog.py); None when it is not installed."""
	got = installed()
	if not got:
		return None
	folder = _package_dir()
	stats = []
	for name in ('catalog.json', 'catalog.py', 'services.py'):
		try:
			st = os.stat(os.path.join(folder or '', name))
			stats.append((name, st.st_mtime_ns, st.st_size))
		except OSError:
			stats.append((name, None, None))
	return (folder, got[0], got[1], tuple(stats))


def _refresh():
	"""Import monitoreo again when the installed package is not the one imported (updated by the job 'gee-update' or
	by a download), keeping the connection to Earth Engine; returns whether it was imported again."""
	with _lock:
		sig = _signature()
		old = _loaded['sig']
		_loaded['sig'] = sig
		if old is None or sig == old or 'monitoreo' not in sys.modules:
			return False
		svc = sys.modules.get('monitoreo.services')
		keep = {k: getattr(svc, k) for k in ('_project', '_shared', '_user') if svc is not None and hasattr(svc, k)}
		for name in [n for n in sys.modules if n == 'monitoreo' or n.startswith('monitoreo.')]:
			del sys.modules[name]
		importlib.invalidate_caches()
		try:
			from monitoreo import services
		except ImportError:
			return True
		for k, v in keep.items():
			if hasattr(services, k):
				setattr(services, k, v)
		return True


def _monitoreo():
	_refresh()
	try:
		with _lock:
			import monitoreo  # noqa: F401
			from monitoreo import catalog, services
			if _loaded['sig'] is None:
				_loaded['sig'] = _signature()
		return catalog, services
	except ImportError as e:
		raise ValueError(f'the package monitoreo is not installed ({e}): press "Update MonitorEO"') from e


def status():
	got = installed()
	commit = got[1] if got else None
	if got and not commit and got[2]:  # an editable copy (development): the commit of its git clone
		import subprocess
		try:
			r = subprocess.run(['git', '-C', got[2], 'rev-parse', 'HEAD'], capture_output=True, text=True, timeout=10,
			                   stdin=subprocess.DEVNULL, creationflags=jobs.NO_WINDOW)
			commit = r.stdout.strip() if r.returncode == 0 else None
		except (OSError, subprocess.SubprocessError):
			pass
	out = {'version': got[0] if got else None, 'commit': commit, 'editable': got[2] if got else None,
	       'credentials': os.path.isfile(os.path.expanduser(CREDENTIALS)),
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
	        'worldcover': WORLDCOVER, 'monitoreo': _label(installed())}


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


def _label(got):
	"""'0.3.33 (b0341eb)' of installed()."""
	if not got:
		return 'not installed'
	return got[0] + (f' ({got[1][:7]})' if got[1] else ' (development copy)' if got[2] else '')


def update(auto=False):
	"""Install the last version of monitoreo of GitHub (as get_vars_gee.py does on every run), in a job. AUTO: the
	update when the app opens. The catalog of the app is read again from the new version (see _refresh)."""
	if jobs.running('gee'):
		raise ValueError('a download is running: update MonitorEO when it finishes')
	if jobs.running('gee-update'):
		raise ValueError('MonitorEO is already being updated')
	code = 'import get_vars_gee as g; g.update_monitoreo(True)'
	before = installed()

	def after(job):
		now = installed()
		changed = bool(now and (not before or now[:2] != before[:2]))
		if changed:
			_refresh()
			job.say(f'\nmonitoreo {_label(before)} -> {_label(now)}: the catalog of the app uses the new version now.')
		return {'before': _label(before), 'now': _label(now), 'version': now[0] if now else None,
		        'commit': now[1] if now else None, 'changed': changed}

	title = 'Automatic update of MonitorEO' if auto else 'Update MonitorEO'
	return jobs.start_commands('gee-update', title, [[sys.executable, '-c', code]], after=after).id


def auto_update():
	"""When the app opens (called by the interface once its window shows): the update of monitoreo in the
	background, once, if the app updates itself (auto_update in its settings and not EHAB_NO_UPDATE). Returns the
	id of the job, or None when it is not started."""
	from . import updater
	with _lock:
		if _auto['started'] or not updater.enabled():
			return None
		if any(j.kind.startswith('gee') for j in jobs.running()):
			return None
		_auto['started'] = True
	try:
		return update(auto=True)
	except ValueError:
		return None


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
	# the area of the ecoregion grid: None = the RESOLVE ecoregions that the study area touches, or a vector file
	rec.setdefault('comparison_area', None)
	rec.setdefault('comparison_area_layer', None)
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
	# a comparison area only with the ecoregion grid (get_vars_gee.py needs its pixel size)
	area = str(rec.get('comparison_area') or '').strip()
	if not area or not (rec.get('ecoregion') or {}).get('resolution'):
		rec.pop('comparison_area', None)
		rec.pop('comparison_area_layer', None)
	else:
		rec['comparison_area'] = area
		if not str(rec.get('comparison_area_layer') or '').strip():
			rec.pop('comparison_area_layer', None)
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


def area_path(path):
	"""Absolute path of the vector file of a comparison area (a relative one is relative to segm_optim, as in the
	recipes of get_vars_gee.py)."""
	path = os.path.expanduser(str(path or '').strip())
	if not path:
		return None
	return os.path.normpath(path if os.path.isabs(path) else os.path.join(paths.SEGM, path))


def area_name(path, layer=None):
	"""Name of a comparison area, as get_vars_gee.py names its polygons (ECO_NAME): the name field of its only
	feature (name, NAME, nombre or Name), else the name of the file."""
	stem = os.path.splitext(os.path.basename(path))[0]
	try:
		from osgeo import ogr
		ds = ogr.Open(path)
		lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
		if lyr is not None and lyr.GetFeatureCount() == 1:
			defn = lyr.GetLayerDefn()
			names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
			field = next((n for n in ('name', 'NAME', 'nombre', 'Name') if n in names), None)
			f = lyr.GetNextFeature()
			if field and f is not None and f.IsFieldSet(field) and str(f.GetField(field)).strip():
				return str(f.GetField(field)).strip()
	except Exception:  # unreadable here: vector_info says why
		pass
	return stem


def area(path, layer=None):
	"""A vector file chosen as the comparison area: its layers, features, CRS, extent and name."""
	p = area_path(path)
	if not p or not os.path.isfile(p):
		raise ValueError(f'{path} not found')
	if not p.lower().endswith(inputs.VECTOR_EXT):
		raise ValueError(f'{os.path.basename(p)} is not a vector file (.gpkg, .shp, .geojson)')
	info = inputs.vector_info(p, layer or None)
	info['name'] = area_name(p, info['layer'])
	return info


def _same_file(a, b):
	try:
		return bool(a and b) and os.path.samefile(a, b)
	except OSError:
		return os.path.normpath(str(a)) == os.path.normpath(str(b))


def area_check(folder, comparison=None, layer=None, only=None):
	"""What a download into inVars/FOLDER does to the ecoregion grid that the folder already has: whether it replaces
	the grid of another area (the RESOLVE ecoregions or another comparison area) with COMPARISON (a vector file; None
	= the RESOLVE ecoregions), and another folder to keep both. ONLY 'study' does not download the ecoregion grid."""
	folder = str(folder or '').strip().strip('/\\')
	out = {'folder': folder, 'existing': None, 'replaces': False, 'suggestion': None, 'new': None}
	if not folder or not re.fullmatch(r'[A-Za-z0-9_-]+', folder):
		return out
	eco_dir = os.path.join(paths.INVARS, folder, paths.ECO_FOLDER)
	comparison = area_path(comparison) if comparison else None
	out['new'] = {'kind': 'comparison area', 'file': comparison, 'layer': layer or None,
	              'name': area_name(comparison, layer) if comparison and os.path.isfile(comparison) else None} \
		if comparison else {'kind': 'ecoregions', 'name': 'the RESOLVE ecoregions that the study area touches'}
	if not os.path.isdir(eco_dir) or not inputs.tifs(eco_dir):
		return out
	old = inputs.ecoregion_area(eco_dir)
	out['existing'] = old
	if only == 'study':
		return out
	if comparison:
		same = old['kind'] == 'comparison area' and _same_file(area_path(old.get('file')), comparison) \
			and (old.get('layer') or None) == (layer or None)
		out['replaces'] = not same
	else:
		out['replaces'] = old['kind'] == 'comparison area'
	if out['replaces']:
		stem = os.path.splitext(os.path.basename(comparison))[0] if comparison else 'resolve'
		base = f"{folder}_{SLUG_RE.sub('_', stem).strip('_')[:40] or 'area'}"
		name, n = base, 2
		while os.path.exists(os.path.join(paths.INVARS, name)):
			name, n = f'{base}_{n}', n + 1
		out['suggestion'] = name
	return out


def run(rec, opts):
	"""Plan (dry_run) or download the recipe REC with get_vars_gee.py. OPTS: dry_run, only, overwrite, use, no_wait,
	drive_dir, no_update."""
	if jobs.running('gee'):
		raise ValueError('a download is already running')
	if jobs.running('gee-update'):
		raise ValueError('MonitorEO is being updated: wait until it finishes (a few seconds)')
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
	comparison = area_path(rec.get('comparison_area'))
	if comparison and not os.path.isfile(comparison):
		raise ValueError(f'the file of the comparison area does not exist: {comparison}')
	if opts.get('only') == 'ecoregion' and not (rec.get('ecoregion') or {}).get('resolution'):
		raise ValueError('tick "Also download the variables of the ecoregion or comparison area" (with its pixel size) '
		                 'to download only that grid')
	save_recipe(rec, paths.APP_RECIPE)
	folder = str(rec['folder']).strip('/\\')
	cmd = [sys.executable, 'get_vars_gee.py', '--recipe', paths.APP_RECIPE, '--folder', folder,
	       '--study-area', sa['file']]
	if sa.get('layer'):
		cmd += ['--layer', sa['layer']]
	if project:
		cmd += ['--project', project]
	if comparison:
		cmd += ['--area', comparison]
		if rec.get('comparison_area_layer'):
			cmd += ['--area-layer', str(rec['comparison_area_layer'])]
	if opts.get('dry_run'):
		cmd.append('--dry-run')
	if opts.get('only') in ('study', 'ecoregion'):
		cmd += ['--only', opts['only']]
	if opts.get('overwrite'):
		cmd.append('--overwrite')
	# only the ecoregion or comparison area: not a set of variables for the HFT model (no GeoTIFFs of the study area)
	if opts.get('use') and not opts.get('dry_run') and opts.get('only') != 'ecoregion':
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
	before = installed()

	def after(job):
		if not opts.get('dry_run'):
			dest = os.path.join(paths.INVARS, folder)
			if os.path.isdir(dest):
				shutil.copy2(paths.APP_RECIPE, os.path.join(dest, 'gee_vars.yaml'))
		now = installed()  # get_vars_gee.py updates monitoreo first: the catalog of the app then uses it too
		changed = bool(now and (not before or now[:2] != before[:2]))
		if changed:
			_refresh()
		return {'folder': folder, 'only': opts.get('only') or None,
		        'monitoreo': {'before': _label(before), 'now': _label(now), 'changed': changed}}

	intro = ('Shared account of MonitorEO' if account == 'shared' else f'Own Earth Engine account, project {project}')
	return jobs.start_commands('gee', title, [cmd], env=env, after=after, intro=intro).id
