"""Functions of the app called by the interface (window.pywebview.api, or POST <token>/api/<function> in a web
browser). Each one returns {'ok': True, 'result': ...} or {'ok': False, 'error': message}."""

import functools
import os
import platform
import subprocess
import sys
import traceback
import webbrowser

from . import __version__, envsetup, gee, inputs, jobs, model, paths, results

import hft_state


def exposed(fn):
	@functools.wraps(fn)
	def wrapper(self, *args, **kwargs):
		try:
			return {'ok': True, 'result': fn(self, *args, **kwargs)}
		except SystemExit as e:
			return {'ok': False, 'error': str(e.code)}
		except Exception as e:
			if not isinstance(e, (ValueError, FileNotFoundError, PermissionError)):
				traceback.print_exc()
			return {'ok': False, 'error': str(e) or type(e).__name__}
	return wrapper


FILE_TYPES = {
	'tif': ('GeoTIFF (*.tif;*.tiff)', 'All files (*.*)'),
	'vector': ('Vector files (*.gpkg;*.shp;*.geojson;*.json)', 'All files (*.*)'),
	'yaml': ('YAML recipe (*.yaml;*.yml)', 'All files (*.*)'),
}


class Api:
	def __init__(self):
		self._window = None
		self._activation = None
		self._tools = None

	# General ------------------------------------------------------------------------------------------------------
	@exposed
	def info(self):
		cfg = model.read_config()
		return {'version': __version__, 'root': paths.ROOT, 'segm': paths.SEGM, 'platform': platform.system(),
		        'python': sys.executable, 'activation': self._activation, 'window': self._window is not None,
		        'study_area': cfg.get('STUDY_AREA'), 'invars_folder': hft_state.invars_folder(cfg)}

	@exposed
	def tools(self):
		if self._tools is None:
			self._tools = envsetup.tools()
		return self._tools

	@exposed
	def job(self, job_id, since=0):
		j = jobs.get(job_id)
		if not j:
			raise ValueError(f'no job {job_id}')
		return j.snapshot(since)

	@exposed
	def cancel(self, job_id):
		return jobs.cancel(job_id)

	@exposed
	def running(self):
		return [{'id': j.id, 'kind': j.kind, 'title': j.title} for j in jobs.running()]

	@exposed
	def open_folder(self, path):
		path = path if os.path.isabs(path) else os.path.join(paths.ROOT, path)
		if not os.path.exists(path):
			raise ValueError(f'{path} not found')
		if sys.platform == 'darwin':
			subprocess.Popen(['open', '-R', path] if os.path.isfile(path) else ['open', path])
		elif os.name == 'nt':
			if os.path.isfile(path):
				subprocess.Popen(['explorer', '/select,', os.path.normpath(path)])
			else:
				os.startfile(path)  # noqa
		else:
			subprocess.Popen(['xdg-open', path if os.path.isdir(path) else os.path.dirname(path)])
		return True

	@exposed
	def open_url(self, url):
		if not str(url).startswith(('https://', 'http://')):
			raise ValueError('only web links')
		webbrowser.open(url)
		return True

	# Dialogs ------------------------------------------------------------------------------------------------------
	def _dialog_type(self, name):
		import webview
		fd = getattr(webview, 'FileDialog', None)
		if fd is not None:
			return getattr(fd, name)
		return getattr(webview, {'OPEN': 'OPEN_DIALOG', 'FOLDER': 'FOLDER_DIALOG', 'SAVE': 'SAVE_DIALOG'}[name])

	@exposed
	def pick_files(self, kind='tif', multiple=True, start=None):
		if self._window is None:
			return None  # the interface asks for the path
		got = self._window.create_file_dialog(self._dialog_type('OPEN'), directory=start or '', allow_multiple=bool(multiple),
		                                      file_types=FILE_TYPES.get(kind, ()))
		return list(got) if got else []

	@exposed
	def pick_folder(self, start=None):
		if self._window is None:
			return None
		got = self._window.create_file_dialog(self._dialog_type('FOLDER'), directory=start or '')
		return got[0] if got else ''

	@exposed
	def pick_save(self, kind='yaml', name='recipe.yaml', start=None):
		if self._window is None:
			return None
		got = self._window.create_file_dialog(self._dialog_type('SAVE'), directory=start or paths.SEGM, save_filename=name,
		                                      file_types=FILE_TYPES.get(kind, ()))
		if not got:
			return ''
		return got if isinstance(got, str) else got[0]

	# Inputs -------------------------------------------------------------------------------------------------------
	@exposed
	def sets(self):
		cfg = model.read_config()
		return {'sets': inputs.list_sets(), 'active': hft_state.invars_folder(cfg)}

	@exposed
	def get_set(self, folder):
		return inputs.get_set(folder)

	@exposed
	def create_set(self, name):
		return inputs.create_set(name)

	@exposed
	def add_rasters(self, folder, scale, files):
		if scale not in ('study', 'ecoregion'):
			raise ValueError('unknown scale')
		return jobs.start_function('inputs', 'Add variables', inputs.add_files, folder, scale, list(files)).id

	@exposed
	def set_included(self, folder, scale, file, included):
		return inputs.set_included(folder, scale, file, bool(included))

	@exposed
	def use_set(self, folder):
		inputs.set_dir(folder)
		hft_state.write_option('INVARS_FOLDER', repr((folder or '').strip('/\\')))
		return True

	@exposed
	def raster_preview(self, path):
		if not paths.inside(path, paths.INVARS):
			raise ValueError('only the variables of inVars')
		return inputs.preview(path)

	@exposed
	def study_area(self):
		sa = model.study_area()
		if not sa:
			return None
		info = inputs.vector_info(sa['file'], sa.get('layer'))
		try:
			all_ids = inputs.vector_ids(sa['file'], info['layer'], sa.get('col_id') or 'cat')
		except ValueError as e:
			all_ids, info['id_error'] = [], str(e)
		return dict(sa, info=info, all_ids=all_ids, ids=sa.get('ids') or all_ids)

	@exposed
	def inspect_vector(self, path, layer=None, col_id='cat'):
		info = inputs.vector_info(path, layer)
		try:
			ids = inputs.vector_ids(path, info['layer'], col_id or 'cat')
			err = None
		except ValueError as e:
			ids, err = [], str(e)
		return {'info': info, 'ids': ids, 'id_error': err}

	@exposed
	def set_study_area(self, file, layer=None, col_id='cat', ids=None):
		info = inputs.vector_info(file, layer)
		all_ids = inputs.vector_ids(file, info['layer'], col_id or 'cat')
		ids = sorted({int(i) for i in (ids or all_ids)})
		if not set(ids) <= set(all_ids):
			raise ValueError('some IDs are not in the study area')
		target = inputs.copy_to_pa(file)
		paths.save_settings(study_area={'file': target, 'layer': info['layer'], 'col_id': col_id or 'cat',
		                                'ids': ids if len(ids) < len(all_ids) else None})
		model.write_palist(ids)
		return self.study_area()['result']

	@exposed
	def study_geojson(self):
		sa = model.study_area()
		if not sa:
			return None
		return inputs.vector_geojson(sa['file'], sa.get('layer'))

	@exposed
	def vector_preview(self, path, layer=None):
		"""Polygons of a vector file chosen by the user, before using it (any folder)."""
		if not str(path).lower().endswith(inputs.VECTOR_EXT):
			raise ValueError('not a vector file')
		return inputs.vector_geojson(path, layer)

	@exposed
	def vector_geojson(self, path):
		if not (paths.inside(path, paths.INVARS) or paths.inside(path, paths.PA_DIR)):
			raise ValueError('only files of inVars or segm_optim/pa')
		return inputs.vector_geojson(path)

	# HFT model ----------------------------------------------------------------------------------------------------
	@exposed
	def params(self):
		return model.params()

	@exposed
	def save_params(self, values, nproc=None):
		return model.save_params(values, nproc)

	@exposed
	def status(self):
		return model.status()

	@exposed
	def scenario_check(self, name):
		return model.scenario_check(name)

	@exposed
	def run_model(self, values, nproc=None, first=None):
		return model.run(values, nproc, first or None)

	# MonitorEO ----------------------------------------------------------------------------------------------------
	@exposed
	def gee_status(self):
		return gee.status()

	@exposed
	def gee_catalog(self):
		return gee.catalog()

	@exposed
	def gee_describe(self, ref):
		return gee.describe(ref)

	@exposed
	def gee_connect(self, account, project=None):
		return gee.connect(account, project)

	@exposed
	def gee_settings(self, account, project=None):
		"""The account of Earth Engine that the downloads use ('shared' or 'own' and its project)."""
		if account not in ('shared', 'own'):
			raise ValueError('unknown account')
		paths.save_settings(ee_account=account, ee_project=(project or '').strip() or None)
		return True

	@exposed
	def gee_auth(self):
		return gee.authenticate()

	@exposed
	def gee_update(self):
		return gee.update()

	@exposed
	def gee_recipe(self):
		return gee.last_recipe()

	@exposed
	def gee_default_recipe(self):
		return {'recipe': gee.default_recipe(), 'source': paths.rel(gee.RECIPE)}

	@exposed
	def gee_load_recipe(self, path):
		return gee.load_recipe(path)

	@exposed
	def gee_save_recipe(self, recipe, path):
		return gee.save_recipe(recipe, path)

	@exposed
	def gee_run(self, recipe, opts):
		return gee.run(recipe, opts or {})

	# Results ------------------------------------------------------------------------------------------------------
	@exposed
	def results_list(self):
		return results.list_results()

	@exposed
	def results(self, study_area, scenario):
		r = results.scenario(study_area, scenario)
		for rec in r['per_id'].values():
			for c in rec['cuts']:
				c['shp'] = paths.rel(c['shp'])
				c['colours'] = paths.rel(c['colours']) if c['colours'] else None
		return r

	@exposed
	def hft_layer(self, study_area, scenario, id_, cut):
		r = results.scenario(study_area, scenario)
		rec = r['per_id'].get(str(id_))
		if not rec:
			raise ValueError(f'no results of the ID {id_}')
		c = next((c for c in rec['cuts'] if c['cut'] == cut), None)
		if not c:
			raise ValueError(f'no results of the cut {cut}')
		return results.hft_layer(c['shp'], c['colours'], r.get('metadata'))

	@exposed
	def study_outline(self, study_area):
		return results.study_outline(study_area)
