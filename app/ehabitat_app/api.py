"""Functions of the app called by the interface (window.pywebview.api, or POST <token>/api/<function> in a web
browser). Each one returns {'ok': True, 'result': ...} or {'ok': False, 'error': message}."""

import functools
import os
import platform
import subprocess
import sys
import threading
import time
import traceback
import webbrowser

from . import __version__, ecosim, envsetup, gee, inputs, jobs, model, paths, results, updater

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

	# Updates ------------------------------------------------------------------------------------------------------
	@exposed
	def update_status(self):
		return updater.status()

	@exposed
	def check_updates(self):
		"""Check GitHub and update now (a job); with a new environment.yml, the environment is updated after it."""
		def work(job):
			res = updater.update(job.say)
			job.say(res['message'])
			if res.get('environment'):
				job.say('environment.yml changed: updating the environment in the background')
				updater.start_environment_update()
			return res
		return jobs.start_function('update', 'Update of eHabitat+', work).id

	@exposed
	def set_auto_update(self, on):
		paths.save_settings(auto_update=bool(on))
		return bool(on)

	@exposed
	def restart(self):
		"""Open the app again (with its new version) and close this one."""
		if jobs.running():
			raise ValueError('there are jobs running: wait until they finish or stop them')
		updater.spawn()

		def close():
			time.sleep(0.5)
			if self._window is not None:
				self._window.destroy()
			else:
				os._exit(0)
		threading.Thread(target=close, daemon=True).start()
		return True

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
			folder = path if os.path.isdir(path) else os.path.dirname(path)
			if not (paths.wsl() and paths.open_in_windows(folder)):  # WSL: the File Explorer of Windows
				subprocess.Popen(['xdg-open', folder])
		return True

	@exposed
	def open_url(self, url):
		if not str(url).startswith(('https://', 'http://')):
			raise ValueError('only web links')
		if not (paths.wsl() and paths.open_in_windows(url)):  # WSL: the web browser of Windows
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
		d = inputs.set_dir(folder)
		if inputs.area_only(d):
			raise ValueError(f'inVars/{folder} has only the variables of an ecoregion or comparison area (no GeoTIFFs of '
			                 'a study area): it cannot be the input of the HFT model; compare the HFTs with it in '
			                 '"Similarity"')
		hft_state.write_option('INVARS_FOLDER', repr((folder or '').strip('/\\')))
		return True

	@exposed
	def raster_preview(self, path):
		if not paths.inside(path, paths.INVARS):
			raise ValueError('only the variables of inVars')
		return inputs.preview(path)

	@exposed
	def study_area(self):
		"""The study areas chosen: their file, layer, ID and name columns, the areas of the file ({id, label, ha}) and
		the IDs to process (palist.csv)."""
		sa = model.study_area()
		if not sa:
			return None
		info = inputs.vector_info(sa['file'], sa.get('layer'))
		label_col = sa.get('label_col', inputs.guess_label_col(info['fields']))
		try:
			areas = inputs.vector_areas(sa['file'], info['layer'], sa.get('col_id') or 'cat', label_col)
		except ValueError as e:
			areas, info['id_error'] = [], str(e)
		all_ids = [a['id'] for a in areas]
		ids = [i for i in (sa.get('ids') or all_ids) if i in set(all_ids)] or all_ids
		return dict(sa, info=info, label_col=label_col, areas=areas, all_ids=all_ids, ids=ids)

	@exposed
	def inspect_vector(self, path, layer=None, col_id='cat', label_col=''):
		"""A vector file before using it: its layers, fields and study areas. LABEL_COL '' guesses the column of the
		names, None leaves them out."""
		info = inputs.vector_info(path, layer)
		if label_col == '':
			label_col = inputs.guess_label_col(info['fields'])
		try:
			areas = inputs.vector_areas(path, info['layer'], col_id or 'cat', label_col)
			err = None
		except ValueError as e:
			areas, err = [], str(e)
		return {'info': info, 'areas': areas, 'ids': [a['id'] for a in areas], 'label_col': label_col, 'id_error': err}

	@exposed
	def set_study_area(self, file, layer=None, col_id='cat', label_col=None):
		"""Use the study areas of FILE. The areas to process are kept if the file, layer and ID column are the same,
		otherwise all of them (chosen then in HFT mapping)."""
		info = inputs.vector_info(file, layer)
		areas = inputs.vector_areas(file, info['layer'], col_id or 'cat', label_col)
		all_ids = [a['id'] for a in areas]
		if not all_ids:
			raise ValueError('the file has no study areas (no IDs)')
		target = inputs.copy_to_pa(file)
		old = model.study_area() or {}
		same = old.get('file') == target and (old.get('layer') or None) == (info['layer'] or None) \
			and (old.get('col_id') or 'cat') == (col_id or 'cat')
		ids = [i for i in (old.get('ids') or []) if i in set(all_ids)] if same else []
		ids = ids or all_ids
		paths.save_settings(study_area={'file': target, 'layer': info['layer'], 'col_id': col_id or 'cat',
		                                'label_col': label_col, 'ids': ids if len(ids) < len(all_ids) else None})
		model.write_palist(ids)
		return self.study_area()['result']

	@exposed
	def set_ids(self, ids):
		"""The study areas (IDs) to process: palist.csv."""
		if jobs.running('model'):
			raise ValueError('the model is running: change the areas when it finishes')
		sa = model.study_area()
		if not sa:
			raise ValueError('choose the study areas first in "Input variables"')
		all_ids = inputs.vector_ids(sa['file'], sa.get('layer'), sa.get('col_id') or 'cat')
		ids = sorted({int(i) for i in ids})
		if not ids:
			raise ValueError('choose at least one study area')
		bad = [i for i in ids if i not in set(all_ids)]
		if bad:
			raise ValueError(f"IDs not in the study areas: {', '.join(map(str, bad))}")
		paths.save_settings(study_area=dict(sa, ids=ids if len(ids) < len(all_ids) else None))
		model.write_palist(ids)
		return ids

	@exposed
	def study_geojson(self):
		sa = model.study_area()
		if not sa:
			return None
		info = inputs.vector_info(sa['file'], sa.get('layer'))
		label_col = sa.get('label_col', inputs.guess_label_col(info['fields']))
		return inputs.vector_geojson(sa['file'], sa.get('layer'), col_id=sa.get('col_id') or 'cat', label_col=label_col)

	@exposed
	def vector_preview(self, path, layer=None, col_id=None, label_col=None):
		"""Polygons of a vector file chosen by the user, before using it (any folder)."""
		if not str(path).lower().endswith(inputs.VECTOR_EXT):
			raise ValueError('not a vector file')
		return inputs.vector_geojson(path, layer, col_id=col_id, label_col=label_col)

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
	def gee_community(self):
		return gee.community()

	@exposed
	def gee_describe_asset(self, asset, account=None, project=None):
		return gee.describe_asset(asset, account, project)

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
	def gee_auto_update(self):
		"""The update of monitoreo when the app opens (once, if the app updates itself): the id of its job or None."""
		return gee.auto_update()

	@exposed
	def gee_area(self, path, layer=None):
		"""A vector file chosen as the comparison area of the ecoregion grid (any folder)."""
		return gee.area(path, layer)

	@exposed
	def gee_area_check(self, folder, comparison=None, layer=None, only=None):
		return gee.area_check(folder, comparison, layer, only)

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

	# Similarity ---------------------------------------------------------------------------------------------------
	@exposed
	def ecosim_options(self):
		return ecosim.options()

	@exposed
	def ecosim_save(self, values):
		return ecosim.save(values)

	@exposed
	def ecosim_sources(self, study_area, scenario):
		return ecosim.sources(study_area, scenario)

	@exposed
	def ecosim_plan(self, study_area, scenario, cut=None, res=None, min_cover=None, ecoregions=None, area=None):
		return ecosim.plan(study_area, scenario, cut, res, min_cover, ecoregions, area)

	@exposed
	def ecosim_areas(self, study_area, scenario):
		return ecosim.areas(study_area, scenario)

	@exposed
	def ecosim_run(self, values, study_area, scenario, cut=''):
		return ecosim.run(values, study_area, scenario, cut)

	@exposed
	def ecosim_runs(self, study_area, scenario):
		return ecosim.runs(study_area, scenario)

	@exposed
	def ecosim_results(self, folder):
		return ecosim.run_results(folder)

	@exposed
	def ecosim_preview(self, tif, kind='sim', bbox=None):
		return ecosim.preview(tif, kind, bbox)

	@exposed
	def ecosim_outline(self, path, field=None, ids=None, folder=None):
		return ecosim.outline(path, field, ids, folder)
