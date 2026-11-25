#### eHabitat+: test of the parameters files (parameters.py, run_hft.py, the app)
#### License: CC BY-SA 3.0
"""
Checks that the parameters of one study area (config.py, the recipe of its set of variables and the order of its
variables) are saved without what belongs to the study area or the computer, and that loading them into another one
changes only the options of the method: the study area, its variables, the scenario and the comparison area stay; the
recipe becomes the last recipe of the app, and the order is applied to the variables of the same names (with their
order prefix). Also that a file with a value that is not valid, of another format or of a newer version is refused
without writing anything, and that the app (presets.py) loads it with the checks of its forms.

    cd segm_optim && python -s tests/test_parameters.py
"""

import os
import shutil
import sys
import tempfile

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SEGM = os.path.dirname(HERE)
sys.path.insert(0, SEGM)
sys.path.insert(0, os.path.join(os.path.dirname(SEGM), 'app'))  # ehabitat_app

import hft_state  # noqa: E402
import parameters  # noqa: E402

VARS = ['SN_NDMI_early.tif', 'SN_WAVI_min.tif', 'SN_slope.tif']
RECIPE = {'folder': 'gee_sn', 'period': ['2017-06-01', '2025-09-30'], 'comparison_area': 'pa/parks.gpkg',
          'comparison_area_layer': 'parks', 'project': 'my-ee-project',
          'variables': [{'index': 'NDMI', 'sensor': 'Sentinel-2', 'metric': 'median'}]}


def config(path, opts):
	shutil.copy(os.path.join(SEGM, 'config.py'), path)
	with open(path) as f:
		text = f.read()
	for k, v in opts.items():
		text = hft_state.set_option(text, k, v)
	with open(path, 'w') as f:
		f.write(text)


def env_vars(folder, files):
	return '{' + ', '.join(f"'{folder}_{os.path.splitext(n)[0]}': {{'file': '{folder}/{n}', 'nodata': -9999}}"
	                      for n in files) + '}'


def main():
	tmp = tempfile.mkdtemp()
	saved = (hft_state.CONFIG, hft_state.INVARS, parameters.APP_RECIPE)
	try:
		inv = os.path.join(tmp, 'inVars')
		for d in ('gee_sn', 'gee_park'):
			os.makedirs(os.path.join(inv, d))
		for n in VARS:  # the set of the first area, in the order SN_slope, SN_NDMI_early, SN_WAVI_min
			open(os.path.join(inv, 'gee_sn', n), 'w').close()
		for i, n in enumerate(['SN_WAVI_min.tif', 'SN_slope.tif', 'SN_NDMI_early.tif', 'other.tif'], 1):  # the other area
			open(os.path.join(inv, 'gee_park', f'{i}_{n}'), 'w').close()
		with open(os.path.join(inv, 'gee_sn', 'gee_vars.yaml'), 'w') as f:
			yaml.safe_dump(RECIPE, f)
		hft_state.INVARS = inv
		parameters.APP_RECIPE = os.path.join(tmp, 'app_recipe.yaml')

		a, b = os.path.join(tmp, 'config_a.py'), os.path.join(tmp, 'config_b.py')
		order_a = ['SN_slope.tif', 'SN_NDMI_early.tif', 'SN_WAVI_min.tif']
		config(a, {'GRASSDB': "'/data/grass_a'", 'STUDY_AREA': "'SierraNevada'", 'INVARS_FOLDER': "'gee_sn'",
		           'SCENARIO': "'borreguil'", 'ENV_VARS': env_vars('gee_sn', order_a), 'SEGM_THRESHOLDS': '[0.05, 0.07]',
		           'THRESHOLD_METHOD': '0.07', 'MIN_SEGM_SIZE': '1', 'HFT_SUBCLASSES': 'True', 'SEGM_PCA': 'False',
		           'SIMIL_CUTS': '[5, 10, 20]', 'SMALL_SEGM_MERGE': "'Similar'", 'CLIP_PERCENTILE': '0.5',
		           'ECOSIM_AREA': "'parques/Ordesa/ecoregion'", 'ECOSIM_CLASSES': "'subclass'", 'ECOSIM_CUT': "'p5'"})
		config(b, {'GRASSDB': "'/data/grass_b'", 'STUDY_AREA': "'Ordesa'", 'INVARS_FOLDER': "'gee_park'",
		           'SCENARIO': "'test'", 'ECOSIM_AREA': "'elsewhere/ecoregion'"})

		# save: the method, without the study area nor the computer
		hft_state.CONFIG = a
		cfg_a = hft_state.read_config()
		path = os.path.join(tmp, 'p.yaml')
		parameters.save(path, parameters.collect())
		data = parameters.read(path)
		for name in parameters.LOCAL:
			assert name not in data['model'] and name not in data['similarity'], name
		assert data['model']['SEGM_THRESHOLDS'] == [0.05, 0.07] and data['model']['THRESHOLD_METHOD'] == 0.07
		assert data['model']['SMALL_SEGM_MERGE'] == 'similar', data['model']['SMALL_SEGM_MERGE']  # as the model uses it
		assert data['similarity']['ECOSIM_CLASSES'] == 'subclass' and 'ECOSIM_AREA' not in data['similarity']
		assert set(data['model']) == set(parameters.MODEL) and set(data['similarity']) == set(parameters.SIMILARITY)
		assert data['variables_order'] == order_a
		assert data['recipe']['period'] == RECIPE['period']
		for k in parameters.RECIPE_LOCAL:
			assert k not in data['recipe'], k
		assert data['from']['study_area'] == 'SierraNevada' and data['from']['recipe'].endswith('gee_vars.yaml')
		print('save: ok')

		# load into the other area: the method changes, the area does not
		hft_state.CONFIG = b
		p = parameters.apply(data)
		cfg_b = hft_state.read_config()
		for name in list(parameters.MODEL) + list(parameters.SIMILARITY):
			assert parameters.model_value(name, cfg_b) == parameters.model_value(name, cfg_a), name
		for name, want in (('GRASSDB', '/data/grass_b'), ('STUDY_AREA', 'Ordesa'), ('INVARS_FOLDER', 'gee_park'),
		                   ('SCENARIO', 'test'), ('ECOSIM_AREA', 'elsewhere/ecoregion')):
			assert cfg_b[name] == want, (name, cfg_b[name])
		with open(parameters.APP_RECIPE) as f:
			rec = yaml.safe_load(f)
		assert rec['variables'] == RECIPE['variables'] and 'comparison_area' not in rec and 'project' not in rec
		assert rec['folder'] == 'gee_Ordesa', rec['folder']  # not the set of the first area
		# the folder of the study area chosen (before its setup.py, config.py has the previous one), a new one
		assert parameters.plan(data, study_area='Picos')['recipe']['folder'] == 'gee_Picos'
		os.makedirs(os.path.join(inv, 'gee_Picos'))
		assert parameters.plan(data, study_area='Picos')['recipe']['folder'] == 'gee_Picos_2'
		import var_order
		got = var_order.read_order(os.path.join(inv, 'gee_park'))
		assert got == ['2_SN_slope.tif', '3_SN_NDMI_early.tif', '1_SN_WAVI_min.tif', '4_other.tif'], got
		assert p['order']['missing'] == []
		# loading it again changes nothing
		assert parameters.plan(data)['options'] == []
		print('load: ok')

		# refused without writing anything
		with open(b) as f:
			before = f.read()
		bad = dict(data, model=dict(data['model'], THRESHOLD_METHOD=0.35))
		try:
			parameters.apply(bad)
			raise AssertionError('a fixed threshold that is not a segmentation threshold was loaded')
		except ValueError as e:
			assert 'THRESHOLD_METHOD' in str(e), e
		for wrong in (dict(data, format='other'), dict(data, version=parameters.VERSION + 1),
		              dict(data, model=dict(data['model'], SIMIL=4)), dict(data, model=dict(data['model'], SEGM_PCA='yes'))):
			wp = os.path.join(tmp, 'wrong.yaml')
			with open(wp, 'w') as f:
				yaml.safe_dump(wrong, f)
			try:
				parameters.apply(parameters.read(wp))
				raise AssertionError(f'loaded: {wrong.get("format")} {wrong.get("version")}')
			except ValueError:
				pass
		with open(b) as f:
			assert f.read() == before, 'config.py changed by a file that is not valid'
		# values that the forms refuse, or that would break config.py, are refused too (also on the command line)
		wp = os.path.join(tmp, 'wrong.yaml')
		for name, value in (('PA_BUFFER', float('inf')), ('PA_BUFFER', float('nan')), ('MIN_SEGM_SIZE', float('inf')),
		                    ('SMALL_SEGM_RADIUS', 0), ('THRESHOLD_ALPHA', 0), ('CLIP_PERCENTILE', 60),
		                    ('ECOSIM_MIN_COVER', 0), ('ECOSIM_MIN_COVER', 1.5)):
			sec = 'model' if name in parameters.MODEL else 'similarity'
			with open(wp, 'w') as f:
				yaml.safe_dump(dict(data, **{sec: dict(data[sec], **{name: value})}), f)
			assert parameters.plan(parameters.read(wp))['errors'], (name, value)
			try:
				parameters.main(['load', wp])
				raise AssertionError(f'{name} = {value} loaded on the command line')
			except SystemExit as e:
				assert str(e).startswith('ERROR'), e
		with open(b) as f:
			assert f.read() == before, 'config.py changed by a value that is not valid'
		# a recipe that the app could not read is refused
		for rec in ({'folder': 'x'}, dict(RECIPE, study=None), dict(RECIPE, period=2020), dict(RECIPE, variables='NDVI'),
		            dict(RECIPE, **{'where': 5}), {**RECIPE, 2020: 'x'}):
			with open(wp, 'w') as f:
				yaml.safe_dump(dict(data, recipe=rec), f)
			try:
				parameters.read(wp)
				raise AssertionError(f'recipe loaded: {rec}')
			except ValueError as e:
				assert 'recipe' in str(e), e
		# whole numbers written as decimals in a choice (SIMIL 1.0): loaded as the whole number
		p = parameters.plan(dict(data, model=dict(data['model'], SIMIL=1.0)))
		assert [o['literal'] for o in p['options'] if o['name'] == 'SIMIL'] == ['1'], p['options']
		# a fixed CRS of the recipe is reported
		assert parameters.plan(dict(data, recipe=dict(RECIPE, study={'crs': 'EPSG:32630'})))['recipe']['crs'] == 'EPSG:32630'
		# an option of a study area in the file is left out, not loaded
		p = parameters.plan(dict(data, model=dict(data['model'], GRASSDB='/x')))
		assert any(x.startswith('GRASSDB') for x in p['ignored']) and not p['errors']
		print('refused: ok')

		# an older config.py without some options: the file has the values the scripts use for them (not those of the
		# other area when loaded)
		old = os.path.join(tmp, 'config_old.py')
		with open(os.path.join(SEGM, 'config.py')) as f:
			text = f.read()
		for name in ('SEGM_THRESHOLDS', 'SMALL_SEGM_RADIUS', 'HFT_SUBCLASSES', 'SIMIL', 'ECOSIM_CLASSES'):
			text = '\n'.join(x for x in text.splitlines() if not x.startswith(name + ' ='))
		with open(old, 'w') as f:
			f.write(text)
		hft_state.CONFIG = old
		d = parameters.collect(recipe=False)
		assert d['model']['SEGM_THRESHOLDS'] == hft_state.DEFAULT_THRESHOLDS and d['model']['SMALL_SEGM_RADIUS'] == 5
		assert d['model']['HFT_SUBCLASSES'] is False and d['model']['SIMIL'] == 2 and d['similarity']['ECOSIM_CLASSES'] == 'hft'
		hft_state.CONFIG = b  # which has the options of the first area now
		parameters.apply(d)
		cfg_b = hft_state.read_config()
		assert cfg_b['SEGM_THRESHOLDS'] == hft_state.DEFAULT_THRESHOLDS and cfg_b['HFT_SUBCLASSES'] is False, cfg_b['SEGM_THRESHOLDS']
		assert cfg_b['SMALL_SEGM_RADIUS'] == 5 and cfg_b['SIMIL'] == 2
		# the command line: a relative file is that of the folder where it is run
		cwd = os.getcwd()
		try:
			os.chdir(tmp)
			parameters.START = tmp
			parameters.main(['save', 'rel.yaml'])
			assert os.path.isfile(os.path.join(tmp, 'rel.yaml'))
		finally:
			os.chdir(cwd)
		# a config.py edited by hand outside the limits of the forms: saved, with a warning
		hft_state.CONFIG = old
		hft_state.write_option('THRESHOLD_ALPHA', '0.02')
		assert any('THRESHOLD_ALPHA' in w for w in parameters.problems(parameters.collect(recipe=False)))
		print('older config.py and command line: ok')

		# the app: the same file into a third area, with the checks of its forms
		c = os.path.join(tmp, 'config_c.py')
		config(c, {'STUDY_AREA': "'Teide'", 'INVARS_FOLDER': "'gee_park'"})
		hft_state.CONFIG = c
		from ehabitat_app import model, presets
		chosen = model.study_area
		model.study_area = lambda: {'file': os.path.join(tmp, 'Teide_area.gpkg')}  # the study area chosen in the app
		pv = presets.preview(path)
		names = {o['name'] for o in pv['options']}
		assert {'SEGM_THRESHOLDS', 'THRESHOLD_METHOD', 'HFT_SUBCLASSES', 'ECOSIM_CLASSES'} <= names, names
		assert all(o['label'] != o['name'] for o in pv['options'])  # the labels of the forms
		r = presets.load(path)
		cfg_c = hft_state.read_config()
		for name in list(parameters.MODEL) + list(parameters.SIMILARITY):
			assert parameters.model_value(name, cfg_c) == parameters.model_value(name, cfg_a), name
		assert cfg_c['STUDY_AREA'] == 'Teide' and r['recipe']
		from ehabitat_app import gee
		rec = gee.load_recipe(parameters.APP_RECIPE)['recipe']  # as the page Download from MonitorEO reads it
		# the folder of the study area chosen in the app, not that of config.py (of the last setup.py)
		assert rec['folder'] == 'gee_Teide_area' and rec['variables'] == RECIPE['variables'] and gee.clean_recipe(rec)
		out = os.path.join(tmp, 'exported')
		e = presets.export(out)
		assert e['path'] == out + '.yaml' and os.path.isfile(e['path'])
		# a value that the forms refuse is refused by the app without writing anything
		with open(c) as f:
			before = f.read()
		with open(wp, 'w') as f:
			yaml.safe_dump(dict(data, model=dict(data['model'], PA_BUFFER=float('nan'))), f)
		try:
			presets.load(wp)
			raise AssertionError('the app loaded PA_BUFFER = nan')
		except ValueError:
			pass
		with open(c) as f:
			assert f.read() == before
		for v in ('inf', 'nan', float('inf')):
			try:
				model._literal(model.BY_NAME['PA_BUFFER'], v)
				raise AssertionError(f'the form took PA_BUFFER = {v}')
			except ValueError:
				pass
		print('app: ok')
		model.study_area = chosen
	finally:
		hft_state.CONFIG, hft_state.INVARS, parameters.APP_RECIPE = saved
		shutil.rmtree(tmp, ignore_errors=True)
	print('all ok')


if __name__ == '__main__':
	main()
