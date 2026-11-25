#### eHabitat+: test of the order of the variables of a set (var_order.py, setup.py, hft_state.py, the app)
#### License: CC BY-SA 3.0
"""
Checks the order of the variables of the model: that of variables_order.txt of the set, then that of the last setup.py
(ENV_VARS of config.py), then by name (sorted, capitals first); that the app lists the variables of the study area in
that order (inputs.get_set) and writes it when they are moved (inputs.set_order, without renaming the files), that
hft_state asks for setup.py again when the order changed, and that the ecoregion pairs (variable_key) are kept.

    cd segm_optim && python -s tests/test_var_order.py
"""

import os
import shutil
import sys
import tempfile

import numpy as np
from osgeo import gdal, osr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'app'))  # ehabitat_app

import hft_state  # noqa: E402
import var_order  # noqa: E402

gdal.UseExceptions()
NAMES = ['SN_NDMI_SD.tif', 'SN_NDMI_early.tif', 'SN_WAVI_min.tif', 'SN_slope.tif']  # sorted(): capitals first


def tif(path):
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(32630)
	ds = gdal.GetDriverByName('GTiff').Create(path, 4, 4, 1, gdal.GDT_Float32)
	ds.SetGeoTransform((450000, 10, 0, 4110000, 0, -10))
	ds.SetProjection(srs.ExportToWkt())
	ds.GetRasterBand(1).WriteArray(np.ones((4, 4), np.float32))
	ds = None


def test_ordered():
	tmp = tempfile.mkdtemp()
	try:
		paths = [os.path.join(tmp, n) for n in NAMES]
		assert var_order.ordered(list(reversed(paths)), tmp) == sorted(paths)  # by name: SD before early, WAVI before slope
		# the order of the last setup.py, the new files after it by name
		prev = ['SN_slope.tif', 'SN_NDMI_early.tif', 'gone.tif']
		got = [os.path.basename(p) for p in var_order.ordered(paths, tmp, prev)]
		assert got == ['SN_slope.tif', 'SN_NDMI_early.tif', 'SN_NDMI_SD.tif', 'SN_WAVI_min.tif'], got
		# variables_order.txt first (its comments, repeated names, paths and missing files left out)
		with open(os.path.join(tmp, var_order.ORDER_FILE), 'w') as f:
			f.write('# comment\nSN_WAVI_min.tif\n\nSN_WAVI_min.tif\nother/SN_slope.tif\nmissing.tif\n')
		assert var_order.read_order(tmp) == ['SN_WAVI_min.tif', 'missing.tif']
		got = [os.path.basename(p) for p in var_order.ordered(paths, tmp, prev)]
		assert got == ['SN_WAVI_min.tif', 'SN_slope.tif', 'SN_NDMI_early.tif', 'SN_NDMI_SD.tif'], got
		var_order.write_order(tmp, ['SN_slope.tif', 'SN_NDMI_SD.tif'])
		assert var_order.read_order(tmp) == ['SN_slope.tif', 'SN_NDMI_SD.tif']
		env = {'a': {'file': 'set/SN_slope.tif'}, 'b': {'file': 'other/x.tif'}, 'c': {'file': 'set/SN_NDMI_SD.tif'},
		       'd': {'file': 'y.tif'}}
		assert var_order.previous_order(env, 'set') == ['SN_slope.tif', 'SN_NDMI_SD.tif']
		assert var_order.previous_order(env, '') == ['y.tif']
		# a BOM (Notepad of Windows) is left out of the first name
		with open(os.path.join(tmp, var_order.ORDER_FILE), 'w', encoding='utf-8-sig') as f:
			f.write('SN_slope.tif\nSN_NDMI_SD.tif\n')
		assert var_order.read_order(tmp) == ['SN_slope.tif', 'SN_NDMI_SD.tif']
		assert var_order.only_order([var_order.ORDER_CHANGED + ' (x)']) and not var_order.only_order([]) and \
			not var_order.only_order([var_order.ORDER_CHANGED, 'new file inVars/x.tif'])
		print('test_ordered: OK')
	finally:
		shutil.rmtree(tmp)


def test_app_and_state():
	from ehabitat_app import inputs, paths
	tmp = tempfile.mkdtemp()
	import setup
	real = (paths.INVARS, hft_state.INVARS, hft_state.read_config, hft_state.get_setup, setup.INVARS)
	try:
		paths.INVARS = hft_state.INVARS = os.path.join(tmp, 'inVars')
		d = os.path.join(paths.INVARS, 'set')
		os.makedirs(os.path.join(d, 'ecoregion'))
		for n in NAMES:
			tif(os.path.join(d, n))
			tif(os.path.join(d, 'ecoregion', n))
		env = {n: {'file': f'set/{n}', 'nodata': None} for n in ['SN_slope.tif', 'SN_NDMI_SD.tif', 'SN_NDMI_early.tif',
		                                                          'SN_WAVI_min.tif']}
		cfg = {'INVARS_FOLDER': 'set', 'STUDY_AREA': 'sa', 'ENV_VARS': env}
		hft_state.read_config = lambda: cfg
		rec = {'invars_folder': 'set', 'variables': {f'set/{n}': hft_state.file_sig(os.path.join(d, n)) for n in NAMES}}
		hft_state.get_setup = lambda name, state=None: rec
		# the set used by the model: in the order of config.py
		st = inputs.get_set('set')
		assert [r['file'] for r in st['study']] == list(env), [r['file'] for r in st['study']]
		assert hft_state.setup_changes(cfg) == []
		# moved: written to variables_order.txt, the files keep their names, setup.py has to run again
		new = ['SN_WAVI_min.tif', 'SN_slope.tif', 'SN_NDMI_SD.tif', 'SN_NDMI_early.tif']
		assert inputs.set_order('set', new) == new
		assert sorted(os.listdir(d)) == sorted(NAMES + ['ecoregion', var_order.ORDER_FILE])
		assert [r['file'] for r in inputs.get_set('set')['study']] == new
		why = hft_state.setup_changes(cfg)
		assert len(why) == 1 and 'order of the variables changed' in why[0], why
		# once setup.py wrote the new order in config.py, nothing to do
		cfg['ENV_VARS'] = {n: env[n] for n in new}
		assert hft_state.setup_changes(cfg) == []
		# the pairs with the ecoregion are those of the names (variable_key), whatever the order
		assert inputs.check_scales(inputs.get_set('set')['study'], inputs.get_set('set')['ecoregion'])[0]['level'] == 'ok'
		# another set (not used by the model): by name, or variables_order.txt
		o = os.path.join(paths.INVARS, 'other')
		os.makedirs(o)
		for n in NAMES:
			tif(os.path.join(o, n))
		assert [r['file'] for r in inputs.get_set('other')['study']] == sorted(NAMES)
		inputs.set_order('other', list(reversed(sorted(NAMES))))
		assert [r['file'] for r in inputs.get_set('other')['study']] == list(reversed(sorted(NAMES)))
		# a list that is not that of the set is refused
		for bad in (NAMES[:3], NAMES + ['x.tif'], NAMES[:3] + NAMES[:1]):
			try:
				inputs.set_order('other', bad)
			except ValueError:
				pass
			else:
				raise AssertionError(f'{bad} was not refused')
		# setup.py takes the same order as the app (read_variables, without GRASS)
		import setup
		setup.INVARS = paths.INVARS
		got = [os.path.basename(v['path']) for v in setup.read_variables(
			'other', var_order.previous_order(cfg['ENV_VARS'], 'other'))[0].values()]
		assert got == [r['file'] for r in inputs.get_set('other')['study']], got
		got = [os.path.basename(v['path']) for v in setup.read_variables('set', var_order.previous_order(cfg['ENV_VARS'], 'set'))[0].values()]
		assert got == [r['file'] for r in inputs.get_set('set')['study']] == new, got
		# a file that cannot be read is in the order too (setup.py would read it): the list of the app includes it
		with open(os.path.join(o, 'broken.tif'), 'w') as f:
			f.write('not a GeoTIFF')
		files = [r['file'] for r in inputs.get_set('other')['study'] if r['included']]
		assert 'broken.tif' in files and inputs.set_order('other', ['broken.tif'] + [x for x in files if x != 'broken.tif'])[0] == 'broken.tif'
		# a config.py that cannot be read does not hide the sets
		def broken():
			raise SyntaxError('invalid syntax')
		hft_state.read_config = broken
		assert len(inputs.get_set('set')['study']) == len(NAMES)
		print('test_app_and_state: OK')
	finally:
		paths.INVARS, hft_state.INVARS, hft_state.read_config, hft_state.get_setup, setup.INVARS = real
		shutil.rmtree(tmp)


if __name__ == '__main__':
	test_ordered()
	test_app_and_state()
	print('All the tests passed')
