#### eHabitat+: test of the axes of the PCA of the segmentation (pca_axes.py, SEGM_PCA_AXES, SEGM_PCA_VARIANCE)
#### License: CC BY-SA 3.0
"""
Checks that the output of i.pca is read, that SEGM_PCA_VARIANCE 0 keeps the first 3 axes and a % takes the fewest
first axes that reach it (all of them if never reached), that SEGM_PCA_AXES takes the first N (overriding the %), that
a value outside 0 / 50-100 (or 0 / 1 or more axes) is refused by config.py, the app and the parameters files, that a %
or a number of axes makes the segmentation out of date (and 0 does not, as before the options), that a parameters file
saved before the options sets them to 0, and, with GRASS, that i.pca of the segmentation gives the table that
pca_axes.parse reads and the loadings of pca_preview.

    cd segm_optim && python -s tests/test_pca_axes.py
"""

import math
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SEGM = os.path.dirname(HERE)
sys.path.insert(0, SEGM)
sys.path.insert(0, os.path.join(os.path.dirname(SEGM), 'app'))

import hft_state  # noqa: E402
import pca_axes  # noqa: E402
import pca_preview  # noqa: E402

OUT = """PC1      4.44 ( 0.4608, 0.0914, 0.4476, 0.4392, 0.3311, 0.3530,-0.0953, 0.3779) [55.45%]
PC2      1.31 ( 0.0922,-0.8027,-0.0030, 0.1494,-0.4730, 0.1012,-0.2500, 0.1686) [16.32%]
PC3      0.98 ( 0.0262,-0.1248, 0.1053, 0.1888,-0.1822,-0.1268, 0.9274, 0.1661) [12.27%]
PC4      0.59 (-0.0196, 0.1403,-0.0190, 0.1635,-0.0474,-0.7819,-0.2242, 0.5376) [ 7.32%]
PC5      0.35 (-0.0491,-0.4518, 0.1996, 0.2317, 0.5720,-0.3947, 0.0125,-0.4657) [ 4.31%]
PC6      0.26 (-0.0734, 0.3026, 0.5043, 0.3593,-0.5353,-0.1206,-0.1288,-0.4494) [ 3.22%]
PC7      0.07 ( 0.7792, 0.0939,-0.4762, 0.1424,-0.1176,-0.1711, 0.0239,-0.3054) [ 0.83%]
PC8      0.02 (-0.4039, 0.0874,-0.5170, 0.7206, 0.0611, 0.1954,-0.0267, 0.0034) [ 0.28%]
"""


def test_choose():
	axes = pca_axes.parse(OUT)
	assert [a[0] for a in axes] == list(range(1, 9)) and axes[3][2] == 7.32 and axes[0][1] == 4.44
	p = [x for _, _, x in axes]
	assert pca_axes.choose(p, 0) == 3  # 84.04 %
	assert pca_axes.choose(p, 84.04) == 3 and pca_axes.choose(p, 84.05) == 4
	assert pca_axes.choose(p, 90) == 4 and pca_axes.choose(p, 95) == 5 and pca_axes.choose(p, 100) == 8
	assert pca_axes.choose(p[:2], 0) == 2  # fewer than 3 axes
	assert pca_axes.choose([60.0, 40.0], 100) == 2
	assert pca_axes.choose([50.0, 34.2, 15.81, 0.0], 100) == 4  # 100.01 % or 99.99 %: all of them
	assert pca_axes.choose([50.0, 34.2, 15.79, 0.0], 100) == 4
	assert len(pca_axes.parse('PC1  1.00 ( 1.0) [100.00%]\r\nPC2  0.00 ( 0.0) [-0.00%]')) == 2
	assert pca_axes.summary([(1, 1, 60.0), (2, 1, 40.0)], 2, 0).endswith('(all the axes, fewer than 3)')
	s = pca_axes.summary(axes, 3, 0)
	assert s == 'PC1-PC3 of 8 axes, 84.0 % of the variance (the first 3 axes, SEGM_PCA_VARIANCE = 0)', s
	assert pca_axes.summary(axes, 4, 90).startswith('PC1-PC4 of 8 axes, 91.4 %')
	# a number of axes: the first N whatever the %, all of them if there are fewer
	assert pca_axes.choose(p, 0, 5) == 5 and pca_axes.choose(p, 90, 2) == 2 and pca_axes.choose(p, 0, 12) == 8
	assert pca_axes.choose(p, 95, 0) == 5 and pca_axes.choose(p, 0, 1) == 1
	assert pca_axes.summary(axes, 5, 0, 5) == 'PC1-PC5 of 8 axes, 95.7 % of the variance (SEGM_PCA_AXES = 5)'
	assert pca_axes.summary(axes, 8, 90, 12).endswith('(SEGM_PCA_AXES = 12, but the area has 8 axes)')
	assert pca_axes.summary(axes, 1, 0, 1).startswith('PC1 of 8 axes, 55.5 %')
	with tempfile.TemporaryDirectory() as t:
		f = os.path.join(t, 'a.csv')
		pca_axes.write(f, axes, 4)
		rows = open(f).read().split('\n')
		assert rows[0] == 'axis eigenvalue variance_pct cumulative_pct used' and rows[4].endswith(' 91.36 1') and rows[5].endswith(' 0')
	print('test_choose: OK')


def test_option():
	for ok in (0, None, 50, 90, 99.5, 100):
		pca_axes.variance_option({'SEGM_PCA_VARIANCE': ok})
	assert pca_axes.variance_option({}) == 0
	for bad in (30, 49.9, 101, -5, True, 'x', float('nan')):
		try:
			pca_axes.variance_option({'SEGM_PCA_VARIANCE': bad})
			raise AssertionError(f'{bad!r} accepted')
		except ValueError:
			pass
	cfg = hft_state.read_config()
	assert hft_state.options_error(dict(cfg, SEGM_PCA_VARIANCE=30)).startswith('SEGM_PCA_VARIANCE')
	# the segmentation is out of date with a %, not with 0 or without the option (as before it), nor without the PCA
	base = hft_state.step_inputs('segm', dict(cfg, SEGM_PCA=True, SEGM_PCA_VARIANCE=0))
	assert 'SEGM_PCA_VARIANCE' not in base
	assert hft_state.step_inputs('segm', {k: v for k, v in dict(cfg, SEGM_PCA=True).items() if k != 'SEGM_PCA_VARIANCE'}) == base
	assert hft_state.step_inputs('segm', dict(cfg, SEGM_PCA=True, SEGM_PCA_VARIANCE=90))['SEGM_PCA_VARIANCE'] == '90'
	assert 'SEGM_PCA_VARIANCE' not in hft_state.step_inputs('segm', dict(cfg, SEGM_PCA=False, SEGM_PCA_VARIANCE=90))
	# a number of axes: recorded instead of the % (which it overrides), not with 0 or without the PCA
	for ok in (0, None, 1, 3, 4.0, 12):
		pca_axes.axes_option({'SEGM_PCA_AXES': ok})
	assert pca_axes.axes_option({}) == 0 and pca_axes.axes_option({'SEGM_PCA_AXES': 4.0}) == 4
	for bad in (-1, 2.5, True, '4', float('nan'), float('inf')):
		try:
			pca_axes.axes_option({'SEGM_PCA_AXES': bad})
			raise AssertionError(f'{bad!r} accepted')
		except ValueError:
			pass
	assert hft_state.options_error(dict(cfg, SEGM_PCA_AXES=2.5)).startswith('SEGM_PCA_AXES')
	assert hft_state.step_inputs('segm', dict(cfg, SEGM_PCA=True, SEGM_PCA_VARIANCE=0, SEGM_PCA_AXES=0)) == base
	ax = hft_state.step_inputs('segm', dict(cfg, SEGM_PCA=True, SEGM_PCA_VARIANCE=90, SEGM_PCA_AXES=4))
	assert ax['SEGM_PCA_AXES'] == '4' and 'SEGM_PCA_VARIANCE' not in ax
	assert hft_state.step_inputs('segm', dict(cfg, SEGM_PCA=True, SEGM_PCA_AXES=5))['SEGM_PCA_AXES'] == '5'
	assert 'SEGM_PCA_AXES' not in hft_state.step_inputs('segm', dict(cfg, SEGM_PCA=False, SEGM_PCA_AXES=4))
	# the app and the parameters files
	from ehabitat_app import model
	import parameters
	p = model.BY_NAME['SEGM_PCA_VARIANCE']
	assert model._literal(p, 0) == '0' and model._literal(p, 90) == '90' and model._literal(p, None) == '0'
	for bad in (30, 101, 'nan', 'inf'):
		try:
			model._literal(p, bad)
			raise AssertionError(f'the form took {bad!r}')
		except ValueError:
			pass
	assert parameters.check('SEGM_PCA_VARIANCE', 95) == '95' and parameters.check('SEGM_PCA_VARIANCE', 90.0) == '90'
	assert 'SEGM_PCA_VARIANCE' in parameters.collect(cfg, recipe=False)['model']
	# the number of axes in the form: a number typed (a text) of 1 or more, or 0 (the other modes)
	q = model.BY_NAME['SEGM_PCA_AXES']
	assert q['within'] == 'SEGM_PCA_VARIANCE'
	assert model._literal(q, 0) == '0' and model._literal(q, None) == '0' and model._literal(q, '4') == '4'
	assert model._literal(q, 6) == '6' and model._literal(q, 4.0) == '4' and model._literal(q, ' 3 ') == '3'
	for bad in ('', '0', '2.5', '-1', 'x', 2.5, -1, 'nan', float('inf')):
		try:
			model._literal(q, bad)
			raise AssertionError(f'the form took {bad!r}')
		except ValueError:
			pass
	assert parameters.check('SEGM_PCA_AXES', 4) == '4' and parameters.check('SEGM_PCA_AXES', 4.0) == '4'
	for bad in (2.5, -1, True, 'x'):
		try:
			parameters.check('SEGM_PCA_AXES', bad)
			raise AssertionError(f'the file took {bad!r}')
		except ValueError:
			pass
	assert parameters.collect(dict(cfg, SEGM_PCA_AXES=None), recipe=False)['model']['SEGM_PCA_AXES'] == 0
	# a file saved before the options (without them) sets them to 0, so that its % is not overridden by the number of
	# axes of config.py; a file with them, as they are
	now = dict(cfg, SEGM_PCA_AXES=4, SEGM_PCA_VARIANCE=0)
	old = {'model': {'SEGM_PCA': True, 'PA_BUFFER': cfg.get('PA_BUFFER', 0)}}
	ch = {o['name']: o['new'] for o in parameters.plan(old, now)['options']}
	assert ch.get('SEGM_PCA_AXES') == 0 and 'SEGM_PCA_VARIANCE' not in ch, ch
	ch = {o['name']: o['new'] for o in parameters.plan({'model': dict(old['model'], SEGM_PCA_VARIANCE=90)}, now)['options']}
	assert ch.get('SEGM_PCA_AXES') == 0 and ch.get('SEGM_PCA_VARIANCE') == 90, ch
	ch = {o['name']: o['new'] for o in parameters.plan({'model': dict(old['model'], SEGM_PCA_AXES=6)}, now)['options']}
	assert ch.get('SEGM_PCA_AXES') == 6, ch
	# not a change when config.py has no option and the file is older (0, as the scripts use it)
	bare = {k: v for k, v in cfg.items() if k not in ('SEGM_PCA_AXES', 'SEGM_PCA_VARIANCE')}
	assert not [o for o in parameters.plan(old, bare)['options'] if o['name'].startswith('SEGM_PCA_')]
	print('test_option: OK')


def test_grass():
	"""i.pca of GRASS (as the segmentation calls it) gives the table that parse reads."""
	code = ("import grass.script as gs\n"
	        "gs.run_command('g.region', n=60, s=0, e=60, w=0, res=1)\n"
	        "for i in range(5):\n"
	        "    gs.mapcalc(f'v{i} = rand(0,100) + {i}*row() + col()*{i % 2}', seed=i + 1, quiet=True)\n"
	        "print(gs.read_command('i.pca', flags='n', input='v0,v1,v2,v3,v4', output='pc', overwrite=True, quiet=True))\n"
	        "gs.mapcalc('c = 0.1', quiet=True)\n"
	        "print('CONSTANT>>' + gs.read_command('i.pca', flags='n', input='v0,v1,c', output='pcc', overwrite=True, quiet=True))\n")
	with tempfile.TemporaryDirectory() as t:
		f = os.path.join(t, 'p.py')
		open(f, 'w').write(code)
		try:
			r = subprocess.run(['grass', '--tmp-project', 'XY', '--exec', sys.executable, f], capture_output=True, text=True,
			                   timeout=300)
		except FileNotFoundError:
			print('test_grass: skipped (no grass)')
			return
	normal, constant = r.stdout.split('CONSTANT>>')
	axes = pca_axes.parse(normal)
	assert len(axes) == 5 and abs(sum(p for _, _, p in axes) - 100) < 0.1, r.stdout + r.stderr
	# with a constant variable i.pca gives nan, which the segmentation stops on (segmentation_pca_par.py)
	import re
	assert re.search(r'\bnan\b', constant, re.I) and len(pca_axes.parse(constant)) < 3, constant
	print('test_grass: OK')


def _data(t):
	"""Four GeoTIFFs (with NoData, not in the same pixels) and two study areas in a shapefile and a GeoPackage."""
	import numpy as np
	from osgeo import gdal, ogr, osr
	gdal.UseExceptions()
	ogr.UseExceptions()
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(32630)
	rng = np.random.default_rng(1)
	rows, cols = np.indices((120, 150))
	base = [rng.normal(size=(120, 150)) + k * rows / 50 + (k % 2) * cols / 40 for k in range(4)]
	base[1][:, :20] = -9999  # NoData in different places in each variable
	base[2][100:, :] = np.nan
	base[3][50:60, 50:60] = -9999
	files = []
	for k, a in enumerate(base):
		f = os.path.join(t, f'v{k}.tif')
		d = gdal.GetDriverByName('GTiff').Create(f, 150, 120, 1, gdal.GDT_Float32)
		d.SetGeoTransform((450000, 10, 0, 4110000, 0, -10))
		d.SetProjection(srs.ExportToWkt())
		b = d.GetRasterBand(1)
		b.WriteArray(a.astype(np.float32))
		if k != 2:
			b.SetNoDataValue(-9999)
		d = None
		files.append(f)
	vec = {}
	for drv, ext in (('ESRI Shapefile', 'shp'), ('GPKG', 'gpkg')):
		f = os.path.join(t, f'areas.{ext}')
		ds = ogr.GetDriverByName(drv).CreateDataSource(f)
		lyr = ds.CreateLayer('areas', srs, ogr.wkbPolygon)
		lyr.CreateField(ogr.FieldDefn('name', ogr.OFTString))
		for n, wkt in enumerate(('POLYGON((450105 4109895,450605 4109895,450605 4109405,450105 4109405,450105 4109895))',
		                         'POLYGON((450700 4109300,451400 4109300,451400 4108900,450700 4108900,450700 4109300))')):
			ft = ogr.Feature(lyr.GetLayerDefn())
			ft.SetField('name', 'x')
			ft.SetGeometry(ogr.CreateGeometryFromWkt(wkt))
			if ext == 'gpkg':
				ft.SetFID(5 + 4 * n)  # FIDs 5 and 9: GRASS numbers them 1 and 2
			lyr.CreateFeature(ft)
		ds = None
		vec[ext] = f
	# a fifth variable that covers only the western 100 columns (area 2 is partly outside it)
	f = os.path.join(t, 'v4.tif')
	d = gdal.GetDriverByName('GTiff').Create(f, 100, 120, 1, gdal.GDT_Float32)
	d.SetGeoTransform((450000, 10, 0, 4110000, 0, -10))
	d.SetProjection(srs.ExportToWkt())
	d.GetRasterBand(1).WriteArray((rng.normal(size=(120, 100)) + cols[:, :100] / 30).astype(np.float32))
	d = None
	files.append(f)
	return files, vec


def test_preview_grass():
	"""preview() gives the % of i.pca -n of GRASS on the same rectangle (the region of the segmentation: the polygons of
	the ID, aligned to the grid, grown by PA_BUFFER), with NoData in some variables; the IDs as v.in.ogr numbers them."""
	import shutil
	if not shutil.which('grass'):
		print('test_preview_grass: skipped (no grass)')
		return
	with tempfile.TemporaryDirectory() as t:
		files, vec = _data(t)
		for ext in ('shp', 'gpkg'):
			for area, buffer in (('1', 0), ('2', 25)):
				got = pca_preview.preview(files, vec[ext], 'cat', [area], buffer)[0]
				assert got['axes'], (ext, area)
				code = ("import grass.script as gs, math\n"
				        f"files = {files!r}\n"
				        "for i, f in enumerate(files):\n"
				        "    gs.run_command('r.in.gdal', input=f, output=f'v{i}', quiet=True, overwrite=True)\n"
				        f"gs.run_command('v.in.ogr', input={vec[ext]!r}, output='pa', quiet=True, overwrite=True)\n"
				        f"gs.run_command('v.extract', input='pa', output='one', where='cat = {area}', quiet=True, overwrite=True)\n"
				        f"gs.run_command('g.region', vector='one', align='v0', grow=int(math.ceil({buffer} / 10)) if {buffer} else 0)\n"
				        "print(gs.read_command('i.pca', flags='n', input=','.join(f'v{i}' for i in range(len(files))), output='pc', "
				        "overwrite=True, quiet=True))\n")
				f = os.path.join(t, 'g.py')
				open(f, 'w').write(code)
				r = subprocess.run(['grass', '--tmp-project', files[0], '--exec', sys.executable, f], capture_output=True,
				                   text=True, timeout=300)
				want = [p for _, _, p in pca_axes.parse(r.stdout)]
				have = [p for _, _, p in got['axes']]
				assert want and max(abs(a - b) for a, b in zip(want, have)) <= 0.011, (ext, area, buffer, want, have, r.stderr[-500:])
				# the loadings: the eigenvectors of i.pca (in parentheses, 4 decimals) times the square root of the
				# eigenvalues, with the sign of the axis that of its largest loading
				vecs = [[float(x) for x in m.split(',')] for m in re.findall(r'^\s*PC\d+\s+\S+\s+\((.*)\)', r.stdout, re.M)]
				assert len(vecs) == len(files) == len(got['loadings']), (vecs, got['loadings'])
				for i, (gv, load) in enumerate(zip(vecs, got['loadings'])):
					eig = got['axes'][i][1]
					sign = 1 if gv[max(range(len(gv)), key=lambda j: abs(gv[j]))] >= 0 else -1
					mine = [x / math.sqrt(eig) for x in load]
					assert max(abs(sign * a - b) for a, b in zip(gv, mine)) <= 0.002, (ext, area, i, gv, mine)
					assert abs(sum(x * x for x in load) - eig) < 1e-3  # the squared loadings of an axis add up to its eigenvalue
				# all the axes take all of every variable; the first ones, less
				assert all(abs(c - 100) < 0.05 for c in pca_preview.captured(got['loadings'], len(files)))
				one = pca_preview.captured(got['loadings'], 1)
				assert all(0 <= c <= 100 for c in one) and abs(sum(one) - 100 * got['axes'][0][1]) < 0.1
		# a constant variable: no PCA (i.pca gives nan), said
		from osgeo import gdal
		c = os.path.join(t, 'const.tif')
		gdal.Translate(c, files[0], scaleParams=[[-100, 100, 0, 0]], outputType=gdal.GDT_Float32)
		r = pca_preview.preview(files[:3] + [c], vec['shp'], 'cat', ['1'])[0]
		assert r['axes'] == [] and r['constant'] == ['const.tif'], r
		# a constant that is not a whole number (0.1 in float32, rounding errors of one pass of sums)
		import numpy as np
		c2 = os.path.join(t, 'const01.tif')
		d = gdal.GetDriverByName('GTiff').Create(c2, 150, 120, 1, gdal.GDT_Float32)
		d.SetGeoTransform((450000, 10, 0, 4110000, 0, -10))
		d.SetProjection(gdal.Open(files[0]).GetProjection())
		d.GetRasterBand(1).WriteArray(np.full((120, 150), 0.1, np.float32))
		d = None
		for area, buf in (('1', 0), ('2', 500)):
			r = pca_preview.preview(files + [c2], vec['shp'], 'cat', [area], buf)[0]
			assert r['axes'] == [] and r['constant'] == ['const01.tif'], (area, r)
		assert len(pca_preview.preview(files, vec['gpkg'], 'cat', None, limit=1)) == 1
		# the variables that load most on an axis
		assert pca_preview.main_variables([0.1, -0.9, 0.5, 0.2], list('abcd')) == [('b', -0.9), ('c', 0.5)]
		assert pca_preview.main_variables([0.1, 0.3, 0.2], list('abc')) == [('b', 0.3)]
	print('test_preview_grass: OK')


def test_app_preview():
	"""The PCA of the app (model.pca_preview, behind "See the PCA in the study area"): the % and the loadings of every
	area chosen, as pca_preview gives them (it failed, the module hidden by the function of the same name)."""
	from ehabitat_app import model
	with tempfile.TemporaryDirectory() as t:
		files, vec = _data(t)
		saved = model.pca_view.study_inputs, model.study_area, model.read_config
		try:
			model.pca_view.study_inputs = lambda cfg: (files, vec['shp'], None, 'cat', None)
			model.study_area = lambda: {'file': vec['shp'], 'layer': None, 'col_id': 'cat', 'ids': ['1', '2']}
			model.read_config = lambda: {'PA_BUFFER': 0}
			r = model.pca_preview(25)
		finally:
			model.pca_view.study_inputs, model.study_area, model.read_config = saved
		want = pca_preview.preview(files, vec['shp'], 'cat', ['1', '2'], 25)
	assert r['variables'] == [os.path.basename(f) for f in files] and r['more'] == 0
	assert [a['id'] for a in r['areas']] == ['1', '2']
	for a, w in zip(r['areas'], want):
		assert a['percents'] == [p for _, _, p in w['axes']] and a['loadings'] == w['loadings'] and a['constant'] == []
		assert len(a['loadings']) == len(files) and all(len(x) == len(files) for x in a['loadings'])
	print('test_app_preview: OK')


if __name__ == '__main__':
	test_choose()
	test_option()
	test_grass()
	test_preview_grass()
	test_app_preview()
	print('All the tests passed')
