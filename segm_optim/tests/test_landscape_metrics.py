#### eHabitat+: test of the landscape metrics of maps of patches (landscape_metrics.py) and of the similar areas
#### License: CC BY-SA 3.0
"""
Checks that landscape_metrics.py measures the same habitat as a raster and as polygons (on the same grid, inside the
same boundary) with the same metrics; writes its tables, the patches as polygons (one per patch) and its record; warns
when the maps of a run are not comparable (other pixel size or landscape); keeps every class with --mode classes (and the
landscape level, the same codes of the classes in every map, a class 0 apart from the cells without polygon, and
RPR relative to the classes of the run); gives a date without habitat its row (no area, no patches); writes valid
polygons for patches of 8 neighbours; keeps the pixel of a raster in a run with polygons; leaves no folder .tmp when a
run fails; and that similarity.py gives the metrics of lsm.py for its similar areas (categories 1 and 2, without
the HFT, 3, and the cells without data, 255).

    cd segm_optim && python -s tests/test_landscape_metrics.py
"""

import csv
import json
import os
import sys
import tempfile

import numpy as np
from osgeo import gdal, ogr, osr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import landscape_metrics as lm  # noqa: E402
import lsm  # noqa: E402
import similarity  # noqa: E402

gdal.UseExceptions()
ogr.UseExceptions()
X0, Y1, RES = 440000.0, 4110000.0, 10.0


def utm():
	s = osr.SpatialReference()
	s.ImportFromEPSG(32630)
	return s


def habitat_map(seed, cuts):
	"""A 60 x 80 map: 3 forest (habitat), 2 scrub, 1 crop; CUTS lines of crop that split the forest."""
	from scipy import ndimage
	rng = np.random.default_rng(seed)
	z = ndimage.gaussian_filter(rng.normal(size=(60, 80)), 5)
	a = np.where(z > np.quantile(z, 0.5), 3, np.where(z > np.quantile(z, 0.2), 2, 1)).astype(np.uint8)
	for i in range(cuts):
		a[:, 10 + 15 * i] = 1
	return a


def write_tif(path, a):
	ds = gdal.GetDriverByName('GTiff').Create(path, a.shape[1], a.shape[0], 1, gdal.GDT_Byte)
	ds.SetGeoTransform((X0, RES, 0, Y1, 0, -RES))
	ds.SetProjection(utm().ExportToWkt())
	ds.GetRasterBand(1).SetNoDataValue(255)
	ds.GetRasterBand(1).WriteArray(a)
	ds = None


def write_polygons(path, a):
	src = gdal.GetDriverByName('MEM').Create('', a.shape[1], a.shape[0], 1, gdal.GDT_Byte)
	src.SetGeoTransform((X0, RES, 0, Y1, 0, -RES))
	src.SetProjection(utm().ExportToWkt())
	src.GetRasterBand(1).WriteArray(a)
	ds = gdal.GetDriverByName('GPKG').Create(path, 0, 0, 0, gdal.GDT_Unknown)
	lyr = ds.CreateLayer('cover', srs=utm(), geom_type=ogr.wkbPolygon)
	lyr.CreateField(ogr.FieldDefn('code', ogr.OFTInteger))
	gdal.Polygonize(src.GetRasterBand(1), None, lyr, 0)
	lyr.CreateField(ogr.FieldDefn('cover', ogr.OFTString))
	names = {1: 'crop', 2: 'scrub', 3: 'forest'}
	lyr.StartTransaction()
	for f in lyr:
		f.SetField('cover', names[f.GetField('code')])
		lyr.SetFeature(f)
	lyr.CommitTransaction()
	ds = None


def write_boundary(path):
	ds = gdal.GetDriverByName('GPKG').Create(path, 0, 0, 0, gdal.GDT_Unknown)
	lyr = ds.CreateLayer('b', srs=utm(), geom_type=ogr.wkbPolygon)
	f = ogr.Feature(lyr.GetLayerDefn())  # a rectangle on the lines of the pixels, inside the maps
	f.SetGeometry(ogr.CreateGeometryFromWkt(f'POLYGON (({X0 + 50} {Y1 - 40}, {X0 + 700} {Y1 - 40}, {X0 + 700} {Y1 - 560}, '
	                                        f'{X0 + 50} {Y1 - 560}, {X0 + 50} {Y1 - 40}))'))
	lyr.CreateFeature(f)
	ds = None


def table(path):
	with open(path, newline='', encoding='utf-8') as f:
		return list(csv.DictReader(f, delimiter=' '))


def close(a, b):
	return abs(float(a) - float(b)) <= 1e-6 * max(1.0, abs(float(a)))


def test_runs():
	with tempfile.TemporaryDirectory() as d:
		lm.OUT = os.path.join(d, 'landscape')
		maps = {2000: habitat_map(1, 0), 2020: habitat_map(1, 3)}
		for y, a in maps.items():
			write_tif(os.path.join(d, f'cover_{y}.tif'), a)
			write_polygons(os.path.join(d, f'cover_{y}.gpkg'), a)
		bnd = os.path.join(d, 'boundary.gpkg')
		write_boundary(bnd)
		common = {'mode': 'habitat', 'boundary': bnd, 'directions': 8, 'edge_depth': 1, 'metrics': None}
		ras = lm.run(dict(common, name='ras', inputs=[{'file': os.path.join(d, f'cover_{y}.tif'), 'label': str(y),
		                                                'habitat': [3]} for y in maps]), plots=False)
		vec = lm.run(dict(common, name='vec', res=RES, inputs=[{'file': os.path.join(d, f'cover_{y}.gpkg'), 'label': str(y),
		                                                         'field': 'cover', 'habitat': ['forest']} for y in maps]),
		             plots=False)
		for f in ('metrics_class.csv', 'metrics_patch.csv', 'patches_1.gpkg', 'patches_2.gpkg', 'classes_1.tif',
		          'landscape_run.json', 'landscape_metadata.txt'):
			assert os.path.isfile(os.path.join(ras, f)), f
		assert not os.path.isfile(os.path.join(ras, 'metrics_landscape.csv'))  # mode habitat: no landscape level
		cr, cv = table(os.path.join(ras, 'metrics_class.csv')), table(os.path.join(vec, 'metrics_class.csv'))
		assert [r['map'] for r in cr] == ['2000', '2020'] and all(r['class_name'] == 'habitat' for r in cr)
		for a, b in zip(cr, cv):  # the same habitat as a raster and as polygons, on the same grid
			for k in ('ca', 'pland', 'np', 'ed', 'mesh', 'enn_mn', 'cohesion'):
				assert close(a[k], b[k]), (k, a[k], b[k])
		assert int(cr[1]['np']) > int(cr[0]['np']) and float(cr[1]['mesh']) < float(cr[0]['mesh'])  # split by the lines
		# the polygons of the patches: one per patch, with their metrics
		rec = json.load(open(os.path.join(ras, 'landscape_run.json')))
		ds = ogr.Open(os.path.join(ras, 'patches_2.gpkg'))
		lyr = ds.GetLayer(0)
		assert lyr.GetFeatureCount() == rec['maps'][1]['patches'] == int(cr[1]['np'])
		total = sum(f.GetField('area') for f in lyr)
		assert close(total, cr[1]['ca'])
		ds = None
		# boundary: the landscape is the rectangle of 65 x 52 pixels
		assert rec['maps'][0]['cells'] == 65 * 52 and not rec['warnings']
		# not comparable: another pixel size and no boundary
		odd = lm.run({'name': 'odd', 'mode': 'habitat', 'res': 20, 'inputs': [
			{'file': os.path.join(d, 'cover_2000.tif'), 'label': 'a', 'habitat': [3]},
			{'file': os.path.join(d, 'cover_2020.gpkg'), 'label': 'b', 'field': 'cover', 'habitat': ['forest']}]}, plots=False)
		w = json.load(open(os.path.join(odd, 'landscape_run.json')))['warnings']
		assert not any('pixel sizes' in x for x in w)  # both at 20 m
		odd = lm.run({'name': 'odd2', 'mode': 'habitat', 'inputs': [
			{'file': os.path.join(d, 'cover_2000.tif'), 'label': 'a', 'habitat': [3]},
			{'file': os.path.join(d, 'cover_2020.tif'), 'label': 'b', 'habitat': [3]}], 'res': None,
			'boundary': None}, plots=False)
		assert not json.load(open(os.path.join(odd, 'landscape_run.json')))['warnings']
		# every class, with the landscape level
		cls = lm.run({'name': 'cls', 'mode': 'classes', 'boundary': bnd, 'metrics': ['c_pland', 'l_shdi', 'l_contag'],
		              'inputs': [{'file': os.path.join(d, 'cover_2000.gpkg'), 'label': '2000', 'field': 'cover'}],
		              'res': RES}, plots=False)
		rows = table(os.path.join(cls, 'metrics_class.csv'))
		assert {r['class_name'] for r in rows} == {'crop', 'scrub', 'forest'}
		assert close(sum(float(r['pland']) for r in rows), 100)
		land = table(os.path.join(cls, 'metrics_landscape.csv'))
		assert len(land) == 1 and float(land[0]['shdi']) > 0
		# errors of the user
		for bad, why in (({'name': 'x y'}, 'name'), ({'name': 'x', 'inputs': []}, 'maps'),
		                 ({'name': 'x', 'inputs': [{'file': os.path.join(d, 'cover_2000.tif'), 'label': 'a'}]}, 'habitat')):
			try:
				lm.run(bad, plots=False)
			except ValueError as e:
				assert why in str(e), str(e)
			else:
				raise AssertionError(f'no error for {bad}')
	print('runs: ok')


def write_coded(path, rows, field_type=ogr.OFTString, column='cover'):
	"""Polygons with a COLUMN of the values ROWS: [(value, column of the 10 m pixels from 0 to 60, row 0 to 52)] as
	squares of one pixel."""
	ds = gdal.GetDriverByName('GPKG').Create(path, 0, 0, 0, gdal.GDT_Unknown)
	lyr = ds.CreateLayer('p', srs=utm(), geom_type=ogr.wkbPolygon)
	lyr.CreateField(ogr.FieldDefn(column, field_type))
	for v, c, r in rows:
		x, y = X0 + 50 + c * RES, Y1 - 40 - r * RES
		f = ogr.Feature(lyr.GetLayerDefn())
		f.SetField(column, v)
		f.SetGeometry(ogr.CreateGeometryFromWkt(f'POLYGON (({x} {y}, {x + RES} {y}, {x + RES} {y - RES}, {x} {y - RES}, {x} {y}))'))
		lyr.CreateFeature(f)
	ds = None


def test_review_cases():
	with tempfile.TemporaryDirectory() as d:
		lm.OUT = os.path.join(d, 'landscape')
		bnd = os.path.join(d, 'boundary.gpkg')
		write_boundary(bnd)
		# text classes: the same code for the same value at every date, whatever values each date has
		write_coded(os.path.join(d, 'a.gpkg'), [('forest', 1, 1), ('grass', 5, 5)])
		write_coded(os.path.join(d, 'b.gpkg'), [('crop', 1, 1), ('forest', 5, 5), ('grass', 9, 9)])
		out = lm.run({'name': 'codes', 'mode': 'classes', 'boundary': bnd, 'res': RES, 'metrics': ['c_ca', 'l_rpr', 'l_pr'],
		              'inputs': [{'file': os.path.join(d, 'a.gpkg'), 'label': 'a', 'field': 'cover'},
		                         {'file': os.path.join(d, 'b.gpkg'), 'label': 'b', 'field': 'cover'}]}, plots=False)
		code = {}
		for r in table(os.path.join(out, 'metrics_class.csv')):
			assert code.setdefault(r['class_name'], r['class']) == r['class'], r
		assert code['(no polygon)'] == '0' and len(set(code.values())) == 4
		land = table(os.path.join(out, 'metrics_landscape.csv'))  # RPR: classes of each map / those of the run (4)
		assert [float(x['rpr']) for x in land] == [75.0, 100.0], land
		# ... also when each date has its own column (as Code_00 and Code_18 of CORINE)
		write_coded(os.path.join(d, 'c00.gpkg'), [('311', 1, 1), ('312', 5, 5)], column='Code_00')
		write_coded(os.path.join(d, 'c18.gpkg'), [('211', 1, 1), ('311', 5, 5)], column='Code_18')
		out = lm.run({'name': 'corine', 'mode': 'classes', 'boundary': bnd, 'res': RES, 'metrics': ['c_ca'],
		              'inputs': [{'file': os.path.join(d, 'c00.gpkg'), 'label': '2000', 'field': 'Code_00'},
		                         {'file': os.path.join(d, 'c18.gpkg'), 'label': '2018', 'field': 'Code_18'}]}, plots=False)
		code = {}
		for r in table(os.path.join(out, 'metrics_class.csv')):
			assert code.setdefault(r['class_name'], r['class']) == r['class'], r
		# whole numbers with a class 0: it is not the cells without polygon
		write_coded(os.path.join(d, 'z.gpkg'), [(0, 1, 1), (2, 5, 5)], ogr.OFTInteger)
		out = lm.run({'name': 'zero', 'mode': 'classes', 'boundary': bnd, 'res': RES, 'metrics': ['c_ca'],
		              'inputs': [{'file': os.path.join(d, 'z.gpkg'), 'label': 'z', 'field': 'cover'}]}, plots=False)
		rows = {r['class_name']: r for r in table(os.path.join(out, 'metrics_class.csv'))}
		assert set(rows) == {'0', '2', '(no polygon)'} and rows['(no polygon)']['class'] == '-1', rows
		assert close(rows['0']['ca'], 0.01)
		# a date without habitat: a row with no area and no patches
		out = lm.run({'name': 'gone', 'mode': 'habitat', 'boundary': bnd, 'res': RES, 'inputs': [
			{'file': os.path.join(d, 'a.gpkg'), 'label': '2000', 'field': 'cover', 'habitat': ['forest']},
			{'file': os.path.join(d, 'a.gpkg'), 'label': '2020', 'field': 'cover', 'habitat': ['crop']}]}, plots=False)
		rows = table(os.path.join(out, 'metrics_class.csv'))
		assert [r['map'] for r in rows] == ['2000', '2020'] and rows[1]['ca'] == '0' and rows[1]['np'] == '0'
		assert rows[1]['enn_mn'] == 'NA'
		# patches of 8 neighbours that touch at corners, with holes: valid multipolygons, one per patch
		a = np.zeros((60, 80), np.uint8)
		a[10:50, 10:50] = 3
		a[20:30, 20:30] = 1  # a hole
		for i in range(0, 20, 2):
			a[52 + (i // 2) % 6, 55 + i // 2] = 3  # a diagonal line of cells
		write_tif(os.path.join(d, 'diag.tif'), a)
		out = lm.run({'name': 'diag', 'inputs': [{'file': os.path.join(d, 'diag.tif'), 'label': 'x', 'habitat': [3]}]},
		             plots=False)
		ds = ogr.Open(os.path.join(out, 'patches_1.gpkg'))
		lyr = ds.GetLayer(0)
		np_ = int(table(os.path.join(out, 'metrics_class.csv'))[0]['np'])
		assert lyr.GetFeatureCount() == np_
		for f in lyr:
			g = f.GetGeometryRef()
			assert g.IsValid() and g.GetGeometryType() == ogr.wkbMultiPolygon
			assert close(g.GetArea() / 1e4, f.GetField('area'))
		ds = None
		# a raster in a run with polygons keeps its pixel unless one is given
		mixed = lm.run({'name': 'mixed', 'inputs': [{'file': os.path.join(d, 'diag.tif'), 'label': 'r', 'habitat': [3]},
		                                            {'file': os.path.join(d, 'a.gpkg'), 'label': 'v'}]}, plots=False)
		assert json.load(open(os.path.join(mixed, 'landscape_run.json')))['maps'][0]['grid']['res_m'] == RES
		# a run that fails leaves no folder .tmp, and the previous results of its name
		try:
			lm.run({'name': 'diag', 'inputs': [{'file': os.path.join(d, 'diag.tif'), 'label': 'x', 'habitat': [3]},
			                                   {'file': bnd, 'label': 'y', 'field': 'nothing', 'habitat': ['1']}]},
			       plots=False)
		except (ValueError, RuntimeError):
			pass
		else:
			raise AssertionError('a column that does not exist ran')
		assert not os.path.exists(os.path.join(lm.OUT, 'diag.tmp'))
		assert os.path.isfile(os.path.join(lm.OUT, 'diag', 'landscape_run.json'))
	print('review cases: ok')


def test_similar_metrics():
	rng = np.random.default_rng(3)
	cat = rng.choice([0, 1, 2], size=(40, 50), p=[0.6, 0.3, 0.1]).astype(np.uint8)
	cat[5:9, 5:9] = 3  # the HFT
	cat[:, :3] = similarity.CAT_NODATA
	got = similarity.similar_metrics(cat, 100.0, 8)
	valid = (cat != similarity.CAT_NODATA) & (cat != 3)
	ref = lsm.compute(((cat == 1) | (cat == 2)).astype(int), valid, 100.0, directions=8)
	c1 = next(r for r in ref['class'] if r['class'] == 1)
	for col, m in similarity.SIM_LSM:
		assert (got[col] is None and c1[m] is None) or close(got[col], c1[m]), (col, got[col], c1[m])
	none = similarity.similar_metrics(np.where(cat == 3, 3, 0).astype(np.uint8), 100.0, 8)
	assert none['PLAND'] == 0 and none['PD'] == 0 and none['MESH'] == 0 and none['AI'] is None
	print('similar areas: ok')


if __name__ == '__main__':
	test_runs()
	test_review_cases()
	test_similar_metrics()
	print('All the tests passed')
