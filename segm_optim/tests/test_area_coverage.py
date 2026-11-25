#### eHabitat+: test of the part of the study areas that the variables cover (area_coverage.py, inputs.coverage)
#### License: CC BY-SA 3.0
"""
With synthetic GeoTIFFs (UTM 30N, 10 m) and study areas, checks that area_coverage.study_coverage gives the part of
each study area inside the extent of all the variables and with data in the first one (a study area fully covered,
partly outside, outside, with NoData inside, smaller than a pixel of the reduced reading, in another CRS, several
features with the same ID, only the IDs to process; NaN as NoData with overviews, a raster south-up, features without
ID, an ID column of text, many parts under one ID), the notes of coverage_checks (an error with no data, a warning
when partly covered, nothing when covered), and inputs.coverage of the app on a set of inVars.

    cd segm_optim && python -s tests/test_area_coverage.py
"""

import os
import shutil
import sys
import tempfile

import numpy as np
from osgeo import gdal, ogr, osr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'app'))  # ehabitat_app

import area_coverage  # noqa: E402

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

X0, Y0 = 450000.0, 4110000.0  # top left corner of the variables
RES = 10.0
SIZE = 2000  # pixels a side: 20 km; read reduced (MAX_SIZE 1000)


def utm():
	s = osr.SpatialReference()
	s.ImportFromEPSG(32630)
	s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return s


def raster(path, data, x0=X0, y0=Y0, nodata=-9999.0):
	ds = gdal.GetDriverByName('GTiff').Create(path, data.shape[1], data.shape[0], 1, gdal.GDT_Float32)
	ds.SetGeoTransform((x0, RES, 0, y0, 0, -RES))
	ds.SetProjection(utm().ExportToWkt())
	b = ds.GetRasterBand(1)
	if nodata is not None:
		b.SetNoDataValue(nodata)
	b.WriteArray(data)
	ds = None
	return path


def box(xmin, ymin, xmax, ymax):
	return ogr.CreateGeometryFromWkt(f'POLYGON(({xmin} {ymin},{xmax} {ymin},{xmax} {ymax},{xmin} {ymax},{xmin} {ymin}))')


def vector(path, features, srs=None):
	"""FEATURES: [(id, name, geometry)] in SRS (UTM 30N by default)."""
	ds = ogr.GetDriverByName('GPKG').CreateDataSource(path)
	lyr = ds.CreateLayer('areas', srs or utm(), ogr.wkbMultiPolygon)
	lyr.CreateField(ogr.FieldDefn('pa_id', ogr.OFTInteger))
	lyr.CreateField(ogr.FieldDefn('name', ogr.OFTString))
	for fid, name, g in features:
		f = ogr.Feature(lyr.GetLayerDefn())
		if fid is not None:
			f.SetField('pa_id', fid)
		f.SetField('name', name)
		f.SetGeometry(ogr.ForceToMultiPolygon(g))
		lyr.CreateFeature(f)
	ds = None
	return path


def close(a, b, tol=0.01):
	return abs(a - b) <= tol


def test_study_coverage():
	tmp = tempfile.mkdtemp()
	try:
		data = np.ones((SIZE, SIZE), np.float32)
		data[:, 1500:] = -9999  # the last 5 km (east) without data in the first variable
		first = raster(os.path.join(tmp, 'a.tif'), data)
		# a second variable covering only the north half (10 km): the extent of all of them is that half
		second = raster(os.path.join(tmp, 'b.tif'), np.ones((SIZE // 2, SIZE), np.float32))
		x1, ymid, ymin = X0 + SIZE * RES, Y0 - SIZE // 2 * RES, Y0 - SIZE * RES
		feats = [
			(1, 'Covered', box(X0 + 1000, ymid + 1000, X0 + 9000, Y0 - 1000)),
			(2, 'Half outside', box(X0 + 1000, ymid - 4000, X0 + 9000, ymid + 4000)),  # half in the south half
			(3, 'Outside', box(X0 + 1000, ymin - 5000, X0 + 9000, ymin - 1000)),
			(4, 'Half NoData', box(x1 - 9000, ymid + 1000, x1 - 1000, Y0 - 1000)),  # 4 of its 8 km in the NoData
			(5, 'Tiny', box(X0 + 5000, Y0 - 3000, X0 + 5005, Y0 - 2995)),  # 5 x 5 m
			(6, 'Two parts', box(X0 + 1000, ymid + 500, X0 + 2000, ymid + 1500)),
			(6, None, box(X0 + 1000, ymin + 500, X0 + 2000, ymin + 1500)),  # its second part, in the south half
			(7, 'Tiny in NoData', box(x1 - 2000, Y0 - 3000, x1 - 1995, Y0 - 2995)),
		]
		vec = vector(os.path.join(tmp, 'areas.gpkg'), feats)
		cov = {a['id']: a for a in area_coverage.study_coverage([first, second], vec, None, 'pa_id', None, 'name')}
		assert sorted(cov) == [1, 2, 3, 4, 5, 6, 7], cov
		assert close(cov[1]['covered'], 1) and close(cov[1]['inside'], 1), cov[1]
		assert close(cov[2]['inside'], 0.5) and close(cov[2]['covered'], 0.5), cov[2]
		assert cov[3]['inside'] == 0 and cov[3]['covered'] == 0, cov[3]
		assert close(cov[4]['inside'], 1) and close(cov[4]['covered'], 0.5), cov[4]
		assert cov[5]['covered'] == 1, cov[5]
		assert close(cov[6]['covered'], 0.5) and cov[6]['label'] == 'Two parts', cov[6]
		assert cov[7]['inside'] == 1 and cov[7]['covered'] == 0, cov[7]
		assert close(cov[1]['km2'], 64, 1e-6), cov[1]
		# only the IDs to process
		assert [a['id'] for a in area_coverage.study_coverage([first], vec, None, 'pa_id', [4, 1])] == [1, 4]
		# without the column of the IDs: 1, 2, ... as GRASS (cat); an unknown column is an error
		assert len(area_coverage.study_coverage([first], vec, None, 'cat')) == len(feats)
		try:
			area_coverage.study_coverage([first], vec, None, 'nope')
		except ValueError as e:
			assert 'nope' in str(e)
		else:
			raise AssertionError('an unknown column of the IDs was not refused')

		# a study area in geographic coordinates (WGS 84) is reprojected to the CRS of the variables
		wgs = osr.SpatialReference()
		wgs.ImportFromEPSG(4326)
		wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		g = feats[1][2].Clone()  # half outside the north half
		g.Transform(osr.CoordinateTransformation(utm(), wgs))
		geo = vector(os.path.join(tmp, 'geo.gpkg'), [(1, 'Geo', g)], wgs)
		a = area_coverage.study_coverage([first, second], geo, None, 'pa_id')[0]
		assert close(a['covered'], 0.5), a

		# the notes
		notes = area_coverage.coverage_checks(list(cov.values()))
		levels = [(n['level'], n['text']) for n in notes]
		assert sum(1 for lv, _ in levels if lv == 'error') == 2, levels  # 3 and 7
		assert any('study area 3 (Outside)' in t and 'outside them' in t for _, t in levels), levels
		assert any(lv == 'warning' and 'about 50 % of study area 2 (Half outside)' in t for lv, t in levels), levels
		assert not any('study area 1 ' in t or 'study area 5 ' in t for _, t in levels), levels
		assert area_coverage.coverage_checks([cov[1], cov[5]]) == []
		many = [dict(cov[2], id=i) for i in range(20)]
		notes = area_coverage.coverage_checks(many)
		assert len(notes) == area_coverage.MAX_NOTES + 1 and '12 more study areas' in notes[-1]['text'], notes[-1]
		assert area_coverage.coverage_checks([dict(cov[2], covered=0.004)])[0]['text'].startswith(
			'The variables cover about less than 1 %'), 'below 1 %'
		print('test_study_coverage: OK')
	finally:
		shutil.rmtree(tmp)


def test_edge_cases():
	"""NaN as NoData with overviews (their averages fill the holes: not read), a raster south-up (not checked), an ID
	NULL (left out, as GRASS), an ID column of text (refused), many parts under one ID (fast)."""
	import time
	tmp = tempfile.mkdtemp()
	try:
		data = np.ones((3000, 3000), np.float32)
		data[np.random.default_rng(0).random(data.shape) < 0.5] = np.nan  # half of the pixels without data
		path = raster(os.path.join(tmp, 'nan.tif'), data, nodata=float('nan'))
		ds = gdal.Open(path, gdal.GA_Update)
		ds.BuildOverviews('AVERAGE', [2, 4, 8])
		ds = None
		feats = [(1, 'Stripes', box(X0 + 1000, Y0 - 20000, X0 + 20000, Y0 - 1000)), (None, 'No ID', box(X0, Y0 - 500, X0 + 500, Y0))]
		vec = vector(os.path.join(tmp, 'a.gpkg'), feats)
		cov = area_coverage.study_coverage([path], vec, None, 'pa_id')
		assert [a['id'] for a in cov] == [1], cov  # the feature without ID is left out
		assert close(cov[0]['covered'], 0.5, 0.02), cov
		# an ID column of text
		ds = ogr.GetDriverByName('GPKG').CreateDataSource(os.path.join(tmp, 't.gpkg'))
		lyr = ds.CreateLayer('t', utm(), ogr.wkbPolygon)
		lyr.CreateField(ogr.FieldDefn('code', ogr.OFTString))
		f = ogr.Feature(lyr.GetLayerDefn())
		f.SetField('code', 'A')
		f.SetGeometry(box(X0 + 100, Y0 - 900, X0 + 900, Y0 - 100))
		lyr.CreateFeature(f)
		ds = None
		try:
			area_coverage.study_coverage([path], os.path.join(tmp, 't.gpkg'), None, 'code')
		except ValueError as e:
			assert 'not of integers' in str(e), e
		else:
			raise AssertionError('an ID column of text was not refused')
		# a raster south-up: not checked (no false error)
		south = os.path.join(tmp, 'south.tif')
		d = gdal.GetDriverByName('GTiff').Create(south, 100, 100, 1, gdal.GDT_Float32)
		d.SetGeoTransform((X0, RES, 0, Y0 - 1000, 0, RES))
		d.SetProjection(utm().ExportToWkt())
		d.GetRasterBand(1).Fill(1)
		d = None
		assert area_coverage.study_coverage([south], vec, None, 'pa_id') == []
		# 2000 parts under one ID (a fragmented study area): joined at once
		parts = [(7, 'Many', box(X0 + 200 * (k % 50) + 10, Y0 - 200 * (k // 50) - 190, X0 + 200 * (k % 50) + 190,
		                         Y0 - 200 * (k // 50) - 10)) for k in range(2000)]
		many = vector(os.path.join(tmp, 'many.gpkg'), parts)
		t = time.time()
		cov = area_coverage.study_coverage([path], many, None, 'pa_id')
		assert time.time() - t < 5, time.time() - t
		assert len(cov) == 1 and close(cov[0]['covered'], 0.5, 0.05), cov
		print('test_edge_cases: OK')
	finally:
		shutil.rmtree(tmp)


def test_app_coverage():
	"""inputs.coverage: the notes of a set of inVars for the study areas of the app; none for variables that setup.py
	would refuse (another CRS)."""
	from ehabitat_app import inputs, paths
	tmp = tempfile.mkdtemp()
	real = paths.INVARS
	try:
		paths.INVARS = os.path.join(tmp, 'inVars')
		os.makedirs(os.path.join(paths.INVARS, 'set'))
		raster(os.path.join(paths.INVARS, 'set', 'a.tif'), np.ones((100, 100), np.float32))
		vec = vector(os.path.join(tmp, 'areas.gpkg'), [(1, 'In', box(X0 + 100, Y0 - 900, X0 + 900, Y0 - 100)),
		                                               (2, 'Half', box(X0 + 500, Y0 - 900, X0 + 1500, Y0 - 100))])
		sa = {'file': vec, 'layer': None, 'col_id': 'pa_id', 'ids': None}
		notes = inputs.coverage('set', sa)
		assert len(notes) == 1 and notes[0]['level'] == 'warning' and 'study area 2 (Half)' in notes[0]['text'], notes
		assert inputs.coverage('set', dict(sa, ids=[1])) == []
		assert inputs.coverage('set', None) == []
		assert inputs.coverage('set', dict(sa, col_id='nope'))[0]['level'] == 'info'
		other = osr.SpatialReference()
		other.ImportFromEPSG(25830)
		ds = gdal.GetDriverByName('GTiff').Create(os.path.join(paths.INVARS, 'set', 'b.tif'), 10, 10, 1, gdal.GDT_Float32)
		ds.SetGeoTransform((X0, RES, 0, Y0, 0, -RES))
		ds.SetProjection(other.ExportToWkt())
		ds = None
		assert inputs.coverage('set', sa) == []  # check_scale gives the error of the CRS
		print('test_app_coverage: OK')
	finally:
		paths.INVARS = real
		shutil.rmtree(tmp)


if __name__ == '__main__':
	test_study_coverage()
	test_edge_cases()
	test_app_coverage()
	print('All the tests passed')
