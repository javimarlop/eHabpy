#### eHabitat+ HFT model: test of the shapefile of the HFTs with the subclasses and of the similarity of the subclasses
#### License: CC BY-SA 3.0
"""
Checks the shapefile of the HFTs (park_segm_*_class.shp) that moranvar_plots.R writes, and that similarity.py reads:

- its fields: segm_id, segm_cat, area_id, hft, area_ha, silhouette and the variables (named by var_fields), with
  park_segm_*_class_fields.csv next to it;
- with HFT_SUBCLASSES, its segments split into the subclasses of their pixels (subclass_pieces of subclasses.R): one
  polygon per segment and subclass, even when its pixels do not touch, with the label hft_<n>_sub_<y>, the code
  100 x HFT + subclass (sub_code) and its area; the pixels without a subclass in a polygon without it;
- that similarity.py draws the HFTs (hft, or hclst_m in the shapefiles of eHabpy 1.6.6 and before) or the subclasses
  (sub_code, ECOSIM_CLASSES = 'subclass') on its grid, and refuses the subclasses of a shapefile without them;
- that the app tells the runs of the subclasses (similarity_<res>m[_<area>]_subclasses) from those of another area.

The part of R needs Rscript with sf and terra (skipped without them).

    cd segm_optim && python -s tests/test_class_shapefile.py
"""

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from osgeo import gdal, ogr, osr

HERE = os.path.dirname(os.path.abspath(__file__))
SEGM = os.path.dirname(HERE)
sys.path.insert(0, SEGM)
sys.path.insert(0, os.path.join(os.path.dirname(SEGM), 'app'))  # ehabitat_app

import similarity  # noqa: E402

gdal.UseExceptions()
ogr.UseExceptions()

# A 6 x 4 grid of 10 m pixels: two segments (columns 0-2 and 3-5) of the HFTs 1 and 2. The subclasses of the pixels:
# 101 and 102 in the first segment (102 in two pixels that do not touch), 201 in the second one but for a pixel
# without data (NA)
R_SCRIPT = r'''
suppressPackageStartupMessages({library(sf); library(terra)})
source('plot_utils.R')
src <- readLines('moranvar_plots.R')
eval(parse(text=src[grep('^dbf_names <- ', src):(grep('^# Colores jer', src) - 1)]))
source('subclasses.R')
out <- commandArgs(trailingOnly=TRUE)[1]
g <- rast(nrows=4, ncols=6, xmin=500000, xmax=500060, ymin=4100000, ymax=4100040, crs='EPSG:32630')
seg <- rast(g); values(seg) <- rep(c(1, 1, 1, 2, 2, 2), 4)
segs <- st_as_sf(as.polygons(seg, dissolve=TRUE))
names(segs)[1] <- 'cat'
segs$segm_id <- c('1100001', '1100002')
segs$hclust_mean <- c(1L, 2L)
vn <- c('gee_x_1_NDVI_S2_p95', 'gee_x_2_LSWI_S2_median')
vars <- setNames(data.frame(c(0, 1), c(1, 0)), var_fields(vn))
layer <- class_layer(segs, '1', vars, c(0.5, 0.25))
write_class_layer(layer, file.path(out, 'park_segm_1_1_p10_class.shp'), setNames(vn, names(vars)))
sc <- matrix(c(101, 101, 102, 201, 201, 201,
               101, 101, 101, 201,  NA, 201,
               102, 101, 101, 201, 201, 201,
               101, 101, 101, 201, 201, 201), nrow=4, byrow=TRUE)
pieces <- subclass_pieces(layer, g, as.vector(t(sc)))
write_class_layer(pieces, file.path(out, 'park_segm_1_1_p25_class.shp'), setNames(vn, names(vars)))
'''


def _features(path):
	ds = ogr.Open(path)
	lyr = ds.GetLayer(0)
	d = lyr.GetLayerDefn()
	fields = [d.GetFieldDefn(i).GetName() for i in range(d.GetFieldCount())]
	rows = []
	for f in lyr:
		r = f.items()
		r['_area'] = f.GetGeometryRef().GetArea()
		r['_parts'] = f.GetGeometryRef().GetGeometryCount() if f.GetGeometryRef().GetGeometryType() in (
			ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D) else 1
		rows.append(r)
	return fields, rows


def _write_shp(path, field, values, field_type=ogr.OFTInteger):
	"""Two polygons, the halves of the 6 x 4 grid of R_SCRIPT, with FIELD = VALUES."""
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(32630)
	ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(path)
	lyr = ds.CreateLayer('x', srs=srs, geom_type=ogr.wkbPolygon)
	lyr.CreateField(ogr.FieldDefn(field, field_type))
	for i, v in enumerate(values):
		x0 = 500000 + 30 * i
		g = ogr.CreateGeometryFromWkt(f'POLYGON (({x0} 4100000, {x0 + 30} 4100000, {x0 + 30} 4100040, {x0} 4100040, '
		                              f'{x0} 4100000))')
		f = ogr.Feature(lyr.GetLayerDefn())
		f.SetField(field, v)
		f.SetGeometry(g)
		lyr.CreateFeature(f)
	ds = None


def _grid():
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(32630)
	native = similarity.Grid((500000, 10, 0, 4100040, 0, -10), 6, 4, srs.ExportToWkt())
	return similarity.CoarseGrid(native, 10)


def test_r_shapefile():
	if not shutil.which('Rscript'):
		print('test_r_shapefile: skipped (no Rscript)')
		return
	tmp = tempfile.mkdtemp()
	try:
		r = subprocess.run(['Rscript', '-e', R_SCRIPT, tmp], cwd=SEGM, capture_output=True, text=True)
		if r.returncode != 0 and ('there is no package' in r.stderr or 'no hay paquete' in r.stderr):
			print('test_r_shapefile: skipped (R without sf or terra)')
			return
		assert r.returncode == 0, r.stderr
		# the segments
		fields, rows = _features(os.path.join(tmp, 'park_segm_1_1_p10_class.shp'))
		assert fields == ['segm_id', 'segm_cat', 'area_id', 'hft', 'area_ha', 'silhouette', 'NDVI_S2_p9', 'LSWI_S2_me'], fields
		assert [r['hft'] for r in rows] == [1, 2] and [r['segm_id'] for r in rows] == ['1100001', '1100002']
		assert all(abs(r['area_ha'] - 0.12) < 1e-9 for r in rows), rows
		with open(os.path.join(tmp, 'park_segm_1_1_p10_class_fields.csv')) as f:
			doc = f.read()
		assert '"hft" "HFT of the segment in this cut"' in doc and 'Mean of gee_x_2_LSWI_S2_median' in doc, doc
		# the segments split into the subclasses
		fields, rows = _features(os.path.join(tmp, 'park_segm_1_1_p25_class.shp'))
		assert fields[:8] == ['segm_id', 'segm_cat', 'area_id', 'hft', 'subclass', 'sub_code', 'area_ha',
		                      'silhouette'], fields
		got = {(r['segm_id'], r['sub_code']): r for r in rows}
		assert set(got) == {('1100001', 101), ('1100001', 102), ('1100002', 201), ('1100002', None)}, set(got)
		pix = {k: round(r['area_ha'] / 0.01) for k, r in got.items()}
		assert pix == {('1100001', 101): 10, ('1100001', 102): 2, ('1100002', 201): 11, ('1100002', None): 1}, pix
		assert all(abs(r['_area'] / 1e4 - r['area_ha']) < 1e-9 for r in rows)
		assert got[('1100001', 102)]['_parts'] == 2  # its two pixels do not touch: one polygon of two parts
		assert got[('1100001', 102)]['subclass'] == 'hft_1_sub_2' and got[('1100002', 201)]['subclass'] == 'hft_2_sub_1'
		assert got[('1100002', None)]['subclass'] is None and got[('1100002', None)]['hft'] == 2
		assert all(r['silhouette'] == (0.5 if r['hft'] == 1 else 0.25) for r in rows)
		assert all(r['NDVI_S2_p9'] == (0 if r['hft'] == 1 else 1) for r in rows)
		# what similarity.py takes from it: the HFTs or the subclasses, on its grid
		shp = os.path.join(tmp, 'park_segm_1_1_p25_class.shp')
		cg = _grid()
		cov = similarity.hft_cover(shp, cg, 10, 0.5, 'hft')
		assert cov['hfts'] == [1, 2] and (cov['assign'] == np.repeat([[1, 2]], 3, axis=1).repeat(4, axis=0)).all()
		cov = similarity.hft_cover(shp, cg, 10, 0.5, 'subclass')
		assert cov['hfts'] == [101, 102, 201], cov['hfts']
		assert int((cov['assign'] == 102).sum()) == 2 and int((cov['assign'] == 0).sum()) == 1
		areas = similarity.hft_areas(shp, 'subclass')
		assert {k: round(v / 0.01) for k, v in areas.items()} == {101: 10, 102: 2, 201: 11}, areas
		assert {k: round(v / 0.01) for k, v in similarity.hft_areas(shp).items()} == {1: 12, 2: 12}
		print('test_r_shapefile: OK')
	finally:
		shutil.rmtree(tmp, ignore_errors=True)


def test_old_shapefiles():
	tmp = tempfile.mkdtemp()
	try:
		old = os.path.join(tmp, 'park_segm_1_1_class.shp')
		_write_shp(old, 'hclst_m', [3, 5])  # as eHabpy 1.6.6 and before
		assert sorted(similarity.hft_areas(old)) == [3, 5]
		assert similarity.hft_cover(old, _grid(), 10, 0.5)['hfts'] == [3, 5]
		t = {'id': '1', 'shp': old, 'classes': 'subclass'}
		inp = {'ref': None, 'ref_note': None, 'area': {'folder': tmp}}
		for f in (lambda: similarity.hft_cover(old, _grid(), 10, 0.5, 'subclass'),
		          lambda: similarity.hft_areas(old, 'subclass'),
		          lambda: similarity.hft_values(t, inp, _grid(), _grid(), 0.5, 10)):
			try:
				f()
				raise AssertionError('a shapefile without subclasses was accepted')
			except ValueError as e:
				assert 'HFT_SUBCLASSES' in str(e), e
		print('test_old_shapefiles: OK')
	finally:
		shutil.rmtree(tmp, ignore_errors=True)


def test_options():
	cfg = {'STUDY_AREA': 'x', 'SCENARIO': 'y'}
	assert similarity._settings(cfg)['classes'] == 'hft'
	assert similarity._settings(cfg, classes='subclass')['classes'] == 'subclass'
	try:
		similarity._settings(cfg, classes='segments')
		raise AssertionError('ECOSIM_CLASSES = segments was accepted')
	except ValueError as e:
		assert 'ECOSIM_CLASSES' in str(e)
	t = {'dir': '/r/p10', 'classes': 'subclass'}
	cg = _grid()
	assert similarity._run_folder(t, cg, '') == '/r/p10/similarity_10m_subclasses'
	assert similarity._run_folder(dict(t, classes='hft'), cg, 'gee_x') == '/r/p10/similarity_10m_gee_x'
	assert similarity._run_folder(t, cg, 'gee_x') == '/r/p10/similarity_10m_gee_x_subclasses'
	assert similarity.RUN_DIR_RE.match('similarity_10m_gee_x_subclasses')
	assert similarity.class_label(1004, 'subclass') == 'subclass hft_10_sub_4' and similarity.class_label(7, 'hft') == 'HFT 7'
	from ehabitat_app import ecosim
	for name, res, tag, sub in (('similarity_100m', '100', None, None), ('similarity_100m_subclasses', '100', None, '_subclasses'),
	                            ('similarity_100m_gee_x', '100', 'gee_x', None),
	                            ('similarity_12_5m_gee_x_subclasses', '12_5', 'gee_x', '_subclasses')):
		m = ecosim.RUN_RE.match(name)
		assert m and (m.group('res'), m.group('tag'), m.group('sub')) == (res, tag, sub), name
	print('test_options: OK')


if __name__ == '__main__':
	test_r_shapefile()
	test_old_shapefiles()
	test_options()
	print('All the tests passed')
