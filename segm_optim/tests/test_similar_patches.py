#### eHabitat+ HFT model: test of the similar patches as polygons (similarity.py)
#### License: CC BY-SA 3.0
"""
Checks that similarity.py writes the similar patches of every HFT or subclass (categories 1 and 2 of
similar_<id>_hft<h>.tif) as polygons, one per patch, in similar_<id>_hft_<h>.gpkg or similar_<id>_hft_<h>_sub_<s>.gpkg,
with their fields, and that the patches follow ECOSIM_NEIGHBOURS (8 or 4 neighbours), valid and in their place also
for a large patch with holes and cells joined only at corners; and that the map All of the app
(ecosim.patches) draws them, from the GeoPackages or, in the runs made before them, from the maps of the similar areas.

    cd segm_optim && python -s tests/test_similar_patches.py
"""

import json
import os
import shutil
import sys
import tempfile

import numpy as np
from osgeo import ogr
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'app'))  # ehabitat_app

import similarity  # noqa: E402

WKT = 'PROJCS["WGS 84 / UTM zone 30N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],' \
      'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],' \
      'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],PARAMETER["scale_factor",0.9996],' \
      'PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1],AUTHORITY["EPSG","32630"]]'


class Grid:
	def __init__(self, w, h, res):
		self.width, self.height, self.res, self.wkt = w, h, res, WKT
		self.gt = (450000.0, res, 0.0, 4100000.0, 0.0, -res)

	srs = similarity.CoarseGrid.srs


# 0 not similar, 1 similar, 2 a large similar patch, 3 the HFT, 255 NoData. The patch of 2 has a pixel that only
# touches it diagonally (row 3, column 3): one patch with 8 neighbours, two with 4.
CAT = np.array([[1, 1, 0, 0, 0, 255],
                [0, 0, 0, 2, 2, 0],
                [3, 3, 0, 2, 2, 0],
                [3, 3, 0, 0, 0, 2],
                [0, 0, 0, 0, 0, 0]], dtype=np.uint8)


def read(path):
	ds = ogr.Open(path)
	lyr = ds.GetLayer(0)
	fields = [lyr.GetLayerDefn().GetFieldDefn(i).GetName() for i in range(lyr.GetLayerDefn().GetFieldCount())]
	feats = []
	for feat in lyr:
		g = feat.GetGeometryRef()
		assert g.GetGeometryType() == ogr.wkbMultiPolygon and g.IsValid(), g.ExportToWkt()
		feats.append({f: feat.GetField(f) for f in fields})
	return fields, sorted(feats, key=lambda f: (f['category'], f['pixels']))


def main():
	cg = Grid(CAT.shape[1], CAT.shape[0], 10.0)
	with tempfile.TemporaryDirectory() as d:
		ctx = {'classes': 'hft', 'structure': ndimage.generate_binary_structure(2, 2)}
		name = similarity.write_patches(d, CAT, cg, ctx, '555', 10)
		assert name == 'similar_555_hft_10.gpkg', name
		fields, feats = read(os.path.join(d, name))
		assert fields == ['category', 'id', 'hft', 'large', 'pixels', 'area_ha'], fields
		got = [(f['category'], f['large'], f['pixels'], f['area_ha'], f['hft'], f['id']) for f in feats]
		assert got == [(1, 0, 2, 0.02, 10, '555'), (2, 1, 5, 0.05, 10, '555')], got

		ctx = {'classes': 'subclass', 'structure': ndimage.generate_binary_structure(2, 1)}
		name = similarity.write_patches(d, CAT, cg, ctx, '555', 1004)
		assert name == 'similar_555_hft_10_sub_4.gpkg', name
		fields, feats = read(os.path.join(d, name))
		assert 'subclass' in fields
		got = [(f['category'], f['pixels'], f['hft'], f['subclass']) for f in feats]
		assert got == [(1, 2, 10, 'hft_10_sub_4'), (2, 1, 10, 'hft_10_sub_4'), (2, 4, 10, 'hft_10_sub_4')], got

		name = similarity.write_patches(d, np.zeros_like(CAT), cg, ctx, '555', 1005)  # nothing similar: no polygon
		assert read(os.path.join(d, name))[1] == []
	large_patches()
	app_patches(cg)
	print('test_similar_patches: OK')


def large_patches():
	"""A large patch with holes (one with a similar island) and cells joined only at corners, among NoData and the HFT:
	one valid multipolygon per patch, with the pixels and category of the patch, for 8 and 4 neighbours, also when the
	similar pixels are far from the corner of the grid (the polygons are drawn in their box)."""
	cat = np.zeros((120, 150), dtype=np.uint8)
	cat[:5, :] = similarity.CAT_NODATA
	cat[30:100, 40:130] = 2  # the large patch
	cat[40:60, 50:70] = 0  # a hole
	cat[45:55, 55:65] = 1  # an island of similar pixels in it
	cat[70:80, 90:100] = 3  # the HFT, inside the patch
	for i in range(12):  # holes that touch at corners, a diagonal line of similar cells and a staircase hole
		cat[85 + i % 6, 45 + 2 * i] = 0
		cat[105 + i % 6 if i < 6 else 110 - i % 6, 20 + i] = 1
	cg = Grid(cat.shape[1], cat.shape[0], 30.0)
	for nb in (8, 4):
		st = ndimage.generate_binary_structure(2, 2 if nb == 8 else 1)
		lab, n = ndimage.label((cat == 1) | (cat == 2), structure=st)
		want = sorted((int(cat[lab == i][0]), int((lab == i).sum())) for i in range(1, n + 1))
		with tempfile.TemporaryDirectory() as d:
			name = similarity.write_patches(d, cat, cg, {'classes': 'hft', 'structure': st}, '7', 3)
			ds = ogr.Open(os.path.join(d, name))
			got = []
			for f in ds.GetLayer(0):
				g = f.GetGeometryRef()
				assert g.GetGeometryType() == ogr.wkbMultiPolygon and g.IsValid(), (nb, g.ExportToWkt()[:200])
				assert abs(g.GetArea() - f.GetField('pixels') * 900) < 1e-6, nb
				assert f.GetField('large') == int(f.GetField('category') == 2)
				got.append((f.GetField('category'), f.GetField('pixels')))
				if f.GetField('category') == 2:  # the large patch where it is: columns 40-130, rows 30-100
					x0, y0 = cg.gt[0], cg.gt[3]
					assert g.GetEnvelope() == (x0 + 40 * 30, x0 + 130 * 30, y0 - 100 * 30, y0 - 30 * 30), g.GetEnvelope()
			ds = None
		assert sorted(got) == want, (nb, sorted(got)[:5], want[:5])


def app_patches(cg):
	from ehabitat_app import ecosim
	from osgeo import gdal
	os.makedirs(similarity.RESULTS, exist_ok=True)
	d = tempfile.mkdtemp(prefix='test_patches_', dir=similarity.RESULTS)
	try:
		ctx = {'classes': 'hft', 'structure': ndimage.generate_binary_structure(2, 2)}
		hfts = []
		for h, gpkg in ((10, True), (11, False)):  # 11: a run made before the GeoPackages
			tif = f'similar_555_hft{h}.tif'
			similarity.write_raster(os.path.join(d, tif), CAT, cg, gdal.GDT_Byte, similarity.CAT_NODATA, 'NEAREST')
			name = similarity.write_patches(d, CAT, cg, ctx, '555', h) if gpkg else None
			hfts.append({'id': '555', 'hft': h, 'status': 'ok', 'similar': tif, 'patches': name})
		hfts.append({'id': '555', 'hft': 12, 'status': 'too_few_pixels', 'similar': None, 'patches': None})
		with open(os.path.join(d, ecosim.RUN_JSON), 'w') as f:
			json.dump({'classes': 'hft', 'options': {'neighbours': 8}, 'hfts': hfts, 'shapefiles': {}}, f)
		r = ecosim.patches(d)
		assert r['summary'] == {'555_10': {'patches': 2, 'large': 1, 'ha': 0.07, 'source': 'gpkg'},
		                        '555_11': {'patches': 2, 'large': 1, 'ha': 0.07, 'source': 'tif'}}, r['summary']
		feats = r['patches']['features']
		assert sorted((f['properties']['key'], f['properties']['category']) for f in feats) == \
			[('555_10', 1), ('555_10', 2), ('555_11', 1), ('555_11', 2)]
		lon, lat = feats[0]['geometry']['coordinates'][0][0][0]  # WGS 84: near -3.5, 37
		assert -4 < lon < -3 and 36.5 < lat < 37.5, (lon, lat)
		assert not r['truncated'] and r['study']['features'] == []
		# only the largest patches beyond MAX_PATCHES; an unreadable file leaves out its class only
		with open(os.path.join(d, 'similar_555_hft_10.gpkg'), 'wb') as f:
			f.write(b'not a GeoPackage')
		os.utime(os.path.join(d, ecosim.RUN_JSON))  # a new run for the cache of patches()
		old, ecosim.MAX_PATCHES = ecosim.MAX_PATCHES, 1
		try:
			r = ecosim.patches(d)
		finally:
			ecosim.MAX_PATCHES = old
		assert r['summary']['555_10']['error'] and r['summary']['555_11']['patches'] == 2, r['summary']
		assert r['truncated'] and [f['properties']['category'] for f in r['patches']['features']] == [2]
	finally:
		shutil.rmtree(d)


if __name__ == '__main__':
	main()
