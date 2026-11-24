#### eHabitat+ HFT model: test of the similar patches as polygons (similarity.py)
#### License: CC BY-SA 3.0
"""
Checks that similarity.py writes the similar patches of every HFT or subclass (categories 1 and 2 of
similar_<id>_hft<h>.tif) as polygons, one per patch, in similar_<id>_hft_<h>.gpkg or similar_<id>_hft_<h>_sub_<s>.gpkg,
with their fields, and that the patches follow ECOSIM_NEIGHBOURS (8 or 4 neighbours).

    cd segm_optim && python -s tests/test_similar_patches.py
"""

import os
import sys
import tempfile

import numpy as np
from osgeo import ogr
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

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
	print('test_similar_patches: OK')


if __name__ == '__main__':
	main()
