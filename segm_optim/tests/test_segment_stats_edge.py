#### eHabitat+ HFT model: more cases of the test of the statistics of the segments (ehab_segm_core.py)
#### License: CC BY-SA 3.0
"""
Cases that tests/test_segment_stats.py does not cover, against the same per-segment reference (old_segment):
- NoData of the variable None (no NoData value) or float('nan'), as setup.py writes them into ENV_VARS, with NaN gaps,
  for Float32 and Float64 values;
- a segment raster without segments, and one whose pixels outside the segments are a negative NoData value;
- a segment raster that goes past the right and bottom edges of a variable: the segments inside it are as before (the
  old code stopped with an error on the others);
- ehabitat() of the wrappers ehabitat1..9 of get_mean_median_var.py, whose threshold ECOR can differ from K;
- a segmentation without its table cat -> segm_id (made by an older eHabpy): ehabitat() stops with a message.

    cd segm_optim && python -s tests/test_segment_stats_edge.py
"""

import os
import sys
import tempfile

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim
sys.path.insert(0, HERE)

from osgeo import gdal  # noqa: E402
import ehab_segm_core as core  # noqa: E402
import test_segment_stats as t  # noqa: E402

gdal.UseExceptions()


def check_against_old(lab, values, nodata, aggfun, skip_box=None):
	"""segment_stats() against old_segment() for every segment of LAB (SKIP_BOX(box) True: not compared)."""
	boxes = ndimage.find_objects(lab)
	n, stat, var, allnan = core.segment_stats(lab, values, aggfun, nodata)
	nnone = 0
	for s in np.unique(lab[lab > 0]):
		box = boxes[s - 1]
		if skip_box is not None and skip_box(box):
			continue
		mask = (lab[box] == s).astype(np.int32)
		ref = t.old_segment(mask, values[box], nodata, aggfun, mask.sum())
		assert n[s] == mask.sum(), s
		if ref[0] == 'None':
			assert allnan[s], s
			nnone += 1
			continue
		scale = abs(ref[0]) + 1
		assert not allnan[s], s
		assert t.close(float(stat[s]), ref[0], 1e-15 * scale), (aggfun.__name__, nodata, s, stat[s], ref[0])
		assert t.close(float(var[s]), ref[1], 1e-24 * scale ** 2), (aggfun.__name__, nodata, s, var[s], ref[1])
		if aggfun is not np.mean:
			assert float(stat[s]) == ref[0], (aggfun.__name__, s)
	return nnone


def test_nodata_none_nan():
	rng = np.random.default_rng(11)
	lab = t.synthetic_labels(rng)
	v = rng.normal(5, 2, lab.shape)
	v[rng.random(lab.shape) < 0.1] = np.nan
	v[lab == 1003] = np.nan  # a segment without data
	for nodata in (None, float('nan')):
		for dtype in (np.float64, np.float32):
			for aggfun in (np.mean, np.median):
				assert check_against_old(lab, v.astype(dtype), nodata, aggfun) > 0
	print("test_nodata_none_nan: OK")


def test_label_rasters():
	n, stat, var, allnan = core.segment_stats(np.zeros((5, 4), np.int32), np.ones((5, 4)), np.mean, None)
	assert n.tolist() == [0] and not allnan.any()
	rng = np.random.default_rng(12)
	lab = t.synthetic_labels(rng)
	v = rng.normal(0, 1, lab.shape)
	lab2 = lab.copy()
	lab2[lab2 == 0] = np.iinfo(np.int32).min  # NoData of an Int32 raster with NULL cells (r.out.gdal)
	for aggfun in (np.mean, np.median):
		a, b = core.segment_stats(lab, v, aggfun), core.segment_stats(lab2, v, aggfun)
		for x, y in zip(a, b):
			assert np.array_equal(x, y, equal_nan=True)
	print("test_label_rasters: OK")


def test_window_past_edges():
	"""_read_pixels on a variable that ends 5 columns and 3 rows before the segment raster: NaN outside, and the
	segments inside are as before."""
	rng = np.random.default_rng(13)
	lab = t.synthetic_labels(rng)
	rows, cols = lab.shape
	full = t.synthetic_variables(rng, lab)[1]  # NoData and NaN
	small = full[:rows - 3, :cols - 5].astype(np.float32)
	with tempfile.TemporaryDirectory() as tmp:
		path = os.path.join(tmp, 'v.tif')
		t.write_tif(path, small, (0.0, t.RES, 0.0, 0.0, 0.0, -t.RES), gdal.GDT_Float32, t.NODATA)
		ds = gdal.Open(path)
		core.res = t.RES
		info = {'band': ds.GetRasterBand(1), 'gt': ds.GetGeoTransform(), 'nodata': t.NODATA}
		seg = core.SegmentPixels(lab)
		v = core._read_pixels(info, (0.0, t.RES, 0.0, 0.0, 0.0, -t.RES), lab.shape, seg.idx)
		r, c = np.divmod(seg.idx, cols)
		assert np.isnan(v[(r >= rows - 3) | (c >= cols - 5)]).all()
		inside = (r < rows - 3) & (c < cols - 5)
		ref = np.where(small == t.NODATA, np.nan, small).astype(np.float64)[r[inside], c[inside]]
		assert np.array_equal(v[inside], ref, equal_nan=True)
		padded = np.full(lab.shape, np.nan)
		padded[:rows - 3, :cols - 5] = np.where(small == t.NODATA, np.nan, small)
		for aggfun in (np.mean, np.median):
			check_against_old(lab, padded, t.NODATA, aggfun,
			                  skip_box=lambda box: box[0].stop > rows - 3 or box[1].stop > cols - 5)
		info = ds = None
	print("test_window_past_edges: OK")


def test_ecor_not_k_and_missing_table():
	rng = np.random.default_rng(14)
	with tempfile.TemporaryDirectory() as tmp:
		old, new, nvars = t._scenario(tmp, rng, 4, 6, (500000.0, 4100000.0), [(500000.0, 4100000.0)] * 3, np.mean)
		work = os.path.join(tmp, 'work')
		# ehabitat2(6, ...): the segments of threshold 6 in the results of k = 2, with 6 in the column ecoregion
		core.gmaps, core.global_maps = 0, {}
		core.ehabitat(6, work, '', 2, np.mean)
		with open(os.path.join(work, 'csv/park_4_hri_results2.csv')) as f:
			rows2 = f.readlines()[1:]
		assert rows2 == new, (rows2[:2], new[:2])
		assert t.compare_rows(old, rows2, nvars, exact_stat=False) < 1e-12
		os.remove(os.path.join(work, 'csv/segm_4_6_ids.csv'))
		try:
			core.ehabitat(6, work, '', 2, np.mean)
		except SystemExit as e:
			assert 'segm_4_6_ids.csv not found' in str(e), e
		else:
			raise AssertionError('ehabitat() did not stop without the table cat -> segm_id')
		for v in core.global_maps.values():
			v['ds'] = v['band'] = None
	print("test_ecor_not_k_and_missing_table: OK")


if __name__ == '__main__':
	test_nodata_none_nan()
	test_label_rasters()
	test_window_past_edges()
	test_ecor_not_k_and_missing_table()
	print("All the tests passed")
