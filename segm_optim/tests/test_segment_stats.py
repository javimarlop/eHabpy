#### eHabitat+ HFT model: test of the statistics of the segments (ehab_segm_core.py)
#### License: CC BY-SA 3.0
"""
Checks that the statistics of all the segments computed at once from one segment raster (ehab_segm_core.SegmentPixels,
segment_stats and ehabitat) are those of the per-segment algorithm of eHabpy up to 1.5 (one mask tiffs/pa_<segm_id>.tif
per segment), copied here as the reference: mean, median, variance, variance * n and number of pixels, the 'None' of a
variable without data in a segment, the filling of the NoData pixels, the order of the rows and the skipped segments.
Synthetic data: segments of odd and even sizes (1 and 2 pixels too), masked pixels, NoData and NaN in the variables
(at the first and last pixel of a segment too), a segment without any data in a variable, variables with other
origins and an integer variable.

    cd segm_optim && python -s tests/test_segment_stats.py
"""

import os
import sys
import tempfile

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim: config.py, hft_state.py, ehab_segm_core.py

from osgeo import gdal  # noqa: E402
import ehab_segm_core as core  # noqa: E402

gdal.UseExceptions()

RTOL = 1e-12
NODATA = -9999
RES = 10.0


# --- Reference: the per-segment algorithm of ehab_segm_core.ehabitat() up to 1.5 (last in commit cdb2a2a) ------------

def old_segment(pa_mask0, var_bb0, nodata, aggfun, sum_pa_mask):
	"""Statistics of one variable in one segment: PA_MASK0 the mask of the segment (its bounding box, > 0 inside),
	VAR_BB0 the window of the variable under it. Returns ('None', 'None', 'None') or (mean, var, var2) as floats."""
	pa_mask = pa_mask0.flatten()
	ind = pa_mask > 0
	var_bb = var_bb0.astype(np.float64).flatten()
	var_pa0 = var_bb[ind]
	var_pa = np.where(var_pa0 == nodata, float('NaN'), var_pa0)
	mask2var = np.isnan(var_pa)
	if mask2var.all() == True:  # noqa: E712 (as in the original)
		return ('None', 'None', 'None')
	var_pa[mask2var] = np.interp(np.flatnonzero(mask2var), np.flatnonzero(~mask2var), var_pa[~mask2var])
	var_mean = float(aggfun(var_pa))
	var_var = float(np.var(var_pa))
	var_var2 = var_var * int(sum_pa_mask)
	return (var_mean, var_var, var_var2)


def old_ehabitat(ecor, nwpath, k, aggfun, env_vars, global_maps, res, park):
	"""ehab_segm_core.ehabitat() up to 1.5, with its globals as arguments: the rows it appended to
	csv/park_<park>_hri_results<k>.csv, one per segment of csv/park_<park>_<ecor>.csv read from tiffs/pa_<segm_id>.tif."""
	lines = []
	ecoparksf = os.path.join(nwpath, f'csv/park_{park}_{ecor}.csv')
	if not os.path.isfile(ecoparksf):
		return lines
	pa_list0 = np.genfromtxt(ecoparksf, dtype=int)
	pa_list = np.unique(pa_list0)
	for px in range(len(pa_list)):
		pa = pa_list[px]
		outfile = os.path.join(nwpath, f'csv/park_{pa}_hri_results{ecor}.csv')
		pa4 = os.path.join(nwpath, f'tiffs/pa_{pa}.tif')
		if not os.path.isfile(outfile) and os.path.isfile(pa4):
			src_ds_pa = gdal.Open(pa4)
			par = src_ds_pa.GetRasterBand(1)
			pa_mask0 = par.ReadAsArray(0, 0, par.XSize, par.YSize).astype(np.int32)
			pa_mask = pa_mask0.flatten()
			ind = pa_mask > 0
			sum_pa_mask = sum(pa_mask[ind])
			if sum_pa_mask <= 0:
				continue
			gt_pa = src_ds_pa.GetGeoTransform()
			ref_var = list(env_vars.keys())[0]
			gt_ref = global_maps[ref_var]['gt']
			xoff_test = int(round((gt_pa[0] - gt_ref[0]) / res))
			yoff_test = int(round((gt_ref[3] - gt_pa[3]) / res))
			if xoff_test >= 0 and yoff_test >= 0:
				results = {}
				for var in env_vars:
					band_info = global_maps[var]
					gt_var = band_info['gt']
					xoff = int(round((gt_pa[0] - gt_var[0]) / res))
					yoff = int(round((gt_var[3] - gt_pa[3]) / res))
					var_bb0 = band_info['band'].ReadAsArray(xoff, yoff, par.XSize, par.YSize)
					r = old_segment(pa_mask0, var_bb0, band_info['nodata'], aggfun, sum_pa_mask)
					results[var] = [x if x == 'None' else repr(x) for x in r]
				means = ' '.join(results[v][0] for v in env_vars)
				vars1 = ' '.join(results[v][1] for v in env_vars)
				vars2 = ' '.join(results[v][2] for v in env_vars)
				lines.append(f"{ecor} {pa} {means} {vars1} {sum_pa_mask} {vars2}\n")
	return lines


# --- Synthetic data ---------------------------------------------------------------------------------------------------

def synthetic_labels(rng, rows=48, cols=64, nseeds=40):
	"""Segments (Voronoi cells of random seeds, labels with gaps), masked pixels (0), a 1-pixel and a 2-pixel segment."""
	r, c = np.mgrid[0:rows, 0:cols]
	sr, sc = rng.integers(0, rows, nseeds), rng.integers(0, cols, nseeds)
	d = (r[..., None] - sr) ** 2 + (c[..., None] - sc) ** 2 + rng.random((rows, cols, nseeds)) * 0.1
	lab = (np.argmin(d, axis=2) * 3 + 2).astype(np.int32)  # labels 2, 5, 8... (not consecutive)
	lab[rng.random((rows, cols)) < 0.08] = 0  # masked pixels (CLIP_TO_PA)
	lab[:, -3:][rng.random((rows, 3)) < 0.5] = 0
	lab[5, 7] = 1000  # 1 pixel
	lab[20, 30:32] = 1001  # 2 pixels
	lab[40:43, 50:53] = 1002  # 9 pixels
	lab[30:32, 10:14] = 1003  # 8 pixels
	return lab


def synthetic_variables(rng, lab):
	"""Three variables (float with a large offset and a small variance, float with NaN, integer), with NoData."""
	shape = lab.shape
	v1 = 1000.0 + rng.normal(0, 0.01, shape)
	v2 = rng.normal(0.4, 0.2, shape)
	v3 = rng.integers(-50, 300, shape).astype(np.float64)
	for v in (v1, v2, v3):
		v[rng.random(shape) < 0.07] = NODATA
	v2[rng.random(shape) < 0.03] = np.nan
	v2[lab == 1003] = NODATA  # a segment without any data in v2
	v1[lab == 1002] = NODATA  # a segment with data only in its middle pixel: gaps at both ends
	v1[41, 51] = 1000.5
	ys, xs = np.nonzero(lab == 1001)
	v3[ys[0], xs[0]] = NODATA  # 2 pixels, one without data
	return v1, v2, v3


def write_tif(path, arr, gt, gtype, nodata=None):
	ds = gdal.GetDriverByName('GTiff').Create(path, arr.shape[1], arr.shape[0], 1, gtype)
	ds.SetGeoTransform(gt)
	band = ds.GetRasterBand(1)
	if nodata is not None:
		band.SetNoDataValue(nodata)
	band.WriteArray(arr)
	ds = None


def close(a, b, atol):
	return abs(a - b) <= RTOL * max(abs(a), abs(b)) + atol


# --- Tests ------------------------------------------------------------------------------------------------------------

def test_arrays():
	"""segment_stats() against the per-segment algorithm, segment by segment, for mean, median and another function."""
	rng = np.random.default_rng(1)
	lab = synthetic_labels(rng)
	boxes = ndimage.find_objects(lab)
	worst = 0.0
	for values in synthetic_variables(rng, lab):
		for aggfun in (np.mean, np.median, np.max):
			n, stat, var, allnan = core.segment_stats(lab, values, aggfun, NODATA)
			for s in np.unique(lab[lab > 0]):
				box = boxes[s - 1]
				mask = (lab[box] == s).astype(np.int32)
				ref = old_segment(mask, values[box], NODATA, aggfun, mask.sum())
				assert n[s] == mask.sum(), (s, n[s], mask.sum())
				if ref[0] == 'None':
					assert allnan[s], s
					continue
				assert not allnan[s], s
				m, vv = float(stat[s]), float(var[s])
				scale = abs(ref[0]) + 1
				assert close(m, ref[0], 1e-15 * scale), (aggfun.__name__, s, m, ref[0])
				assert close(vv, ref[1], 1e-24 * scale ** 2), (aggfun.__name__, s, vv, ref[1])
				assert close(vv * int(n[s]), ref[2], 1e-24 * scale ** 2 * n[s]), (aggfun.__name__, s)
				if aggfun is not np.mean:
					assert m == ref[0], (aggfun.__name__, s, m, ref[0])  # the same values: exactly the same result
				worst = max(worst, abs(m - ref[0]) / max(abs(ref[0]), 1e-300), abs(vv - ref[1]) / max(abs(ref[1]), 1e-300))
	sizes = np.bincount(lab[lab > 0])
	sizes = sizes[sizes > 0]
	assert (sizes % 2 == 0).any() and (sizes % 2 == 1).any() and 1 in sizes and 2 in sizes
	print(f"test_arrays: OK ({len(sizes)} segments, largest relative difference {worst:.2e})")


def _scenario(tmp, rng, park, k, seg_origin, var_origins, aggfun):
	"""Files of a run (segment raster + its table and lists, the per-segment masks of the old algorithm, variables) in
	TMP/work and TMP/inVars; returns the rows of the old algorithm and of ehab_segm_core.ehabitat()."""
	work, invars = os.path.join(tmp, 'work'), os.path.join(tmp, 'inVars')
	for d in (os.path.join(work, 'csv'), os.path.join(work, 'tiffs'), invars):
		os.makedirs(d, exist_ok=True)
	lab = synthetic_labels(rng)
	rows, cols = lab.shape
	gt_seg = (seg_origin[0], RES, 0.0, seg_origin[1], 0.0, -RES)
	cats = sorted(set(np.unique(lab[lab > 0]).tolist()) | {78})  # cat 78: a segment without pixels (masked)
	sid = {c: int(f'{park}{k}{c:05d}') for c in cats}

	# Segmentation (new): segment raster, cat -> segm_id, segm_id sorted as text
	write_tif(os.path.join(work, f'tiffs/segm_{park}_{k}.tif'), lab, gt_seg, gdal.GDT_Int32)
	with open(os.path.join(work, f'csv/segm_{park}_{k}_ids.csv'), 'w') as f:
		f.write('cat segm_id\n' + ''.join(f'{c} {sid[c]}\n' for c in cats))
	with open(os.path.join(work, f'csv/park_{park}_{k}.csv'), 'w') as f:
		f.write(''.join(s + '\n' for s in np.unique([str(sid[c]) for c in cats])))
	# Segmentation (old): one mask per segment on its bounding box, 1 inside and 0 outside
	boxes = ndimage.find_objects(lab)
	for c in cats:
		if c < len(boxes) + 1 and boxes[c - 1] is not None:
			box = boxes[c - 1]
			mask = (lab[box] == c).astype(np.uint8)
			gt = (gt_seg[0] + box[1].start * RES, RES, 0.0, gt_seg[3] - box[0].start * RES, 0.0, -RES)
		else:
			mask, gt = np.zeros((1, 1), np.uint8), gt_seg
		write_tif(os.path.join(work, f'tiffs/pa_{sid[c]}.tif'), mask, gt, gdal.GDT_Byte)

	# Variables: on the same grid, each with its own origin and size (covering the segment raster to its right and
	# bottom edges), NoData -9999; v3 integer
	env_vars, global_maps = {}, {}
	full = synthetic_variables(rng, lab)
	for i, (values, (ox, oy)) in enumerate(zip(full, var_origins), 1):
		dx, dy = int(round((seg_origin[0] - ox) / RES)), int(round((oy - seg_origin[1]) / RES))
		big = np.full((rows + dy + 4, cols + dx + 3), float(NODATA))
		y0, x0 = max(dy, 0), max(dx, 0)
		big[y0:dy + rows, x0:dx + cols] = values[y0 - dy:, x0 - dx:]
		name = f'v{i}.tif'
		gtype = gdal.GDT_Int16 if i == 3 else gdal.GDT_Float32
		write_tif(os.path.join(invars, name), big.astype(np.int16) if i == 3 else big,
		          (ox, RES, 0.0, oy, 0.0, -RES), gtype, NODATA)
		env_vars[f'v{i}'] = {'file': name, 'nodata': NODATA}
	for var, props in env_vars.items():
		ds = gdal.Open(os.path.join(invars, props['file']))
		global_maps[var] = {'ds': ds, 'band': ds.GetRasterBand(1), 'gt': ds.GetGeoTransform(), 'nodata': props['nodata']}

	old = old_ehabitat(k, work, k, aggfun, env_vars, global_maps, RES, park)

	# New: ehabitat() of ehab_segm_core with this run as its globals
	core.ENV_VARS, core.res, core.gmaps, core.global_maps, core.park = env_vars, RES, 0, {}, park
	out = os.path.join(work, f'csv/park_{park}_hri_results{k}.csv')
	for f in (out, os.path.join(work, f'csv/ecoregs_done{k}.csv')):
		if os.path.exists(f):
			os.remove(f)
	core.ehabitat(k, work, '', k, aggfun)
	with open(out) as f:
		header, *new = f.readlines()
	expected = ('ecoregion segm_id ' + ' '.join(f'{v}pamean' for v in env_vars) + ' '
	            + ' '.join(f'{v}pavar' for v in env_vars) + ' sumpamask ' + ' '.join(f'{v}pavar2' for v in env_vars) + '\n')
	assert header == expected, header
	for maps in (global_maps, core.global_maps):  # close the files
		for v in maps.values():
			v['ds'] = v['band'] = None
	return old, new, len(env_vars)


def compare_rows(old, new, nvars, exact_stat):
	assert len(old) == len(new), (len(old), len(new))
	worst = 0.0
	for a, b in zip(old, new):
		a, b = a.split(), b.split()
		assert len(a) == len(b) == 3 + 3 * nvars, (a, b)
		assert a[:2] == b[:2] and a[2 + 2 * nvars] == b[2 + 2 * nvars], (a[:2], b[:2])  # ecoregion, segm_id, n
		for i, (x, y) in enumerate(zip(a, b)):
			if i < 2 or i == 2 + 2 * nvars:
				continue
			assert (x == 'None') == (y == 'None'), (a[:2], i, x, y)
			if x == 'None':
				continue
			x, y = float(x), float(y)
			if exact_stat and i < 2 + nvars:
				assert x == y, (a[:2], i, x, y)
			assert close(x, y, 1e-20), (a[:2], i, x, y)
			worst = max(worst, abs(x - y) / max(abs(x), 1e-300))
	return worst


def test_files():
	"""ehabitat() on a segment raster against the old per-segment algorithm on the masks of every segment."""
	rng = np.random.default_rng(7)
	cases = [
		# park, k, origin of the segment raster, origins of the three variables (on the same grid)
		(7, 3, (500000.0, 4100000.0), [(499980.0, 4100030.0), (499900.0, 4100100.0), (500000.0, 4100000.0)]),
		# the segment raster starts 2 columns left of and 1 row above the first variable: the segments that start
		# there are skipped, as before; a negative area ID: the rows sorted by segm_id are not in the order of cat
		(-2, 5, (500000.0, 4100000.0), [(500020.0, 4099990.0), (499950.0, 4100050.0), (499990.0, 4100010.0)]),
	]
	for park, k, seg_origin, var_origins in cases:
		for aggfun in (np.mean, np.median):
			with tempfile.TemporaryDirectory() as tmp:
				old, new, nvars = _scenario(tmp, rng, park, k, seg_origin, var_origins, aggfun)
				nlisted = sum(1 for _ in open(os.path.join(tmp, 'work', f'csv/park_{park}_{k}.csv')))
				worst = compare_rows(old, new, nvars, exact_stat=aggfun is np.median)
				nnone = sum(line.split().count('None') for line in new)
				print(f"test_files: OK (area {park}, threshold {k}, {aggfun.__name__}: {len(new)} rows of "
				      f"{nlisted} segments listed, {nnone} 'None', largest relative difference {worst:.2e})")
				assert nnone > 0 and len(new) < nlisted


def test_wrappers():
	"""The backward-compatible functions of get_mean_median_var.py are still there."""
	import get_mean_median_var as g
	for name in ['run_batch_all', 'run_park', 'run_batch', 'ehabitat'] + [f'{p}{k}' for p in ('ehabitat', 'run_batch')
	                                                                    for k in range(1, 10)]:
		assert callable(getattr(g, name)), name
	print("test_wrappers: OK")


if __name__ == '__main__':
	test_arrays()
	test_files()
	test_wrappers()
	print("All the tests passed")
