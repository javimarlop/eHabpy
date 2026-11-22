#### eHabitat+ HFT model: more cases of the test of the statistics of the segments (ehab_segm_core.py)
#### License: CC BY-SA 3.0
"""
Cases that tests/test_segment_stats.py and tests/test_segment_stats_edge.py do not cover, against the same per-segment
reference (old_segment, old_ehabitat: ehab_segm_core.py up to 1.5):
- segments made of several separate patches (in the same rows, with a hole taken by another segment), whose gaps are
  filled by interpolation from one patch to the other in row-major order, as the per-segment masks did;
- constant segments of whole numbers: the same mean and a variance of exactly 0 (so that a constant variable keeps
  giving Moran's I = nan downstream);
- run_batch_all() as get_mean_median_var.py runs it: two study areas of csv/segm_done.csv (one with a negative ID) on
  the same variables, the nine thresholds, csv/ecoregs<k>.csv shared by both, a threshold without segments in one
  area (its results file keeps only the header, as before), and the progress (##UNITS / ##DONE with EHAB_PROGRESS) adding up to the same total.

    cd segm_optim && python -s tests/test_segment_stats_areas.py
"""

import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
SEGM_OPTIM = os.path.dirname(HERE)
sys.path.insert(0, SEGM_OPTIM)
sys.path.insert(0, HERE)

from osgeo import gdal  # noqa: E402
import ehab_segm_core as core  # noqa: E402
import test_segment_stats as t  # noqa: E402

gdal.UseExceptions()


def check(lab, values, nodata, aggfuns=(np.mean, np.median)):
	boxes = ndimage.find_objects(lab)
	for aggfun in aggfuns:
		n, stat, var, allnan = core.segment_stats(lab, values, aggfun, nodata)
		for s in np.unique(lab[lab > 0]):
			box = boxes[s - 1]
			mask = (lab[box] == s).astype(np.int32)
			ref = t.old_segment(mask, values[box], nodata, aggfun, mask.sum())
			assert n[s] == mask.sum(), s
			if ref[0] == 'None':
				assert allnan[s], s
				continue
			scale = abs(ref[0]) + 1
			assert t.close(float(stat[s]), ref[0], 1e-15 * scale), (aggfun.__name__, s, stat[s], ref[0])
			assert t.close(float(var[s]), ref[1], 1e-24 * scale ** 2), (aggfun.__name__, s, var[s], ref[1])
			if aggfun is np.median:
				assert float(stat[s]) == ref[0], (s, stat[s], ref[0])


def test_patches():
	"""A segment of three separate patches (two in the same rows, one with a hole of another segment)."""
	lab = np.full((12, 20), 3, np.int32)
	lab[:, -2:] = 0
	lab[2:5, 1:4] = 7
	lab[2:5, 10:14] = 7
	lab[8:10, 5:8] = 7
	lab[3, 11] = 9  # hole of segment 7 taken by segment 9
	lab[6, 15:17] = 11  # 2 pixels: even median
	rng = np.random.default_rng(21)
	v = rng.normal(10, 3, lab.shape)
	v[2, 3] = v[2, 10] = v[3, 1] = np.nan  # gaps at the end of a patch and the start of the next one, in a row
	v[4, 13] = v[8, 5] = np.nan  # last pixel of the second patch and first one of the third
	v[9, 7] = t.NODATA  # last pixel of the segment
	check(lab, v, t.NODATA)
	print("test_patches: OK")


def test_constant_whole_numbers():
	"""Whole numbers (integer rasters) are summed exactly: the same mean as before and a variance of exactly 0."""
	lab = np.repeat(np.arange(1, 7, dtype=np.int32), [1, 2, 9, 128, 1000, 20000]).reshape(1, -1)
	for c in (0.0, 7.0, -3.0, 2500.0):
		v = np.full(lab.shape, c)
		for aggfun in (np.mean, np.median):
			n, stat, var, allnan = core.segment_stats(lab, v, aggfun, t.NODATA)
			assert (stat[1:] == c).all() and (var[1:] == 0).all(), (c, stat, var)
	print("test_constant_whole_numbers: OK")


DRIVER = r'''
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import hft_state
hft_state.fixed_level = lambda cfg=None: None  # the nine thresholds, whatever THRESHOLD_METHOD config.py has
import ehab_segm_core as core
core.ENV_VARS, core.res = json.loads(sys.argv[2]), float(sys.argv[3])
core.run_batch_all(np.median if sys.argv[4] == 'median' else np.mean)
'''


def build_area(work, rng, park, origin, grid_origin, full_vars, empty_k=None):
	"""Segment rasters, tables and lists of the nine thresholds of area PARK (and the per-segment masks of the old
	algorithm) on the grid of the variables; threshold EMPTY_K without segments. Returns the number of segments of each k."""
	counts = {}
	for k in range(1, 10):
		lab = t.synthetic_labels(rng, rows=30, cols=40, nseeds=6 + 3 * k)
		if k == empty_k:
			lab[:] = 0
		rows, cols = lab.shape
		gt_seg = (origin[0], t.RES, 0.0, origin[1], 0.0, -t.RES)
		cats = np.unique(lab[lab > 0]).tolist()
		sid = {c: int(f'{park}{k}{c:05d}') for c in cats}
		t.write_tif(os.path.join(work, f'tiffs/segm_{park}_{k}.tif'), lab, gt_seg, gdal.GDT_Int32)
		with open(os.path.join(work, f'csv/segm_{park}_{k}_ids.csv'), 'w') as f:
			f.write('cat segm_id\n' + ''.join(f'{c} {sid[c]}\n' for c in cats))
		if cats:  # as segmentation_pca_par.py: no list without segments
			with open(os.path.join(work, f'csv/park_{park}_{k}.csv'), 'w') as f:
				f.write(''.join(s + '\n' for s in np.unique([str(sid[c]) for c in cats])))
		boxes = ndimage.find_objects(lab)
		for c in cats:
			box = boxes[c - 1]
			gt = (gt_seg[0] + box[1].start * t.RES, t.RES, 0.0, gt_seg[3] - box[0].start * t.RES, 0.0, -t.RES)
			t.write_tif(os.path.join(work, f'tiffs/pa_{sid[c]}.tif'), (lab[box] == c).astype(np.uint8), gt, gdal.GDT_Byte)
		counts[k] = len(cats)
		# some gaps of every variable in the pixels of this area
		dy, dx = int(round((grid_origin[1] - origin[1]) / t.RES)), int(round((origin[0] - grid_origin[0]) / t.RES))
		for v in full_vars:
			win = v[dy:dy + rows, dx:dx + cols]
			win[(rng.random(win.shape) < 0.05) & (lab > 0)] = t.NODATA
	return counts


def test_run_batch_all_two_areas():
	rng = np.random.default_rng(31)
	grid_origin = (400000.0, 4200000.0)
	gshape = (70, 90)
	full_vars = [1000.0 + rng.normal(0, 0.01, gshape), rng.normal(0.4, 0.2, gshape),
	             rng.integers(-50, 300, gshape).astype(np.float64)]
	areas = {12: (400100.0, 4199850.0), -3: (400350.0, 4199700.0)}  # overlapping windows of the same grid
	for aggname, aggfun in (('mean', np.mean), ('median', np.median)):
		with tempfile.TemporaryDirectory() as tmp:
			work, invars = os.path.join(tmp, 'work'), os.path.join(tmp, 'inVars')
			for d in (os.path.join(work, 'csv'), os.path.join(work, 'tiffs'), invars):
				os.makedirs(d)
			vars_now = [v.copy() for v in full_vars]
			counts = {park: build_area(work, rng, park, origin, grid_origin, vars_now, 4 if park == 12 else None)
			          for park, origin in areas.items()}
			# csv/ecoregs<k>.csv: shared by the areas, k once per segment of each one
			for k in range(1, 10):
				with open(os.path.join(work, f'csv/ecoregs{k}.csv'), 'w') as f:
					f.write(f'{k}\n' * sum(c[k] for c in counts.values()))
			with open(os.path.join(work, 'csv/segm_done.csv'), 'w') as f:
				f.write(''.join(f'{p}\n' for p in areas))
			env_vars = {}
			for i, v in enumerate(vars_now, 1):
				name = f'v{i}.tif'
				gtype = gdal.GDT_Int16 if i == 3 else gdal.GDT_Float32
				t.write_tif(os.path.join(invars, name), v.astype(np.int16) if i == 3 else v,
				            (grid_origin[0], t.RES, 0.0, grid_origin[1], 0.0, -t.RES), gtype, t.NODATA)
				env_vars[f'v{i}'] = {'file': name, 'nodata': t.NODATA}

			env = dict(os.environ, EHAB_PROGRESS='1', PYTHONNOUSERSITE='1')
			p = subprocess.run([sys.executable, '-s', '-c', DRIVER, SEGM_OPTIM, json.dumps(env_vars), str(t.RES), aggname],
			                   cwd=work, env=env, capture_output=True, text=True)
			assert p.returncode == 0, p.stdout[-2000:] + p.stderr[-2000:]
			units = sum(float(line.split()[1]) for line in p.stdout.splitlines() if line.startswith('##UNITS '))
			done = sum(float(line.split()[1]) for line in p.stdout.splitlines() if line.startswith('##DONE '))
			total = sum(sum(c.values()) for c in counts.values())
			assert units == total and abs(done - units) < 1e-3 * 9 * len(areas), (units, done, total)

			global_maps = {}
			for var, props in env_vars.items():
				ds = gdal.Open(os.path.join(invars, props['file']))
				global_maps[var] = {'ds': ds, 'band': ds.GetRasterBand(1), 'gt': ds.GetGeoTransform(),
				                    'nodata': props['nodata']}
			worst, nrows = 0.0, 0
			for park in areas:
				for k in range(1, 10):
					old = t.old_ehabitat(k, work, k, aggfun, env_vars, global_maps, t.RES, park)
					with open(os.path.join(work, f'csv/park_{park}_hri_results{k}.csv')) as f:
						header, *new = f.readlines()
					assert header.startswith('ecoregion segm_id v1pamean'), header
					assert (park == 12 and k == 4) == (len(new) == 0), (park, k, len(new))
					worst = max(worst, t.compare_rows(old, new, len(env_vars), exact_stat=aggfun is np.median))
					nrows += len(new)
			for v in global_maps.values():
				v['ds'] = v['band'] = None
			print(f"test_run_batch_all_two_areas: OK ({aggname}: areas {list(areas)}, {nrows} rows, progress "
			      f"{done:.4f} of {units:.0f} units, largest relative difference {worst:.2e})")


if __name__ == '__main__':
	test_patches()
	test_constant_whole_numbers()
	test_run_batch_all_two_areas()
	print("All the tests passed")
