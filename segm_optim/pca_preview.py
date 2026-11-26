#### eHabitat+: the axes of the PCA in each study area before the segmentation (SEGM_PCA_AXES, SEGM_PCA_VARIANCE)
#### License: CC BY-SA 3.0

"""
The standardized PCA of the variables of the model in each study area, as i.pca computes it in the segmentation (the
pixels with data in every variable of the rectangle of the polygons of each area, on the grid of the first variable,
grown by PA_BUFFER), from the GeoTIFFs and without GRASS: how much of the variance the first 3 axes explain, how many
axes a % needs (pca_axes.choose), how much each variable loads on each axis (its correlation with the axis) and how
much of each variable the first axes take, to choose SEGM_PCA_AXES or SEGM_PCA_VARIANCE before segmenting. The app uses
it in HFT mapping (See the PCA in the study area); from segm_optim:

    python3 pca_preview.py [--variance 90 | --axes 4] [--ids 1 2 ...]

It is not part of the code of the segmentation (hft_state.STEPS): changing it does not make the segmentations out of
date.
"""

import math
import os

import pca_axes


def _window(gt, nx, ny, box, grow):
	"""(column, row, columns, rows) of the grid GT (NX x NY) that covers BOX (xmin, ymin, xmax, ymax), its edges on the
	lines of the grid (as g.region align), grown by GROW cells, clipped to the grid."""
	x0 = math.floor((box[0] - gt[0]) / gt[1] + 1e-9) - grow
	x1 = math.ceil((box[2] - gt[0]) / gt[1] - 1e-9) + grow
	y0 = math.floor((box[3] - gt[3]) / gt[5] + 1e-9) - grow
	y1 = math.ceil((box[1] - gt[3]) / gt[5] - 1e-9) + grow
	x0, y0, x1, y1 = max(x0, 0), max(y0, 0), min(x1, nx), min(y1, ny)
	return x0, y0, max(x1 - x0, 0), max(y1 - y0, 0)


def area_boxes(vector, col_id, ids=None, layer=None, srs_wkt=None):
	"""{id: (xmin, ymin, xmax, ymax)} of the polygons of every ID (COL_ID) of the vector file, in the CRS SRS_WKT."""
	from osgeo import ogr, osr
	ds = ogr.Open(vector)
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	tr = None
	if srs_wkt and lyr.GetSpatialRef() is not None:
		dst = osr.SpatialReference(wkt=srs_wkt)
		src = lyr.GetSpatialRef()
		for s in (src, dst):
			s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		if not src.IsSame(dst):
			tr = osr.CoordinateTransformation(src, dst)
	want = None if ids is None else {str(i) for i in ids}
	# without the column COL_ID (cat, made when setup.py imports it): the number of the feature from 1 in the order they
	# are read, as v.import and v.in.ogr number them (and inputs.vector_ids), whatever the FIDs of the file
	use_n = not (col_id and lyr.GetLayerDefn().GetFieldIndex(col_id) >= 0)
	out = {}
	for i, f in enumerate(lyr, 1):
		key = str(i) if use_n else str(f.GetField(col_id))
		if want is not None and key not in want:
			continue
		g = f.GetGeometryRef()
		if g is None:
			continue
		g = g.Clone()
		if tr is not None:
			g.Transform(tr)
		x0, x1, y0, y1 = g.GetEnvelope()
		b = out.get(key)
		out[key] = (x0, y0, x1, y1) if b is None else (min(b[0], x0), min(b[1], y0), max(b[2], x1), max(b[3], y1))
	return out


def preview(files, vector, col_id, ids=None, buffer=0, layer=None, rows=512, limit=None):
	"""The standardized PCA of the GeoTIFFs FILES (in the order of the model) in the rectangle of the polygons of every ID
	of the vector file, grown by BUFFER (map units, rounded up to cells as the segmentation does), with the pixels that
	have data in every variable: [{'id', 'pixels', 'axes': [(axis, eigenvalue, % of the variance)]}], as pca_axes.parse
	reads them from i.pca, 'loadings' (for each axis, the loading of each file, in their order: its correlation with the
	axis, from -1 to 1, the eigenvector times the square root of the eigenvalue; the sign of an axis is arbitrary, and is
	that of its largest loading) and 'constant' (the files with the same value in all those pixels: i.pca, and so the
	segmentation, cannot compute the PCA with them, and 'axes' and 'loadings' are then empty). The grid is that of the
	first file; the others must have the same pixel size and CRS, and may cover less (their pixels outside have no data).
	LIMIT: at most so many areas (the first ones)."""
	import numpy as np
	from osgeo import gdal
	gdal.UseExceptions()
	dss = [gdal.Open(f) for f in files]
	if len(dss) < 2:
		raise ValueError('the PCA needs at least 2 variables')
	g0 = dss[0].GetGeoTransform()
	for f, d in zip(files, dss):
		gt = d.GetGeoTransform()
		if abs(gt[1] - g0[1]) > 1e-9 * abs(g0[1]) or abs(gt[5] - g0[5]) > 1e-9 * abs(g0[5]):
			raise ValueError(f'{os.path.basename(f)} has another pixel size than {os.path.basename(files[0])}')
	res = abs(g0[1])
	grow = int(math.ceil(float(buffer or 0) / res - 1e-9)) if buffer else 0
	boxes = area_boxes(vector, col_id, ids, layer, dss[0].GetProjection())
	nds = []
	for d in dss:
		b = d.GetRasterBand(1)
		nds.append(b.GetNoDataValue())
	out = []
	keys = sorted(boxes, key=lambda k: (len(k), k)) if ids is None else [str(i) for i in ids if str(i) in boxes]
	for key in keys[:limit] if limit else keys:
		x0, y0, w, h = _window(g0, dss[0].RasterXSize, dss[0].RasterYSize, boxes[key], grow)
		p = len(dss)
		# sums of the values minus those of the first pixel with data (fewer rounding errors than of the raw values), and
		# the minimum and maximum of each variable, which tell a constant one
		n, s1, s2, ref = 0, np.zeros(p), np.zeros((p, p)), None
		lo, hi = np.full(p, np.inf), np.full(p, -np.inf)
		for r in range(y0, y0 + h, rows):
			hh = min(rows, y0 + h - r)
			cols = []
			for d, nd in zip(dss, nds):
				gt = d.GetGeoTransform()
				# the same pixel of every variable: its offset from the grid of the first one; outside the raster, NoData
				dx = int(round((g0[0] - gt[0]) / gt[1])) + x0
				dy = int(round((g0[3] - gt[3]) / gt[5])) + r
				a = np.full((hh, w), np.nan)
				cx0, cy0 = max(dx, 0), max(dy, 0)
				cx1, cy1 = min(dx + w, d.RasterXSize), min(dy + hh, d.RasterYSize)
				if cx1 > cx0 and cy1 > cy0:
					b = d.GetRasterBand(1).ReadAsArray(cx0, cy0, cx1 - cx0, cy1 - cy0).astype(np.float64)
					if nd is not None and not math.isnan(nd):
						b[b == nd] = np.nan
					a[cy0 - dy:cy1 - dy, cx0 - dx:cx1 - dx] = b
				cols.append(a.ravel())
			x = np.column_stack(cols)
			x = x[np.all(np.isfinite(x), axis=1)]
			if len(x):
				lo, hi = np.minimum(lo, x.min(axis=0)), np.maximum(hi, x.max(axis=0))
				if ref is None:
					ref = x[0].copy()
				x = x - ref
				n += len(x)
				s1 += x.sum(axis=0)
				s2 += x.T @ x
		if n < p + 1:
			out.append({'id': key, 'pixels': int(n), 'axes': [], 'loadings': []})
			continue
		mu = s1 / n
		cov = (s2 - n * np.outer(mu, mu)) / (n - 1)
		sd = np.sqrt(np.clip(np.diag(cov), 0, None))
		constant = [os.path.basename(f) for f, a, b in zip(files, lo, hi) if a == b]
		if constant:  # i.pca -n divides by its standard deviation: nan in every axis, as the segmentation would get
			out.append({'id': key, 'pixels': int(n), 'axes': [], 'loadings': [], 'constant': constant})
			continue
		corr = cov / np.outer(sd, sd)
		eig, vec = np.linalg.eigh(corr)
		order = np.argsort(eig)[::-1]
		eig, vec = eig[order].clip(0, None), vec[:, order]
		pct = 100 * eig / eig.sum()
		load = vec * np.sqrt(eig)  # column i: the correlation of each variable with axis i
		load *= np.where(load[np.abs(load).argmax(axis=0), np.arange(p)] < 0, -1, 1)
		out.append({'id': key, 'pixels': int(n), 'axes': [(i + 1, float(e), float(round(q, 2))) for i, (e, q) in
		                                                 enumerate(zip(eig, pct))],
		            'loadings': [[float(round(x, 4)) for x in load[:, i]] for i in range(p)], 'constant': []})
	return out


def captured(loadings, k):
	"""How much of the variance of each variable (in % of it, in the order of LOADINGS, from preview()) the first K axes
	take: the sum of its squared loadings on them (100 % with all the axes)."""
	if not loadings:
		return []
	return [100 * sum(loadings[i][j] ** 2 for i in range(min(k, len(loadings)))) for j in range(len(loadings[0]))]


def main_variables(loadings_axis, names, top=3, least=0.4):
	"""The variables that load most on an axis (LOADINGS_AXIS, one per name of NAMES): [(name, loading)], the TOP ones
	with |loading| of at least LEAST, from the largest (the largest one also when smaller)."""
	pairs = sorted(zip(names, loadings_axis), key=lambda x: -abs(x[1]))
	return [x for i, x in enumerate(pairs[:top]) if i == 0 or abs(x[1]) >= least]


def study_inputs(cfg):
	"""(files of the variables of the model, in its order; vector file of the study area; its layer; COL_ID; IDs of
	palist.csv) of config.py CFG, as setup.py and the segmentation use them."""
	import glob
	import hft_state
	import var_order
	folder = os.path.join(hft_state.INVARS, hft_state.invars_folder(cfg))
	paths = [p for p in glob.glob(os.path.join(folder, '*')) if p.lower().endswith(('.tif', '.tiff'))]
	ev = cfg.get('ENV_VARS') or {}
	previous = [os.path.basename(v['file']) for v in ev.values() if isinstance(v, dict) and v.get('file')] \
		if isinstance(ev, dict) else []
	files = var_order.ordered(paths, folder, previous)
	rec = hft_state.get_setup(cfg.get('STUDY_AREA')) or {}
	vector = rec.get('study_area_file')
	if not (vector and os.path.isfile(vector)):
		pa = os.path.join(hft_state.HERE, 'pa')
		cands = [os.path.join(pa, f) for f in sorted(os.listdir(pa)) if f.lower().endswith(('.gpkg', '.shp', '.geojson'))] \
			if os.path.isdir(pa) else []
		same = [c for c in cands if os.path.splitext(os.path.basename(c))[0] == cfg.get('STUDY_AREA')]
		vector = (same or cands or [None])[0]
	ids = open(hft_state.PALIST).read().split() if os.path.isfile(hft_state.PALIST) else None
	return files, vector, rec.get('study_area_layer'), cfg.get('COL_ID') or 'cat', ids


def main(argv=None):
	import argparse
	import hft_state
	ap = argparse.ArgumentParser(description='The axes of the PCA of the variables in each study area, before the '
	                                         'segmentation: how much of the variance the first 3 explain, and how many '
	                                         'are needed for a %.')
	how = ap.add_mutually_exclusive_group()
	how.add_argument('--variance', type=float, help='%% of the variance (default: SEGM_PCA_VARIANCE of config.py, or 90)')
	how.add_argument('--axes', type=int, help='a number of axes (default: SEGM_PCA_AXES of config.py, when it is not 0)')
	ap.add_argument('--ids', nargs='*', help='IDs of the study area (default: those of palist.csv)')
	a = ap.parse_args(argv)
	cfg = hft_state.read_config()
	try:  # SEGM_PCA_VARIANCE of config.py only when it is used (not with a number of axes)
		na = a.axes if a.axes is not None else (0 if a.variance is not None else pca_axes.axes_option(cfg))
		v = a.variance if a.variance is not None else 90 if na else (pca_axes.variance_option(cfg) or 90)
	except ValueError as e:
		raise SystemExit(f'ERROR: {e}')
	if a.axes is not None and a.axes < 1:
		raise SystemExit('ERROR: --axes must be 1 or more')
	if not na and not pca_axes.MIN_VARIANCE <= v <= 100:
		raise SystemExit(f'ERROR: --variance must be from {pca_axes.MIN_VARIANCE} to 100')
	files, vector, layer, col_id, ids = study_inputs(cfg)
	if not files or not vector:
		raise SystemExit('ERROR: no variables or no study area: copy them to ../inVars and pa/ (see setup.py)')
	print(f"Variables ({len(files)}): {', '.join(os.path.basename(f) for f in files)}")
	res = preview(files, vector, col_id, a.ids or ids, cfg.get('PA_BUFFER', 0) or 0, layer)
	if not res:
		raise SystemExit(f'ERROR: no study area of {os.path.basename(vector)} with those IDs (column {col_id})')
	for r in res:
		if r.get('constant'):
			print(f"ID {r['id']}: {', '.join(r['constant'])} has the same value in every pixel: the PCA (and the "
			      'segmentation with SEGM_PCA) cannot be computed with it; take it out of the set')
			continue
		if not r['axes']:
			print(f"ID {r['id']}: too few pixels with data in every variable ({r['pixels']})")
			continue
		p = [x for _, _, x in r['axes']]
		k = pca_axes.choose(p, v, na)
		print(f"ID {r['id']} ({r['pixels']:,} pixels): the first 3 axes explain {sum(p[:3]):.1f} % of the variance; "
		      + (f"the first {k} axes, {sum(p[:k]):.1f} %" if na else f"{v:g} % needs {k} axes ({sum(p[:k]):.1f} %)"))
		names = [os.path.splitext(os.path.basename(f))[0] for f in files]
		cum = 0.0
		for n, eig, x in r['axes']:
			cum += x
			main = ', '.join(f'{nm} {ld:+.2f}' for nm, ld in main_variables(r['loadings'][n - 1], names))
			print(f'  PC{n:<3} {x:6.2f} %  cumulative {cum:6.2f} %' + ('  *' if n <= k else '   ') + f'  mainly {main}')
		cap = captured(r['loadings'], k)
		print(f'  Each variable in the {k} axes marked * (the % of its variance they take; under 50 %, hardly segmented):')
		for nm, c in sorted(zip(names, cap), key=lambda x: x[1]):
			print(f'    {c:5.1f} %  {nm}' + ('  <' if c < 50 else ''))


if __name__ == '__main__':
	os.chdir(os.path.dirname(os.path.abspath(__file__)))
	main()
