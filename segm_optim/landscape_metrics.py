#### eHabitat+: landscape metrics of maps of patches of the user (rasters or polygons), e.g. a habitat at several dates
#### Author: Javier Martinez-Lopez (UTF-8) 2026
#### License: CC BY-SA 3.0

"""
The landscape metrics (lsm.py: those of the R package landscapemetrics, which follows FRAGSTATS) of one or more maps of
patches given by the user, without the similarity of the HFTs: rasters of classes (GeoTIFF) or polygons (GeoPackage,
shapefile, GeoJSON), e.g. the same habitat mapped at several dates, to follow its fragmentation. Run from segm_optim:

    python3 landscape_metrics.py --spec landscape/<run>/landscape_spec.json    # as the app runs it
    python3 landscape_metrics.py habitat_2005.tif habitat_2020.tif --habitat 1 --name habitat   # habitat: value 1
    python3 landscape_metrics.py lc_2000.tif lc_2020.tif --mode classes --boundary park.gpkg      # every class
    python3 landscape_metrics.py --inspect map.gpkg       # what a file holds (JSON): classes, fields, CRS, pixel
    python3 landscape_metrics.py --catalogue              # the metrics (JSON)

Each map is a raster of classes (whole numbers; NoData outside the landscape) or polygons, rasterized at the pixel size
(--res; every polygon map of a run on the same grid: that of the boundary, or the rectangle around all of them, so that
the dates are comparable). With --mode habitat (the default) the classes or polygons chosen as habitat (--habitat:
values of the raster or of the column --field; all the polygons without --field) are one class, 1, and the rest of the
landscape 0: the metrics are those of the habitat (patch and class levels). With --mode classes every class is kept (the
values of the raster or of --field, with the same codes in every map; the cells inside the landscape without a polygon
are one more class, '(no polygon)') and the landscape level is computed too. --boundary: the polygons of the landscape
(e.g. the study area): cells outside are left out. Maps in geographic coordinates are projected to their UTM zone.

The results go to landscape/<name> (--name): metrics_landscape.csv, metrics_class.csv and metrics_patch.csv (one row
per map, per map and class, per patch; space-separated, NA), patches_<n>.gpkg (the patches of the map n as polygons,
with their patch metrics), classes_<n>.tif (the map as it was measured), metrics_trend.png/.csv (the class metrics of
every map, by landscape_plots.R), landscape_metadata.txt and landscape_run.json (written last). Patches of 8 neighbours
(--neighbours 4: only those that share a side), core areas farther than --edge-depth cells from their edge.
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime

import numpy as np
from osgeo import gdal, ogr, osr

import lsm

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'landscape')
RUN_JSON = 'landscape_run.json'
SPEC_JSON = 'landscape_spec.json'
RUN_TXT = 'landscape_metadata.txt'
PLOTS_R = os.path.join(HERE, 'landscape_plots.R')
RASTER_EXT = ('.tif', '.tiff', '.img', '.vrt', '.asc')
NAME_RE = re.compile(r'^[A-Za-z0-9_-]{1,80}$')
NODATA = 255  # of the maps classes_<n>.tif of a run with fewer than 255 classes (else -2147483648)
MAX_CELLS = 400e6  # of a map (rows x columns): beyond it, a coarser pixel size is needed
# peak memory of a run per cell of a map (measured with fragmented maps): one class (mode habitat), every class
BYTES_PER_CELL = {'habitat': 170, 'classes': 290}
MAX_CLASSES = 1000  # distinct values of a raster of classes
# The metrics of a run unless others are chosen: those that say most about the amount, subdivision, shape, isolation
# and connectedness of a habitat (FRAGSTATS names; the app shows them first)
FRAGMENTATION = ['c_ca', 'c_pland', 'c_np', 'c_pd', 'c_lpi', 'c_area_mn', 'c_area_cv', 'c_te', 'c_ed', 'c_lsi',
                 'c_shape_mn', 'c_frac_mn', 'c_gyrate_mn', 'c_enn_mn', 'c_enn_cv', 'c_tca', 'c_cpland', 'c_cai_mn',
                 'c_ai', 'c_clumpy', 'c_pladj', 'c_cohesion', 'c_division', 'c_mesh', 'c_split',
                 'l_np', 'l_pd', 'l_lpi', 'l_ed', 'l_lsi', 'l_contag', 'l_shdi', 'l_shei', 'l_sidi', 'l_pr',
                 'l_mesh', 'l_division', 'l_ai', 'l_iji', 'l_ent', 'l_condent', 'l_mutinf', 'l_relmutinf']
PATCH_METRICS = lsm.keys('patch')  # always, for the polygons of the patches and metrics_patch.csv
GTIFF = ['TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=2', 'BIGTIFF=IF_SAFER']


def rel(path):
	if not path:
		return None
	try:
		r = os.path.relpath(path, ROOT)
	except ValueError:
		return os.path.abspath(path)
	return os.path.abspath(path) if r.startswith('..') else r.replace(os.sep, '/')


def user_path(path):
	path = os.path.expanduser(str(path or '').strip())
	return path if not path or os.path.isabs(path) else os.path.normpath(os.path.join(ROOT, path))


def is_raster(path):
	return str(path).lower().endswith(RASTER_EXT)


def _srs(wkt_or_srs):
	if wkt_or_srs is None:
		return None
	s = wkt_or_srs.Clone() if isinstance(wkt_or_srs, osr.SpatialReference) else osr.SpatialReference()
	if not isinstance(wkt_or_srs, osr.SpatialReference):
		s.ImportFromWkt(wkt_or_srs)
	s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return s


def crs_label(srs):
	if srs is None:
		return 'no CRS'
	s = srs.Clone()
	try:
		s.AutoIdentifyEPSG()
	except RuntimeError:
		pass
	code = s.GetAuthorityCode(None)
	return f'{s.GetAuthorityName(None)}:{code}' if code else (s.GetName() or 'unknown CRS')


def utm_for(lon, lat):
	u = osr.SpatialReference()
	u.SetWellKnownGeogCS('WGS84')
	u.SetUTM(int((lon + 180) // 6) % 60 + 1, lat >= 0)
	u.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return u


def nice(v):
	"""A round number near V (1, 2, 2.5 or 5 x a power of ten), at least V."""
	if v <= 0:
		return 1.0
	p = 10 ** math.floor(math.log10(v))
	for m in (1, 2, 2.5, 5, 10):
		if m * p >= v - 1e-9:
			return float(m * p)
	return float(10 * p)


def memory_total():
	"""The physical memory of this computer (bytes), or None if it cannot be known."""
	try:
		return os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
	except (ValueError, OSError, AttributeError):
		pass
	try:  # Windows
		import ctypes

		class Status(ctypes.Structure):
			_fields_ = [('dwLength', ctypes.c_ulong), ('dwMemoryLoad', ctypes.c_ulong)] + [
				(n, ctypes.c_ulonglong) for n in ('ullTotalPhys', 'ullAvailPhys', 'ullTotalPageFile', 'ullAvailPageFile',
				                                  'ullTotalVirtual', 'ullAvailVirtual', 'ullAvailExtendedVirtual')]
		st = Status()
		st.dwLength = ctypes.sizeof(Status)
		if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
			return int(st.ullTotalPhys)
	except Exception:  # noqa: BLE001 - not Windows, or no answer
		pass
	return None


def check_size(width, height, what, mode='classes'):
	"""ValueError when a grid of WIDTH x HEIGHT cells is too large: more than MAX_CELLS, or more memory than this
	computer has (about BYTES_PER_CELL of the MODE per cell, for the metrics and the polygons of the patches)."""
	cells = float(width) * float(height)
	if cells > MAX_CELLS:
		raise ValueError(f'{what} has {width} x {height} pixels: too many; choose a larger pixel size')
	need, total = cells * BYTES_PER_CELL.get(mode, max(BYTES_PER_CELL.values())) + 4e8, memory_total()
	if total and need > 0.8 * total:
		raise ValueError(f'{what} has {width} x {height} pixels, which need about {need / 2 ** 30:.1f} GB of memory, and '
		                 f'this computer has {total / 2 ** 30:.1f} GB: choose a larger pixel size, or a boundary that '
		                 f'leaves out the rest of the map')


def _num(v, digits=10):
	if v is None or (isinstance(v, float) and not math.isfinite(v)):
		return 'NA'
	if isinstance(v, (int, np.integer)):
		return str(int(v))
	return f'{float(v):.{digits}g}'


def _text(v):
	"""A text value of the space-separated tables: quoted when it has spaces or quotes (read.table reads it)."""
	s = str(v)
	return '"' + s.replace('"', "'") + '"' if re.search(r'[\s"#\']', s) or not s else s


# Inspecting a file ----------------------------------------------------------------------------------------------------

def _layer(ds, layer):
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	if lyr is None:
		raise ValueError(f'no layer {layer!r} in the file')
	return lyr


def raster_classes(path, limit=MAX_CLASSES):
	"""{value: cells} of the first band of PATH (without its NoData), read in blocks; ValueError beyond LIMIT values or
	with values that are not whole numbers."""
	ds = gdal.Open(path)
	b = ds.GetRasterBand(1)
	nd = b.GetNoDataValue()
	counts = {}
	bw, bh = b.GetBlockSize()
	step = max(bh, int(4e6 // max(ds.RasterXSize, 1)))
	for r0 in range(0, ds.RasterYSize, step):
		a = b.ReadAsArray(0, r0, ds.RasterXSize, min(step, ds.RasterYSize - r0))
		ok = np.ones(a.shape, bool) if a.dtype.kind in 'iub' else np.isfinite(a)
		if nd is not None and not (isinstance(nd, float) and math.isnan(nd)):
			ok &= a != nd
		v = a[ok]
		if v.dtype.kind == 'f' and v.size and np.any(v != np.round(v)):
			raise ValueError(f'{os.path.basename(path)} has values that are not whole numbers: it is not a map of '
			                 f'classes (patches need classes, e.g. 1 = habitat, 0 = not)')
		u, c = np.unique(v, return_counts=True)
		for x, n in zip(u.tolist(), c.tolist()):
			counts[int(x)] = counts.get(int(x), 0) + n
		if len(counts) > limit:
			raise ValueError(f'{os.path.basename(path)} has more than {limit} distinct values: it is not a map of '
			                 f'classes')
	return dict(sorted(counts.items()))


def inspect(path, layer=None):
	"""What a map of patches holds, for the app: a raster (CRS, pixel, size, its classes and cells) or a vector file
	(layers, fields with their values when they have few, number of polygons, CRS, extent)."""
	path = user_path(path)
	if not os.path.isfile(path):
		raise ValueError(f'{path} not found')
	out = {'path': path, 'file': os.path.basename(path), 'name': os.path.splitext(os.path.basename(path))[0]}
	m = re.search(r'(?<!\d)(1[89]\d\d|20\d\d|21\d\d)(?!\d)', out['name'])
	out['year'] = m.group(1) if m else None
	if is_raster(path):
		ds = gdal.Open(path)
		gt = ds.GetGeoTransform()
		srs = _srs(ds.GetSpatialRef())
		out.update({'kind': 'raster', 'crs': crs_label(srs), 'projected': bool(srs and srs.IsProjected()),
		            'res': [abs(gt[1]), abs(gt[5])], 'width': ds.RasterXSize, 'height': ds.RasterYSize,
		            'nodata': ds.GetRasterBand(1).GetNoDataValue()})
		if out['nodata'] is not None and not math.isfinite(out['nodata']):
			out['nodata'] = str(out['nodata'])
		try:
			counts = raster_classes(path)
			cell = abs(gt[1] * gt[5]) * (srs.GetLinearUnits() ** 2 if srs and srs.IsProjected() else float('nan'))
			out['classes'] = [{'value': k, 'cells': n, 'ha': round(n * cell / 1e4, 2) if math.isfinite(cell) else None}
			                  for k, n in counts.items()]
		except ValueError as e:
			out['classes'], out['error'] = [], str(e)
		return out
	ds = ogr.Open(path)
	if ds is None:
		raise ValueError(f'{os.path.basename(path)} is neither a raster nor a vector file that GDAL can read')
	layers = [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]
	lyr = _layer(ds, layer if layer in layers else None)
	srs = _srs(lyr.GetSpatialRef())
	defn = lyr.GetLayerDefn()
	fields = []
	for i in range(defn.GetFieldCount()):
		fd = defn.GetFieldDefn(i)
		if fd.GetType() not in (ogr.OFTInteger, ogr.OFTInteger64, ogr.OFTString):
			continue
		fields.append({'name': fd.GetName(), 'integer': fd.GetType() != ogr.OFTString, 'values': {}})
	n = 0
	gtypes = set()
	for f in lyr:
		n += 1
		g = f.GetGeometryRef()
		if g is not None:
			gtypes.add(ogr.GeometryTypeToName(g.GetGeometryType()))
		for fd in fields:
			if fd['values'] is None:
				continue
			v = f.GetField(fd['name'])
			if v is None:
				continue
			fd['values'][str(v)] = fd['values'].get(str(v), 0) + 1
			if len(fd['values']) > 200:
				fd['values'] = None  # too many: not a field of classes
	xmin, xmax, ymin, ymax = lyr.GetExtent()
	for fd in fields:
		if fd['values'] is not None:
			fd['values'] = [{'value': k, 'polygons': c} for k, c in sorted(
				fd['values'].items(), key=lambda kv: (0, int(kv[0])) if fd['integer'] else (1, kv[0]))]
	out.update({'kind': 'vector', 'layers': layers, 'layer': lyr.GetName(), 'fields': fields, 'count': n,
	            'geometry': sorted(gtypes), 'crs': crs_label(srs), 'projected': bool(srs and srs.IsProjected()),
	            'extent': [xmin, ymin, xmax, ymax]})
	if gtypes and not any('POLYGON' in t.upper() for t in gtypes):
		out['error'] = f'it has no polygons ({", ".join(sorted(gtypes))}): the patches need polygons'
	return out


# The grids ------------------------------------------------------------------------------------------------------------

class Grid:
	def __init__(self, x0, y1, res, width, height, srs):
		self.x0, self.y1, self.res, self.width, self.height = float(x0), float(y1), float(res), int(width), int(height)
		self.srs = srs

	@property
	def gt(self):
		return (self.x0, self.res, 0.0, self.y1, 0.0, -self.res)

	@property
	def extent(self):
		return [self.x0, self.y1 - self.height * self.res, self.x0 + self.width * self.res, self.y1]

	def record(self):
		return {'crs': crs_label(self.srs), 'res': self.res, 'width': self.width, 'height': self.height,
		        'extent': self.extent}


def _mem_raster(grid, dtype=gdal.GDT_Int32, fill=0):
	ds = gdal.GetDriverByName('MEM').Create('', grid.width, grid.height, 1, dtype)
	ds.SetGeoTransform(grid.gt)
	ds.SetProjection(grid.srs.ExportToWkt())
	ds.GetRasterBand(1).Fill(fill)
	return ds


def _layer_extent(path, layer, srs):
	"""Extent of a vector layer in SRS (its corners densified)."""
	ds = ogr.Open(path)
	lyr = _layer(ds, layer)
	src = _srs(lyr.GetSpatialRef())
	xmin, xmax, ymin, ymax = lyr.GetExtent()
	if src is None or srs is None or src.IsSame(srs):
		return [xmin, ymin, xmax, ymax]
	return list(osr.CoordinateTransformation(src, srs).TransformBounds(xmin, ymin, xmax, ymax, 21))


def _layer_srs(path, layer):
	ds = ogr.Open(path)
	return _srs(_layer(ds, layer).GetSpatialRef())


def vector_grid(items, boundary, res, mode='classes'):
	"""The grid of the polygon maps of a run: in the CRS of the boundary, else of the first map, if projected (else the
	UTM zone of their centre), over the boundary or the rectangle around all the maps, with pixels of RES metres
	(None: a round size near 1/2000 of the longer side)."""
	srs = None
	for path, layer in ([boundary] if boundary else []) + items:
		s = _layer_srs(path, layer)
		if s is not None:
			srs = s
			break
	if srs is None:
		raise ValueError('the polygons have no coordinate reference system (CRS)')
	if not srs.IsProjected():
		ext = _layer_extent(*(boundary or items[0]), srs)
		lon, lat = (ext[0] + ext[2]) / 2, (ext[1] + ext[3]) / 2
		srs = utm_for(lon, lat)
	exts = [_layer_extent(p, l, srs) for p, l in ([boundary] if boundary else items)]
	xmin, ymin = min(e[0] for e in exts), min(e[1] for e in exts)
	xmax, ymax = max(e[2] for e in exts), max(e[3] for e in exts)
	if not (xmax > xmin and ymax > ymin):
		raise ValueError('the polygons have an empty extent')
	unit = srs.GetLinearUnits() or 1.0
	if not res:
		res = nice(max(xmax - xmin, ymax - ymin) * unit / 2000)
	r = float(res) / unit  # in the units of the CRS
	x0, y1 = math.floor(xmin / r) * r, math.ceil(ymax / r) * r
	w, h = int(math.ceil((xmax - x0) / r - 1e-9)), int(math.ceil((y1 - ymin) / r - 1e-9))
	check_size(w, h, f'the grid of the polygons ({res:g} m)', mode)
	return Grid(x0, y1, r, w, h, srs), float(res)


def rasterize(path, layer, grid, burn=None, field=None, where_values=None, codes=None, fill=0):
	"""Polygons of PATH on GRID (int32): BURN (1) for every polygon, or the code of its value of FIELD (CODES: {text:
	code}, else the integer value); with WHERE_VALUES, only the polygons whose FIELD is one of them; FILL where there is
	no polygon."""
	ds = ogr.Open(path)
	lyr = _layer(ds, layer)
	if field and lyr.GetLayerDefn().GetFieldIndex(field) < 0:
		raise ValueError(f'{os.path.basename(path)} has no column {field}')
	src = _srs(lyr.GetSpatialRef())
	mem = ogr.GetDriverByName('Memory').CreateDataSource('') if gdal.GetDriverByName('MEM').GetMetadataItem(
		gdal.DCAP_VECTOR) != 'YES' else gdal.GetDriverByName('MEM').Create('', 0, 0, 0, gdal.GDT_Unknown)
	out = mem.CreateLayer('p', srs=grid.srs, geom_type=ogr.wkbMultiPolygon)
	out.CreateField(ogr.FieldDefn('v', ogr.OFTInteger))
	ct = osr.CoordinateTransformation(src, grid.srs) if src is not None and not src.IsSame(grid.srs) else None
	want = {str(v) for v in where_values} if where_values is not None else None
	n = 0
	for f in lyr:
		g = f.GetGeometryRef()
		if g is None:
			continue
		val = f.GetField(field) if field else None
		if want is not None and (val is None or str(val) not in want):
			continue
		if field and burn is None:
			if val is None:
				continue
			code = codes[str(val)] if codes is not None else int(val)
		else:
			code = burn if burn is not None else 1
		g = g.Clone()
		if ct is not None:
			g.Transform(ct)
		nf = ogr.Feature(out.GetLayerDefn())
		nf.SetGeometry(g)
		nf.SetField('v', int(code))
		out.CreateFeature(nf)
		n += 1
	r = _mem_raster(grid, fill=fill)
	gdal.RasterizeLayer(r, [1], out, options=['ATTRIBUTE=v'])
	return r.GetRasterBand(1).ReadAsArray().astype(np.int64), n


def read_raster(path, res=None, mode='classes'):
	"""(values, valid, grid, notes) of the first band of a raster of classes: in its grid or, when it is in geographic
	coordinates or another pixel size RES is asked for, warped with the most common class of each new pixel (mode) to
	that size (its UTM zone for geographic coordinates)."""
	ds = gdal.Open(path)
	gt = ds.GetGeoTransform()
	srs = _srs(ds.GetSpatialRef())
	notes = []
	if srs is None:
		raise ValueError(f'{os.path.basename(path)} has no coordinate reference system (CRS)')
	if abs(gt[2]) > 1e-12 or abs(gt[4]) > 1e-12 or gt[5] >= 0:
		raise ValueError(f'{os.path.basename(path)} is not a north-up grid')
	native = abs(gt[1])
	square = abs(abs(gt[1]) - abs(gt[5])) <= 1e-6 * abs(gt[1])
	warp_to = None
	if not srs.IsProjected():
		cx = gt[0] + gt[1] * ds.RasterXSize / 2
		cy = gt[3] + gt[5] * ds.RasterYSize / 2
		warp_to = utm_for(cx, cy)
		res = res or nice(abs(gt[1]) * 111320 * math.cos(math.radians(cy)) * 0.999)
		notes.append(f'in geographic coordinates: projected to {crs_label(warp_to)} with pixels of {res:g} m (the most '
		             f'common class of each)')
	elif not square or (res and abs(res - native * (srs.GetLinearUnits() or 1)) > 1e-6 * res):
		warp_to = srs
		res = res or min(abs(gt[1]), abs(gt[5])) * (srs.GetLinearUnits() or 1)
		notes.append(f'resampled to pixels of {res:g} m (the most common class of each)')
	if warp_to is not None:
		unit = warp_to.GetLinearUnits() or 1.0
		# its size before warping it (gdal.Warp makes the whole raster in memory)
		x0, x1 = gt[0], gt[0] + gt[1] * ds.RasterXSize
		y0, y1 = gt[3] + gt[5] * ds.RasterYSize, gt[3]
		b = osr.CoordinateTransformation(srs, warp_to).TransformBounds(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1),
		                                                                21) if not srs.IsSame(warp_to) else (x0, y0, x1, y1)
		check_size(int((b[2] - b[0]) * unit / res) + 2, int((b[3] - b[1]) * unit / res) + 2,
		           f'{os.path.basename(path)} at {res:g} m', mode)
		nd = ds.GetRasterBand(1).GetNoDataValue()
		w = gdal.Warp('', ds, format='MEM', dstSRS=warp_to.ExportToWkt(), xRes=res / unit, yRes=res / unit,
		              resampleAlg='mode', targetAlignedPixels=True, srcNodata=nd, dstNodata=nd if nd is not None else -9999,
		              outputType=gdal.GDT_Float64 if ds.GetRasterBand(1).DataType in (gdal.GDT_Float32, gdal.GDT_Float64)
		              else gdal.GDT_Int32)
		ds = w
		gt = ds.GetGeoTransform()
		srs = warp_to
	check_size(ds.RasterXSize, ds.RasterYSize, os.path.basename(path), mode)
	b = ds.GetRasterBand(1)
	a = b.ReadAsArray()
	nd = b.GetNoDataValue()
	valid = np.isfinite(a) if a.dtype.kind == 'f' else np.ones(a.shape, bool)
	if nd is not None and not (isinstance(nd, float) and math.isnan(nd)):
		valid &= a != nd
	if a.dtype.kind == 'f':
		if np.any(a[valid] != np.round(a[valid])):
			raise ValueError(f'{os.path.basename(path)} has values that are not whole numbers: not a map of classes')
		a = np.where(valid, a, 0).astype(np.int64)
	unit = srs.GetLinearUnits() or 1.0
	grid = Grid(gt[0], gt[3], abs(gt[1]), ds.RasterXSize, ds.RasterYSize, srs)
	grid.metres = abs(gt[1]) * unit
	return a.astype(np.int64), valid, grid, notes


def boundary_mask(boundary, grid):
	"""The cells of GRID inside the polygons of the boundary (their centres)."""
	a, n = rasterize(boundary[0], boundary[1], grid, burn=1)
	if not n:
		raise ValueError(f'{os.path.basename(boundary[0])} has no polygons')
	return a == 1


# A run ----------------------------------------------------------------------------------------------------------------

def _check_spec(spec):
	name = str(spec.get('name') or '').strip()
	if not NAME_RE.match(name):
		raise ValueError(f'the name of the run must have letters, digits, - or _ only (up to 80): {name!r}')
	items = spec.get('inputs') or []
	if not items:
		raise ValueError('no maps of patches')
	labels = [str(i.get('label') or '').strip() for i in items]
	if any(not x for x in labels):
		raise ValueError('every map needs a label (e.g. its year)')
	if len(set(labels)) != len(labels):
		raise ValueError('the labels of the maps must be different: ' + ', '.join(sorted({x for x in labels if
		                                                                                  labels.count(x) > 1})))
	mode = spec.get('mode') or 'habitat'
	if mode not in ('habitat', 'classes'):
		raise ValueError(f"mode must be 'habitat' or 'classes', not {mode!r}")
	for it in items:
		p = user_path(it.get('file'))
		if not os.path.isfile(p):
			raise ValueError(f'{p} not found')
		if mode == 'habitat' and (is_raster(p) or it.get('field')) and not it.get('habitat'):
			raise ValueError(f"{os.path.basename(p)}: choose which {'classes' if is_raster(p) else 'values'} are habitat")
	d = int(spec.get('directions') or 8)
	if d not in (4, 8):
		raise ValueError('the neighbours of the patches must be 8 or 4')
	e = int(spec.get('edge_depth') or 1)
	if e < 1:
		raise ValueError('the edge depth must be at least 1 cell')
	metrics = spec.get('metrics')
	if metrics:
		bad = [m for m in metrics if m not in lsm.KEYS]
		if bad:
			raise ValueError('unknown metrics: ' + ', '.join(bad))
	res = spec.get('res')
	if res not in (None, '', 0) and float(res) <= 0:
		raise ValueError('the pixel size must be positive')
	return name, mode, d, e


def field_codes(items):
	"""The codes of the classes of the columns of the polygons of a run in mode classes, the same in every map, whatever
	the name of its column (e.g. Code_00 and Code_18 of CORINE): {column: (codes {value as text: code}, code of the cells
	without polygon, legend)}, the same for every column. When a column has text, every value is numbered from 1 in
	alphabetical order; when all are whole numbers, they keep their value, and the cells without polygon get 0 or, when 0
	is a class, a code below every class."""
	values, text = {}, {}
	for it in items:
		path, field = user_path(it['file']), it.get('field') or None
		if is_raster(path) or not field:
			continue
		ds = ogr.Open(path)
		lyr = _layer(ds, it.get('layer') or None)
		i = lyr.GetLayerDefn().GetFieldIndex(field)
		if i < 0:
			raise ValueError(f'{os.path.basename(path)} has no column {field}')
		is_text = lyr.GetLayerDefn().GetFieldDefn(i).GetType() not in (ogr.OFTInteger, ogr.OFTInteger64)
		text[field] = text.get(field, False) or is_text
		vals = values.setdefault(field, set())
		for f in lyr:
			v = f.GetField(i)
			if v is not None:
				vals.add(str(v))
	every = set().union(*values.values()) if values else set()
	if any(text.values()):
		codes = {v: k + 1 for k, v in enumerate(sorted(every))}
		background = 0
	else:
		codes = {v: int(v) for v in every}
		background = 0 if not codes or min(codes.values()) > 0 else min(codes.values()) - 1
	legend = {background: '(no polygon)', **{c: v for v, c in codes.items()}}
	return {field: (codes, background, legend) for field in values}


def _classes_map(it, mode, grid_v, res_v, res, boundary, codes=None):
	"""(classes, valid, grid, legend {code: name}, notes) of a map of the run: a raster (with the pixel size RES asked
	for, else its own) or polygons (on the grid GRID_V of RES_V metres; CODES: those of field_codes)."""
	path = user_path(it['file'])
	notes = []
	if is_raster(path):
		a, valid, grid, notes = read_raster(path, res, mode)
		legend = None
		if mode == 'habitat':
			hab = np.isin(a, np.asarray([int(v) for v in it['habitat']], dtype=np.int64))
			a = np.where(hab, 1, 0)
	else:
		grid = grid_v
		field = it.get('field') or None
		layer = it.get('layer') or None
		if mode == 'habitat':
			a, n = rasterize(path, layer, grid, burn=1, field=field, where_values=it.get('habitat') if field else None)
			legend = None
		elif field:
			cmap, background, legend = (codes or field_codes([it]))[field]
			a, n = rasterize(path, layer, grid, field=field, codes=cmap, fill=background)
		else:
			a, n = rasterize(path, layer, grid, burn=1)
			legend = {0: '(no polygon)', 1: 'polygons'}
		if not n:
			notes.append('no polygon was drawn (none is habitat, or they are outside the grid)')
		valid = np.ones(a.shape, bool)
		grid.metres = res_v
	if boundary:
		inside = boundary_mask(boundary, grid)
		if not inside.any():
			raise ValueError(f'{os.path.basename(path)} does not overlap the boundary')
		valid &= inside
	if mode == 'habitat':
		legend = {0: 'not habitat', 1: 'habitat'}
	return a, valid, grid, legend, notes


def write_classes(path, a, valid, grid):
	big = a.size and (a[valid].max(initial=0) >= NODATA or a[valid].min(initial=0) < 0)
	dtype, nd = (gdal.GDT_Int32, -2147483648) if big else (gdal.GDT_Byte, NODATA)
	ds = gdal.GetDriverByName('GTiff').Create(path, grid.width, grid.height, 1, dtype, options=GTIFF)
	ds.SetGeoTransform(grid.gt)
	ds.SetProjection(grid.srs.ExportToWkt())
	b = ds.GetRasterBand(1)
	b.SetNoDataValue(nd)
	b.WriteArray(np.where(valid, a, nd))
	ds.FlushCache()
	ds = None


def write_patches(path, labels, grid, rows, legend, map_label):
	"""The patches (LABELS: their numbers on the grid, 0 = none) as polygons in a GeoPackage, one multipolygon per patch
	with its class, its number and its patch metrics (ROWS)."""
	src = gdal.GetDriverByName('MEM').Create('', grid.width, grid.height, 1, gdal.GDT_Int32)
	src.SetGeoTransform(grid.gt)
	src.SetProjection(grid.srs.ExportToWkt())
	band = src.GetRasterBand(1)
	band.WriteArray(labels)
	mask = gdal.GetDriverByName('MEM').Create('', grid.width, grid.height, 1, gdal.GDT_Byte)
	mask.GetRasterBand(1).WriteArray((labels > 0).astype(np.uint8))
	mem = ogr.GetDriverByName('Memory').CreateDataSource('') if gdal.GetDriverByName('MEM').GetMetadataItem(
		gdal.DCAP_VECTOR) != 'YES' else gdal.GetDriverByName('MEM').Create('', 0, 0, 0, gdal.GDT_Unknown)
	tmp = mem.CreateLayer('t', srs=grid.srs, geom_type=ogr.wkbPolygon)
	tmp.CreateField(ogr.FieldDefn('patch', ogr.OFTInteger))
	# with 4 neighbours the polygons are valid (with 8, the rings of large patches touch themselves and repairing them
	# takes minutes); the parts of a patch of 8 neighbours are joined below by its number
	gdal.Polygonize(band, mask.GetRasterBand(1), tmp, 0, [])
	parts = {}
	for f in tmp:
		parts.setdefault(f.GetField('patch'), []).append(f.GetGeometryRef().Clone())
	ds = gdal.GetDriverByName('GPKG').Create(path, 0, 0, 0, gdal.GDT_Unknown)
	lyr = ds.CreateLayer('patches', srs=grid.srs, geom_type=ogr.wkbMultiPolygon)
	fields = [('map', ogr.OFTString), ('patch', ogr.OFTInteger), ('class', ogr.OFTInteger64),
	          ('class_name', ogr.OFTString)] + [(m.split('_', 1)[1], ogr.OFTReal) for m in PATCH_METRICS]
	for name, t in fields:
		lyr.CreateField(ogr.FieldDefn(name, t))
	by_id = {r['id']: r for r in rows}
	lyr.StartTransaction()
	gdal.PushErrorHandler('CPLQuietErrorHandler')
	try:
		for pid, geoms in sorted(parts.items()):
			g = ogr.Geometry(ogr.wkbMultiPolygon)
			for p in geoms:
				if not p.IsValid():
					p = p.MakeValid()
				p = ogr.ForceToMultiPolygon(p)
				for i in range(p.GetGeometryCount()):
					g.AddGeometry(p.GetGeometryRef(i))
			feat = ogr.Feature(lyr.GetLayerDefn())
			feat.SetGeometry(g)
			r = by_id.get(pid, {})
			feat.SetField('map', str(map_label))
			feat.SetField('patch', int(pid))
			if 'class' in r:
				feat.SetField('class', int(r['class']))
				feat.SetField('class_name', str((legend or {}).get(r['class'], r['class'])))
			for m in PATCH_METRICS:
				v = r.get(m.split('_', 1)[1])
				if v is not None:
					feat.SetField(m.split('_', 1)[1], float(v))
			lyr.CreateFeature(feat)
	finally:
		gdal.PopErrorHandler()
	lyr.CommitTransaction()
	ds = mem = src = mask = None


def _write_table(path, header, rows):
	with open(path, 'w', encoding='utf-8') as f:
		f.write(' '.join(header) + '\n')
		for r in rows:
			f.write(' '.join(_text(r.get(h)) if h in ('map', 'class_name', 'file') else _num(r.get(h)) for h in header)
			        + '\n')


def run(spec, plots=True):
	"""Compute the metrics of the maps of SPEC (see main) into landscape/<name>; returns its folder. A run that fails
	leaves the previous results of that name as they were, and no folder <name>.tmp."""
	name = _check_spec(spec)[0]
	try:
		return _run(spec, plots)
	except BaseException:
		shutil.rmtree(os.path.join(OUT, name) + '.tmp', ignore_errors=True)
		raise


def _run(spec, plots):
	t0 = time.time()
	name, mode, directions, edge_depth = _check_spec(spec)
	items = spec['inputs']
	metrics = list(spec.get('metrics') or FRAGMENTATION)
	if mode == 'habitat':
		metrics = [m for m in metrics if not m.startswith('l_')]
	metrics = list(dict.fromkeys(PATCH_METRICS + metrics))
	res = float(spec['res']) if spec.get('res') not in (None, '', 0) else None
	b = spec.get('boundary') or None
	boundary = (user_path(b['file'] if isinstance(b, dict) else b), (b.get('layer') if isinstance(b, dict) else None)
	            or None) if b else None
	if boundary and not os.path.isfile(boundary[0]):
		raise ValueError(f'the boundary {boundary[0]} does not exist')
	classes_max = spec.get('classes_max') or None
	vectors = [(user_path(i['file']), i.get('layer') or None) for i in items if not is_raster(user_path(i['file']))]
	grid_v, res_v = (vector_grid(vectors, boundary, res, mode) if vectors else (None, res))
	final = os.path.join(OUT, name)
	outdir = final + '.tmp'  # written aside: a run that fails keeps the previous results
	if os.path.isdir(outdir):
		shutil.rmtree(outdir)
	os.makedirs(outdir)
	print(f"Landscape metrics of {len(items)} map(s), mode {mode} ({'the habitat as one class' if mode == 'habitat' else 'every class'}), "
	      f'patches of {directions} neighbours, edge depth {edge_depth} cell(s)', flush=True)
	if grid_v is not None:
		print(f'Polygons on a grid of {res_v:g} m, {grid_v.width} x {grid_v.height} pixels, {crs_label(grid_v.srs)}'
		      f"{' (the extent of the boundary)' if boundary else ' (the rectangle around all the polygon maps)'}", flush=True)
	if boundary:
		print(f'Boundary of the landscape: {rel(boundary[0])}', flush=True)
	rec_maps, land_rows, class_rows, patch_rows = [], [], [], []
	codes = field_codes(items) if mode == 'classes' else None  # the same codes of the classes in every map
	for n, it in enumerate(items, 1):
		t1 = time.time()
		label = str(it['label']).strip()
		path = user_path(it['file'])
		print(f'{n}/{len(items)} {label}: {os.path.basename(path)}', flush=True)
		a, valid, grid, legend, notes = _classes_map(it, mode, grid_v, res_v, res, boundary, codes)
		for x in notes:
			print(f'  {x}', flush=True)
		cell_m = getattr(grid, 'metres', None) or grid.res
		if not valid.any():
			raise ValueError(f'{label}: the map has no cells with data inside the landscape')
		out = lsm.compute(a, valid, cell_m, metrics, directions=directions, edge_depth=edge_depth,
		                  classes_max=classes_max, only_classes=[1] if mode == 'habitat' else None, keep=True)
		L = out.pop('landscape_object')
		if mode == 'habitat' and not out['class']:  # no habitat at this date: its amounts are 0, the rest NA
			out['class'] = [lsm.absent_class(1, metrics)]
		if legend is None:
			legend = {c: str(c) for c in out['classes']}
		cls_file, gp = f'classes_{n}.tif', f'patches_{n}.gpkg'
		write_classes(os.path.join(outdir, cls_file), a, valid, grid)
		write_patches(os.path.join(outdir, gp), L.patches[0][1:-1, 1:-1], grid, out['patch'], legend, label)
		for r in out['patch']:
			patch_rows.append(dict(r, map=label, class_name=legend.get(r['class'], r['class'])))
		for r in out['class']:
			class_rows.append(dict(r, map=label, class_name=legend.get(r['class'], r['class'])))
		if out['landscape']:
			land_rows.append(dict(out['landscape'], map=label))
		ha = out['cells'] * cell_m * cell_m / 1e4
		summary = {'map': label, 'file': path, 'kind': 'raster' if is_raster(path) else 'vector',
		           'layer': it.get('layer') or None, 'field': it.get('field') or None, 'habitat': it.get('habitat'),
		           'grid': dict(grid.record(), res_m=cell_m), 'cells': out['cells'], 'landscape_ha': round(ha, 4),
		           'classes': [{'class': c, 'name': legend.get(c, str(c))} for c in out['classes']],
		           'patches': len(out['patch']), 'classes_tif': cls_file, 'patches_gpkg': gp, 'notes': notes,
		           'seconds': round(time.time() - t1, 1)}
		hab = next((r for r in out['class'] if r['class'] == 1), None) if mode == 'habitat' else None
		if hab:
			summary['habitat_metrics'] = {'ca': hab.get('ca'), 'pland': hab.get('pland'), 'np': hab.get('np')}
		rec_maps.append(summary)
		print(f"  {out['cells']:,} cells of {cell_m:g} m ({ha:,.1f} ha), {len(out['classes'])} classes, "
		      f"{len(out['patch']):,} patches ({summary['seconds']} s)", flush=True)
	# RPR: the classes of each map / those of all the maps of the run (unless the spec gives the maximum)
	if land_rows and 'l_rpr' in metrics and not classes_max:
		every = {c['name'] for m in rec_maps for c in m['classes']}
		for row, m in zip(land_rows, rec_maps):
			row['rpr'] = 100.0 * len(m['classes']) / len(every)
	# the maps of a run are comparable when they cover the same landscape at the same pixel size
	warnings = []
	sizes = {m['grid']['res_m'] for m in rec_maps}
	areas = [m['landscape_ha'] for m in rec_maps]
	if len(sizes) > 1:
		warnings.append(f"the maps have different pixel sizes ({', '.join(f'{s:g}' for s in sorted(sizes))} m): the "
		                f"metrics of edges, shapes and patches change with the pixel, so their differences are not only "
		                f"changes of the landscape (choose one pixel size)")
	if areas and max(areas) > min(areas) * 1.01:
		warnings.append(f'the landscapes of the maps have different areas ({min(areas):,.1f} to {max(areas):,.1f} ha): '
		                f'PLAND, PD, ED, LPI, MESH... are relative to the landscape; give a boundary to compare the same '
		                f'area')
	for w in warnings:
		print('WARNING: ' + w, flush=True)
	pm = [m.split('_', 1)[1] for m in metrics if m.startswith('p_')]
	cm = [m.split('_', 1)[1] for m in metrics if m.startswith('c_')]
	lmx = [m.split('_', 1)[1] for m in metrics if m.startswith('l_')]
	_write_table(os.path.join(outdir, 'metrics_patch.csv'), ['map', 'class', 'class_name', 'id', 'cells'] + pm,
	             patch_rows)
	_write_table(os.path.join(outdir, 'metrics_class.csv'), ['map', 'class', 'class_name'] + cm, class_rows)
	tables = {'patch': 'metrics_patch.csv', 'class': 'metrics_class.csv'}
	if land_rows:
		_write_table(os.path.join(outdir, 'metrics_landscape.csv'), ['map'] + lmx, land_rows)
		tables['landscape'] = 'metrics_landscape.csv'
	plot_files = []
	if plots and len(items) > 1 and os.path.isfile(PLOTS_R):
		cmd = ['Rscript', PLOTS_R, f'--dir={outdir}', f'--mode={mode}']
		print('Plots (landscape_plots.R)...', flush=True)
		try:
			r = subprocess.run(cmd, cwd=HERE)
			if r.returncode == 0:
				plot_files = sorted(f for f in os.listdir(outdir) if f.lower().endswith('.png'))
			else:
				print(f'WARNING: landscape_plots.R failed (exit code {r.returncode}): no plots', flush=True)
		except OSError as e:
			print(f'WARNING: landscape_plots.R could not be run ({e}): is R (Rscript) installed?', flush=True)
	record = {'version': 1, 'date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'name': name, 'mode': mode,
	          'options': {'directions': directions, 'edge_depth': edge_depth, 'res': res_v,
	                      'classes_max': classes_max or (len({c['name'] for m in rec_maps for c in m['classes']})
	                                                     if mode == 'classes' else None),
	                      'count_boundary': False, 'consider_boundary': False, 'neighbourhood': 4, 'base': 'log2'},
	          'boundary': {'file': rel(boundary[0]), 'layer': boundary[1]} if boundary else None,
	          'vector_grid': grid_v.record() if grid_v is not None else None,
	          'metrics': metrics, 'maps': rec_maps,
	          'tables': tables, 'plots': plot_files, 'warnings': warnings,
	          'references': [lsm.REFERENCE, lsm.FRAGSTATS], 'seconds': round(time.time() - t0, 1)}
	write_metadata(os.path.join(outdir, RUN_TXT), record)
	with open(os.path.join(outdir, SPEC_JSON), 'w', encoding='utf-8') as f:
		json.dump(spec, f, indent=1, ensure_ascii=False)
	with open(os.path.join(outdir, RUN_JSON), 'w', encoding='utf-8') as f:
		json.dump(record, f, indent=1, ensure_ascii=False)
	if os.path.isdir(final):
		shutil.rmtree(final)
	os.replace(outdir, final)
	return final


def write_metadata(path, rec):
	o = rec['options']
	lines = ['eHabitat+ landscape metrics: options and inputs of the files in this folder (landscape_metrics.py)',
	         f"Date: {rec['date']}", f"Run: {rec['name']}",
	         'Mode: ' + ('habitat: the classes or polygons chosen as habitat are class 1, the rest of the landscape 0; '
	                     'the metrics are those of the habitat (patches and class 1)' if rec['mode'] == 'habitat' else
	                     'classes: every class of the maps (the same codes in every map; for polygons, the cells inside '
	                     'the landscape without a polygon are the class "(no polygon)"), and the landscape level; RPR '
	                     f"relative to {o.get('classes_max')} classes"), '',
	         '## Options (as the R package landscapemetrics, calculate_lsm)',
	         f"Patches of {o['directions']} neighbours (directions = {o['directions']})",
	         f"Edge depth of the core areas: {o['edge_depth']} cell(s) (edge_depth); cells next to cells without data or "
	         f"to the edge of the map are edge too (consider_boundary = FALSE)",
	         'Total edge (TE, ED): only the edges between classes, not those with cells without data nor with the edge '
	         'of the map (count_boundary = FALSE)',
	         'Adjacencies of CONTAG, IJI, PLADJ and the entropies: the 4 neighbours of each cell (neighbourhood = 4), '
	         'entropies in bits (log2)',
	         f"Pixel size of the polygons: {o['res']:g} m" if o.get('res') else 'Pixel size: that of each raster',
	         f"Boundary of the landscape: {rec['boundary']['file']}" if rec['boundary'] else
	         'Boundary of the landscape: none (every cell with data of each map; for the polygons, the rectangle around '
	         'all the polygon maps)', '', '## Maps (map label, file, pixel, landscape, patches)']
	for m in rec['maps']:
		hab = f", habitat: {', '.join(str(x) for x in m['habitat'])}" if m.get('habitat') else ''
		lines.append(f"{m['map']}: {rel(m['file'])}{' (layer ' + m['layer'] + ')' if m.get('layer') else ''}"
		             f"{', column ' + m['field'] if m.get('field') else ''}{hab}; {m['grid']['res_m']:g} m, "
		             f"{m['grid']['width']} x {m['grid']['height']} pixels, {m['grid']['crs']}; {m['cells']} cells "
		             f"({m['landscape_ha']} ha); {m['patches']} patches")
		lines += [f'  {n}' for n in m.get('notes') or []]
	if rec['warnings']:
		lines += ['', '## Warnings'] + rec['warnings']
	lines += ['', '## Tables', 'metrics_patch.csv: one row per patch of each map (map, class, id, cells, patch metrics)',
	          'metrics_class.csv: one row per class of each map', 'metrics_landscape.csv: one row per map (mode classes)',
	          'patches_<n>.gpkg: the patches of the map n as polygons, with their patch metrics',
	          'classes_<n>.tif: the map n as it was measured (habitat 1, not habitat 0 in mode habitat)', '',
	          '## Metrics (key, name, unit)']
	for m in rec['metrics']:
		c = lsm.KEYS[m]
		lines.append(f"{m}: {c['name']} ({c['unit']})")
	lines += ['', '## References', 'The metrics are computed as in the R package landscapemetrics 2.2 (checked with it, '
	          'segm_optim/tests/test_lsm.py), which follows FRAGSTATS; CIRCLE takes the smallest circle that contains '
	          'every corner of the cells of a patch.'] + rec['references'] + ['', f"Time: {rec['seconds']} s"]
	with open(path, 'w', encoding='utf-8') as f:
		f.write('\n'.join(lines) + '\n')


def runs():
	"""The finished runs (folders of landscape/ with landscape_run.json), newest first."""
	out = []
	if not os.path.isdir(OUT):
		return out
	for d in os.listdir(OUT):
		p = os.path.join(OUT, d, RUN_JSON)
		if NAME_RE.match(d) and os.path.isfile(p):
			try:
				with open(p, encoding='utf-8') as f:
					rec = json.load(f)
			except (OSError, ValueError):
				continue
			out.append({'name': d, 'folder': os.path.join(OUT, d), 'date': rec.get('date'), 'mode': rec.get('mode'),
			            'maps': [m.get('map') for m in rec.get('maps') or []], 'mtime': os.path.getmtime(p)})
	return sorted(out, key=lambda x: x['mtime'], reverse=True)


def main():
	ap = argparse.ArgumentParser(description='Landscape metrics (FRAGSTATS / landscapemetrics) of maps of patches.')
	ap.add_argument('files', nargs='*', help='rasters of classes or polygon files, one per date (in order)')
	ap.add_argument('--spec', help='a JSON file with the inputs and options of a run (as the app writes it)')
	ap.add_argument('--name', help='name of the run: its folder in landscape/')
	ap.add_argument('--labels', nargs='+', help='label of each file (default: the year in its name, else its name)')
	ap.add_argument('--mode', choices=('habitat', 'classes'), default='habitat')
	ap.add_argument('--habitat', nargs='+', help='values (classes of the rasters, or of --field) that are habitat')
	ap.add_argument('--field', help='column of the classes of the polygons')
	ap.add_argument('--res', type=float, help='pixel size (m) of the polygons, or to resample the rasters')
	ap.add_argument('--boundary', help='polygons of the landscape (cells outside are left out)')
	ap.add_argument('--neighbours', type=int, choices=(4, 8), default=8)
	ap.add_argument('--edge-depth', type=int, default=1, help='edge depth of the core areas, in cells')
	ap.add_argument('--metrics', nargs='+', help="keys of the metrics (c_pland, l_shdi...; 'all' = every one)")
	ap.add_argument('--no-plots', action='store_true')
	ap.add_argument('--inspect', help='print what a file holds (JSON)')
	ap.add_argument('--layer', help='layer of --inspect')
	ap.add_argument('--catalogue', action='store_true', help='print the metrics (JSON)')
	a = ap.parse_args()
	try:
		if a.catalogue:
			print(json.dumps(lsm.CATALOGUE, indent=1, ensure_ascii=False))
			return
		if a.inspect:
			print(json.dumps(inspect(a.inspect, a.layer), indent=1, ensure_ascii=False))
			return
		if a.spec:
			with open(a.spec, encoding='utf-8') as f:
				spec = json.load(f)
		else:
			if not a.files:
				ap.error('give the maps (files) or --spec')
			labels = a.labels or []
			if labels and len(labels) != len(a.files):
				ap.error('give one label per file')
			items = []
			for i, f in enumerate(a.files):
				info = {'name': os.path.splitext(os.path.basename(f))[0]}
				m = re.search(r'(?<!\d)(1[89]\d\d|20\d\d|21\d\d)(?!\d)', info['name'])
				items.append({'file': f, 'label': labels[i] if labels else (m.group(1) if m else info['name']),
				              'field': a.field, 'habitat': a.habitat})
			metrics = None if not a.metrics else lsm.keys() if a.metrics == ['all'] else a.metrics
			spec = {'name': a.name or 'landscape', 'inputs': items, 'mode': a.mode, 'res': a.res,
			        'boundary': a.boundary, 'directions': a.neighbours, 'edge_depth': a.edge_depth, 'metrics': metrics}
		folder = run(spec, plots=not a.no_plots)
	except (ValueError, RuntimeError) as e:
		sys.exit(f'ERROR: {e}')
	print(f'Results in {rel(folder)}')


if __name__ == '__main__':
	os.chdir(HERE)
	main()
