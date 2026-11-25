#### eHabitat+: how much of each study area the variables cover
#### License: CC BY-SA 3.0

"""
The part of each study area that has data in the variables: inside the extent of all of them and with data (not NoData)
in the first one, which defines the pixels analysed. A study area that the variables only partly cover gets HFTs only
where they have data, so the app ("Input variables") and setup.py warn about it.

To stay fast with large variables, the first variable is read over the rectangle of each study area at a reduced
resolution (at most MAX_SIZE pixels a side; not from its overviews, whose averages fill the NoData), and the study
area is rasterized on that same reduced grid.
"""

import math
import os

import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

MAX_SIZE = 1000  # pixels a side of the reduced reading of the first variable over a study area
FULL = 0.98  # from this fraction on, a study area is taken as covered (the pixels of its edge, read reduced)
MAX_NOTES = 8  # notes on single study areas; the others in one line


def _srs(s):
	s = s.Clone()
	s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return s


def _box(xmin, ymin, xmax, ymax):
	ring = ogr.Geometry(ogr.wkbLinearRing)
	for x, y in ((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)):
		ring.AddPoint_2D(x, y)
	g = ogr.Geometry(ogr.wkbPolygon)
	g.AddGeometry(ring)
	return g


def _extent(ds):
	gt = ds.GetGeoTransform()
	x0, x1 = gt[0], gt[0] + gt[1] * ds.RasterXSize
	y0, y1 = gt[3], gt[3] + gt[5] * ds.RasterYSize
	return min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)


def _areas(vector, layer, col_id, ids, label_col, srs):
	"""{id: (label, geometry in SRS)} of the study areas, as setup.py takes them: the integer column COL_ID or, with
	'cat' when the file has no such column, 1, 2, ... in the order of the features. Only IDS, if given."""
	ds = ogr.Open(vector)
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	defn = lyr.GetLayerDefn()
	i = defn.GetFieldIndex(col_id or 'cat')
	if i < 0 and (col_id or 'cat') != 'cat':
		raise ValueError(f'{os.path.basename(vector)} has no column {col_id}')
	if i >= 0 and defn.GetFieldDefn(i).GetType() not in (ogr.OFTInteger, ogr.OFTInteger64):
		raise ValueError(f'the column {col_id} is not of integers')
	j = defn.GetFieldIndex(label_col) if label_col else -1
	src = lyr.GetSpatialRef()
	ct = osr.CoordinateTransformation(_srs(src), _srs(srs)) if src is not None and not src.IsSame(srs) else None
	want = {int(x) for x in ids} if ids else None
	labels, parts = {}, {}
	for n, f in enumerate(lyr, 1):
		if i >= 0 and not f.IsFieldSetAndNotNull(i):  # GRASS leaves out the features without ID
			continue
		fid = f.GetFieldAsInteger64(i) if i >= 0 else n
		g = f.GetGeometryRef()
		if (want is not None and fid not in want) or g is None or g.IsEmpty():
			continue
		g = ogr.ForceToMultiPolygon(g.GetLinearGeometry() if g.HasCurveGeometry() else g.Clone())
		if g is None or g.IsEmpty():
			continue
		if ct is not None:
			g.Transform(ct)
		if not g.IsValid():
			g = g.MakeValid() if hasattr(g, 'MakeValid') else g.Buffer(0)  # as v.import, both halves of a bow tie
		if not labels.get(fid) and j >= 0 and f.IsFieldSetAndNotNull(j):
			labels[fid] = f.GetFieldAsString(j).strip() or None
		parts.setdefault(fid, []).append(g)
	out = {}
	for fid, gs in parts.items():
		m = ogr.Geometry(ogr.wkbGeometryCollection)
		for g in gs:
			m.AddGeometry(g)
		u = m.UnaryUnion() if hasattr(m, 'UnaryUnion') else m.UnionCascaded()  # once, not part after part
		out[fid] = (labels.get(fid), u)
	return out


def _valid(a, nodata):
	ok = np.ones(a.shape, bool) if nodata is None else a != nodata
	if a.dtype.kind == 'f':
		ok &= ~np.isnan(a)
	return ok


def _data_fraction(ds, g, max_size):
	"""Fraction of the geometry G (inside the extent of DS) with data in the first band of DS, read reduced."""
	gt = ds.GetGeoTransform()
	band = ds.GetRasterBand(1)
	nodata = band.GetNoDataValue()
	xmin, xmax, ymin, ymax = g.GetEnvelope()
	c0 = max(0, int(math.floor((xmin - gt[0]) / gt[1])))
	c1 = min(ds.RasterXSize, int(math.ceil((xmax - gt[0]) / gt[1])))
	r0 = max(0, int(math.floor((ymax - gt[3]) / gt[5])))
	r1 = min(ds.RasterYSize, int(math.ceil((ymin - gt[3]) / gt[5])))
	w, h = c1 - c0, r1 - r0
	if w <= 0 or h <= 0:
		return 0.0
	step = max(1, math.ceil(max(w, h) / max_size))
	bx, by = math.ceil(w / step), math.ceil(h / step)
	a = band.ReadAsArray(c0, r0, w, h, buf_xsize=bx, buf_ysize=by)
	ok = _valid(a, nodata)
	# the study area on the same reduced grid
	mem = gdal.GetDriverByName('MEM').Create('', bx, by, 1, gdal.GDT_Byte)
	mem.SetGeoTransform((gt[0] + c0 * gt[1], gt[1] * w / bx, 0, gt[3] + r0 * gt[5], 0, gt[5] * h / by))
	mem.SetProjection(ds.GetProjection())
	vds = (ogr.GetDriverByName('MEM') or ogr.GetDriverByName('Memory')).CreateDataSource('')  # MEM since GDAL 3.11
	vl = vds.CreateLayer('a', ds.GetSpatialRef(), ogr.wkbUnknown)
	feat = ogr.Feature(vl.GetLayerDefn())
	feat.SetGeometry(g)
	vl.CreateFeature(feat)
	gdal.RasterizeLayer(mem, [1], vl, burn_values=[1])
	inside = mem.ReadAsArray().astype(bool)
	n = int(inside.sum())
	if n == 0:  # smaller than a reduced pixel: the pixel of a point of it, at full resolution
		p = g.PointOnSurface()
		c = min(ds.RasterXSize - 1, max(0, int((p.GetX() - gt[0]) / gt[1])))
		r = min(ds.RasterYSize - 1, max(0, int((p.GetY() - gt[3]) / gt[5])))
		return float(_valid(band.ReadAsArray(c, r, 1, 1), nodata)[0, 0])
	return int((ok & inside).sum()) / n


def study_coverage(rasters, vector, layer=None, col_id='cat', ids=None, label_col=None, max_size=MAX_SIZE):
	"""How much of each study area of VECTOR (the IDS of its column COL_ID, all of them without IDS) the GeoTIFFs
	RASTERS cover: [{id, label, km2, inside (fraction inside the extent of all of them), covered (fraction with data in
	the first one too)}]. RASTERS must have the same CRS (the first one gives it), projected; [] when they are not
	north-up or have no CRS."""
	if not rasters:
		return []
	dss = [gdal.OpenEx(r, gdal.OF_RASTER, open_options=['OVERVIEW_LEVEL=NONE']) for r in rasters]
	first = dss[0]
	srs = first.GetSpatialRef()
	gt = first.GetGeoTransform()
	if srs is None or gt[2] or gt[4] or gt[1] <= 0 or gt[5] >= 0:  # only north-up
		return []
	xmin, ymin, xmax, ymax = _extent(first)
	for d in dss[1:]:
		x0, y0, x1, y1 = _extent(d)
		xmin, ymin, xmax, ymax = max(xmin, x0), max(ymin, y0), min(xmax, x1), min(ymax, y1)
	ext = _box(xmin, ymin, xmax, ymax) if xmin < xmax and ymin < ymax else None
	factor = (srs.GetLinearUnits() ** 2) / 1e6 if srs.IsProjected() else None
	out = []
	for fid, (label, g) in sorted(_areas(vector, layer, col_id, ids, label_col, srs).items()):
		area = g.GetArea()
		if area <= 0:
			continue
		part = g.Intersection(ext) if ext is not None else None
		inside = part.GetArea() / area if part is not None and not part.IsEmpty() else 0.0
		data = _data_fraction(first, part, max_size) if inside > 0 else 0.0
		out.append({'id': fid, 'label': label, 'km2': round(area * factor, 3) if factor else None,
		            'inside': min(1.0, inside), 'covered': min(1.0, inside * data)})
	return out


def _name(a):
	return f"study area {a['id']}" + (f" ({a['label']})" if a.get('label') else '')


def _pct(f):
	p = f * 100
	return 'less than 1' if 0 < p < 1 else f'{p:.0f}'


def coverage_checks(cov, full=FULL, max_notes=MAX_NOTES):
	"""The notes of study_coverage as those of the checks of the app ({level, text}): an error when a study area has
	no data at all, a warning when it is only partly covered, nothing when the variables cover all of them."""
	none = [a for a in cov if a['covered'] <= 0]
	part = [a for a in cov if 0 < a['covered'] < full]
	out = []
	for a in none[:max_notes]:
		out.append({'level': 'error', 'text': f"The variables have no data in {_name(a)}"
		            + (' (it is outside them)' if a['inside'] <= 0 else '') + ': it will have no HFTs.'})
	for a in part[:max(0, max_notes - len(none))]:
		out.append({'level': 'warning', 'text': f"The variables cover about {_pct(a['covered'])} % of {_name(a)}: "
		                                        'the rest will have no HFTs.'})
	rest = len(none) + len(part) - len(out)
	if rest > 0:
		empty = len(none) > max_notes  # some of those not named have no data at all
		out.append({'level': 'error' if empty else 'warning',
		            'text': f"And {rest} more study area{'s' if rest > 1 else ''} not fully covered by the variables"
		                    + (' (some of them without any data)' if empty else '') + '.'})
	return out
