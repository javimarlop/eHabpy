"""Input variables and study area.

A set of variables is a subfolder of inVars (INVARS_FOLDER in config.py; '' is inVars itself) with the GeoTIFFs of
the study area scale (parcela, the fine pixel size that setup.py imports into GRASS) and, in its subfolder
'ecoregion', those of the ecoregion scale (a coarser pixel size, for the similarity of the HFTs in the ecoregion), as
get_vars_gee.py leaves them. The ecoregion grid covers the RESOLVE ecoregions that the study area touches or a
comparison area of the user (its gee_vars.json says which); when that area does not contain the study area,
'ecoregion/study_area' has the study area at the pixel of the ecoregion grid (for the statistics of the HFTs in the
similarity). A set with only 'ecoregion' (no GeoTIFFs of the study area) is a comparison area for the similarity, not
a set for the HFT model. The variables left out of a set go to its subfolder '_excluded', which setup.py does not
read. The study area is a vector file that the app copies to segm_optim/pa (where moranvar_plots.R looks for it)."""

import base64
import glob
import json
import math
import os
import re
import shutil
import uuid

from osgeo import gdal, ogr, osr

from . import paths

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

RASTER_EXT = ('.tif', '.tiff')
VECTOR_EXT = ('.gpkg', '.shp', '.geojson', '.json')
SHP_PARTS = ('.shp', '.shx', '.dbf', '.prj', '.cpg', '.qix', '.sbn', '.sbx', '.qmd', '.shp.xml')
NAME_RE = re.compile(r'[A-Za-z0-9_-]+(/[A-Za-z0-9_-]+)?')
SAME_CRS = ['IGNORE_DATA_AXIS_TO_SRS_AXIS_MAPPING=YES', 'CRITERION=EQUIVALENT']

# Colours of the quick look of a raster (viridis)
VIRIDIS = ['440154', '482878', '3e4989', '31688e', '26828e', '1f9e89', '35b779', '6ece58', 'b5de2b', 'fde725']


def _wgs84():
	s = osr.SpatialReference()
	s.ImportFromEPSG(4326)
	s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return s


def _crs_label(srs):
	if srs is None:
		return None
	srs = srs.Clone()
	try:
		srs.AutoIdentifyEPSG()
	except RuntimeError:
		pass
	code = srs.GetAuthorityCode(None)
	auth = srs.GetAuthorityName(None)
	name = srs.GetName()
	return f'{auth}:{code} ({name})' if code else name


def _bbox_wgs84(srs, xmin, ymin, xmax, ymax):
	"""[west, south, east, north] of a rectangle of SRS (its edges densified)."""
	if srs is None:
		return None
	try:
		src = srs.Clone()
		src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		ct = osr.CoordinateTransformation(src, _wgs84())
		return list(ct.TransformBounds(xmin, ymin, xmax, ymax, 21))
	except Exception:
		return None


# Rasters -------------------------------------------------------------------------------------------------------------

def raster_info(path):
	"""What the model needs to know of a GeoTIFF: CRS, pixel size, size, NoData, bands, extent."""
	out = {'file': os.path.basename(path), 'path': path, 'name': os.path.splitext(os.path.basename(path))[0]}
	try:
		ds = gdal.Open(path)
	except RuntimeError as e:
		out['error'] = str(e)
		return out
	gt = ds.GetGeoTransform()
	srs = ds.GetSpatialRef()
	band = ds.GetRasterBand(1)
	xmin, ymax = gt[0], gt[3]
	xmax, ymin = gt[0] + gt[1] * ds.RasterXSize, gt[3] + gt[5] * ds.RasterYSize
	out.update({
		'crs': _crs_label(srs), 'projected': bool(srs and srs.IsProjected()), 'wkt': srs.ExportToWkt() if srs else None,
		'res': [abs(gt[1]), abs(gt[5])], 'width': ds.RasterXSize, 'height': ds.RasterYSize, 'bands': ds.RasterCount,
		'nodata': band.GetNoDataValue(), 'dtype': gdal.GetDataTypeName(band.DataType),
		'extent': [min(xmin, xmax), min(ymin, ymax), max(xmin, xmax), max(ymin, ymax)],
		'bbox': _bbox_wgs84(srs, min(xmin, xmax), min(ymin, ymax), max(xmin, xmax), max(ymin, ymax)),
		'size_mb': round(os.path.getsize(path) / 2 ** 20, 2),
	})
	if out['nodata'] is not None and math.isnan(out['nodata']):
		out['nodata'] = 'nan'
	return out


def _fmt(v):
	return f'{v:g}' if isinstance(v, (int, float)) else str(v)


def check_scale(rasters, scale):
	"""Problems of the GeoTIFFs of a scale: [{level: error|warning|ok, text}]. The study area scale needs the same CRS,
	a projected one, and the same square pixel size in all of them (setup.py stops otherwise)."""
	out = []
	good = [r for r in rasters if not r.get('error')]
	for r in rasters:
		if r.get('error'):
			out.append({'level': 'error', 'text': f"{r['file']}: cannot be read ({r['error']})"})
	if not good:
		if not rasters:
			out.append({'level': 'warning' if scale == 'ecoregion' else 'error',
			            'text': 'No variables at this scale.'})
		return out
	first = good[0]
	srs0 = osr.SpatialReference()
	srs0.ImportFromWkt(first['wkt']) if first.get('wkt') else None
	for r in good:
		if not r.get('wkt'):
			out.append({'level': 'error', 'text': f"{r['file']} has no coordinate reference system (CRS)."})
			continue
		srs = osr.SpatialReference()
		srs.ImportFromWkt(r['wkt'])
		if first.get('wkt') and not (srs.IsSame(srs0, SAME_CRS) or srs.ExportToProj4() == srs0.ExportToProj4()):
			out.append({'level': 'error', 'text': f"{r['file']} has another CRS ({r['crs']}) than {first['file']} ({first['crs']})."})
		if not r['projected']:
			out.append({'level': 'error', 'text': f"{r['file']} has geographic coordinates ({r['crs']}): "
			                                      'reproject it to a projected CRS in metres (e.g. UTM).'})
		if not math.isclose(r['res'][0], r['res'][1], rel_tol=1e-6):
			out.append({'level': 'error', 'text': f"{r['file']}: the pixels are not square ({_fmt(r['res'][0])} x {_fmt(r['res'][1])})."})
		elif not math.isclose(r['res'][0], first['res'][0], rel_tol=1e-6):
			out.append({'level': 'error', 'text': f"{r['file']} has another pixel size ({_fmt(r['res'][0])}) than "
			                                      f"{first['file']} ({_fmt(first['res'][0])})."})
		if r['nodata'] is None:
			out.append({'level': 'warning', 'text': f"{r['file']} has no NoData value: all its pixels are taken as data."})
		if r['bands'] > 1:
			out.append({'level': 'warning', 'text': f"{r['file']} has {r['bands']} bands: only the first one is used."})
		if r['extent'] != first['extent']:
			out.append({'level': 'info', 'text': f"{r['file']} covers another extent than {first['file']} "
			                                     '(the region of the study area is used).'})
	if scale == 'study' and len(good) < 3:
		out.append({'level': 'warning', 'text': 'With fewer than 3 variables the segmentation cannot use the PCA (SEGM_PCA).'})
	if not any(o['level'] == 'error' for o in out):
		out.insert(0, {'level': 'ok', 'text': f"{len(good)} variables, {first['crs']}, {_fmt(first['res'][0])} m pixels."})
	# the extent notes only once
	seen, uniq = set(), []
	for o in out:
		if o['level'] == 'info':
			if 'info' in seen:
				continue
			seen.add('info')
		uniq.append(o)
	return uniq


def variable_key(file):
	"""Name of a variable without the number of order of get_vars_gee.py (1_, 01_) and its extension."""
	return re.sub(r'^\d+_', '', os.path.splitext(os.path.basename(file))[0]).lower()


def check_scales(study, eco):
	"""Notes on the two scales together: the similarity in the ecoregion needs the same variables at both."""
	if not eco:
		return [{'level': 'info', 'text': 'No ecoregion variables: the HFT mapping only uses those of the study area; '
		                                  'the similarity in the ecoregion will need them.'}]
	s = {variable_key(r['file']) for r in study}
	e = {variable_key(r['file']) for r in eco}
	out = []
	if s - e:
		out.append({'level': 'warning', 'text': 'Study area variables without their pair in the ecoregion: ' + ', '.join(sorted(s - e))})
	if e - s:
		out.append({'level': 'warning', 'text': 'Ecoregion variables without their pair in the study area: ' + ', '.join(sorted(e - s))})
	if not out:
		out.append({'level': 'ok', 'text': 'Both scales have the same variables.'})
	return out


def _tifs(folder):
	return sorted(f for f in glob.glob(os.path.join(folder, '*')) if f.lower().endswith(RASTER_EXT) and os.path.isfile(f))


tifs = _tifs

STUDY_AREA_FOLDER = 'study_area'  # ecoregion/study_area: the study grid at the pixel of the ecoregion grid
RECORD = 'gee_vars.json'  # what get_vars_gee.py downloaded in a folder
COMPARISON = 'Comparison area'  # BIOME_NAME of the polygons of a comparison area (ECO_ID -1, -2, ...)


def _read_record(folder):
	try:
		with open(os.path.join(folder, RECORD), encoding='utf-8') as f:
			rec = json.load(f)
		return rec if isinstance(rec, dict) else {}
	except (OSError, ValueError):
		return {}


def _polygon_props(gpkg):
	"""ECO_ID, ECO_NAME and BIOME_NAME of the polygons of ecoregions.gpkg."""
	out = []
	try:
		ds = ogr.Open(gpkg)
		lyr = ds.GetLayer(0)
		defn = lyr.GetLayerDefn()
		names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
		for f in lyr:
			out.append({k: f.GetField(k) for k in ('ECO_ID', 'ECO_NAME', 'BIOME_NAME') if k in names})
	except (RuntimeError, AttributeError):
		pass
	return out


def ecoregion_area(eco_dir):
	"""What the ecoregion grid of a set (its folder ECO_DIR) covers: {kind: 'ecoregions' (the RESOLVE ecoregions that
	the study area touches), 'comparison area' (a vector file of the user) or None (unknown: GeoTIFFs of the user),
	name, names (of its polygons), features, file and layer (of a comparison area, as downloaded), polygons
	(ecoregions.gpkg or None)}. A comparison area is known by the key comparison_area of gee_vars.json or, else, by its
	polygons (BIOME_NAME 'Comparison area' or a negative ECO_ID)."""
	rec = _read_record(eco_dir)
	gpkg = os.path.join(eco_dir, 'ecoregions.gpkg')
	out = {'kind': None, 'name': None, 'names': [], 'features': None, 'file': None, 'layer': None,
	       'polygons': gpkg if os.path.isfile(gpkg) else None}
	feats = [e for e in rec.get('ecoregions') or [] if isinstance(e, dict)]
	if not feats and out['polygons']:
		feats = _polygon_props(gpkg)
	names = [str(f['ECO_NAME']) for f in feats if f.get('ECO_NAME')]
	ca = rec.get('comparison_area')
	negative = any(isinstance(f.get('ECO_ID'), (int, float)) and f['ECO_ID'] < 0 for f in feats)
	if isinstance(ca, dict) or negative or any(f.get('BIOME_NAME') == COMPARISON for f in feats):
		ca = ca if isinstance(ca, dict) else {}
		name = ca.get('name') or (names[0] if len(set(names)) == 1 else None) or \
			(os.path.splitext(os.path.basename(ca['file']))[0] if ca.get('file') else 'comparison area')
		out.update(kind='comparison area', name=name, names=names, features=ca.get('features') or len(feats) or None,
		           file=ca.get('file'), layer=ca.get('layer'))
	elif feats:
		out.update(kind='ecoregions', name=', '.join(names) if names else 'RESOLVE ecoregions', names=names,
		           features=len(feats))
	return out


def area_only(d):
	"""Whether the set of the folder D has only the variables of an ecoregion or comparison area (no GeoTIFF of the
	study area, not even left out): not a set for the HFT model."""
	return not _tifs(d) and not _tifs(os.path.join(d, paths.EXCLUDED)) and bool(_tifs(os.path.join(d, paths.ECO_FOLDER)))


def set_dir(folder):
	folder = (folder or '').strip().strip('/\\')
	path = os.path.normpath(os.path.join(paths.INVARS, folder))
	if not paths.inside(path, paths.INVARS):
		raise ValueError(f'{folder} is not a subfolder of inVars')
	return path


def list_sets():
	"""Sets of variables: inVars and its subfolders (two levels) with GeoTIFFs, without 'ecoregion' and '_excluidas'."""
	out = []
	if not os.path.isdir(paths.INVARS):
		return out
	skip = {paths.ECO_FOLDER, paths.EXCLUDED}
	for d, dirs, _ in os.walk(paths.INVARS):
		rel = os.path.relpath(d, paths.INVARS).replace(os.sep, '/')
		depth = 0 if rel == '.' else rel.count('/') + 1
		dirs[:] = sorted(x for x in dirs if x not in skip and not x.startswith('.')) if depth < 2 else []
		folder = '' if rel == '.' else rel
		n = len(_tifs(d))
		eco_dir = os.path.join(d, paths.ECO_FOLDER)
		ne = len(_tifs(eco_dir))
		nx = len(_tifs(os.path.join(d, paths.EXCLUDED))) + len(_tifs(os.path.join(eco_dir, paths.EXCLUDED)))
		if n or ne or nx or (depth == 1 and not os.listdir(d)):
			a = ecoregion_area(eco_dir) if ne else None
			out.append({'folder': folder, 'label': 'inVars' if not folder else folder, 'study': n, 'ecoregion': ne,
			            'excluded': nx, 'area': a['kind'] if a else None, 'area_name': a['name'] if a else None,
			            'area_only': area_only(d),
			            'study_eco': len(_tifs(os.path.join(eco_dir, STUDY_AREA_FOLDER))) if ne else 0})
	return out


def get_set(folder):
	"""The variables of a set at both scales, with the left out ones and the checks."""
	d = set_dir(folder)
	if not os.path.isdir(d):
		raise ValueError(f'folder {paths.rel(d)} not found')
	study = [raster_info(f) for f in _tifs(d)]
	eco = [raster_info(f) for f in _tifs(os.path.join(d, paths.ECO_FOLDER))]
	for r in study + eco:
		r['included'] = True
	ex_study = [dict(raster_info(f), included=False) for f in _tifs(os.path.join(d, paths.EXCLUDED))]
	ex_eco = [dict(raster_info(f), included=False) for f in _tifs(os.path.join(d, paths.ECO_FOLDER, paths.EXCLUDED))]
	eco_dir = os.path.join(d, paths.ECO_FOLDER)
	eco_poly = os.path.join(eco_dir, 'ecoregions.gpkg')
	sa_dir = os.path.join(eco_dir, STUDY_AREA_FOLDER)
	eco_study = [dict(raster_info(f), included=True) for f in _tifs(sa_dir)] if eco else []
	only = area_only(d)
	area = ecoregion_area(eco_dir) if eco else None
	both = check_scales(study, eco)
	if only:
		both = [{'level': 'info', 'text': 'Only the variables of ' + ('a comparison area' if area and area['kind'] ==
		         'comparison area' else 'an ecoregion') + ', without those of a study area: this set is not for the HFT '
		         'model. Compare the HFTs of a scenario with it in "Similarity" (Compare with), when its variables have '
		         'the same names as those of the scenario.'}]
	return {
		'folder': folder or '', 'path': d, 'rel': paths.rel(d),
		'study': study + ex_study, 'ecoregion': eco + ex_eco,
		'checks': {'study': check_scale(study, 'study') if not only else
		           [{'level': 'info', 'text': 'No variables of a study area: a set with only an ecoregion or comparison '
		                                      'area.'}],
		           'ecoregion': check_scale(eco, 'ecoregion') if eco else
		           [{'level': 'info', 'text': 'No variables at this scale (folder ecoregion).'}],
		           'both': both, 'eco_study': check_eco_study(eco, eco_study) if eco_study else []},
		'ecoregions_file': eco_poly if os.path.isfile(eco_poly) else None,
		'area': area, 'area_only': only,
		'eco_study': eco_study, 'eco_study_rel': paths.rel(sa_dir) if eco_study else None,
	}


def check_eco_study(eco, eco_study):
	"""Notes on ecoregion/study_area (the study area at the pixel of the ecoregion grid): the statistics of the HFTs
	in the similarity come from it when it has all the variables of the ecoregion grid, on its CRS and pixel."""
	out = []
	keys = {variable_key(r['file']) for r in eco}
	have = {variable_key(r['file']) for r in eco_study if not r.get('error')}
	good = [r for r in eco if not r.get('error') and r.get('res')]
	first = good[0] if good else None
	for r in eco_study:
		if r.get('error'):
			out.append({'level': 'error', 'text': f"{STUDY_AREA_FOLDER}/{r['file']}: cannot be read ({r['error']})"})
		elif first and not math.isclose(r['res'][0], first['res'][0], rel_tol=1e-6):
			out.append({'level': 'warning', 'text': f"{STUDY_AREA_FOLDER}/{r['file']} has another pixel size "
			                                        f"({_fmt(r['res'][0])} m) than the ecoregion grid ({_fmt(first['res'][0])} m)."})
	if keys - have:
		out.append({'level': 'warning', 'text': f'{STUDY_AREA_FOLDER} lacks variables of the ecoregion grid: '
		                                        + ', '.join(sorted(keys - have)) + ' (the similarity takes the HFT '
		                                        'statistics from the ecoregion grid itself then).'})
	if not any(o['level'] in ('error', 'warning') for o in out):
		res = _fmt(eco_study[0]['res'][0]) if eco_study and eco_study[0].get('res') else '?'
		out.insert(0, {'level': 'ok', 'text': f'The study area at the pixel of the ecoregion grid ({len(eco_study)} '
		                                      f'variables, {res} m): the similarity takes the statistics of the HFTs '
		                                      'from it.'})
	return out


def create_set(name):
	name = (name or '').strip().strip('/\\')
	if not NAME_RE.fullmatch(name) or name in (paths.ECO_FOLDER, paths.EXCLUDED):
		raise ValueError('the name can only have letters (without accents), digits, - and _')
	d = set_dir(name)
	if os.path.exists(d):
		raise ValueError(f'inVars/{name} already exists')
	os.makedirs(os.path.join(d, paths.ECO_FOLDER))
	return name


def _scale_dir(folder, scale):
	d = set_dir(folder)
	return os.path.join(d, paths.ECO_FOLDER) if scale == 'ecoregion' else d


def add_files(job, folder, scale, files):
	"""Put FILES (GeoTIFFs) in the set FOLDER at SCALE: a hard link when they are on the same disk (no space used),
	else a copy."""
	dest = _scale_dir(folder, scale)
	os.makedirs(dest, exist_ok=True)
	added = []
	for f in files:
		job.check()
		if not f.lower().endswith(RASTER_EXT):
			job.say(f'{os.path.basename(f)}: not a GeoTIFF (.tif), skipped')
			continue
		target = os.path.join(dest, os.path.basename(f))
		if os.path.exists(target):
			if os.path.samefile(f, target):
				job.say(f'{os.path.basename(f)}: already in the set')
				continue
			raise ValueError(f'there is already a file {os.path.basename(f)} in {paths.rel(dest)}: rename one of them')
		try:
			os.link(f, target)
			job.say(f'{os.path.basename(f)}: linked in {paths.rel(dest)}')
		except OSError:
			job.say(f'{os.path.basename(f)}: copying to {paths.rel(dest)} ({os.path.getsize(f) / 2 ** 20:,.0f} MB)...')
			shutil.copy2(f, target + '.part')
			os.replace(target + '.part', target)
		added.append(os.path.basename(f))
	return added


def set_included(folder, scale, file, included):
	"""Move a GeoTIFF of a set to its subfolder _excluded (INCLUDED false) or back."""
	base = _scale_dir(folder, scale)
	file = os.path.basename(file)
	src, dst = (os.path.join(base, paths.EXCLUDED, file), os.path.join(base, file)) if included else \
		(os.path.join(base, file), os.path.join(base, paths.EXCLUDED, file))
	if not os.path.isfile(src):
		raise ValueError(f'{file} not found')
	if os.path.exists(dst):
		raise ValueError(f'{paths.rel(dst)} already exists')
	os.makedirs(os.path.dirname(dst), exist_ok=True)
	os.replace(src, dst)
	return True


def warp_preview(path, max_size=900, resample='near', overviews=True, bbox=None):
	"""A raster in Web Mercator, at most MAX_SIZE pixels wide, for the map: (its values as float64, the mask of the
	pixels with data, [west, south, east, north] in WGS 84). RESAMPLE: the method of gdal.Warp (near, average, max...);
	OVERVIEWS False reads the full resolution (for max: the overviews of a map of classes may have lost small patches).
	BBOX ([west, south, east, north] in WGS 84): only that rectangle, with pixels not smaller than those of the raster."""
	import numpy as np
	ds = gdal.Open(path)
	nodata = ds.GetRasterBand(1).GetNoDataValue()
	w, h = ds.RasterXSize, ds.RasterYSize
	merc = osr.SpatialReference()
	merc.ImportFromEPSG(3857)
	extra = {}
	if bbox and ds.GetSpatialRef() is not None:
		xmin, ymin, xmax, ymax = osr.CoordinateTransformation(_wgs84(), merc).TransformBounds(*map(float, bbox), 21)
		# the size of a pixel of the raster in Web Mercator, at the centre of BBOX
		src = ds.GetSpatialRef().Clone()
		src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
		x0, y0, _ = osr.CoordinateTransformation(merc, src).TransformPoint(cx, cy)
		x1, _, _ = osr.CoordinateTransformation(src, merc).TransformPoint(x0 + abs(ds.GetGeoTransform()[1]), y0)
		width = max(1, min(max_size, int(math.ceil((xmax - xmin) / (abs(x1 - cx) or 1)))))
		extra['outputBounds'] = (xmin, ymin, xmax, ymax)
	else:
		width = max(1, int(round(w / max(max(w, h) / max_size, 1))))
	tmp = f'/vsimem/preview_{uuid.uuid4().hex}.tif'
	opts = gdal.WarpOptions(dstSRS='EPSG:3857', format='GTiff', outputType=gdal.GDT_Float32, width=width, height=0,
	                        resampleAlg=resample, srcNodata=nodata, dstNodata=-3.4e38, srcBands=[1], dstBands=[1],
	                        overviewLevel='AUTO' if overviews else 'NONE', **extra)
	try:
		out = gdal.Warp(tmp, ds, options=opts)
		a = out.GetRasterBand(1).ReadAsArray().astype('float64')
		gt = out.GetGeoTransform()
		ow, oh = out.RasterXSize, out.RasterYSize
		out = None
	finally:
		gdal.Unlink(tmp)
	valid = np.isfinite(a) & (a > -3e38)
	if nodata is not None and not (isinstance(nodata, float) and math.isnan(nodata)):
		valid &= a != nodata
	xmin, ymax = gt[0], gt[3]
	xmax, ymin = gt[0] + gt[1] * ow, gt[3] + gt[5] * oh
	return a, valid, _bbox_wgs84(merc, xmin, ymin, xmax, ymax)


def png_uri(rgba):
	"""PNG (data URI) of an RGBA image: an array of bytes (4, rows, columns)."""
	_, oh, ow = rgba.shape
	mem = gdal.GetDriverByName('MEM').Create('', ow, oh, 4, gdal.GDT_Byte)
	for i in range(4):
		mem.GetRasterBand(i + 1).WriteArray(rgba[i])
	png = f'/vsimem/preview_{uuid.uuid4().hex}.png'
	gdal.GetDriverByName('PNG').CreateCopy(png, mem)
	try:
		f = gdal.VSIFOpenL(png, 'rb')
		gdal.VSIFSeekL(f, 0, 2)
		size = gdal.VSIFTellL(f)
		gdal.VSIFSeekL(f, 0, 0)
		data = gdal.VSIFReadL(1, size, f)
		gdal.VSIFCloseL(f)
	finally:
		gdal.Unlink(png)
	return 'data:image/png;base64,' + base64.b64encode(data).decode()


def preview(path, max_size=900):
	"""Quick look of a raster in Web Mercator for the map: PNG (data URI), its corners and the range of values."""
	import numpy as np
	a, valid, bbox = warp_preview(path, max_size)
	if not valid.any():
		raise ValueError('the raster has no data')
	lo, hi = np.percentile(a[valid], [2, 98])
	if hi <= lo:
		hi = lo + 1
	idx = np.clip(((a - lo) / (hi - lo) * 255), 0, 255).astype('uint8')
	lut = _ramp(VIRIDIS)
	rgba = np.zeros((4,) + a.shape, dtype='uint8')
	for i in range(3):
		rgba[i] = lut[idx, i]
	rgba[3] = np.where(valid, 230, 0)
	return {'image': png_uri(rgba), 'bbox': bbox, 'min': float(lo), 'max': float(hi), 'ramp': ['#' + c for c in VIRIDIS]}


def _ramp(hexes):
	import numpy as np
	cols = np.array([[int(h[i:i + 2], 16) for i in (0, 2, 4)] for h in hexes], dtype='float64')
	x = np.linspace(0, 1, len(cols))
	t = np.linspace(0, 1, 256)
	return np.stack([np.interp(t, x, cols[:, i]) for i in range(3)], axis=1).astype('uint8')


# Study area ----------------------------------------------------------------------------------------------------------

def _int_type(t):
	return t in (ogr.OFTInteger, ogr.OFTInteger64)


def vector_info(path, layer=None):
	"""Layers of a vector file and, of LAYER (or the only one), its fields, number of features, CRS and extent."""
	if not os.path.isfile(path):
		raise ValueError(f'{path} not found')
	ds = ogr.Open(path)
	if ds is None:
		raise ValueError(f'{path} is not a vector file that GDAL can read')
	layers = [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]
	if layer and layer not in layers:
		layer = None
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	defn = lyr.GetLayerDefn()
	fields = [{'name': defn.GetFieldDefn(i).GetName(), 'type': defn.GetFieldDefn(i).GetTypeName(),
	           'integer': _int_type(defn.GetFieldDefn(i).GetType())} for i in range(defn.GetFieldCount())]
	srs = lyr.GetSpatialRef()
	xmin, xmax, ymin, ymax = lyr.GetExtent()
	area = None
	if srs is not None and srs.IsProjected():
		area = 0.0
		for f in lyr:
			g = f.GetGeometryRef()
			if g is not None:
				area += g.GetArea()
		lyr.ResetReading()
		area = area * (srs.GetLinearUnits() ** 2) / 1e6
	return {'path': path, 'file': os.path.basename(path), 'layers': layers, 'layer': lyr.GetName() if len(layers) > 1 else None,
	        'fields': fields, 'count': lyr.GetFeatureCount(), 'crs': _crs_label(srs),
	        'projected': bool(srs and srs.IsProjected()), 'bbox': _bbox_wgs84(srs, xmin, ymin, xmax, ymax),
	        'area_km2': round(area, 3) if area is not None else None,
	        'in_pa': paths.inside(path, paths.PA_DIR)}


def vector_ids(path, layer, col_id):
	"""IDs of the study area as setup.py takes them: the integer column COL_ID or, with 'cat' when the file has no
	such column, 1, 2, ... in the order of the features (the categories that GRASS gives them)."""
	ds = ogr.Open(path)
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	defn = lyr.GetLayerDefn()
	i = defn.GetFieldIndex(col_id)
	if i < 0:
		if col_id == 'cat':
			return list(range(1, lyr.GetFeatureCount() + 1))
		raise ValueError(f'{os.path.basename(path)} has no column {col_id}')
	if not _int_type(defn.GetFieldDefn(i).GetType()):
		raise ValueError(f'the column {col_id} is not of integers')
	return sorted({f.GetFieldAsInteger64(i) for f in lyr if f.IsFieldSet(i)})


NAME_RE_FIELD = re.compile(r'(^|_)(name|nombre|nom|label|site|desig|title|titulo)', re.I)


def guess_label_col(fields):
	"""A text column with the names of the areas (name, NOMBRE, SITE_NAME...), or None."""
	text = [f['name'] for f in fields if f['type'] in ('String', 'WideString')]
	for f in text:
		if NAME_RE_FIELD.search(f):
			return f
	return None


def _feature_ids(lyr, col_id):
	"""(index of COL_ID or -1 for the categories of GRASS, 1, 2, ... in the order of the features)."""
	defn = lyr.GetLayerDefn()
	i = defn.GetFieldIndex(col_id)
	if i < 0 and col_id != 'cat':
		raise ValueError(f'the layer has no column {col_id}')
	if i >= 0 and not _int_type(defn.GetFieldDefn(i).GetType()):
		raise ValueError(f'the column {col_id} is not of integers')
	return i


def vector_areas(path, layer, col_id, label_col=None):
	"""The study areas of a vector file (one per ID, as setup.py and palist.csv take them): [{id, label, ha}],
	with the name of LABEL_COL (if any) and the area in hectares (projected CRS)."""
	ds = ogr.Open(path)
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	i = _feature_ids(lyr, col_id)
	j = lyr.GetLayerDefn().GetFieldIndex(label_col) if label_col else -1
	srs = lyr.GetSpatialRef()
	factor = (srs.GetLinearUnits() ** 2) / 1e4 if srs is not None and srs.IsProjected() else None
	areas = {}
	for n, f in enumerate(lyr, 1):
		if i >= 0 and not f.IsFieldSet(i):
			continue
		fid = f.GetFieldAsInteger64(i) if i >= 0 else n
		a = areas.setdefault(fid, {'id': fid, 'label': None, 'ha': 0.0 if factor else None})
		if j >= 0 and a['label'] is None and f.IsFieldSet(j):
			a['label'] = f.GetFieldAsString(j).strip() or None
		g = f.GetGeometryRef()
		if factor and g is not None:
			a['ha'] += g.GetArea() * factor
	for a in areas.values():
		if a['ha'] is not None:
			a['ha'] = round(a['ha'], 1)
	return [areas[k] for k in sorted(areas)]


def vector_geojson(path, layer=None, max_points=60000, col_id=None, label_col=None):
	"""Polygons of a vector file in WGS 84 (GeoJSON), simplified to about MAX_POINTS vertices for the map. With
	COL_ID, each feature has the ID of its study area (and the name of LABEL_COL)."""
	ds = ogr.Open(path)
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	try:
		id_index = _feature_ids(lyr, col_id) if col_id else None
	except ValueError:
		id_index = None
	label_index = lyr.GetLayerDefn().GetFieldIndex(label_col) if label_col else -1
	srs = lyr.GetSpatialRef()
	ct = None
	if srs is not None:
		src = srs.Clone()
		src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		ct = osr.CoordinateTransformation(src, _wgs84())
	xmin, xmax, ymin, ymax = lyr.GetExtent()
	npts = 0
	geoms = []
	for n, f in enumerate(lyr, 1):
		g = f.GetGeometryRef()
		if g is None:
			continue
		g = g.Clone()
		npts += _count_points(g)
		props = {'fid': f.GetFID()}
		if id_index is not None:
			props['id'] = (f.GetFieldAsInteger64(id_index) if f.IsFieldSet(id_index) else None) if id_index >= 0 else n
		if label_index >= 0 and f.IsFieldSet(label_index):
			props['label'] = f.GetFieldAsString(label_index)
		geoms.append((props, g))
	tol = 0
	if npts > max_points:
		tol = max(xmax - xmin, ymax - ymin) / 2000
	feats = []
	for props, g in geoms:
		if tol:
			g = g.SimplifyPreserveTopology(tol)
		if ct is not None:
			g.Transform(ct)
		feats.append({'type': 'Feature', 'properties': props, 'geometry': _json(g)})
	return {'type': 'FeatureCollection', 'features': feats}


def _json(g):
	import json
	return json.loads(g.ExportToJson(['COORDINATE_PRECISION=6']))


def _count_points(g):
	n = g.GetPointCount()
	for i in range(g.GetGeometryCount()):
		n += _count_points(g.GetGeometryRef(i))
	return n


def vector_parts(path):
	"""Files of a vector file (a shapefile has several)."""
	stem, ext = os.path.splitext(path)
	if ext.lower() != '.shp':
		return [path]
	return [stem + e for e in SHP_PARTS if os.path.exists(stem + e)] + \
	       [p for p in glob.glob(stem + '.*') if p.lower().endswith('.shp.xml') and p not in (stem + '.shp.xml',)]


def copy_to_pa(path):
	"""Copy the study area file to segm_optim/pa (unless it is there) and return its path there."""
	if paths.inside(path, paths.PA_DIR):
		return path
	os.makedirs(paths.PA_DIR, exist_ok=True)
	target = os.path.join(paths.PA_DIR, os.path.basename(path))
	if os.path.exists(target):
		if os.path.samefile(path, target):
			return target
		raise ValueError(f'there is already a file {os.path.basename(path)} in segm_optim/pa: rename one of them')
	for p in vector_parts(path):
		shutil.copy2(p, os.path.join(paths.PA_DIR, os.path.basename(p)))
	return target


def grass_name(path):
	"""Name of the GRASS project and map of a study area file (as setup.py)."""
	stem = os.path.splitext(os.path.basename(path))[0]
	name = re.sub(r'\W', '_', stem, flags=re.ASCII)
	return name if name[:1].isalpha() else 'v_' + name
