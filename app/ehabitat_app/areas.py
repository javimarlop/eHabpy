"""Areas from MonitorEO: vector files for the study areas or a comparison area, made with the package monitoreo in
Earth Engine (with the account of the downloads):

- the burned area of a fire: the perimeter of the pixels burned between a period before and a period after it, from the
  dNBR, the RBR or the RdNBR (services.burn_severity with perimeter, monitoreo >= 0.3.47), and the unburned ring around
  it (a band of a given width outside the perimeter, starting at a given distance from it), e.g. the area to compare
  the HFTs of the burned area with;
- protected areas of WDPA (also only its Ramsar sites) and administrative units of FAO GAUL 2025 (levels 1 and 2),
  found by their name.

Each area is a folder of segm_optim/areas with a GeoPackage of polygons in the WGS 84 / UTM zone of their centre,
whose fields cat (1, 2, ...: the ID of each feature, a study area for setup.py) and name make it a file of study areas
as it is, its record (area.json: how it was made) and a README with its references; the burned area keeps there the
files of monitoreo too (the perimeter in EPSG:4326, the table of the burn index and its README)."""

import datetime as dt
import inspect
import json
import math
import os
import re

from osgeo import ogr, osr

from . import gee, jobs, model, paths

ogr.UseExceptions()
osr.UseExceptions()

AREAS_DIR = os.path.join(paths.SEGM, 'areas')
RECORD = 'area.json'
README = 'README.txt'
MONITOREO_DIR = 'monitoreo'  # subfolder of a burned area: the files of monitoreo (perimeter, table, README)
SLUG_RE = re.compile(r'[^A-Za-z0-9_-]+')
SIMPLIFY_M = 10  # metres: the polygons of WDPA and GAUL fetched from Earth Engine (as monitoreo.areas.WDPA_MAX_ERROR)
SEARCH_BUFFER = 1000  # metres around the study area where a fire is searched, by default
RING_WIDTH = 500  # metres: the default width of the unburned ring
KIND_TEXT = {'burned': 'burned area of a fire', 'protected': 'protected areas of WDPA',
             'admin': 'administrative units of FAO GAUL 2025'}


# Options and areas made before -----------------------------------------------------------------------------------

def options():
	"""What the dialog of the areas offers with the installed monitoreo: {'burn': {indices, sensors (of each index),
	burned_min (the value of the dNBR and the RdNBR from which a pixel is burned; the RBR has none), min_ha, buffer,
	ring} or None (a version of monitoreo without the perimeter of the burned area), 'protected', 'admin' (whether it
	can find them), 'monitoreo' (its version), 'made' (the areas made before), 'study' (the file of the study area)}."""
	_, services = gee._monitoreo()
	burn = None
	fn = getattr(services, 'burn_severity', None)
	if fn is not None and 'perimeter' in inspect.signature(fn).parameters:
		from monitoreo import fire
		burn = {'indices': list(getattr(services, 'BURN_INDICES', ('dNBR', 'RBR', 'RdNBR'))),
		        'sensors': {k: list(v) for k, v in getattr(services, 'BURN_SENSORS', {}).items()},
		        'burned_min': {k: c['burned'] for k, c in getattr(fire, 'CLASSES', {}).items()},
		        'min_ha': getattr(fire, 'PERIMETER_MIN_HA', 1.0), 'buffer': SEARCH_BUFFER, 'ring': RING_WIDTH}
	sa = model.study_area()
	return {'burn': burn, 'protected': hasattr(services, 'find_protected_areas'),
	        'admin': hasattr(services, 'find_admin_units'), 'monitoreo': gee._label(gee.installed()),
	        'made': made(), 'study': paths.rel(sa['file']) if sa else None}


def _record(folder):
	try:
		with open(os.path.join(folder, RECORD), encoding='utf-8') as f:
			rec = json.load(f)
	except (OSError, ValueError):
		return None
	files = []
	for x in rec.get('files') or []:
		p = os.path.join(folder, x['file'])
		if os.path.isfile(p):
			files.append(dict(x, path=p, rel=paths.rel(p), comparison=comparison_path(p)))
	return dict(rec, folder=os.path.basename(folder), files=files) if files else None


def made():
	"""The areas made before (folders of segm_optim/areas with their record), the last first."""
	if not os.path.isdir(AREAS_DIR):
		return []
	out = [r for r in (_record(os.path.join(AREAS_DIR, d)) for d in os.listdir(AREAS_DIR)
	                   if os.path.isdir(os.path.join(AREAS_DIR, d))) if r]
	return sorted(out, key=lambda r: r.get('made') or '', reverse=True)


def comparison_path(path):
	"""PATH as the comparison area of a recipe takes it: relative to segm_optim when it is inside (as get_vars_gee.py
	reads it), else absolute."""
	path = os.path.abspath(path)
	return os.path.relpath(path, paths.SEGM).replace(os.sep, '/') if paths.inside(path, paths.SEGM) else path


# Finding protected areas and administrative units ----------------------------------------------------------------

def _connected():
	_, services = gee._monitoreo()
	gee._earth_engine(services)
	return services


def find_protected(text, ramsar=False):
	"""Protected areas of WDPA whose name (or English name) contains TEXT (with RAMSAR, only the Ramsar sites):
	[{id (WDPA ID), name, desig, iucn, iso3, km2, year}]."""
	text = str(text or '').strip()
	if len(text) < 2:
		raise ValueError('write at least two letters of the name of the protected area')
	services = _connected()
	rows = services.find_protected_areas(text, 50, 'ramsar' if ramsar else None)
	return [{'id': r.get('SITE_ID') or r.get('WDPAID'), 'name': r.get('NAME') or r.get('NAME_ENG'),
	         'desig': r.get('DESIG_ENG') or r.get('DESIG'), 'iucn': r.get('IUCN_CAT'), 'iso3': r.get('ISO3'),
	         'km2': r.get('REP_AREA'), 'year': r.get('STATUS_YR')} for r in rows]


def find_admin(text, level=2, country=None):
	"""Administrative units of FAO GAUL 2025 of LEVEL (1 or 2) whose name contains TEXT, in COUNTRY (ISO3 code) if
	given: [{id (GAUL1_CODE or GAUL2_CODE), name (that of the unit), full (country, region...), iso3, level}]."""
	text = str(text or '').strip()
	if len(text) < 2:
		raise ValueError('write at least two letters of the name of the administrative unit')
	level = int(level)
	if level not in (1, 2):
		raise ValueError('the level is 1 (regions, states, provinces...) or 2 (provinces, districts, municipalities...)')
	country = str(country or '').strip().upper() or None
	if country and not re.fullmatch(r'[A-Z]{3}', country):
		raise ValueError('write the country as its ISO3 code (e.g. ESP), or leave it empty')
	services = _connected()
	rows = services.find_admin_units(text, level, country, 50)
	return [{'id': r.get(f'GAUL{level}_CODE'), 'name': r.get(f'GAUL{level}_NAME'), 'full': r.get('DISP_EN'),
	         'iso3': r.get('ISO3_CODE'), 'level': level} for r in rows]


# Geometries ------------------------------------------------------------------------------------------------------

def _wgs84():
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(4326)
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return srs


def _utm(lon, lat):
	"""WGS 84 / UTM zone of a point (as monitoreo.areas.utm_of)."""
	zone = min(int(math.floor((lon + 180) / 6)) + 1, 60)
	srs = osr.SpatialReference()
	srs.ImportFromEPSG((32600 if lat >= 0 else 32700) + zone)
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return srs


def _polygons(g):
	"""MultiPolygon of the polygons of G (valid; a simplified geometry can be a collection with lines or points)."""
	if g is None:
		return None
	if not g.IsValid():
		g = g.MakeValid()
	out = ogr.Geometry(ogr.wkbMultiPolygon)
	flat = ogr.GT_Flatten(g.GetGeometryType())
	if flat in (ogr.wkbPolygon, ogr.wkbMultiPolygon):
		g = ogr.ForceToMultiPolygon(g)
		for i in range(g.GetGeometryCount()):
			out.AddGeometry(g.GetGeometryRef(i))
	elif flat == ogr.wkbGeometryCollection:
		for i in range(g.GetGeometryCount()):
			part = _polygons(g.GetGeometryRef(i))
			for j in range(part.GetGeometryCount() if part else 0):
				out.AddGeometry(part.GetGeometryRef(j))
	return out if out.GetGeometryCount() else None


def _from_geojson(gj):
	"""MultiPolygon of the polygons of the GeoJSON geometry GJ (None if it has none)."""
	if not gj or (gj.get('coordinates') is None and not gj.get('geometries')):
		return None
	return _polygons(ogr.CreateGeometryFromJson(json.dumps(gj)))


def _union(geoms):
	geoms = [g for g in geoms if g is not None]
	if not geoms:
		return None
	acc = ogr.Geometry(ogr.wkbMultiPolygon)
	for g in geoms:
		for i in range(g.GetGeometryCount()):
			acc.AddGeometry(g.GetGeometryRef(i))
	return _polygons(acc.UnionCascaded())


def _filled(g):
	"""G without its holes: the union of the outer rings of its polygons (by a buffer of 0: a polygon that was an island
	in a hole of another one is inside it now, which makes a union of GEOS crawl)."""
	out = ogr.Geometry(ogr.wkbMultiPolygon)
	for i in range(g.GetGeometryCount()):
		p = ogr.Geometry(ogr.wkbPolygon)
		p.AddGeometry(g.GetGeometryRef(i).GetGeometryRef(0).Clone())
		out.AddGeometry(p)
	return _polygons(out.Buffer(0))


RING_FIRST_STEP = 10.0  # metres: the first buffer of the ring (then 4 times bigger each time)
RING_TOLERANCE = 0.5  # metres: the simplification of the polygons after each buffer (the vertices of their arcs)


def grow(g, distance):
	"""G (in metres) grown by DISTANCE metres, in buffers of growing size (10, 40, 160 m...): the buffer of a buffer is
	the buffer of the sum of both distances, and the first small buffers join the patches that are close, so that the
	big ones are fast (one buffer of 500 m of the 2000 patches of a burned area took minutes, so 2 seconds)."""
	done, step = 0.0, RING_FIRST_STEP
	while done < distance:
		d = min(step, distance - done)
		g = g.Buffer(d, 8).SimplifyPreserveTopology(RING_TOLERANCE)
		done += d
		step *= 4
	return g


def ring(g, width, gap=0.0):
	"""The ring of WIDTH metres around the polygons G (in metres) that starts GAP metres from them: what is between GAP
	and GAP + WIDTH metres outside G with its holes filled (its unburned islands are inside, not around it)."""
	filled = _filled(g)
	inner = grow(filled, float(gap)) if gap > 0 else filled
	return _polygons(grow(filled, float(gap) + float(width)).Difference(inner))


def _centre(geoms, srs):
	"""(lon, lat) of the centre of the extent of GEOMS (in SRS)."""
	env = [g.GetEnvelope() for g in geoms]
	x = (min(e[0] for e in env) + max(e[1] for e in env)) / 2
	y = (min(e[2] for e in env) + max(e[3] for e in env)) / 2
	pt = ogr.Geometry(ogr.wkbPoint)
	pt.AddPoint_2D(x, y)
	if not srs.IsSame(_wgs84()):
		pt.Transform(osr.CoordinateTransformation(srs, _wgs84()))
	return pt.GetX(), pt.GetY()


def _write(path, layer, srs, fields, rows):
	"""Write ROWS ([(geometry, {field: value})]) to the GeoPackage PATH (layer LAYER, CRS SRS) with the FIELDS
	([(name, ogr type)]); returns PATH."""
	drv = ogr.GetDriverByName('GPKG')
	if os.path.exists(path):
		drv.DeleteDataSource(path)
	ds = drv.CreateDataSource(path)
	lyr = ds.CreateLayer(layer, srs, ogr.wkbMultiPolygon)
	for name, kind in fields:
		lyr.CreateField(ogr.FieldDefn(name, kind))
	for g, props in rows:
		f = ogr.Feature(lyr.GetLayerDefn())
		for name, _ in fields:
			if props.get(name) is not None:
				f.SetField(name, props[name])
		f.SetGeometry(g)
		lyr.CreateFeature(f)
	ds = None
	return path


# Making an area --------------------------------------------------------------------------------------------------

def _slug(text):
	return SLUG_RE.sub('_', str(text or '').strip()).strip('_')[:60] or 'area'


def _new_folder(name):
	"""A new folder of segm_optim/areas for the area NAME: its slug, or slug_2, slug_3... if it exists (a file of the
	study areas is copied to segm_optim/pa by its name, so that two areas never have the same file name)."""
	base = _slug(name)
	slug, n = base, 2
	taken = {os.path.splitext(x)[0].lower() for x in os.listdir(paths.PA_DIR)} if os.path.isdir(paths.PA_DIR) else set()
	while os.path.exists(os.path.join(AREAS_DIR, slug)) or slug.lower() in taken:
		slug, n = f'{base}_{n}', n + 1
	folder = os.path.join(AREAS_DIR, slug)
	os.makedirs(folder)
	return folder, slug


def _date(text, what):
	try:
		return dt.date.fromisoformat(str(text)).isoformat()
	except ValueError as e:
		raise ValueError(f'write the dates of {what} (YYYY-MM-DD)') from e


def _number(v, what, low=0.0, default=None):
	if v in (None, ''):
		return default
	try:
		x = float(v)
	except (TypeError, ValueError) as e:
		raise ValueError(f'{what} is a number') from e
	if not math.isfinite(x) or x < low:
		raise ValueError(f'{what} cannot be less than {low:g}')
	return x


def check(kind, p):
	"""The parameters P of an area of KIND ('burned', 'protected' or 'admin'), checked (ValueError if not)."""
	if kind not in KIND_TEXT:
		raise ValueError(f'unknown kind of area {kind}')
	p = dict(p or {})
	p['name'] = str(p.get('name') or '').strip()
	if not p['name']:
		raise ValueError('write a name for the area')
	if kind == 'burned':
		before = [_date(x, 'the period before the fire') for x in (p.get('before') or [None, None])]
		after = [_date(x, 'the period after the fire') for x in (p.get('after') or [None, None])]
		if before[0] > before[1] or after[0] > after[1]:
			raise ValueError('each period goes from its first day to its last day')
		if before[1] >= after[0]:
			raise ValueError('the period before the fire has to end before the period after it starts')
		p.update(before=before, after=after)
		p['index'] = {'dnbr': 'dNBR', 'rbr': 'RBR', 'rdnbr': 'RdNBR'}.get(str(p.get('index') or 'dNBR').lower())
		if not p['index']:
			raise ValueError('the burn index is dNBR, RBR or RdNBR')
		p['burned_min'] = _number(p.get('burned_min'), 'the value from which a pixel is burned', -10)
		if p['index'] == 'RBR' and p['burned_min'] is None:
			raise ValueError('the RBR has no classes: write the RBR from which a pixel is burned (e.g. 0.1)')
		p['min_ha'] = _number(p.get('min_ha'), 'the smallest patch', 0)
		p['ring'] = _number(p.get('ring'), 'the width of the ring', 0, 0)
		p['gap'] = _number(p.get('gap'), 'the distance of the ring from the burned area', 0, 0)
		p['resolution'] = _number(p.get('resolution'), 'the pixel size', 1)
		p['buffer'] = _number(p.get('buffer'), 'the buffer of the search area', 0, SEARCH_BUFFER)
		p['search'] = p.get('search') or 'study'
		if p['search'] == 'study':
			if not model.study_area():
				raise ValueError('choose the study area first, or search the fire in a vector file')
		elif p['search'] == 'file':
			f = gee.area_path(p.get('file'))
			if not f or not os.path.isfile(f):
				raise ValueError('choose the vector file of the area where the fire is searched')
			p['file'] = f
		else:
			raise ValueError('search the fire in the study area or in a vector file')
	else:
		items = [x for x in p.get('items') or [] if x.get('id') is not None]
		if not items:
			raise ValueError('choose at least one ' + ('protected area' if kind == 'protected' else
			                                           'administrative unit'))
		p['items'] = items
		if kind == 'admin':
			p['level'] = int(p.get('level') or items[0].get('level') or 2)
	return p


def make(kind, params):
	"""Make an area of KIND with PARAMS in a job (Earth Engine may take a while): its result is the record of the area
	(see made)."""
	if jobs.running('gee-area'):
		raise ValueError('an area is already being made: wait until it finishes')
	if jobs.running('gee-update'):
		raise ValueError('MonitorEO is being updated: wait until it finishes (a few seconds)')
	p = check(kind, params)
	title = {'burned': 'Burned area from MonitorEO', 'protected': 'Protected areas from MonitorEO',
	         'admin': 'Administrative units from MonitorEO'}[kind]

	def work(job):
		services = _connected()
		os.makedirs(AREAS_DIR, exist_ok=True)
		folder, slug = _new_folder(p['name'])
		try:
			fn = {'burned': _burned, 'protected': _regions, 'admin': _regions}[kind]
			files, cited, extra = fn(job, kind, p, services, folder, slug)
		except BaseException:
			import shutil
			shutil.rmtree(folder, ignore_errors=True)
			raise
		rec = {'kind': kind, 'name': p['name'], 'made': dt.datetime.now().isoformat(timespec='seconds'),
		       'monitoreo': gee._label(gee.installed()), 'params': p, 'files': files, **extra}
		with open(os.path.join(folder, RECORD), 'w', encoding='utf-8') as f:
			json.dump(rec, f, indent=1, ensure_ascii=False, default=str)
		_readme(folder, rec, cited)
		job.say(f'\nDone: {paths.rel(folder)}')
		for x in files:
			job.say(f"  {x['file']}: {x['label']} ({x['features']} {'area' if x['features'] == 1 else 'areas'}, "
			        f"{x['ha']:,.1f} ha)")
		return _record(folder)

	return jobs.start_function('gee-area', title, work).id


def _search_region(p):
	"""(region of monitoreo, its description) where the fire is searched: the rectangle of the study area grown by
	BUFFER metres, or the polygons of a vector file (each one a region with its own perimeter)."""
	if p['search'] == 'file':
		return {'file': p['file']}, f"the polygons of {paths.rel(p['file'])}"
	sa = model.study_area()
	ds = ogr.Open(sa['file'])
	lyr = ds.GetLayerByName(sa['layer']) if sa.get('layer') else ds.GetLayer(0)
	srs = lyr.GetSpatialRef()
	if srs is None:
		raise ValueError(f"the study area {os.path.basename(sa['file'])} has no CRS")
	srs = srs.Clone()
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	xmin, xmax, ymin, ymax = lyr.GetExtent()
	ds = None
	if not srs.IsProjected():  # a rectangle in metres: in the UTM zone of its centre
		utm = _utm((xmin + xmax) / 2, (ymin + ymax) / 2)
		box = ogr.CreateGeometryFromWkt(f'POLYGON(({xmin} {ymin},{xmax} {ymin},{xmax} {ymax},{xmin} {ymax},{xmin} {ymin}))')
		box.Transform(osr.CoordinateTransformation(srs, utm))
		xmin, xmax, ymin, ymax = box.GetEnvelope()
		srs = utm
	b = p['buffer'] / (srs.GetLinearUnits() or 1)
	xmin, ymin, xmax, ymax = xmin - b, ymin - b, xmax + b, ymax + b
	box = ogr.Geometry(ogr.wkbLinearRing)
	n = 16  # points on each side, so that the rectangle keeps its shape in EPSG:4326
	for i in range(n):
		box.AddPoint_2D(xmin + (xmax - xmin) * i / n, ymin)
	for i in range(n):
		box.AddPoint_2D(xmax, ymin + (ymax - ymin) * i / n)
	for i in range(n):
		box.AddPoint_2D(xmax - (xmax - xmin) * i / n, ymax)
	for i in range(n + 1):
		box.AddPoint_2D(xmin, ymax - (ymax - ymin) * i / n)
	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(box)
	poly.Transform(osr.CoordinateTransformation(srs, _wgs84()))
	gj = {'type': 'FeatureCollection', 'features': [
		{'type': 'Feature', 'properties': {'name': 'search area'}, 'geometry': json.loads(poly.ExportToJson())}]}
	return {'geojson': gj, 'name_field': 'name'}, \
		f"the rectangle of the study area {os.path.basename(sa['file'])} grown by {p['buffer']:g} m"


def _burned(job, kind, p, services, folder, slug):
	"""The burned area (and its ring) of the parameters P: burn_severity of monitoreo with perimeter, then its polygons
	in UTM with cat and name."""
	region, where = _search_region(p)
	raw = os.path.join(folder, MONITOREO_DIR)
	os.makedirs(raw)
	index, sensor = p['index'], p.get('sensor') or 'Sentinel-2'
	perimeter = os.path.join(raw, f'perimeter_{index}.gpkg')
	csv = os.path.join(raw, f'burn_{index}.csv')
	tif = os.path.join(raw, f'burn_{index}.tif') if p.get('geotiff') else None
	job.say(f"Burned area from the {index} of {sensor} (median of the images of each period): before the fire "
	        f"{p['before'][0]} to {p['before'][1]}, after it {p['after'][0]} to {p['after'][1]}, in {where}.")
	job.say('Earth Engine computes the burn index, its table and the perimeter of the burned pixels (this may take a '
	        'few minutes)...')
	kw = {}
	if tif:
		kw['progress'] = lambda done, total: job.log(f'\r  GeoTIFF of the burn index: {done}/{total} tiles')
	out = services.burn_severity(region, sensor, p['before'][0], p['before'][1], p['after'][0], p['after'][1],
	                             index=index, resolution=p.get('resolution'), csv=csv, path=tif, perimeter=perimeter,
	                             burned_min=p.get('burned_min'), min_ha=p.get('min_ha'), **kw)
	if tif:
		job.say('')
	job.check()
	# the perimeter of monitoreo (EPSG:4326, a feature for each region) in the UTM zone of its centre
	src = ogr.Open(perimeter)
	lyr = src.GetLayer(0)
	feats = []
	for f in lyr:
		g = _polygons(f.GetGeometryRef().Clone()) if f.GetGeometryRef() is not None else None
		props = {k: f.GetField(k) for k in ('region', 'ha', 'patches', 'index', 'burned_min', 'min_ha', 'sensor',
		                                    'before', 'after') if lyr.GetLayerDefn().GetFieldIndex(k) >= 0}
		if g is None or g.IsEmpty():
			job.say(f"WARNING: no burned pixels in {props.get('region') or 'the search area'}")
			continue
		feats.append((g, props))
	src = None
	if not feats:
		raise ValueError(f"no burned area: no pixel has a {index} of at least {out.get('burned_min')} in patches of "
		                 f"{out.get('min_ha')} ha or more. Check the dates (the fire between both periods), the index "
		                 'and the value from which a pixel is burned.')
	wgs84 = _wgs84()
	utm = _utm(*_centre([g for g, _ in feats], wgs84))
	to_utm = osr.CoordinateTransformation(wgs84, utm)
	for g, _ in feats:
		g.Transform(to_utm)
	several = len(feats) > 1
	burned_min, min_ha = out.get('burned_min'), out.get('min_ha')
	rows = []
	for i, (g, props) in enumerate(feats, 1):
		name = f"{p['name']} ({props.get('region')})" if several else p['name']
		rows.append((g, dict(props, cat=i, name=name, ha=round(props.get('ha') or g.GetArea() / 1e4, 2))))
	fields = [('cat', ogr.OFTInteger), ('name', ogr.OFTString), ('ha', ogr.OFTReal), ('patches', ogr.OFTInteger),
	          ('index', ogr.OFTString), ('burned_min', ogr.OFTReal), ('min_ha', ogr.OFTReal),
	          ('sensor', ogr.OFTString), ('before', ogr.OFTString), ('after', ogr.OFTString),
	          ('region', ogr.OFTString)]
	area_file = os.path.join(folder, f'{slug}.gpkg')
	_write(area_file, 'burned_area', utm, fields, rows)
	total = sum(r[1]['ha'] for r in rows)
	patches = sum(r[1].get('patches') or 0 for r in rows)
	label = f'Burned area ({index} ≥ {burned_min:g}, patches of {min_ha:g} ha or more)'
	job.say(f'Burned area: {total:,.1f} ha in {patches} patches ({index} ≥ {burned_min:g}).')
	files = [{'file': os.path.basename(area_file), 'role': 'burned', 'label': label, 'features': len(rows),
	          'ha': round(total, 2)}]
	if p['ring'] > 0:
		ring_rows = []
		for g, props in rows:
			r = ring(g, p['ring'], p['gap'])
			if r is None:
				continue
			ring_rows.append((r, {'cat': props['cat'], 'name': f"{props['name']}: unburned ring of {p['ring']:g} m"
			                      + (f" from {p['gap']:g} m" if p['gap'] else ''), 'ha': round(r.GetArea() / 1e4, 2),
			                      'width_m': p['ring'], 'gap_m': p['gap'], 'region': props.get('region')}))
		suffix = f"_ring{p['ring']:g}m" + (f"_gap{p['gap']:g}m" if p['gap'] else '')
		ring_file = os.path.join(folder, f'{slug}{suffix}.gpkg')
		_write(ring_file, 'unburned_ring', utm, [('cat', ogr.OFTInteger), ('name', ogr.OFTString), ('ha', ogr.OFTReal),
		                                         ('width_m', ogr.OFTReal), ('gap_m', ogr.OFTReal),
		                                         ('region', ogr.OFTString)], ring_rows)
		ring_ha = sum(r[1]['ha'] for r in ring_rows)
		job.say(f"Unburned ring of {p['ring']:g} m" + (f" starting {p['gap']:g} m from the burned area" if p['gap']
		                                                else '') + f': {ring_ha:,.1f} ha.')
		files.append({'file': os.path.basename(ring_file), 'role': 'ring',
		              'label': f"Unburned ring of {p['ring']:g} m around the burned area"
		                       + (f", from {p['gap']:g} m" if p['gap'] else ''), 'features': len(ring_rows),
		              'ha': round(ring_ha, 2)})
	cited = {'references': out.get('references') or [], 'credits': out.get('credits') or []}
	extra = {'crs': utm.GetAuthorityCode(None) and f"EPSG:{utm.GetAuthorityCode(None)}", 'search': where,
	         'burned_min': burned_min, 'min_ha': min_ha, 'scale_m': out.get('resolution'), 'title': out.get('title'),
	         'monitoreo_files': sorted(os.listdir(raw))}
	return files, cited, extra


def _regions(job, kind, p, services, folder, slug):
	"""Protected areas of WDPA or administrative units of GAUL (the items of P), each one a feature (or their union
	with P['join'])."""
	import ee
	from monitoreo import areas as mo_areas
	from monitoreo import references
	feats = []
	for item in p['items']:
		job.check()
		code = int(item['id'])
		if kind == 'protected':
			job.say(f"WDPA {code}: {item.get('name')}...")
			fc = mo_areas.region({'wdpa_id': code})  # simplified by 10 m in monitoreo
			props = {'name': item.get('name'), 'wdpa_id': code, 'desig': item.get('desig'), 'iucn_cat': item.get('iucn'),
			         'iso3': item.get('iso3')}
		else:
			job.say(f"GAUL {p['level']} {code}: {item.get('full') or item.get('name')}...")
			fc = mo_areas.region({f"admin{p['level']}": code}).map(lambda f: f.simplify(SIMPLIFY_M))
			props = {'name': item.get('name'), 'gaul_code': code, 'gaul_level': p['level'], 'full_name': item.get('full'),
			         'iso3': item.get('iso3')}
		try:
			info = fc.getInfo()
		except ee.EEException as e:
			raise ValueError(f'Earth Engine: {e}') from e
		g = _union([_from_geojson(f.get('geometry')) for f in info.get('features') or []])
		if g is None:
			job.say(f"WARNING: {item.get('name')} has no polygons in Earth Engine")
			continue
		feats.append((g, props))
	if not feats:
		raise ValueError('none of them has polygons in Earth Engine')
	wgs84 = _wgs84()
	utm = _utm(*_centre([g for g, _ in feats], wgs84))
	to_utm = osr.CoordinateTransformation(wgs84, utm)
	for g, _ in feats:
		g.Transform(to_utm)
	if p.get('join') and len(feats) > 1:
		g = _union([g for g, _ in feats])
		names = [x['name'] for _, x in feats]
		feats = [(g, {'name': p['name'], 'iso3': ','.join(sorted({x.get('iso3') or '' for _, x in feats} - {''})),
		              'parts': '; '.join(names)})]
	rows = [(g, dict(props, cat=i, ha=round(g.GetArea() / 1e4, 2))) for i, (g, props) in enumerate(feats, 1)]
	if kind == 'protected':
		fields = [('cat', ogr.OFTInteger), ('name', ogr.OFTString), ('ha', ogr.OFTReal), ('wdpa_id', ogr.OFTInteger64),
		          ('desig', ogr.OFTString), ('iucn_cat', ogr.OFTString), ('iso3', ogr.OFTString)]
		region = {'wdpa_id': [int(x['id']) for x in p['items']]}
	else:
		fields = [('cat', ogr.OFTInteger), ('name', ogr.OFTString), ('ha', ogr.OFTReal), ('gaul_code', ogr.OFTInteger64),
		          ('gaul_level', ogr.OFTInteger), ('full_name', ogr.OFTString), ('iso3', ogr.OFTString)]
		region = {f"admin{p['level']}": [int(x['id']) for x in p['items']]}
	if p.get('join') and len(rows) > 1:
		fields.append(('parts', ogr.OFTString))
	area_file = os.path.join(folder, f'{slug}.gpkg')
	_write(area_file, 'protected_areas' if kind == 'protected' else 'admin_units', utm, fields, rows)
	total = sum(r[1]['ha'] for r in rows)
	label = (('Ramsar sites' if p.get('ramsar') else 'Protected areas') + ' of WDPA') if kind == 'protected' else \
		f"Administrative units of FAO GAUL 2025 (level {p['level']})"
	if p.get('join') and len(p['items']) > 1:
		label += ', joined'
	files = [{'file': os.path.basename(area_file), 'role': kind, 'label': label, 'features': len(rows),
	          'ha': round(total, 2)}]
	cited = references.cite('find_protected_areas' if kind == 'protected' else 'find_admin_units', {'region': region})
	extra = {'crs': f"EPSG:{utm.GetAuthorityCode(None)}", 'simplify_m': SIMPLIFY_M}
	return files, cited, extra


def _readme(folder, rec, cited):
	"""README.txt of the area: its files, how it was made and what to cite with it."""
	p = rec['params']
	lines = [f"{rec['name']}: {KIND_TEXT[rec['kind']]}", '',
	         f"Made by the eHabitat+ app with MonitorEO (package monitoreo {rec['monitoreo']}) in Google Earth Engine "
	         f"on {rec['made'][:10]}. Polygons in {rec.get('crs') or 'UTM'} (WGS 84 / UTM zone of their centre); the "
	         'field cat is the ID of each feature (a study area for setup.py) and name its name.', '', 'Files:']
	for x in rec['files']:
		lines.append(f"  {x['file']}: {x['label']}, {x['features']} {'feature' if x['features'] == 1 else 'features'}, "
		             f"{x['ha']:,.1f} ha")
	lines += ['', 'How it was made:']
	if rec['kind'] == 'burned':
		lines += [f"  Burn index {p['index']} of {p.get('sensor') or 'Sentinel-2'}: the median NBR of the images of the "
		          f"period before the fire ({p['before'][0]} to {p['before'][1]}) and of the period after it "
		          f"({p['after'][0]} to {p['after'][1]}), with the pixel filters of MonitorEO (clouds, shadows, water, "
		          f"snow), on pixels of {rec.get('scale_m') or 0:g} m (services.burn_severity).",
		          f"  Searched in {rec.get('search')}.",
		          f"  Burned: the pixels with {p['index']} of at least {rec.get('burned_min')}, joined when they touch "
		          f"(also by a corner), without the patches of less than {rec.get('min_ha')} ha and with the unburned "
		          f"islands of less than {rec.get('min_ha')} ha inside them filled (the perimeter of monitoreo).",
		          f"  The files of monitoreo are in {MONITOREO_DIR}/: the perimeter in EPSG:4326, the table of the burn "
		          'index and its README.']
		if p.get('ring'):
			lines.append(f"  Unburned ring: between {p.get('gap') or 0:g} and {(p.get('gap') or 0) + p['ring']:g} m "
			             'outside the burned area with its holes filled (computed by the app with GDAL).')
	else:
		what = 'WDPA ID' if rec['kind'] == 'protected' else f"GAUL{p.get('level')}_CODE"
		lines.append(f"  {what}: " + ', '.join(f"{x['id']} ({x.get('full') or x.get('name')})" for x in p['items']))
		lines.append(f"  Polygons fetched from Earth Engine, simplified by {SIMPLIFY_M} m"
		             + (', joined into one feature.' if p.get('join') and len(p['items']) > 1 else
		                ', one feature each (its parts joined).'))
	refs = cited.get('references') or []
	if refs:
		lines += ['', 'References:'] + [f"  {r.get('text') or r.get('citation')}" for r in refs]
	if cited.get('credits'):
		lines += ['', 'Credits:'] + [f'  {c}' for c in cited['credits']]
	with open(os.path.join(folder, README), 'w', encoding='utf-8') as f:
		f.write('\n'.join(lines) + '\n')
