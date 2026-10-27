"""Earth Engine image of a variable described by a dictionary (a 'spec'), of one of three kinds:

catalog     a variable of MonitorEO: index (NDVI, LST, ...), sensor (Sentinel-2, Landsat, Modis) and, optionally,
            collection (when the catalog has several for that variable and sensor) or ref (ID of the entry).
collection  any other image collection of Earth Engine: collection and band (e.g. IDAHO_EPSCOR/TERRACLIMATE, pr),
            with optional scale and offset of the band.
image       a fixed image (or a collection of tiles, mosaicked): image and band (e.g. COPERNICUS/DEM/GLO30, DEM),
            with optional scale, offset and terrain (slope or aspect, in degrees, computed from the band).

The variables of the catalog and of other collections are aggregated in time with metric, temporal and
interannual over period and season (see composite). Other keys:
  quality            pixel filters (clouds, water, snow, aerosols) and scene filter of the catalog variables
                     (see sensors.QUALITY; by default all the pixel filters of MonitorEO)
  resample           'auto' (bilinear when the output pixels are smaller than those of the data), 'nearest',
                     'bilinear' or 'bicubic'
  aggregate          fixed images on pixels more than twice as big as theirs: 'auto' or 'mean' (mean of the
                     pixels of the data inside each output pixel) or 'none' (Earth Engine pyramids)
  exclude_landcover  classes of ESA WorldCover 2021 (10 m) set to NoData, e.g. [40, 50, 80]: crops, built-up, water
"""

import re

import ee

from . import catalog, composite, sensors

WORLDCOVER = 'ESA/WorldCover/v200'
SENSOR_ABBR = {'Sentinel-2': 'S2', 'Landsat': 'L', 'Modis': 'MODIS'}
KEYS = {'name', 'index', 'sensor', 'collection', 'ref', 'image', 'band', 'scale', 'offset', 'terrain', 'metric',
        'metrics', 'temporal', 'interannual', 'period', 'season', 'quality', 'resample', 'aggregate',
        'exclude_landcover'}
DEFAULTS = {'temporal': 'period', 'interannual': 'mean', 'resample': 'auto', 'aggregate': 'auto', 'quality': {},
            'season': None, 'exclude_landcover': None}

_native = {}
_types = {}


def kind(spec):
	if spec.get('index'):
		return 'catalog'
	if spec.get('image'):
		return 'image'
	if spec.get('collection'):
		return 'collection'
	raise ValueError(f'variable {spec.get("name", spec)}: give index (and sensor), collection and band, or image')


def _metric_label(spec):
	m = str(spec['metric']).lower()
	if spec['temporal'] == 'period':
		return m
	return f"{spec['temporal']}-{m}_{str(spec['interannual']).lower()}"


def _default_name(spec):
	k = kind(spec)
	if k == 'catalog':
		base = f"{spec['index']}_{SENSOR_ABBR.get(spec.get('sensor'), spec.get('sensor') or '')}".rstrip('_')
	else:
		base = f"{spec.get(k).split('/')[-1]}_{spec.get('band') or ''}".rstrip('_')
		if spec.get('terrain'):
			base += '_' + spec['terrain']
	return base


def specs(var, defaults=None):
	"""Specs of a variable of a recipe: one per metric of 'metric' (a metric or a list of them), with the keys
	missing in VAR taken from DEFAULTS and DEFAULTS of this module, and 'name' (a name for files: letters,
	digits, '-' and '_')."""
	unknown = set(var) - KEYS
	if unknown:
		raise ValueError(f"variable {var.get('name', var)}: unknown keys {', '.join(sorted(unknown))}")
	base = {**DEFAULTS, **{k: v for k, v in (defaults or {}).items() if k in KEYS}, **var}
	# all the quality options, so that the record of a file changes when a default changes
	base['quality'] = {**sensors.QUALITY, **((defaults or {}).get('quality') or {}), **(var.get('quality') or {})}
	unknown = set(base['quality']) - set(sensors.QUALITY)
	if unknown:
		raise ValueError(f"unknown quality options {', '.join(sorted(unknown))}: use {', '.join(sensors.QUALITY)}")
	k = kind(base)
	if k == 'image':
		metrics = [None]
	else:
		metrics = base.pop('metrics', None) or base.get('metric')
		if metrics is None:
			raise ValueError(f"variable {base.get('name') or _default_name(base)}: give metric (e.g. median or [p5, p95])")
		metrics = metrics if isinstance(metrics, list) else [metrics]
		for m in metrics + ([base['interannual']] if base['temporal'] != 'period' else []):
			composite.check_metric(m)
		if base['temporal'] not in composite.TEMPORAL:
			raise ValueError(f"unknown temporal {base['temporal']}: use {', '.join(composite.TEMPORAL)}")
		if not base.get('period'):
			raise ValueError(f'variable {base.get("name") or _default_name(base)}: give period [first day, last day]')
		composite.windows(base['period'], base['season'])
	out = []
	for m in metrics:
		s = dict(base, metric=m)
		s.pop('metrics', None)
		name = var.get('name') or _default_name(s)
		if m is not None and (len(metrics) > 1 or not var.get('name')):
			name = f'{name}_{_metric_label(s)}'
		s['name'] = re.sub(r'[^A-Za-z0-9_-]', '_', str(name))
		out.append(s)
	return out


def entry_of(spec):
	return catalog.find(spec['index'], spec.get('sensor'), spec.get('collection'), spec.get('ref'))


def _span(spec):
	wins = composite.windows(spec['period'], spec['season'])
	return wins[0][1].isoformat(), wins[-1][2].isoformat()


def native_scale(spec):
	"""Pixel size (m) of the data of SPEC."""
	k = kind(spec)
	if k == 'catalog':
		return float(entry_of(spec)['m'])
	key = (spec.get(k), spec.get('band'))
	if key not in _native:
		src = ee.Image(spec['image']) if k == 'image' and _is_image(spec['image']) else ee.ImageCollection(spec.get(k)).first()
		_native[key] = src.select(spec['band']).projection().nominalScale().getInfo()
	return _native[key]


def _is_image(asset):
	if asset not in _types:
		_types[asset] = ee.data.getAsset(asset)['type']
	return _types[asset] == 'IMAGE'


def _resample(spec, res, native):
	method = str(spec.get('resample') or 'auto').lower()
	if method == 'auto':
		return 'bilinear' if res < native else None
	return None if method == 'nearest' else method


def _fixed(spec, res):
	"""Image of a spec of kind 'image'."""
	if _is_image(spec['image']):
		img = ee.Image(spec['image']).select(spec['band'])
	else:
		col = ee.ImageCollection(spec['image']).select(spec['band'])
		img = col.mosaic().setDefaultProjection(col.first().projection())
	native = native_scale(spec)
	img = img.multiply(float(spec.get('scale') or 1)).add(float(spec.get('offset') or 0))
	method = _resample(spec, res, native)
	if method:
		img = img.resample(method)
	terrain = spec.get('terrain')
	if terrain:
		ops = {'slope': ee.Terrain.slope, 'aspect': ee.Terrain.aspect}
		if terrain not in ops:
			raise ValueError(f"unknown terrain {terrain}: use {', '.join(ops)}")
		img = ops[terrain](img)
	agg = str(spec.get('aggregate') or 'auto').lower()
	if agg in ('auto', 'mean') and res > 2 * native:
		img = img.reduceResolution(ee.Reducer.mean(), True, 65535)
	return img.toFloat().rename('value')


def _from_collection(spec, region, res):
	"""Collection of single band images ('value') of a spec of kind 'collection'."""
	start, end = _span(spec)
	scale, offset = float(spec.get('scale') or 1), float(spec.get('offset') or 0)
	method = _resample(spec, res, native_scale(spec))

	def prep(img):
		img = ee.Image(img)
		out = img.select(spec['band']).multiply(scale).add(offset)
		if method:
			out = out.resample(method)
		return out.toFloat().rename('value').copyProperties(img, ['system:time_start'])

	return ee.ImageCollection(spec['collection']).filterBounds(region).filterDate(start, end).map(prep)


def collection(spec, region, res):
	"""Collection of single band images ('value') of a spec of kind 'catalog' or 'collection', before the temporal
	aggregation."""
	if kind(spec) == 'collection':
		return _from_collection(spec, region, res)
	entry = entry_of(spec)
	start, end = _span(spec)
	return sensors.collection(entry, region, start, end, spec['quality'], _resample(spec, res, float(entry['m'])))


def image(spec, region, res):
	"""Image ('value') of the variable SPEC over REGION (ee.Geometry) for output pixels of RES metres."""
	if kind(spec) == 'image':
		img = _fixed(spec, res)
	else:
		img = composite.composite(collection(spec, region, res), spec['metric'], spec['temporal'],
		                          spec['interannual'], spec['period'], spec['season'])
	classes = spec.get('exclude_landcover')
	if classes:
		wc = ee.ImageCollection(WORLDCOVER).first()
		img = img.updateMask(wc.remap(list(classes), [0] * len(classes), 1))
	return img


def n_images(spec, region):
	"""Number of images of the period (and season) of SPEC over REGION (None for a fixed image)."""
	if kind(spec) == 'image':
		return None
	wins = composite.windows(spec['period'], spec['season'])
	col = collection(dict(spec, resample='nearest'), region, 1e9)
	return col.filter(composite.date_filter(wins)).size().getInfo()


def describe(spec):
	"""What SPEC computes, for the record of the downloaded files."""
	k = kind(spec)
	d = {'kind': k}
	if k == 'catalog':
		e = entry_of(spec)
		d.update({'monitoreo_ref': e['ref'], 'variable': e['vname_en'], 'collection': e['col'],
		          'formula': e['ft'], 'scale': e.get('scale'), 'offset': e.get('offset'), 'unit': e['unit'],
		          'native_resolution': e['m'], 'pixel_filters': sensors.filters(e, spec['quality'])})
		if e['ref'] in catalog.FIXES:
			d['fixed'] = catalog.FIXES[e['ref']]
		if e.get('notes', {}).get('ft'):
			d['note'] = e['notes']['ft']
	return d
