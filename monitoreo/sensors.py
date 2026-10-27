"""Collections of the sensors of MonitorEO: quality masks, scale and offset of the bands and formula of the variable.

Port of sensors/*_cut, modules/sensorList and processFt (modules/auxCalcIndex) of MonitorEO, with these changes:
  - the pixel filters of the app ('remove pixels with clouds, water, snow, aerosols'), which are switched off
    there, with the QA bands of the catalog (*_bandQA): SCL for Sentinel-2 (Cloud Score+ for its clouds),
    QA_PIXEL and SR_QA_AEROSOL for Landsat and the QA bands of each MODIS product;
  - the scale and offset are applied to the bands before the formula (the app evaluates the formulas with spaces,
    e.g. EVI or ALBEDO, on the raw bands, and the + 1 of EVI or the - 0.018 of the Landsat albedo then mean nothing);
  - band names of the formulas are matched without case (LST_day_1km is the band LST_Day_1km);
  - flag values of the MODIS products (e.g. 250 = cloud in NDSI_Snow_Cover) are masked;
  - Landsat 7 is only used until 2003-05-31, as the catalog says (failure of its scan line corrector).
"""

import re

import ee

# Bands of Landsat 5 and 7, and of Landsat 8 and 9, with the names used by the formulas of the catalog
LANDSAT_57 = {'SR_B1': 'blue', 'SR_B2': 'green', 'SR_B3': 'red', 'SR_B4': 'nir', 'SR_B5': 'swir1', 'SR_B7': 'swir2',
              'ST_B6': 'st', 'QA_PIXEL': 'QA_PIXEL'}
LANDSAT_89 = {'SR_B2': 'blue', 'SR_B3': 'green', 'SR_B4': 'red', 'SR_B5': 'nir', 'SR_B6': 'swir1', 'SR_B7': 'swir2',
              'ST_B10': 'st', 'QA_PIXEL': 'QA_PIXEL', 'SR_QA_AEROSOL': 'SR_QA_AEROSOL'}
LANDSAT = [  # collection, bands, last date
	('LANDSAT/LT05/C02/T1_L2', LANDSAT_57, None),
	('LANDSAT/LE07/C02/T1_L2', LANDSAT_57, '2003-06-01'),
	('LANDSAT/LC08/C02/T1_L2', LANDSAT_89, None),
	('LANDSAT/LC09/C02/T1_L2', LANDSAT_89, None),
]
S2 = 'COPERNICUS/S2_SR_HARMONIZED'
CLOUD_SCORE = 'GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED'

# Scene cloud cover property of each sensor, for quality['scene_clouds']
SCENE_CLOUDS = {'Sentinel-2': 'CLOUDY_PIXEL_PERCENTAGE', 'Landsat': 'CLOUD_COVER_LAND'}

# Quality options and their defaults. The pixel filters of MonitorEO (clouds, water, snow, aerosols) remove, in
# every image, the pixels that its QA bands flag; all of them are on by default. 'auto' = on, except for the
# variables that measure what the filter removes (MEASURES).
QUALITY = {
	'clouds': True,        # clouds and cloud shadows
	'water': True,         # water
	'snow': True,          # snow and ice
	'aerosols': True,      # high aerosol content (Landsat 8 and 9, MODIS NDVI, EVI and LAI)
	'cloud_score': 0.6,    # Sentinel-2 clouds: minimum Cloud Score+ (cs_cdf) of a clear pixel; None = SCL classes,
	                       # as MonitorEO (in Guajares Cloud Score+ gave an NDVI p5 closer to the reference)
	'scene_clouds': None,  # Sentinel-2 and Landsat: maximum cloud cover of a scene (%), None = all scenes
}
FILTERS = ('clouds', 'water', 'snow', 'aerosols')
MEASURES = {  # variables whose signal a filter removes
	'water': {'NDWI', 'Chl-a', 'NDCI'},  # open water and water quality
	'snow': {'NDSI'},
	'aerosols': {'ARSL'},
}
S2_SCL = {'clouds': [3, 8, 9, 10], 'water': [6], 'snow': [11]}  # cloud shadow, clouds, cirrus; water; snow


def bits(band, first, last=None):
	"""Value of the bits FIRST to LAST of a QA band."""
	last = first if last is None else last
	return band.rightShift(first).bitwiseAnd((1 << (last - first + 1)) - 1)


def _modis_cloud_state(qa):
	# bits 3-4: 0 clear, 1 cloudy, 2 mixed, 3 not set (assumed clear)
	state = bits(qa, 3, 4)
	return state.eq(0).Or(state.eq(3))


def _vi_clouds(i):
	# MOD13Q1: SummaryQA 3 = cloudy; DetailedQA bit 8 adjacent cloud, 10 mixed clouds, 15 shadow
	qa = i.select('DetailedQA')
	return i.select('SummaryQA').neq(3).And(bits(qa, 8).eq(0)).And(bits(qa, 10).eq(0)).And(bits(qa, 15).eq(0))


# MODIS (and other) products: band of valid raw values (flags above max) and the mask of each pixel filter
# {collection: {'valid': (band, max raw value), filter: f(img) -> 1 = keep}}
PRODUCTS = {
	'MODIS/061/MOD13Q1': {  # DetailedQA: bits 6-7 aerosols (3 high), 11-13 land/water (1 land), 14 snow/ice
		'clouds': _vi_clouds,
		'water': lambda i: bits(i.select('DetailedQA'), 11, 13).eq(1),
		'snow': lambda i: i.select('SummaryQA').neq(2).And(bits(i.select('DetailedQA'), 14).eq(0)),
		'aerosols': lambda i: bits(i.select('DetailedQA'), 6, 7).neq(3),
	},
	'MODIS/061/MOD11A1': {'clouds': lambda i: bits(i.select('QC_Day'), 0, 1).eq(0)},  # LST of good quality
	'MODIS/061/MOD16A2': {'valid': ('ET', 32700), 'clouds': lambda i: _modis_cloud_state(i.select('ET_QC'))},
	'MODIS/061/MCD15A3H': {  # FparExtra_QC: bits 0-1 land/sea (0 land), 2 snow, 3 aerosol, 4 cirrus, 5 cloud, 6 shadow
		'valid': ('Lai', 100),
		'clouds': lambda i: _modis_cloud_state(i.select('FparLai_QC')).And(bits(i.select('FparExtra_QC'), 4, 6).eq(0)),
		'water': lambda i: bits(i.select('FparExtra_QC'), 0, 1).eq(0),
		'snow': lambda i: bits(i.select('FparExtra_QC'), 2).eq(0),
		'aerosols': lambda i: bits(i.select('FparExtra_QC'), 3).eq(0),
	},
	# NDSI_Snow_Cover 200-255 are flags (the classes of NDSI_Snow_Cover_Class: 237 inland water, 239 ocean,
	# 250 cloud, ...), not values: masked always
	'MODIS/061/MOD10A1': {
		'valid': ('NDSI_Snow_Cover', 100),
		'clouds': lambda i: i.select('NDSI_Snow_Cover_Class').unmask(0).neq(250),
		'water': lambda i: i.select('NDSI_Snow_Cover_Class').unmask(0).remap([237, 239], [0, 0], 1),
	},
	'MODIS/061/MCD19A2_GRANULES': {  # AOD_QA: bits 0-2 cloud mask (1 clear), 3-4 land (0), water (1), snow, ice
		'clouds': lambda i: bits(i.select('AOD_QA'), 0, 2).eq(1),
		'water': lambda i: bits(i.select('AOD_QA'), 3, 4).neq(1),
		'snow': lambda i: bits(i.select('AOD_QA'), 3, 4).lt(2),
	},
	# 0 = full BRDF inversion, 1 = magnitude inversion (too few clear observations)
	'MODIS/061/MCD43A3': {'clouds': lambda i: i.select('BRDF_Albedo_Band_Mandatory_Quality_shortwave').eq(0)},
}
PRODUCTS['MODIS/061/MYD13Q1'] = PRODUCTS['MODIS/061/MOD13Q1']
PRODUCTS['MODIS/061/MYD11A1'] = PRODUCTS['MODIS/061/MOD11A1']
PRODUCTS['MODIS/061/MYD10A1'] = PRODUCTS['MODIS/061/MOD10A1']


def filters(entry, quality=None):
	"""Pixel filters (clouds, water, snow, aerosols) that apply to ENTRY with QUALITY: {filter: True/False}, with
	'auto' resolved and False for the filters its sensor has no QA band for (e.g. aerosols of Sentinel-2)."""
	quality = dict(QUALITY, **(quality or {}))
	sensor, col = entry['sensor'], entry['col']
	available = {'Sentinel-2': {'clouds', 'water', 'snow'}, 'Landsat': set(FILTERS)}.get(sensor)
	if available is None:
		available = {f for f in FILTERS if f in PRODUCTS.get(col, {})}
	out = {}
	for f in FILTERS:
		on = quality[f]
		if str(on).lower() == 'auto':
			on = entry['index'] not in MEASURES.get(f, ())
		out[f] = bool(on) and f in available
	return out


def conflicts(entry, quality=None):
	"""Filters on for ENTRY that remove what its variable measures (e.g. snow for NDSI)."""
	on = filters(entry, quality)
	return [f for f in FILTERS if on[f] and entry['index'] in MEASURES.get(f, ())]

_band_names = {}


def number(v, default):
	"""Number of a field of the catalog ('#' or missing = DEFAULT)."""
	try:
		return float(v)
	except (TypeError, ValueError):
		return default


def band_names(entry):
	"""Band names of the images of the collection of ENTRY."""
	if entry['sensor'] == 'Landsat':
		return list(LANDSAT_89.values())
	col = entry['col']
	if col not in _band_names:
		_band_names[col] = ee.ImageCollection(col).first().bandNames().getInfo()
	return _band_names[col]


def parse_formula(entry):
	"""(kind, bands) of the formula of ENTRY: ('nd', [b1, b2]) for N(b1,b2), else ('expr', {name: band}) with the
	band of every name of the expression (matched without case)."""
	ft = entry['ft'].strip()
	names = band_names(entry)
	lower = {b.lower(): b for b in names}

	def band(name):
		b = lower.get(name.strip().lower())
		if b is None:
			raise ValueError(f"band {name.strip()} of the formula {ft} of {entry['index']} ({entry['sensor']}) "
			                 f"not found; bands: {', '.join(names)}")
		return b

	if ft.startswith('N(') and ft.endswith(')'):
		return 'nd', [band(b) for b in ft[2:-1].split(',')]
	return 'expr', {name: band(name) for name in set(re.findall(r'[A-Za-z_][A-Za-z0-9_]*', ft))}


def quality_mask(img, entry, quality):
	"""Mask of the pixels of an image of the collection of ENTRY that pass the pixel filters of QUALITY
	(None = no mask)."""
	sensor, col = entry['sensor'], entry['col']
	on = filters(entry, quality)
	masks = []
	if sensor == 'Sentinel-2':
		scl = img.select('SCL')
		for f, classes in S2_SCL.items():
			if not on[f]:
				continue
			if f == 'clouds' and quality.get('cloud_score') is not None:
				masks.append(img.select('cs_cdf').gte(quality['cloud_score']))
			else:
				masks.append(scl.remap(classes, [0] * len(classes), 1))
	elif sensor == 'Landsat':
		qa = img.select('QA_PIXEL')
		if on['clouds']:
			masks.append(bits(qa, 1, 4).eq(0))  # dilated cloud, cirrus, cloud, cloud shadow
		if on['snow']:
			masks.append(bits(qa, 5).eq(0))
		if on['water']:
			masks.append(bits(qa, 7).eq(0))
		if on['aerosols']:
			masks.append(bits(img.select('SR_QA_AEROSOL'), 6, 7).neq(3))  # aerosol level 3 = high
	else:
		product = PRODUCTS.get(col, {})
		if 'valid' in product:
			b, vmax = product['valid']
			masks.append(img.select(b).lte(vmax))
		masks += [product[f](img) for f in FILTERS if on[f]]
	if not masks:
		return None
	mask = masks[0]
	for m in masks[1:]:
		mask = mask.And(m)
	return mask


def raw_collection(entry, region, start, end, quality):
	"""Images of the collection of ENTRY over REGION between START and END (excluded), with the band names of the
	formulas, filtered by the scene cloud cover if asked."""
	sensor = entry['sensor']
	prop = SCENE_CLOUDS.get(sensor)
	scene = quality.get('scene_clouds')

	def scene_filter(col):
		return col.filter(ee.Filter.lte(prop, scene)) if prop and scene is not None else col

	if sensor == 'Landsat':
		cols = []
		for col_id, aliases, last in LANDSAT:
			end_l = min(end, last) if last else end
			if end_l <= start:
				continue
			c = scene_filter(ee.ImageCollection(col_id).filterBounds(region).filterDate(start, end_l))
			c = c.select(list(aliases), list(aliases.values()))
			if 'SR_QA_AEROSOL' not in aliases:  # Landsat 5 and 7 have no aerosol QA: nothing is removed
				c = c.map(lambda i: i.addBands(ee.Image.constant(0).toUint8().rename('SR_QA_AEROSOL')))
			cols.append(c)
		col = cols[0]
		for c in cols[1:]:
			col = col.merge(c)
		return col
	col = scene_filter(ee.ImageCollection(entry['col']).filterBounds(region).filterDate(start, end))
	if sensor == 'Sentinel-2' and filters(entry, quality)['clouds'] and quality.get('cloud_score') is not None:
		col = col.linkCollection(ee.ImageCollection(CLOUD_SCORE), ['cs_cdf'])
	return col


def collection(entry, region, start, end, quality=None, resample=None):
	"""Collection of single band images ('value', float) of the variable of ENTRY over REGION between START and END
	(dates 'YYYY-MM-DD', END excluded), with the quality masks of QUALITY (see QUALITY) and, if given, the
	resampling method RESAMPLE ('bilinear' or 'bicubic') of the bands when the output pixels are smaller."""
	quality = dict(QUALITY, **(quality or {}))
	scale = number(entry.get('scale'), 1.0)
	offset = number(entry.get('offset'), 0.0)
	kind, bands = parse_formula(entry)
	used = list(bands) if kind == 'nd' else sorted(set(bands.values()))

	def compute(img):
		img = ee.Image(img)
		data = img.select(used).multiply(scale).add(offset)
		if resample:
			data = data.resample(resample)
		mask = quality_mask(img, entry, quality)
		if mask is not None:
			data = data.updateMask(mask)
		if kind == 'nd':
			out = data.normalizedDifference(bands)
		else:
			out = data.expression(entry['ft'], {name: data.select(b) for name, b in bands.items()})
		return out.toFloat().rename('value').copyProperties(img, ['system:time_start'])

	return raw_collection(entry, region, start, end, quality).map(compute)
