"""Study areas and ecoregions (RESOLVE Ecoregions 2017, the Ecoregions2017 layer, in Earth Engine)."""

import json
import math

import ee
from osgeo import ogr, osr

ogr.UseExceptions()
osr.UseExceptions()

ECOREGIONS = 'RESOLVE/ECOREGIONS/2017'
ECO_FIELDS = ['ECO_ID', 'ECO_NAME', 'BIOME_NAME', 'REALM']


def wgs84():
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(4326)
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return srs


def srs_of(crs):
	"""OSR spatial reference of a CRS given as 'EPSG:code' (or any other text that OSR understands, e.g. WKT)."""
	srs = osr.SpatialReference()
	srs.SetFromUserInput(crs)
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	return srs


def read_vector(path, layer=None):
	"""Union of the polygons of a vector file (OGR geometry, in the CRS of the file) and that CRS."""
	ds = ogr.Open(path)
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	srs = lyr.GetSpatialRef()
	if srs is None:
		raise ValueError(f'{path} has no CRS')
	srs = srs.Clone()
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	union = ogr.Geometry(ogr.wkbMultiPolygon)
	for feat in lyr:
		g = feat.GetGeometryRef()
		if g is not None:
			union = union.Union(g)
	if union.IsEmpty():
		raise ValueError(f'{path} has no polygons')
	return union, srs


def transform(geom, src, dst):
	g = geom.Clone()
	g.Transform(osr.CoordinateTransformation(src, dst))
	return g


def utm_crs(geom, srs):
	"""'EPSG:326zz' or 'EPSG:327zz': WGS 84 / UTM zone of the centre of GEOM."""
	c = transform(geom, srs, wgs84()).Centroid()
	lon, lat = c.GetX(), c.GetY()
	zone = min(int(math.floor((lon + 180) / 6)) + 1, 60)
	return f'EPSG:{(32600 if lat >= 0 else 32700) + zone}'


def bounds(geom, srs, crs):
	"""(xmin, ymin, xmax, ymax) of GEOM in CRS."""
	xmin, xmax, ymin, ymax = transform(geom, srs, srs_of(crs)).GetEnvelope()
	return xmin, ymin, xmax, ymax


def to_ee(geom, srs, simplify=None):
	"""Earth Engine geometry (EPSG:4326) of an OGR geometry, simplified by SIMPLIFY (in the units of SRS)."""
	if simplify:
		geom = geom.SimplifyPreserveTopology(simplify)
	g = transform(geom, srs, wgs84())
	return ee.Geometry(json.loads(g.ExportToJson()), opt_geodesic=False)


def ecoregions(geom):
	"""Ecoregions that intersect GEOM (ee.Geometry), as an ee.FeatureCollection."""
	return ee.FeatureCollection(ECOREGIONS).filterBounds(geom)


def ecoregion_info(eco):
	"""[{ECO_ID, ECO_NAME, BIOME_NAME, REALM}] of the ecoregions of the collection ECO."""
	return eco.select(ECO_FIELDS).map(lambda f: ee.Feature(None, f.toDictionary())).getInfo()['features']


def ecoregion_bounds(eco, crs, max_error):
	"""(xmin, ymin, xmax, ymax) of the union of the ecoregions ECO in CRS."""
	ring = eco.geometry(max_error).bounds(max_error, crs).getInfo()['coordinates'][0]
	xs, ys = [p[0] for p in ring], [p[1] for p in ring]
	return min(xs), min(ys), max(xs), max(ys)


def save_ecoregions(eco, path, max_error):
	"""Write the ecoregions ECO, simplified by MAX_ERROR metres, to a GeoPackage (EPSG:4326)."""
	fc = eco.select(ECO_FIELDS).map(lambda f: f.simplify(max_error)).getInfo()
	drv = ogr.GetDriverByName('GPKG')
	ds = drv.CreateDataSource(path)
	lyr = ds.CreateLayer('ecoregions', wgs84(), ogr.wkbMultiPolygon)
	lyr.CreateField(ogr.FieldDefn('ECO_ID', ogr.OFTInteger))
	for name in ECO_FIELDS[1:]:
		lyr.CreateField(ogr.FieldDefn(name, ogr.OFTString))
	for f in fc['features']:
		feat = ogr.Feature(lyr.GetLayerDefn())
		for name in ECO_FIELDS:
			feat.SetField(name, f['properties'].get(name))
		g = ogr.CreateGeometryFromJson(json.dumps(f['geometry']))
		feat.SetGeometry(ogr.ForceToMultiPolygon(g))
		lyr.CreateFeature(feat)
	ds = None
