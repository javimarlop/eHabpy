#### eHabitat+ desktop app: test of the search of protected areas by their WDPA ID (areas.wdpa_ids, find_protected)
#### License: CC BY-SA 3.0
"""
Checks which texts of the search of protected areas of the app are WDPA IDs (only numbers, separated by commas,
semicolons or spaces) and which are names, and that a name of one letter is refused while an ID of one digit is not;
and the area of each one found: the reported one (REP_AREA), else that of its polygons (GIS_AREA), never 0.
With EE_PROJECT set, it also finds in Earth Engine the Sierra Nevada National Park (WDPA ID 389011) and Doñana
(389012) by their IDs.

    cd segm_optim && python -s tests/test_wdpa_id.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'app'))  # ehabitat_app

from ehabitat_app import areas  # noqa: E402


def test_wdpa_ids():
	assert areas.wdpa_ids('389011') == [389011]
	assert areas.wdpa_ids(' 389011, 389012 ') == [389011, 389012]
	assert areas.wdpa_ids('389011;389012 555588816') == [389011, 389012, 555588816]
	assert areas.wdpa_ids('7') == [7]
	for name in ('Sierra Nevada', 'Doñana', '', None, '389011a', 'Parque 389011', '389011,'):
		assert areas.wdpa_ids(name) is None, name
	try:
		areas.find_protected('S')
	except ValueError as e:
		assert 'WDPA ID' in str(e)
	else:
		raise AssertionError('a name of one letter was not refused')
	# an ID of one digit reaches the search by ID (without Earth Engine: both replaced), and a name the search by name
	calls = []

	class Services:
		def find_protected_areas(self, text, limit, designation):
			calls.append(('name', text, designation))
			return [{'SITE_ID': 1, 'NAME': 'Sierra'}]

	real = areas._connected, areas._protected_by_id, areas._gis_areas
	areas._connected = Services
	areas._protected_by_id = lambda ids: calls.append(('id', ids)) or [{'SITE_ID': i, 'NAME': f'site {i}', 'REP_AREA': 1}
	                                                                   for i in ids]
	areas._gis_areas = lambda ids: calls.append(('gis', ids)) or {1: 1723.7}
	try:
		assert [r['id'] for r in areas.find_protected('7', ramsar=True)] == [7]
		row = areas.find_protected('Sierra', ramsar=True)[0]
		assert row['name'] == 'Sierra' and row['km2'] == 1723.7 and row['km2_from'] == 'polygons', row
	finally:
		areas._connected, areas._protected_by_id, areas._gis_areas = real
	# the area of the polygons is asked only for the rows without a reported area
	assert calls == [('id', [7]), ('name', 'Sierra', 'ramsar'), ('gis', [1])], calls
	print('test_wdpa_ids: OK')


def test_area_km2():
	"""The area of a protected area in the search: the reported one (REP_AREA), else that of its polygons (GIS_AREA, in
	the row or asked for), else none (not 0, which WDPA gives when it was not reported)."""
	assert areas.km2(858.83) == 858.83 and areas.km2('12.5') == 12.5
	for v in (0, 0.0, '0', None, '', 'n/a', -3, float('nan'), float('inf')):
		assert areas.km2(v) is None, v
	rows = [{'SITE_ID': 389011, 'NAME': 'Sierra Nevada', 'REP_AREA': 858.8289465, 'GIS_AREA': 859.1},
	        {'SITE_ID': 555722086, 'NAME': 'Sierra Nevada', 'REP_AREA': 0},
	        {'WDPAID': 12, 'NAME_ENG': 'Old', 'REP_AREA': 0, 'GIS_AREA': 3.2},
	        {'SITE_ID': 99, 'NAME': 'Nothing', 'REP_AREA': 0}]
	got = areas._protected_rows(rows, {555722086: 1723.73})
	assert [(r['id'], r['km2'], r['km2_from']) for r in got] == [
		(389011, 858.8289465, 'reported'), (555722086, 1723.73, 'polygons'), (12, 3.2, 'polygons'), (99, None, None)], got
	assert got[2]['name'] == 'Old'
	print('test_area_km2: OK')


def test_preview_geometry():
	"""The outlines of the preview: coordinates rounded to 5 decimals, and none when a geometry is empty (a small polygon
	simplified away), also inside a GeometryCollection."""
	poly = {'type': 'Polygon', 'coordinates': [[[-3.123456789, 37.1], [-3.2, 37.123456789], [-3.1, 37.2], [-3.123456789, 37.1]]]}
	g = areas.preview_geometry(poly)
	assert g['coordinates'][0][0] == [-3.12346, 37.1] and g['coordinates'][0][1] == [-3.2, 37.12346], g
	for empty in (None, {}, {'type': 'Polygon', 'coordinates': []}, {'type': 'Polygon', 'coordinates': [[]]},
	              {'type': 'GeometryCollection', 'geometries': []},
	              {'type': 'GeometryCollection', 'geometries': [{'type': 'Polygon', 'coordinates': []}]}):
		assert areas.preview_geometry(empty) is None, empty
	gc = areas.preview_geometry({'type': 'GeometryCollection', 'geometries': [poly, {'type': 'Polygon', 'coordinates': []},
	                                                                         {'type': 'Point', 'coordinates': [1.0, 2.0]}]})
	assert [x['type'] for x in gc['geometries']] == ['Polygon', 'Point'], gc
	multi = {'type': 'MultiPolygon', 'coordinates': [poly['coordinates']]}
	assert areas.preview_geometry(multi)['coordinates'][0][0][0] == [-3.12346, 37.1]
	print('test_preview_geometry: OK')


def test_find_by_id():
	"""In Earth Engine (EE_PROJECT): the two national parks by their WDPA IDs."""
	if not os.environ.get('EE_PROJECT'):
		print('test_find_by_id: skipped (no EE_PROJECT)')
		return
	from ehabitat_app import gee
	gee._monitoreo()[1].init(os.environ['EE_PROJECT'])
	rows = areas.find_protected('389011, 389012')
	got = {r['id']: r['name'] for r in rows}
	assert set(got) == {389011, 389012}, got
	assert 'Sierra Nevada' in got[389011] and 'Doñana' in got[389012], got
	assert all(r['iso3'] == 'ESP' and r['km2'] for r in rows), rows
	print('test_find_by_id: OK', got)


if __name__ == '__main__':
	test_wdpa_ids()
	test_area_km2()
	test_preview_geometry()
	test_find_by_id()
	print('All the tests passed')
