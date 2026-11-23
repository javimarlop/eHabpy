#### eHabitat+ desktop app: test of the search of protected areas by their WDPA ID (areas.wdpa_ids, find_protected)
#### License: CC BY-SA 3.0
"""
Checks which texts of the search of protected areas of the app are WDPA IDs (only numbers, separated by commas,
semicolons or spaces) and which are names, and that a name of one letter is refused while an ID of one digit is not.
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

	real = areas._connected, areas._protected_by_id
	areas._connected = Services
	areas._protected_by_id = lambda ids: calls.append(('id', ids)) or [{'SITE_ID': i, 'NAME': f'site {i}'} for i in ids]
	try:
		assert [r['id'] for r in areas.find_protected('7', ramsar=True)] == [7]
		assert areas.find_protected('Sierra', ramsar=True)[0]['name'] == 'Sierra'
	finally:
		areas._connected, areas._protected_by_id = real
	assert calls == [('id', [7]), ('name', 'Sierra', 'ramsar')], calls
	print('test_wdpa_ids: OK')


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
	test_find_by_id()
	print('All the tests passed')
