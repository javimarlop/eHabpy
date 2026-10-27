"""Catalog of variables of MonitorEO (catalog.json, made by tools/import_catalog.py)."""

import json
import os

CATALOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'catalog.json')

# Corrections of entries of MonitorEO, applied when an entry is read (catalog.json keeps the values of the app):
# ref: {field: value}
FIXES = {
	32: {'scale': 0.1},  # Lai of MCD15A3H has a scale of 0.1 (raw 58 = 5.8 m2/m2 in the Amazon rain forest), not 0.01
}

_catalog = None


def load():
	"""The whole catalog: source, fields, collections (date ranges) and variables."""
	global _catalog
	if _catalog is None:
		with open(CATALOG, encoding='utf-8') as f:
			_catalog = json.load(f)
	return _catalog


def variables():
	"""Entries of the catalog, with FIXES applied."""
	return [dict(e, **FIXES.get(e['ref'], {})) for e in load()['variables']]


def find(index, sensor=None, collection=None, ref=None):
	"""Entry of the catalog of a variable: by REF, or by INDEX (e.g. 'NDVI') and SENSOR ('Sentinel-2', 'Landsat',
	'Modis') and, when several collections have them (e.g. MOD13Q1 and MYD13Q1), COLLECTION. Without COLLECTION
	the first entry of the catalog is used. Names are case insensitive."""
	entries = variables()
	if ref is not None:
		for e in entries:
			if e['ref'] == int(ref):
				return e
		raise KeyError(f'no entry with ref {ref} in the catalog of MonitorEO')
	low = lambda s: str(s).lower().replace(' ', '')
	matches = [e for e in entries if low(e['index']) == low(index)
	           and (sensor is None or low(e['sensor']) == low(sensor))
	           and (collection is None or e['col'] == collection)]
	if not matches:
		options = sorted({f"{e['index']} ({e['sensor']})" for e in entries
		                  if low(e['index']) == low(index) or sensor is None})
		raise KeyError(f"{index} with sensor {sensor}{f' and collection {collection}' if collection else ''} "
		               f"is not in the catalog of MonitorEO. Options: {', '.join(options) or 'see --list'}")
	return matches[0]


def alternatives(entry):
	"""Other entries with the same variable and sensor (other collections)."""
	return [e for e in variables() if e['index'] == entry['index'] and e['sensor'] == entry['sensor']
	        and e['ref'] != entry['ref']]


def table():
	"""Rows (ref, variable, sensor, collection, resolution, dates, unit, formula) to list the catalog."""
	return [(e['ref'], e['index'], e['sensor'], e['col'], e['m'], e['dr'], e['unit'], e['ft'])
	        for e in variables()]
