#### Imports the variable catalog of the MonitorEO app (OBSNEV, Google Earth Engine) into catalog.json
#### License: CC BY-SA 3.0

"""
MonitorEO is a Google Earth Engine app written in JavaScript. Its catalog of
variables (one entry per variable and sensor, with the collection, formula,
scale, offset, unit, ...) is a JavaScript array in its translations modules.
This script reads the source of the app (public, as the JSON the app loads),
parses those arrays as data (the JavaScript is never run) and writes
../catalog.json, so that the Python code uses the same catalog as the app:

    python3 tools/import_catalog.py                 # from the published app
    python3 tools/import_catalog.py --modules FILE  # from a saved modules JSON

Run it again when the catalog of MonitorEO changes.
"""

import argparse
import datetime
import json
import os
import re
import sys
import urllib.request

APP_URL = 'https://ee-ugr-default-project.projects.earthengine.app/view/monitoreo-obsnev'
MODULES_URL = 'https://ee-ugr-default-project.projects.earthengine.app/javascript/monitoreo-obsnev-modules.json'
REPO = 'users/monitoreo-obsnev/monitoreo-obsnev'
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'catalog.json')

# Fields of the entries that are numbers, sometimes written as strings in the app (e.g. "scale": "10000")
NUMERIC = ('ref', 'm', 'scale', 'offset', 'scalet', 'offsett', 'min', 'max', 'minp', 'maxp', 'fe', 'aux',
           'clouds', 'water', 'snow', 'aerosols', 'pclouds', 'pwater', 'psnow', 'paerosols')


def module(modules, name):
	"""Source of the module NAME (e.g. 'translations/es.js') of the app."""
	for key, src in modules['dependencies'].items():
		if key.split(':', 1)[-1] == name:
			return src
	sys.exit(f'ERROR: module {name} not found in the source of MonitorEO')


def split_comment(line):
	"""(code, comment) of a line of JavaScript: the comment starts at // outside a string."""
	in_str = False
	i = 0
	while i < len(line):
		c = line[i]
		if in_str:
			if c == '\\':
				i += 1
			elif c == '"':
				in_str = False
		elif c == '"':
			in_str = True
		elif line.startswith('//', i):
			return line[:i], line[i + 2:].strip()
		i += 1
	return line, ''


def parse_data(src):
	"""Entries of the array 'var data = [...]' of a translations module, and the comments of their fields
	({ref: {field: comment}})."""
	start = src.index('[', src.index('var data'))
	end = src.index('];', start) + 1
	code, notes, ref = [], {}, None
	for line in src[start:end].splitlines():
		line, comment = split_comment(line)
		m = re.search(r'"ref"\s*:\s*(\d+)', line)
		if m:
			ref = int(m.group(1))
		key = re.match(r'\s*"(\w+)"\s*:', line)
		if comment and key and ref is not None:
			notes.setdefault(ref, {})[key.group(1)] = comment
		code.append(line)
	text = re.sub(r',(\s*[\]}])', r'\1', '\n'.join(code))  # trailing commas are not valid JSON
	return json.loads(text), notes


def number(v):
	if isinstance(v, str):
		try:
			f = float(v)
		except ValueError:
			return v
		return int(f) if f.is_integer() else f
	return v


def date_ranges(src):
	"""First and last date of each collection, from predefinedRanges of modules/sensorList."""
	ranges = {}
	for m in re.finditer(r'"([\w/\-]+)"\s*:\s*\{\s*start:\s*"([\d-]+)"(?:\s*,\s*end:\s*"([\d-]+)")?', src):
		ranges[m.group(1)] = {'start': m.group(2), 'end': m.group(3)}
	return ranges


def main():
	parser = argparse.ArgumentParser(description='Import the variable catalog of MonitorEO into catalog.json.')
	parser.add_argument('--modules', metavar='FILE', help=f'saved source of the app (default: download {MODULES_URL})')
	parser.add_argument('--out', default=OUT, help='output file (default: %(default)s)')
	args = parser.parse_args()

	if args.modules:
		with open(args.modules, encoding='utf-8') as f:
			modules = json.load(f)
	else:
		with urllib.request.urlopen(MODULES_URL, timeout=60) as r:
			modules = json.loads(r.read().decode('utf-8'))

	es, notes = parse_data(module(modules, 'translations/es.js'))
	en, _ = parse_data(module(modules, 'translations/en.js'))
	large, _ = parse_data(module(modules, 'translations/es_largeareas.js'))
	en = {e['ref']: e for e in en}
	large = {e['ref'] for e in large}

	entries = []
	for e in es:
		e = {k: number(v) if k in NUMERIC else v for k, v in e.items()}
		ref = e['ref']
		e['vname_en'] = en.get(ref, {}).get('vname', e['vname'])
		e['type_en'] = en.get(ref, {}).get('type', e['type'])
		e['large_areas'] = ref in large
		if ref in notes:
			e['notes'] = notes[ref]
		entries.append(e)

	catalog = {
		'source': {
			'app': APP_URL,
			'script': REPO,
			'modules': MODULES_URL if not args.modules else os.path.basename(args.modules),
			'imported': datetime.date.today().isoformat(),
		},
		'fields': {
			'ref': 'ID of the entry', 'type': 'group of variables', 'index': 'variable', 'vname': 'name',
			'sensor': 'Sentinel-2, Landsat or Modis', 'col': 'Earth Engine collection (LANDSAT = Landsat 5, 7, 8 and 9)',
			'm': 'native pixel size (m)', 'tr': 'revisit time', 'dr': 'dates', 'unit': 'unit',
			'ft': "formula: a band, an expression of the bands, or N(b1,b2) for the normalized difference",
			'scale': 'scale factor of the bands', 'offset': 'offset of the bands',
			'clouds_bandQA': 'QA band for clouds', 'water_bandQA': 'QA band for water', 'snow_bandQA': 'QA band for snow',
			'aerosols_bandQA': 'QA band for aerosols', 'minp': 'minimum of the colour ramp', 'maxp': 'maximum of the colour ramp',
			'pm': 'colour ramp', 'large_areas': 'offered for large areas', 'notes': 'comments of the fields in the app',
		},
		'collections': date_ranges(module(modules, 'modules/sensorList')),
		'variables': entries,
	}
	with open(args.out, 'w', encoding='utf-8') as f:
		json.dump(catalog, f, ensure_ascii=False, indent=1)
		f.write('\n')
	print(f'{len(entries)} variables written to {args.out}')


if __name__ == '__main__':
	main()
