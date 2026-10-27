#### eHabitat+: input variables from Google Earth Engine with the catalog of MonitorEO (OBSNEV)
#### License: CC BY-SA 3.0

"""
Run from the segm_optim folder, with the study area in pa/ (as for setup.py):

    python3 get_vars_gee.py --list      # variables of the catalog of MonitorEO
    python3 get_vars_gee.py --dry-run   # grids, sizes and number of images, without downloading
    python3 get_vars_gee.py             # download the variables of gee_vars.yaml

It reads the recipe gee_vars.yaml (--recipe for another one): the variables, period, season, metrics and quality
masks, and computes each variable with Google Earth Engine (the Python core of MonitorEO, ../monitoreo) on two grids:

  - the study area: the rectangle of the polygons of pa/ grown by PA_BUFFER (and 2 pixels more), with a fine pixel
    size (study: resolution), saved to ../inVars/<folder>, ready for setup.py (INVARS_FOLDER = '<folder>');
  - the ecoregions of RESOLVE Ecoregions 2017 (Ecoregions2017) that the study area touches (their union if it is in
    several), with a coarser pixel size (ecoregion: resolution), saved to ../inVars/<folder>/ecoregion, with the
    polygons of the ecoregions in ecoregions.gpkg. setup.py ignores this subfolder.

Both grids have the same CRS (unless ecoregion: crs says otherwise), edges at multiples of the pixel size, NoData
-9999 and the files are numbered in the order of the recipe (1_, 2_, ... or 01_, 02_, ...). The recipe of every file
is recorded in gee_vars.json of its folder: running it again only downloads the files whose recipe or grid changed
(--overwrite downloads all of them again).

A grid is downloaded directly, in tiles computed by Earth Engine on demand, or, with drive: true in the recipe (or
--drive), by export tasks to the folder drive_folder of Google Drive, without limits of size or time per tile (better
for big grids). The exported files are imported from the local copy of Google Drive for desktop, or from the folder
given with --drive-dir where they were downloaded. Running it again while the tasks run (or with --no-wait) does not
start them again: it waits for them and imports their files.

It needs an Earth Engine project (project: in the recipe, --project or the variable EE_PROJECT) and the credentials of
'earthengine authenticate'.
"""

import argparse
import datetime
import json
import os
import re
import sys
import time

import yaml
from tqdm import tqdm

import hft_state
import setup

sys.path.insert(0, os.path.normpath(os.path.join(hft_state.HERE, '..')))
from monitoreo import areas, catalog, drive, sensors, variable  # noqa: E402
from monitoreo.grid import Grid, download, NODATA  # noqa: E402

RECIPE = os.path.join(hft_state.HERE, 'gee_vars.yaml')
RECORD = 'gee_vars.json'
ECO_FOLDER = 'ecoregion'
MARGIN_PIXELS = 2
RECIPE_KEYS = {'project', 'folder', 'crs', 'study', 'ecoregion', 'variables', 'drive_folder', 'drive_dir'} | (variable.KEYS - {
	'name', 'index', 'sensor', 'collection', 'ref', 'image', 'band', 'scale', 'offset', 'terrain', 'metric', 'metrics'})

fail = setup.fail


def read_recipe(path):
	try:
		with open(path, encoding='utf-8') as f:
			recipe = yaml.safe_load(f) or {}
	except OSError as e:
		fail(f'cannot read the recipe {path}: {e}')
	except yaml.YAMLError as e:
		fail(f'{path} is not valid YAML: {e}')
	unknown = set(recipe) - RECIPE_KEYS
	if unknown:
		fail(f"unknown keys in {path}: {', '.join(sorted(unknown))}")
	if not recipe.get('variables'):
		fail(f'{path} has no variables')
	return recipe


def expand(recipe):
	"""Specs of the files to download, in order, with the file name (numbered) in 'file'."""
	defaults = {k: v for k, v in recipe.items() if k in variable.KEYS}
	out = []
	for i, var in enumerate(recipe['variables'], 1):
		if not isinstance(var, dict):
			fail(f'variable {i} of the recipe is not a dictionary of options')
		try:
			new = variable.specs(var, defaults)
			for s in new:
				if variable.kind(s) == 'catalog':
					entry = variable.entry_of(s)
					for f in sensors.conflicts(entry, s['quality']):
						print(f"WARNING: the {f} filter removes what {entry['index']} measures: "
						      f"give this variable quality: {{{f}: false}}")
		except (ValueError, KeyError) as e:
			fail(f'variable {i} of the recipe: ' + str(e).strip('"'))
		out += new
	width = len(str(len(out))) if len(out) >= 10 else 1
	names = set()
	for i, s in enumerate(out, 1):
		if s['name'] in names:
			fail(f"two variables give the file name {s['name']}: give them different names (name:)")
		names.add(s['name'])
		s['file'] = f"{i:0{width}d}_{s['name']}.tif"
	return out


def print_catalog():
	rows = catalog.table()
	print(f"Catalog of MonitorEO ({catalog.load()['source']['imported']}, {len(rows)} variables). In the recipe: "
	      'index and sensor (and collection when there are several).\n')
	print(f"{'ref':>3}  {'index':<7}{'sensor':<11}{'collection':<35}{'m':>6}  {'unit':<18}formula")
	for ref, index, sensor, col, m, _, unit, ft in rows:
		print(f'{ref:>3}  {index:<7}{sensor:<11}{col:<35}{m:>6}  {unit:<18}{ft}')
	print('\nAlso any other collection (collection: and band:) or image (image: and band:, with terrain: slope) '
	      'of Earth Engine.')


def size_text(pixels):
	mb = pixels * 4 / 2 ** 20
	return f'{pixels:,} pixels, {mb:,.0f} MB per variable' if mb >= 1 else f'{pixels:,} pixels, {mb:.2f} MB per variable'


def load_record(folder):
	try:
		with open(os.path.join(folder, RECORD), encoding='utf-8') as f:
			return json.load(f)
	except (OSError, ValueError):
		return {}


def save_record(folder, rec):
	path = os.path.join(folder, RECORD)
	with open(path + '.tmp', 'w', encoding='utf-8') as f:
		json.dump(rec, f, indent=1, ensure_ascii=False, default=str)
	os.replace(path + '.tmp', path)


def spec_key(spec):
	"""What defines a file: its spec (without the file name) as JSON."""
	return json.loads(json.dumps({k: v for k, v in spec.items() if k != 'file'}, sort_keys=True, default=str))


def wait_and_import(pending, recs, grids, outdirs, local, no_wait):
	"""Wait for the export tasks of PENDING ([(grid name, file)]) and import their files from the local folder LOCAL
	(Google Drive for desktop, or where they were downloaded). True when all of them were imported."""
	if no_wait:
		print('\nThe export tasks go on in Earth Engine: run the same command again to import their files when done.')
		return False
	where_files = local or 'your Google Drive'
	print(f'\nWaiting for the export tasks (also in the Tasks tab of the Code Editor); files read from {where_files}.\n'
	      'Ctrl+C stops waiting, not the tasks: run the same command again to go on.')
	shown, sizes, left = {}, {}, []
	pending = list(pending)
	while pending:
		for name, f in list(pending):
			rec, outdir = recs[name], outdirs[name]
			t = rec['tasks'][f]
			state, progress, error = drive.status(t['operation'])
			msg = None
			if state in drive.FAILED:
				msg = f'{state.lower()}: {error}'
				pending.remove((name, f))
			elif state != drive.DONE:
				msg = f'{state.lower()} {progress:.0%}' if progress else state.lower()
			elif not local or not os.path.isdir(local):
				msg = f"exported to Google Drive/{t['drive_folder']}/{t['drive_name']}*.tif"
				left.append(f)
				pending.remove((name, f))
			else:
				found = drive.find(local, t['drive_name'])
				now = [os.path.getsize(p) for p in found]
				if len(found) < drive.pieces(grids[name]) or sizes.get(f) != now:
					sizes[f] = now
					msg = 'exported, waiting for Google Drive to copy it to this computer'
				else:
					path = os.path.join(outdir, f)
					drive.mosaic(found, grids[name], path, metadata={
						'source': 'Google Earth Engine', 'monitoreo': json.dumps(t['source'], default=str)})
					rec['files'][f] = {k: t[k] for k in ('spec', 'grid', 'source', 'images')}
					rec['files'][f].update({'downloaded': datetime.datetime.now().isoformat(timespec='seconds'),
					                        'via': f"Google Drive/{t['drive_folder']}"})
					del rec['tasks'][f]
					save_record(outdir, rec)
					msg = f"imported ({', '.join(os.path.basename(p) for p in found)})"
					pending.remove((name, f))
			if msg and shown.get(f) != msg:
				print(f'  {f}: {msg}')
				shown[f] = msg
		if pending:
			try:
				time.sleep(30)
			except KeyboardInterrupt:
				print('\nStopped waiting. The export tasks go on in Earth Engine: run the same command again to import '
				      'their files.')
				return False
	if left:
		print('\nDownload the folder of the exported files from your Google Drive and run the same command with '
		      '--drive-dir <that folder> to import them.')
	return not left and all(not r['tasks'] for r in recs.values())


def main():
	sys.stdout.reconfigure(line_buffering=True)
	parser = argparse.ArgumentParser(description='Download the input variables of gee_vars.yaml from Google Earth '
	                                             'Engine (catalog of MonitorEO) for the study area and its ecoregions.')
	parser.add_argument('--recipe', default=RECIPE, help='recipe of the variables (default: gee_vars.yaml)')
	parser.add_argument('--list', action='store_true', help='list the variables of the catalog of MonitorEO and exit')
	parser.add_argument('--dry-run', action='store_true', help='show the grids, sizes and number of images, without downloading')
	parser.add_argument('--project', help='Earth Engine project (default: project of the recipe, or EE_PROJECT)')
	parser.add_argument('--folder', help='subfolder of ../inVars (default: folder of the recipe)')
	parser.add_argument('--study-area', metavar='FILE', help='vector file of the study area (default: the only one in pa/)')
	parser.add_argument('--layer', help='layer of the study area file, if it has several')
	parser.add_argument('--only', choices=('study', 'ecoregion'), help='download only the variables of one grid')
	parser.add_argument('--overwrite', action='store_true', help='download again the files that are up to date')
	parser.add_argument('--tile', type=int, default=512, help='size of the tiles computed by Earth Engine (pixels, default 512)')
	parser.add_argument('--workers', type=int, default=8, help='tiles computed at the same time (default 8)')
	parser.add_argument('--drive', action='store_true', help='export every grid to Google Drive (default: the grids '
	                                                         'with drive: true in the recipe)')
	parser.add_argument('--drive-dir', metavar='FOLDER', help='local folder with the files exported to Google Drive '
	                    '(default: drive_dir of the recipe, or the folder of Google Drive for desktop)')
	parser.add_argument('--no-wait', action='store_true', help='start the export tasks and exit; run again to import')
	parser.add_argument('--use', action='store_true', help="set INVARS_FOLDER to the folder in config.py when done")
	args = parser.parse_args()

	if args.list:
		print_catalog()
		return

	recipe = read_recipe(args.recipe)
	specs = expand(recipe)
	folder = str(args.folder or recipe.get('folder') or 'gee').strip('/\\')
	study_opts = recipe.get('study') or {}
	eco_opts = recipe.get('ecoregion') or {}
	if not study_opts.get('resolution'):
		fail('give the pixel size of the study area (study: resolution) in the recipe')
	cfg = hft_state.read_config()
	buffer = study_opts.get('buffer')
	buffer = float(cfg.get('PA_BUFFER', 0) if buffer is None else buffer)

	study_file = setup.find_study_area(args.study_area)
	if not study_file:
		fail(f'no vector file of the study area in {setup.PA_DIR}: copy it there or use --study-area')
	layer = setup.check_layer(study_file, args.layer)

	project = args.project or recipe.get('project') or os.environ.get('EE_PROJECT')
	if not project:
		fail('give the Earth Engine project: project: in the recipe, --project or the variable EE_PROJECT')
	import ee
	try:
		ee.Initialize(project=project)
	except Exception as e:
		fail(f'cannot connect to Earth Engine with the project {project}: {e}\n'
		     'Run earthengine authenticate (once) and check that the project is registered in Earth Engine.')

	# Grids
	geom, srs = areas.read_vector(study_file, layer)
	crs = recipe.get('crs') or 'auto'
	crs = areas.utm_crs(geom, srs) if str(crs).lower() == 'auto' else str(crs)
	res = float(study_opts['resolution'])
	grids = {}
	grids['study'] = Grid.around(crs, res, *areas.bounds(geom, srs, crs), margin=buffer + MARGIN_PIXELS * res)
	eco = None
	if eco_opts.get('resolution') and args.only != 'study':
		eco_res = float(eco_opts['resolution'])
		eco_crs = str(eco_opts.get('crs') or crs)
		eco = areas.ecoregions(areas.to_ee(geom, srs, simplify=res if srs.IsProjected() else None))
		info = areas.ecoregion_info(eco)
		if not info:
			fail('the study area is not in any ecoregion of RESOLVE Ecoregions 2017')
		grids['ecoregion'] = Grid.around(eco_crs, eco_res, *areas.ecoregion_bounds(eco, eco_crs, eco_res / 2))
	if args.only == 'ecoregion':
		if 'ecoregion' not in grids:
			fail('the recipe has no ecoregion: resolution')
		del grids['study']

	outdirs = {'study': os.path.join(hft_state.INVARS, folder),
	           'ecoregion': os.path.join(hft_state.INVARS, folder, ECO_FOLDER)}
	print(f'Study area: {study_file}')
	for name, g in grids.items():
		print(f"{'Study area' if name == 'study' else 'Ecoregions'} grid: {g.crs}, {g.res:g} m, {g.width} x {g.height}: "
		      f'{size_text(g.pixels)} -> {hft_state.where(os.path.relpath(outdirs[name], hft_state.INVARS))}')
	if eco is not None:
		for e in info:
			p = e['properties']
			print(f"  ecoregion {p['ECO_ID']}: {p['ECO_NAME']} ({p['BIOME_NAME']})")
	print(f'{len(specs)} variables:')
	for s in specs:
		extra = ''
		if variable.kind(s) == 'catalog':
			on = [f for f, v in sensors.filters(variable.entry_of(s), s['quality']).items() if v]
			extra = f"  (pixel filters: {', '.join(on) or 'none'})"
		print(f"  {s['file']}{extra}")

	if args.dry_run:
		print('\nImages of each variable over the study area:')
		for s in specs:
			try:
				n = variable.n_images(s, grids.get('study', next(iter(grids.values()))).region())
			except Exception as e:
				n = f'ERROR: {e}'
			print(f"  {s['file']}: {'fixed image' if n is None else n}")
		return

	# Downloads: directly, in tiles, or by export tasks to Google Drive for the grids with drive: true (or --drive)
	via_drive = {name: args.drive or bool((study_opts if name == 'study' else eco_opts).get('drive')) for name in grids}
	drive_folder = re.sub(r'[^A-Za-z0-9_-]', '_', str(recipe.get('drive_folder') or 'eHabitat_GEE'))
	local = args.drive_dir or recipe.get('drive_dir')
	if any(via_drive.values()) and not local:
		found = drive.local_drives()
		if len(found) > 1:
			fail(f"several Google Drive folders on this computer ({', '.join(found)}): give the folder of the exported "
			     'files with --drive-dir')
		local = os.path.join(found[0], drive_folder) if found else None
	recs, pending = {}, []
	for name, g in grids.items():
		outdir = outdirs[name]
		os.makedirs(outdir, exist_ok=True)
		rec = recs[name] = load_record(outdir)
		rec['grid'] = g.describe()
		rec['study_area'] = os.path.basename(study_file)
		rec['catalog'] = catalog.load()['source']
		if name == 'ecoregion':
			rec['ecoregions'] = [e['properties'] for e in info]
			areas.save_ecoregions(eco, os.path.join(outdir, 'ecoregions.gpkg'), g.res / 2)
		files = rec.setdefault('files', {})
		tasks = rec.setdefault('tasks', {})
		region = g.region()
		print(f"\n{'Study area' if name == 'study' else 'Ecoregions'} ({g.res:g} m"
		      f"{', export to Google Drive/' + drive_folder if via_drive[name] else ''}):")
		for s in specs:
			path = os.path.join(outdir, s['file'])
			key = spec_key(s)
			old = files.get(s['file'], {})
			if (not args.overwrite and os.path.exists(path) and old.get('spec') == key
			        and old.get('grid') == rec['grid']):
				print(f"  {s['file']}: up to date")
				continue
			task = tasks.get(s['file'])
			if (via_drive[name] and not args.overwrite and task and task['spec'] == key and task['grid'] == rec['grid']
			        and drive.status(task['operation'])[0] not in drive.FAILED):
				print(f"  {s['file']}: export task already started")
				pending.append((name, s['file']))
				continue
			img = variable.image(s, region, g.res)
			if name == 'ecoregion':
				img = img.clipToCollection(eco)
			desc = variable.describe(s)
			n = variable.n_images(s, region)
			if n == 0:
				print(f"  WARNING: {s['file']}: no images in the period; the file will be all NoData")
			if via_drive[name]:
				drive_name = f"{folder}_{name}_{os.path.splitext(s['file'])[0]}"
				tasks[s['file']] = {'spec': key, 'grid': rec['grid'], 'source': desc, 'images': n,
				                    'operation': drive.export(img, g, drive_name, drive_folder, project),
				                    'drive_folder': drive_folder, 'drive_name': drive_name,
				                    'started': datetime.datetime.now().isoformat(timespec='seconds')}
				save_record(outdir, rec)
				print(f"  {s['file']}: export task started (Google Drive/{drive_folder}/{drive_name}.tif)")
				pending.append((name, s['file']))
				continue
			bar = tqdm(total=len(g.tiles(args.tile)), desc=f"  {s['file']}", unit='tile', leave=True)

			def progress(done, total, bar=bar):
				bar.update(done - bar.n)

			try:
				download(img, g, path, tile=args.tile, workers=args.workers, progress=progress,
				         metadata={'source': 'Google Earth Engine', 'monitoreo': json.dumps(desc, default=str)})
			except Exception as e:
				bar.close()
				save_record(outdir, rec)
				fail(f"{s['file']}: {e}")
			bar.close()
			files[s['file']] = {'spec': key, 'grid': rec['grid'], 'source': desc, 'images': n,
			                    'downloaded': datetime.datetime.now().isoformat(timespec='seconds')}
			save_record(outdir, rec)
		# files no longer in the recipe stay on disk; say so, since setup.py would import them
		stale = [f for f in files if f not in {s['file'] for s in specs}]
		for f in stale:
			if os.path.exists(os.path.join(outdir, f)):
				print(f'  WARNING: {f} is not in the recipe any more: remove it if it should not be used')

	if pending and not wait_and_import(pending, recs, grids, outdirs, local, args.no_wait):
		return
	print(f"\nDone: variables in {hft_state.where(folder)} (NoData {NODATA}).")
	if args.use:
		hft_state.write_option('INVARS_FOLDER', repr(folder))
		print(f"config.py: INVARS_FOLDER = {folder!r}. Next: python3 run_hft.py")
	elif hft_state.invars_folder(cfg) != folder:
		print(f"To use them, set INVARS_FOLDER = {folder!r} in config.py (or run with --use), then python3 run_hft.py")


if __name__ == '__main__':
	main()
