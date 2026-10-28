#### eHabitat+: input variables from Google Earth Engine with the catalog of MonitorEO (OBSNEV)
#### License: CC BY-SA 3.0

"""
Run from the segm_optim folder, with the study area in pa/ (as for setup.py):

    python3 get_vars_gee.py --list      # variables of the catalog of MonitorEO
    python3 get_vars_gee.py --dry-run   # grids, sizes and number of images, without downloading
    python3 get_vars_gee.py             # download the variables of gee_vars.yaml

It reads the recipe gee_vars.yaml (--recipe for another one): the variables, period, season, metrics and quality
masks, and computes each variable with Google Earth Engine on two grids:

  - the study area: the rectangle of the polygons of pa/ grown by PA_BUFFER (and 2 pixels more), with a fine pixel
    size (study: resolution), saved to ../inVars/<folder>, ready for setup.py (INVARS_FOLDER = '<folder>');
  - the ecoregions of RESOLVE Ecoregions 2017 (Ecoregions2017) that the study area touches (their union if it is in
    several), with a coarser pixel size (ecoregion: resolution), saved to ../inVars/<folder>/ecoregion, with the
    polygons of the ecoregions in ecoregions.gpkg. setup.py ignores this subfolder.

The work is done by the Python package of MonitorEO, monitoreo (https://github.com/monitoreo-obsnev/monitoreo; its
module monitoreo.recipe), so that eHabitat+ computes the variables as the MonitorEO tools do (QGIS, R, MCP). Every
run checks that the installed monitoreo is the last version of GitHub and installs it if not (--no-update skips
this, e.g. without internet; an editable copy, for development, is only compared with GitHub).

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
import importlib
import importlib.metadata
import json
import os
import re
import subprocess
import sys
import urllib.parse
import urllib.request

from tqdm import tqdm

import hft_state
import setup

MONITOREO_REPO = 'monitoreo-obsnev/monitoreo'
MONITOREO_BRANCH = 'main'
RECIPE = os.path.join(hft_state.HERE, 'gee_vars.yaml')
RECORD = 'gee_vars.json'
DRIVE_FOLDER = 'eHabitat_GEE'

fail = setup.fail


# The package monitoreo, kept at the last version of GitHub ------------------------------------------------------

def installed_monitoreo():
	"""(version, commit of GitHub or None, folder of an editable copy or None) of the installed package monitoreo,
	or None if it is not installed."""
	try:
		dist = importlib.metadata.distribution('monitoreo')
	except importlib.metadata.PackageNotFoundError:
		return None
	try:
		url = json.loads(dist.read_text('direct_url.json') or '{}')
	except ValueError:
		url = {}
	commit = url.get('vcs_info', {}).get('commit_id')
	m = re.search(r'/archive/([0-9a-f]{40})\.zip$', url.get('url', ''))
	if m:
		commit = m.group(1)
	folder = None
	if url.get('dir_info', {}).get('editable') and url.get('url', '').startswith('file:'):
		folder = urllib.request.url2pathname(urllib.parse.urlparse(url['url']).path)
	return dist.version, commit, folder


def latest_monitoreo():
	"""Commit of the last version of monitoreo in GitHub."""
	req = urllib.request.Request(f'https://api.github.com/repos/{MONITOREO_REPO}/commits/{MONITOREO_BRANCH}',
	                             headers={'Accept': 'application/vnd.github.sha', 'User-Agent': 'eHabitat'})
	with urllib.request.urlopen(req, timeout=15) as r:
		sha = r.read().decode().strip()
	if not re.fullmatch(r'[0-9a-f]{40}', sha):
		raise ValueError(f'unexpected answer of GitHub: {sha[:80]}')
	return sha


def git(folder, *args):
	"""Output of a git command in FOLDER, or None if it fails."""
	try:
		r = subprocess.run(['git', '-C', folder, *args], capture_output=True, text=True)
	except OSError:
		return None
	return r.stdout.strip() if r.returncode == 0 else None


def pip_install(*args):
	subprocess.run([sys.executable, '-m', 'pip', 'install', '--disable-pip-version-check', '--quiet', *args], check=True)


def update_monitoreo(check=True):
	"""Make sure that the package monitoreo is the last version of GitHub, installing it if it is not. An editable
	copy (development) is only compared with GitHub. CHECK=False only checks that it is installed."""
	installed = installed_monitoreo()
	if not check:
		if not installed:
			fail('the package monitoreo is not installed: run without --no-update to install it from GitHub')
		print(f'monitoreo {installed[0]} (not compared with GitHub: --no-update)')
		return
	try:
		latest = latest_monitoreo()
	except (OSError, ValueError) as e:
		if not installed:
			fail(f'the package monitoreo is not installed and GitHub cannot be reached to install it ({e})')
		print(f'WARNING: GitHub cannot be reached to check the version of monitoreo ({e}): using the installed one, '
		      f'{installed[0]}')
		return
	version, commit, folder = installed or (None, None, None)
	if folder:
		head = git(folder, 'rev-parse', 'HEAD')
		if head == latest:
			print(f'monitoreo {version}: development copy in {folder}, at the last version of GitHub ({latest[:7]})')
		elif head and git(folder, 'merge-base', '--is-ancestor', latest, 'HEAD') is not None:
			print(f'monitoreo {version}: development copy in {folder}, ahead of GitHub ({latest[:7]})')
		else:
			print(f'WARNING: the development copy of monitoreo in {folder} is not the last version of GitHub '
			      f'({latest[:7]}): update it (git pull) or install the package from GitHub')
		return
	if commit == latest:
		print(f'monitoreo {version} ({latest[:7]}): the last version of GitHub')
		return
	print(f"{'Updating' if installed else 'Installing'} monitoreo from GitHub ({MONITOREO_REPO}, {latest[:7]})...")
	url = f'monitoreo @ https://github.com/{MONITOREO_REPO}/archive/{latest}.zip'
	try:
		pip_install('--force-reinstall', '--no-deps', url)  # the version number may be the same as the installed one
		pip_install(url)  # the dependencies that are missing, if any
	except (OSError, subprocess.CalledProcessError) as e:
		if not installed:
			fail(f'cannot install monitoreo from GitHub: {e}')
		print(f'WARNING: cannot update monitoreo from GitHub ({e}): using the installed one, {version}')
		return
	importlib.invalidate_caches()
	print(f'monitoreo {installed_monitoreo()[0]} ({latest[:7]}) installed')


# Variables ----------------------------------------------------------------------------------------------------------

def print_catalog(catalog):
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
	parser.add_argument('--no-update', action='store_true', help='do not compare the package monitoreo with its last '
	                                                             'version in GitHub (e.g. without internet)')
	args = parser.parse_args()

	update_monitoreo(not args.no_update)
	try:
		from monitoreo import catalog, recipe, services
	except ImportError as e:
		fail(f'the package monitoreo is too old or broken ({e}): run without --no-update to install its last version')

	if args.list:
		print_catalog(catalog)
		return

	try:
		rec = recipe.read(args.recipe, extra={'folder'})
	except ValueError as e:
		fail(str(e))
	rec['drive_folder'] = rec.get('drive_folder') or DRIVE_FOLDER
	folder = str(args.folder or rec.get('folder') or 'gee').strip('/\\')
	cfg = hft_state.read_config()

	study_file = setup.find_study_area(args.study_area)
	if not study_file:
		fail(f'no vector file of the study area in {setup.PA_DIR}: copy it there or use --study-area')
	layer = setup.check_layer(study_file, args.layer)

	project = args.project or rec.get('project') or os.environ.get('EE_PROJECT')
	if not project:
		fail('give the Earth Engine project: project: in the recipe, --project or the variable EE_PROJECT')
	try:
		services.init(project)
		plan = recipe.plan(rec, study_file, layer, buffer=float(cfg.get('PA_BUFFER', 0)), only=args.only,
		                   extra={'folder'}, warn=print)
	except (ValueError, RuntimeError) as e:
		fail(str(e))

	outdir = os.path.join(hft_state.INVARS, folder)
	dirs = recipe.folders(outdir)
	summary = plan.summary()
	print(f'Study area: {study_file}')
	for name, g in plan.grids.items():
		print(f"{'Study area' if name == 'study' else 'Ecoregions'} grid: {g.crs}, {g.res:g} m, {g.width} x {g.height}: "
		      f'{size_text(g.pixels)} -> {hft_state.where(os.path.relpath(dirs[name], hft_state.INVARS))}')
	for p in plan.ecoregions:
		print(f"  ecoregion {p['ECO_ID']}: {p['ECO_NAME']} ({p['BIOME_NAME']})")
	print(f'{len(plan.specs)} variables:')
	for f in summary['files']:
		filters = f.get('pixel_filters')
		print(f"  {f['file']}" + (f"  (pixel filters: {', '.join(filters) or 'none'})" if filters is not None else ''))

	if args.dry_run:
		print('\nImages of each variable over the study area:')
		for f, n in recipe.n_images(plan).items():
			print(f"  {f}: {'fixed image' if n is None else n}")
		return

	# Downloads: directly, in tiles, or by export tasks to Google Drive for the grids with drive: true (or --drive)
	bars = {}

	def progress(file, done, total):
		bar = bars.get(file)
		if bar is None:
			bar = bars[file] = tqdm(total=total, desc=f'  {file}', unit='tile', leave=True)
		bar.update(done - bar.n)
		if done == total:
			bar.close()

	try:
		res = recipe.run(plan, outdir, RECORD, args.overwrite, args.drive, args.drive_dir, wait=not args.no_wait,
		                 tile=args.tile, workers=args.workers, project=project, log=print, progress=progress)
	except (ValueError, RuntimeError) as e:
		for bar in bars.values():
			bar.close()
		fail(str(e))
	if not res['done']:
		return
	print(f"\nDone: variables in {hft_state.where(folder)} (NoData -9999).")
	if args.use:
		hft_state.write_option('INVARS_FOLDER', repr(folder))
		print(f"config.py: INVARS_FOLDER = {folder!r}. Next: python3 run_hft.py")
	elif hft_state.invars_folder(cfg) != folder:
		print(f"To use them, set INVARS_FOLDER = {folder!r} in config.py (or run with --use), then python3 run_hft.py")


if __name__ == '__main__':
	main()
