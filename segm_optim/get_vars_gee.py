#### eHabitat+: input variables from Google Earth Engine with the catalog of MonitorEO (OBSNEV)
#### License: CC BY-SA 3.0

"""
Run from the segm_optim folder, with the study area in pa/ (as for setup.py):

    python3 get_vars_gee.py --list        # variables of the catalog of MonitorEO
    python3 get_vars_gee.py --dry-run     # grids, sizes and number of images, without downloading
    python3 get_vars_gee.py               # download the variables of gee_vars.yaml
    python3 get_vars_gee.py --area FILE   # the same with a comparison area instead of the ecoregions

It reads the recipe gee_vars.yaml (--recipe for another one): the variables, period, season, metrics and quality
masks, and computes each variable with Google Earth Engine on two grids:

  - the study area: the rectangle of the polygons of pa/ grown by PA_BUFFER (and 2 pixels more), with a fine pixel
    size (study: resolution), saved to ../inVars/<folder>, ready for setup.py (INVARS_FOLDER = '<folder>');
  - the ecoregion grid, with a coarser pixel size (ecoregion: resolution), saved to ../inVars/<folder>/ecoregion with
    its polygons in ecoregions.gpkg (setup.py ignores this subfolder). Its area is
      - by default, the ecoregions of RESOLVE Ecoregions 2017 (Ecoregions2017) that the study area touches (their
        union if it is in several);
      - or a comparison area: the polygons of a vector file of your own (comparison_area: and comparison_area_layer:
        in the recipe, or --area FILE and --area-layer LAYER), inside or outside the study area, e.g. another
        mountain range to compare the HFTs with. Each polygon (feature) of the file is an area of ecoregions.gpkg,
        with ECO_ID -1, -2, ..., its name (field name, NAME, nombre or Name, else the name of the file) in ECO_NAME
        and 'Comparison area' in BIOME_NAME; ecoregion/gee_vars.json records the file in comparison_area.
    The pixels outside the polygons are NoData.

When the polygons of the ecoregion grid do not contain the whole rectangle of the study grid (always with a
comparison area outside the study area), the study area is also downloaded at the pixel and in the CRS of the
ecoregion grid, to ../inVars/<folder>/ecoregion/study_area: the rectangle of the study grid with edges at multiples of
the ecoregion pixel, not clipped, the same variables, downloaded directly (it is small). The similarity of the HFTs
(similarity.py) takes the statistics of each HFT from it, so that they come from the same product and pixel as the
area they are compared with.

--only study or --only ecoregion downloads one of the grids (the ecoregion grid with its study_area when needed). A
download with --only ecoregion into a new folder is only a comparison area for the similarity (ECOSIM_AREA =
'<folder>/ecoregion' in config.py), not a set of variables for the HFT model.

The work is done by the Python package of MonitorEO, monitoreo (https://github.com/monitoreo-obsnev/monitoreo; its
module monitoreo.recipe), so that eHabitat+ computes the variables as the MonitorEO tools do (QGIS, R, MCP). Every
run checks that the installed monitoreo is the last version of GitHub and installs it if not (--no-update skips
this, e.g. without internet; an editable copy, for development, is only compared with GitHub).

All the grids have the same CRS (unless ecoregion: crs says otherwise), edges at multiples of the pixel size, NoData
-9999 and the files are numbered in the order of the recipe (1_, 2_, ... or 01_, 02_, ...). The recipe of every file
is recorded in gee_vars.json of its folder: running it again only downloads the files whose recipe or grid changed
(--overwrite downloads all of them again), and all those of the ecoregion grid when its area changed (another
comparison area, or the ecoregions instead of a comparison area or the reverse). Until all the files of the new area
are there, ecoregion/ stays that of the old area (comparison_area of its gee_vars.json, ecoregions.gpkg) and records
the new one in comparison_area_pending: running the same command again goes on with the download, and going back to
the old area keeps its files that were not replaced yet.

A grid is downloaded directly, in tiles computed by Earth Engine on demand, or, with drive: true in the recipe (or
--drive), by export tasks to the folder drive_folder of Google Drive, without limits of size or time per tile (better
for big grids). The exported files are imported from the local copy of Google Drive for desktop, or from the folder
given with --drive-dir where they were downloaded. Running it again while the tasks run (or with --no-wait) does not
start them again: it waits for them and imports their files.

With an Earth Engine project (project: in the recipe, --project or the variable EE_PROJECT) it computes with the
credentials of 'earthengine authenticate'; without one, with the shared account of MonitorEO (the research EECUs of
each month of each user; it cannot export to Google Drive, so its grids need drive: false).
"""

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
import re
import shutil
import signal
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
EXTRA_KEYS = {'folder', 'comparison_area', 'comparison_area_layer'}  # keys of the recipe that monitoreo does not use
STUDY_ECO = 'study_area'  # subfolder of ecoregion/: the study grid at the pixel of the ecoregion grid
COMPARISON = 'Comparison area'  # BIOME_NAME of the polygons of a comparison area in ecoregions.gpkg
NAME_FIELDS = ('name', 'NAME', 'nombre', 'Name')  # fields with the name of each polygon of a comparison area
PENDING = 'comparison_area_pending'  # key of ecoregion/gee_vars.json: the area being downloaded (see start_area)
PREVIOUS = 'previous_area'  # key of ecoregion/gee_vars.json meanwhile: the files of the published area
POLYGONS = 'ecoregions.gpkg'  # the polygons of the ecoregion grid (written by monitoreo)
POLYGONS_BEFORE = '.ecoregions_before.gpkg'  # meanwhile, a copy of those of the published area
POLYGONS_PENDING = '.ecoregions_pending.gpkg'  # meanwhile, those of the area being downloaded

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
	imported = (catalog.load().get('source') or {}).get('imported')  # older versions of monitoreo
	got = installed_monitoreo()
	source = f'imported {imported}' if imported else f'monitoreo {got[0]}' if got else 'monitoreo'
	print(f'Catalog of MonitorEO ({source}, {len(rows)} variables). In the recipe: index and sensor (and collection '
	      'when there are several).\n')
	print(f"{'ref':>3}  {'index':<9}{'sensor':<11}{'collection':<35}{'m':>6}  {'unit':<18}formula")
	for ref, index, sensor, col, m, _, unit, ft in rows:
		print(f"{ref:>3}  {index:<9}{sensor:<11}{str(col or ''):<35}{'' if m is None else m:>6}  {str(unit or ''):<18}"
		      f"{ft or ''}")
	print('\nAlso any other collection (collection: and band:) or image (image: and band:, with terrain: slope) '
	      'of Earth Engine.')


def size_text(pixels):
	mb = pixels * 4 / 2 ** 20
	return f'{pixels:,} pixels, {mb:,.0f} MB per variable' if mb >= 1 else f'{pixels:,} pixels, {mb:.2f} MB per variable'


def where(path):
	"""PATH as shown to the user: inVars/<path relative to inVars>."""
	return hft_state.where(os.path.relpath(path, hft_state.INVARS).replace(os.sep, '/'))


def grid_text(g):
	return f'{g.crs}, {g.res:g} m, {g.width} x {g.height}: {size_text(g.pixels)}'


def upgrade_record(folder, specs):
	"""Files recorded in FOLDER by an older monitoreo whose specs lack keys that the installed one adds with the value
	None (an option that is off, e.g. where) were made in the same way: record them with the new specs, so that
	recipe.run does not download them again. Returns the number of files."""
	from monitoreo import recipe
	rec = recipe.load_record(folder, RECORD)
	n = 0
	for s in specs:
		new = recipe.spec_key(s)
		for entry in ((rec.get('files') or {}).get(s['file']), (rec.get('tasks') or {}).get(s['file'])):
			old = (entry or {}).get('spec')
			same = isinstance(old, dict) and {k: v for k, v in new.items() if k in old or v is not None} == old
			if same and old != new:
				entry['spec'] = new
				n += 1
	if n:
		recipe.save_record(folder, rec, RECORD)
	return n


# Comparison area ----------------------------------------------------------------------------------------------------

def _gdal():
	from osgeo import ogr, osr
	ogr.UseExceptions()
	osr.UseExceptions()
	return ogr, osr


def _polygons(g):
	"""MultiPolygon of the polygons of the OGR geometry G (2D, linear, valid), or None if it has none."""
	ogr, _ = _gdal()
	g = g.Clone()
	g.FlattenTo2D()
	if g.HasCurveGeometry():
		g = g.GetLinearGeometry()
	if not g.IsValid():
		g = g.MakeValid()
	t = ogr.GT_Flatten(g.GetGeometryType())
	if t in (ogr.wkbPolygon, ogr.wkbMultiPolygon):
		out = ogr.ForceToMultiPolygon(g)
	elif t == ogr.wkbGeometryCollection:
		out = ogr.Geometry(ogr.wkbMultiPolygon)
		for i in range(g.GetGeometryCount()):
			part = g.GetGeometryRef(i)
			if ogr.GT_Flatten(part.GetGeometryType()) in (ogr.wkbPolygon, ogr.wkbMultiPolygon):
				part = ogr.ForceToMultiPolygon(part)
				for j in range(part.GetGeometryCount()):
					out.AddGeometry(part.GetGeometryRef(j))
	else:
		return None
	return None if out.IsEmpty() else out


def area_source(args, rec):
	"""(absolute path, layer or None) of the comparison area: --area (relative to the current folder; '' = none) or
	comparison_area of the recipe (relative to segm_optim), with --area-layer or comparison_area_layer; None for the
	ecoregions of RESOLVE."""
	if args.area is not None:
		path, layer, base = args.area, args.area_layer, os.getcwd()
	else:
		path, base = rec.get('comparison_area'), hft_state.HERE
		layer = args.area_layer if args.area_layer is not None else rec.get('comparison_area_layer')
	path = str(path or '').strip()
	layer = str(layer).strip() or None if layer is not None else None
	if not path:
		if layer and args.area_layer:
			fail('--area-layer needs a comparison area (--area FILE or comparison_area: in the recipe)')
		return None
	return os.path.normpath(os.path.join(base, os.path.expanduser(path))), layer


def read_area(path, layer):
	"""The comparison area of the vector file PATH (its LAYER): {'file', 'layer', 'name', 'srs', 'geoms' (MultiPolygons
	of its features, in SRS), 'props' (ECO_ID, ECO_NAME, BIOME_NAME, REALM of each), 'union', 'sig' (of the
	geometries: it changes when they change, not when the file is copied or only its names change)}."""
	ogr, osr = _gdal()
	from monitoreo import areas
	if not os.path.exists(path):
		fail(f'comparison area not found: {path}')
	try:
		ds = ogr.Open(path)
	except RuntimeError as e:
		fail(f'cannot open the comparison area {path}: {e}')
	if ds is None:
		fail(f'cannot open the comparison area {path}')
	layers = [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]
	if layer and layer not in layers:
		fail(f"layer {layer} not found in {path}; its layers are: {', '.join(layers)}")
	if not layer and len(layers) > 1:
		fail(f"the comparison area {path} has several layers ({', '.join(layers)}): choose one with --area-layer "
		     '(comparison_area_layer: in the recipe)')
	lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
	srs = lyr.GetSpatialRef()
	if srs is None:
		fail(f'the comparison area {path} has no CRS')
	srs = srs.Clone()
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	defn = lyr.GetLayerDefn()
	fields = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
	name_field = next((f for f in NAME_FIELDS if f in fields), None)
	stem = os.path.splitext(os.path.basename(path))[0]
	geoms, props, skipped = [], [], 0
	for feat in lyr:
		g = feat.GetGeometryRef()
		g = _polygons(g) if g is not None and not g.IsEmpty() else None
		if g is None:
			skipped += 1
			continue
		name = feat.GetField(name_field) if name_field else None
		name = str(name).strip() if name is not None and str(name).strip() else stem
		geoms.append(g)
		props.append({'ECO_ID': -len(geoms), 'ECO_NAME': name, 'BIOME_NAME': COMPARISON, 'REALM': None})
	ds = None
	if not geoms:
		fail(f'the comparison area {path} has no polygons')
	if skipped:
		print(f'WARNING: {skipped} feature(s) of {os.path.basename(path)} without polygons are ignored')
	union = ogr.Geometry(ogr.wkbMultiPolygon)
	for g in geoms:
		union = union.Union(g)
	wgs84 = areas.wgs84()
	sig = hashlib.sha1('\n'.join(sorted(areas.transform(g, srs, wgs84).ExportToJson(['COORDINATE_PRECISION=7'])
	                                    for g in geoms)).encode()).hexdigest()[:16]
	name = props[0]['ECO_NAME'] if len(props) == 1 else stem + (f' ({layer})' if layer else '')
	return {'file': path, 'layer': layer, 'name': name, 'srs': srs, 'geoms': geoms, 'props': props, 'union': union,
	        'sig': sig}


def area_key(area):
	"""The key comparison_area of the record of the ecoregion grid for AREA (None: the ecoregions of RESOLVE)."""
	if area is None:
		return None
	return {'file': area['file'], 'layer': area['layer'], 'name': area['name'], 'features': len(area['props']),
	        'sig': area['sig']}


def reproject(geom, srs, crs, step):
	"""The OGR geometry GEOM (in SRS) in CRS, its segments split to STEP metres first so that they keep their shape."""
	from monitoreo import areas
	g = geom.Clone()
	g.Segmentize(step / (111320.0 if srs.IsGeographic() else srs.GetLinearUnits() or 1.0))
	return areas.transform(g, srs, areas.srs_of(crs))


def envelope(geom):
	xmin, xmax, ymin, ymax = geom.GetEnvelope()
	return xmin, ymin, xmax, ymax


def rectangle(grid):
	"""OGR polygon of the extent of GRID (in its CRS)."""
	ogr, _ = _gdal()
	x0, y0, x1, y1 = grid.extent
	return ogr.CreateGeometryFromWkt(f'POLYGON(({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))')


def use_area(plan, area, eco_crs, eco_res):
	"""Make the ecoregion grid of PLAN that of the comparison AREA: its polygons (simplified by half a pixel, as
	monitoreo saves the ecoregions), their properties and the grid around them (no margin, as for the ecoregions)."""
	import ee
	from monitoreo import areas
	from monitoreo.grid import Grid
	eco_srs = areas.srs_of(eco_crs)
	feats = []
	for g, p in zip(area['geoms'], area['props']):
		g = areas.transform(g, area['srs'], eco_srs)
		feats.append(ee.Feature(areas.to_ee(g, eco_srs, simplify=eco_res / 2 if eco_srs.IsProjected() else None), p))
	plan.eco = ee.FeatureCollection(feats)
	plan.ecoregions = [dict(p) for p in area['props']]
	plan.grids['ecoregion'] = Grid.around(eco_crs, eco_res, *envelope(reproject(area['union'], area['srs'], eco_crs,
	                                                                            eco_res)))


def contains_study(plan, area, study):
	"""Whether the polygons of the ecoregion grid of PLAN (the comparison AREA, or the ecoregions of RESOLVE) contain
	the rectangle of the STUDY grid: with OGR for a comparison area, with Earth Engine for the ecoregions."""
	res = plan.grids['ecoregion'].res
	if area is not None:
		return bool(reproject(area['union'], area['srs'], study.crs, res).Contains(rectangle(study)))
	return bool(plan.eco.geometry(res / 2).contains(study.region(), res / 2).getInfo())


def study_at_eco(study, eco):
	"""Grid of the rectangle of the STUDY grid at the pixel and in the CRS of the ECO grid (edges at multiples of its
	pixel)."""
	from monitoreo import areas
	from monitoreo.grid import Grid
	if areas.srs_of(study.crs).IsSame(areas.srs_of(eco.crs)):
		box = study.extent
	else:
		box = envelope(reproject(rectangle(study), areas.srs_of(study.crs), eco.crs, study.res))
	return Grid.around(eco.crs, eco.res, *box)


# Area of the files of the ecoregion grid ----------------------------------------------------------------------------
#
# comparison_area (absent for the RESOLVE ecoregions), ecoregions and ecoregions.gpkg of ecoregion/ describe the area
# whose download finished: the published area. A download of another area records it in comparison_area_pending,
# keeps the record of the files of the published area in previous_area (with their size and time, to know which are
# still intact) and its polygons in ecoregions.gpkg (a copy in POLYGONS_BEFORE), and publishes the new area only when
# all its files are there. So an interrupted switch leaves the folder labelled as the published area, running it
# again goes on with it (the files and export tasks of the new area are kept), and going back to the published area
# keeps its files that were not replaced yet.

def _signature(path):
	"""[size, modification time in ns] of the file PATH, or None if it does not exist."""
	try:
		st = os.stat(path)
	except OSError:
		return None
	return [st.st_size, st.st_mtime_ns]


def _file_id(path):
	try:
		st = os.stat(path)
	except OSError:
		return None
	return st.st_ino, st.st_size, st.st_mtime_ns


def _remove(path):
	if os.path.exists(path):
		os.remove(path)


def pending_key(plan, area):
	"""The key comparison_area_pending of the record of the ecoregion grid while the files of the area of PLAN (AREA,
	or the ecoregions of RESOLVE: {'resolve': True, 'name', 'ecoregions': ECO_IDs}) are downloaded."""
	if area is not None:
		return area_key(area)
	names = [str(e['ECO_NAME']) for e in plan.ecoregions if e.get('ECO_NAME')]
	return {'resolve': True, 'name': '; '.join(names) or None, 'ecoregions': [e.get('ECO_ID') for e in plan.ecoregions]}


def key_area(key):
	"""(kind, signature) of the area of a key comparison_area or comparison_area_pending."""
	if key.get('resolve'):
		return 'resolve', sorted(str(i) for i in key.get('ecoregions') or [])
	return 'comparison', key.get('sig')


def published_area(rec):
	"""(kind, signature) of the published area of REC, the record of an ecoregion grid (its comparison_area, else its
	ecoregions): ('comparison', sig, None when unknown) or ('resolve', ECO_IDs)."""
	if isinstance(rec.get('comparison_area'), dict):
		return key_area(rec['comparison_area'])
	ecos = [e for e in rec.get('ecoregions') or [] if isinstance(e, dict)]
	if any(e.get('BIOME_NAME') == COMPARISON or (isinstance(e.get('ECO_ID'), (int, float)) and e['ECO_ID'] < 0)
	       for e in ecos):
		return 'comparison', None
	return 'resolve', sorted(str(e.get('ECO_ID')) for e in ecos)


def recorded_area(rec):
	"""(kind, signature) of the area of the files that REC (the record of an ecoregion grid) has: that of
	comparison_area_pending while it is downloaded, else the published one; None if it has no files."""
	if isinstance(rec.get(PENDING), dict):
		return key_area(rec[PENDING])
	if not (rec.get('files') or rec.get('tasks')):
		return None
	return published_area(rec)


def new_area(plan, area):
	return ('comparison', area['sig']) if area else ('resolve', sorted(str(e.get('ECO_ID')) for e in plan.ecoregions))


def key_title(key):
	if key.get('resolve'):
		return f"the RESOLVE ecoregions ({key['name']})" if key.get('name') else 'the RESOLVE ecoregions'
	return f"the comparison area {key.get('name')}"


def published_title(rec):
	if isinstance(rec.get('comparison_area'), dict):
		return key_title(rec['comparison_area'])
	names = '; '.join(str(e['ECO_NAME']) for e in rec.get('ecoregions') or []
	                  if isinstance(e, dict) and e.get('ECO_NAME'))
	if published_area(rec)[0] == 'comparison':
		return f'the comparison area {names}'
	return f'the RESOLVE ecoregions ({names})' if names else 'the RESOLVE ecoregions'


def switch_of(rec, new):
	"""What a download of the area NEW (new_area) does with the files recorded in REC (the record of an ecoregion
	grid): 'same' (they are of that area), 'switch' (they are of another area, or there are none), 'resume' (they are
	of that area, whose download did not finish), 'back' (that area is the published one, after an unfinished download
	of another one) or 'other' (neither the pending nor the published area)."""
	if isinstance(rec.get(PENDING), dict):
		if key_area(rec[PENDING]) == new:
			return 'resume'
		return 'back' if published_area(rec) == new else 'other'
	old = recorded_area(rec)
	return 'same' if old is not None and old[1] is not None and old == new else 'switch'


def intact(folder, previous):
	"""The records of the files of PREVIOUS (previous_area of the record of the ecoregion grid in FOLDER) that did not
	change since it was written."""
	stat = previous.get('stat') or {}
	return {f: e for f, e in (previous.get('files') or {}).items()
	        if stat.get(f) is not None and _signature(os.path.join(folder, f)) == stat[f]}


def area_note(folder, plan, area):
	"""What the download of the area of PLAN (AREA) does with the files of the ecoregion grid in FOLDER when they are
	of another area, or None."""
	from monitoreo import recipe
	rec = recipe.load_record(folder, RECORD)
	what = switch_of(rec, new_area(plan, area))
	shown = where(folder)
	if what == 'switch' and (rec.get('files') or rec.get('tasks')):
		return (f'{shown} has the files of {published_title(rec)}: those of this area replace them, and the folder '
		        'stays that of the other area until all of them are there')
	if what == 'resume':
		return (f"{shown} goes on with the unfinished download of this area ({len(rec.get('files') or {})} of "
		        f'{len(plan.specs)} files there); it keeps {published_title(rec)} until all of them are there')
	if what == 'back':
		prev = rec.get(PREVIOUS) or {}
		return (f'the download of {key_title(rec[PENDING])} into {shown} did not finish: back to this area, with '
		        f"{len(intact(folder, prev))} of its {len(prev.get('files') or {})} files intact")
	if what == 'other':
		return (f'the unfinished download of {key_title(rec[PENDING])} into {shown} is left for this area; the folder '
		        f'keeps {published_title(rec)} until all the files of this one are there')
	return None


def start_area(folder, plan, area, log=print):
	"""Prepare the record of the ecoregion grid in FOLDER for the download of the area of PLAN (AREA, or the RESOLVE
	ecoregions), as explained above, before recipe.run downloads it. Files of another area are forgotten in 'files'
	(recipe.run would keep them when the grid and the variables did not change). Returns (whether another area than
	the published one is being downloaded, the ecoregions that recipe.run records meanwhile: those of the published
	area)."""
	from monitoreo import recipe
	os.makedirs(folder, exist_ok=True)
	rec = recipe.load_record(folder, RECORD)
	before = json.dumps(rec, sort_keys=True, default=str)
	what = switch_of(rec, new_area(plan, area))
	poly, copy, aside = (os.path.join(folder, f) for f in (POLYGONS, POLYGONS_BEFORE, POLYGONS_PENDING))
	going = ('the export tasks of {} go on in Earth Engine: cancel them in the Tasks tab of the Code Editor if they '
	         'are still running')
	if what in ('back', 'other') and rec.get('tasks'):
		log(f'  ({going.format(key_title(rec[PENDING]))})')
	if what == 'back':
		prev = rec.pop(PREVIOUS, None) or {}
		rec.pop(PENDING)
		rec['files'], rec['tasks'] = intact(folder, prev), prev.get('tasks') or {}
		if prev.get('grid'):
			rec['grid'] = prev['grid']
	elif what == 'other':
		rec['files'], rec['tasks'] = {}, {}
		rec[PENDING] = pending_key(plan, area)
		_remove(aside)
	elif what == 'switch':
		files = rec.get('files') or {}
		if rec.get('tasks'):
			log(f'  ({going.format(published_title(rec))})')
		if os.path.exists(poly):
			shutil.copy2(poly, copy)
		else:
			_remove(copy)
		rec[PREVIOUS] = {'files': files, 'tasks': rec.get('tasks') or {}, 'grid': rec.get('grid'),
		                 'stat': {f: s for f in files if (s := _signature(os.path.join(folder, f)))},
		                 'polygons': os.path.exists(poly)}
		rec['files'], rec['tasks'] = {}, {}
		rec[PENDING] = pending_key(plan, area)
		_remove(aside)
	if json.dumps(rec, sort_keys=True, default=str) != before:
		recipe.save_record(folder, rec, RECORD)
	# the polygons of the published area in place (a run killed while those of the new one were there)
	if os.path.exists(copy) and (PENDING in rec or what == 'back'):
		shutil.copy2(copy, poly)
	elif PENDING in rec and not (rec.get(PREVIOUS) or {}).get('polygons'):
		_remove(poly)
	if PENDING not in rec:
		for f in (copy, aside):
			_remove(f)
		return False, plan.ecoregions
	return True, rec.get('ecoregions') or []


def polygon_guard(folder, pending):
	"""A function to call during the download of the ecoregion grid in FOLDER by recipe.run (from its log), which
	writes the polygons of the grid to ecoregions.gpkg when it starts: while another area than the published one is
	PENDING, it moves them aside (end_area puts them in place when all the files are there) and puts back those of the
	published area, so that ecoregions.gpkg always matches comparison_area."""
	poly = os.path.join(folder, POLYGONS)
	first = _file_id(poly)
	moved = [not pending]

	def check():
		if moved[0] or _file_id(poly) == first:
			return
		os.replace(poly, os.path.join(folder, POLYGONS_PENDING))
		if os.path.exists(os.path.join(folder, POLYGONS_BEFORE)):
			shutil.copy2(os.path.join(folder, POLYGONS_BEFORE), poly)
		moved[0] = True
	return check


def end_area(folder, plan, area, done, log=print):
	"""After the download of the ecoregion grid of PLAN (AREA) to FOLDER: when all its files are there (DONE), publish
	its area (comparison_area, removed for the RESOLVE ecoregions; ecoregions; ecoregions.gpkg); else the record keeps
	describing the published area (its grid too). Then remove RESOLVE from the references of the files of a comparison
	area (uncite_resolve)."""
	from monitoreo import recipe
	rec = recipe.load_record(folder, RECORD)
	before = json.dumps(rec, sort_keys=True, default=str)
	pending = rec.get(PENDING) if isinstance(rec.get(PENDING), dict) else None
	if done:
		rec.pop(PENDING, None)
		rec.pop(PREVIOUS, None)
		rec['ecoregions'] = [dict(e) for e in plan.ecoregions]
		key = area_key(area)
		if key:
			rec['comparison_area'] = key
		else:
			rec.pop('comparison_area', None)
	elif pending and 'grid' in (rec.get(PREVIOUS) or {}):
		if rec[PREVIOUS]['grid']:
			rec['grid'] = rec[PREVIOUS]['grid']
		else:
			rec.pop('grid', None)
	if json.dumps(rec, sort_keys=True, default=str) != before:
		recipe.save_record(folder, rec, RECORD)
	if done:
		if os.path.exists(os.path.join(folder, POLYGONS_PENDING)):
			os.replace(os.path.join(folder, POLYGONS_PENDING), os.path.join(folder, POLYGONS))
		_remove(os.path.join(folder, POLYGONS_BEFORE))
		if pending:
			log(f'{where(folder)}: all the files of {key_title(pending)} are there; it is now the folder of that area')
	elif pending:
		log(f'{where(folder)} keeps {published_title(rec)} until all the files of {key_title(pending)} are there: run '
		    'the same command again to go on')
	try:
		uncite_resolve(folder, log)
	except Exception as e:  # it runs in a finally: never hide the outcome of the download
		log(f'WARNING: cannot remove RESOLVE from the references of the comparison area in {where(folder)}: {e}')


# References of a comparison area ------------------------------------------------------------------------------------

def without_resolve(text, spec=None):
	"""The references TEXT (one per line, as recipe.run writes them) without those of RESOLVE Ecoregions 2017 (but
	those that the variable of SPEC cites itself, if any)."""
	from monitoreo import references
	if not text:
		return text
	ids, credits = references.of_datasets([references.ECOREGIONS])
	marks = set(credits)
	for i in ids:
		r = references.get(i)
		marks.update(x for x in (r.get('text'), r.get('url')) if x)
	keep = set()
	if isinstance(spec, dict):
		period = spec.get('period') or [None, None]
		try:
			keep = set(references.texts(references.cite('download_image', dict(spec, start=period[0], end=period[-1]))))
		except (KeyError, ValueError, TypeError, IndexError):
			pass
	return '\n'.join(x for x in text.split('\n') if x in keep or not any(m in x for m in marks))


def uncite_resolve(folder, log=print):
	"""recipe.run cites RESOLVE Ecoregions 2017 as a dataset of every file of the ecoregion grid: remove it from the
	files recorded in the ecoregion grid FOLDER when they are of a comparison area (the metadata item references of
	their GeoTIFFs, and the references of their export tasks, written to the GeoTIFFs when they are imported). It can
	be run alone on a folder downloaded before. Returns the number of files and tasks corrected."""
	from osgeo import gdal
	from monitoreo import recipe
	gdal.UseExceptions()
	rec = recipe.load_record(folder, RECORD)
	if (recorded_area(rec) or (None,))[0] != 'comparison':
		return 0
	n = 0
	for f, e in (rec.get('files') or {}).items():
		path = os.path.join(folder, f)
		if not os.path.isfile(path):
			continue
		try:
			ds = gdal.Open(path)
			text = ds.GetMetadataItem('references')
			ds = None
			new = without_resolve(text, (e or {}).get('spec'))
			if text and new != text:
				ds = gdal.Open(path, gdal.GA_Update)
				ds.SetMetadataItem('references', new)
				ds = None
				n += 1
		except RuntimeError as err:
			log(f'WARNING: cannot remove RESOLVE from the references of {where(path)}: {err}')
	tasks = 0
	for f, t in (rec.get('tasks') or {}).items():
		new = without_resolve(t.get('references'), t.get('spec'))
		if t.get('references') and new != t['references']:
			t['references'] = new
			tasks += 1
	if tasks:
		recipe.save_record(folder, rec, RECORD)
	if n + tasks:
		log(f'{where(folder)}: RESOLVE Ecoregions 2017 removed from the references of {n + tasks} file(s) of the '
		    'comparison area (it is not used)')
	return n + tasks


# Main ---------------------------------------------------------------------------------------------------------------

def main():
	sys.stdout.reconfigure(line_buffering=True)
	parser = argparse.ArgumentParser(description='Download the input variables of gee_vars.yaml from Google Earth '
	                                             'Engine (catalog of MonitorEO) for the study area and its ecoregions '
	                                             'or a comparison area.')
	parser.add_argument('--recipe', default=RECIPE, help='recipe of the variables (default: gee_vars.yaml)')
	parser.add_argument('--list', action='store_true', help='list the variables of the catalog of MonitorEO and exit')
	parser.add_argument('--dry-run', action='store_true', help='show the grids, sizes and number of images, without downloading')
	parser.add_argument('--project', help='Earth Engine project (default: project of the recipe, or EE_PROJECT)')
	parser.add_argument('--folder', help='subfolder of ../inVars (default: folder of the recipe)')
	parser.add_argument('--study-area', metavar='FILE', help='vector file of the study area (default: the only one in pa/)')
	parser.add_argument('--layer', help='layer of the study area file, if it has several')
	parser.add_argument('--area', metavar='FILE', help='vector file of a comparison area for the ecoregion grid, '
	                    'instead of the ecoregions of RESOLVE that the study area touches (default: comparison_area of '
	                    "the recipe; '' = the ecoregions)")
	parser.add_argument('--area-layer', metavar='LAYER', help='layer of the comparison area file, if it has several '
	                                                          '(default: comparison_area_layer of the recipe)')
	parser.add_argument('--only', choices=('study', 'ecoregion'), help='download only the study grid, or only the '
	                    'ecoregion grid (with the study area at its pixel when its polygons do not contain it)')
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
		import ee
		from monitoreo import catalog, recipe, services
	except ImportError as e:
		fail(f'the package monitoreo is too old or broken ({e}): run without --no-update to install its last version')

	if args.list:
		print_catalog(catalog)
		return

	try:
		rec = recipe.read(args.recipe, extra=EXTRA_KEYS)
	except ValueError as e:
		fail(str(e))
	rec['drive_folder'] = rec.get('drive_folder') or DRIVE_FOLDER
	folder = str(args.folder or rec.get('folder') or 'gee').strip('/\\')
	cfg = hft_state.read_config()
	outdir = os.path.join(hft_state.INVARS, folder)
	dirs = recipe.folders(outdir)
	dirs_eco = os.path.join(dirs['ecoregion'], STUDY_ECO)

	study_file = setup.find_study_area(args.study_area)
	if not study_file:
		fail(f'no vector file of the study area in {setup.PA_DIR}: copy it there or use --study-area')
	layer = setup.check_layer(study_file, args.layer)

	# Comparison area (read locally, before connecting to Earth Engine): only for the ecoregion grid
	source = area_source(args, rec)
	eco_opts = rec.get('ecoregion') or {}
	area = None
	if source and args.only == 'study':
		print(f'Comparison area {source[0]}: not used with --only study')
	elif source:
		if not eco_opts.get('resolution'):
			fail('give the pixel size of the ecoregion grid (ecoregion: resolution) in the recipe')
		area = read_area(*source)

	project = args.project or rec.get('project') or os.environ.get('EE_PROJECT')
	if not project:
		if not hasattr(services, 'shared'):
			fail('give the Earth Engine project: project: in the recipe, --project or the variable EE_PROJECT')
		print('No Earth Engine project: using the shared account of MonitorEO (research EECUs of each month, '
		      'only direct downloads, not exports to Google Drive)')
	try:
		services.init(project)  # without a project, the shared account of MonitorEO (monitoreo >= 0.3.19)
		# the study grid is always planned (it decides whether the ecoregion grid needs the study area at its pixel),
		# and the ecoregions of RESOLVE only when there is no comparison area
		plan = recipe.plan(dict(rec, ecoregion=None) if area else rec, study_file, layer,
		                   buffer=float(cfg.get('PA_BUFFER', 0)), only='study' if args.only == 'study' else None,
		                   extra=EXTRA_KEYS, warn=print)
		plan.recipe = rec
		study = plan.grids['study']
		if area:
			use_area(plan, area, str(eco_opts.get('crs') or study.crs), float(eco_opts['resolution']))
		if args.only == 'ecoregion':
			if 'ecoregion' not in plan.grids:
				fail('the recipe has no ecoregion: resolution')
			del plan.grids['study']
		study_eco = None  # the study grid at the pixel of the ecoregion grid, and why
		if 'ecoregion' in plan.grids:
			if not contains_study(plan, area, study):
				study_eco = study_at_eco(study, plan.grids['ecoregion']), (
					f"the {'comparison area' if area else 'ecoregions'} do{'es' if area else ''} not contain the whole "
					'rectangle of the study grid')
			elif os.path.exists(os.path.join(dirs_eco, RECORD)):
				study_eco = study_at_eco(study, plan.grids['ecoregion']), 'it was downloaded before: kept up to date'
	except (ValueError, RuntimeError) as e:
		fail(str(e))
	except ee.EEException as e:
		fail(f'Earth Engine: {e}')
	if getattr(services, 'shared', lambda: False)():
		drive_grids = [n for n in plan.grids if args.drive or (rec.get(n) or {}).get('drive')]
		if drive_grids and not args.dry_run:
			fail(f"the shared account of MonitorEO cannot export to Google Drive ({' and '.join(drive_grids)} grid with "
			     'drive: true): set drive: false in the recipe to download the grid directly in tiles, or give your '
			     'own Earth Engine project (--project, project: in the recipe or EE_PROJECT)')

	rel = where
	summary = plan.summary()
	print(f'Study area: {study_file}')
	for name, g in plan.grids.items():
		label = 'Study area' if name == 'study' else 'Comparison area' if area else 'Ecoregions'
		print(f'{label} grid: {grid_text(g)} -> {rel(dirs[name])}')
	if study_eco:
		print(f'Study area at the pixel of the ecoregion grid: {grid_text(study_eco[0])} -> {rel(dirs_eco)}\n'
		      f'  ({study_eco[1]}: for the statistics of the HFTs in the similarity)')
	if 'ecoregion' in plan.grids:
		if area:
			n = len(area['props'])
			print(f"Comparison area: {area['file']}{' (layer ' + area['layer'] + ')' if area['layer'] else ''} "
			      f"({n} polygon{'s' if n > 1 else ''})")
			for p in plan.ecoregions:
				print(f"  {p['ECO_ID']}: {p['ECO_NAME']}")
		else:
			for p in plan.ecoregions:
				print(f"  ecoregion {p['ECO_ID']}: {p['ECO_NAME']} ({p['BIOME_NAME']})")
		note = area_note(dirs['ecoregion'], plan, area)
		if note:
			print(f'  {note}')
	print(f'{len(plan.specs)} variables:')
	for f in summary['files']:
		filters = f.get('pixel_filters')
		print(f"  {f['file']}" + (f"  (pixel filters: {', '.join(filters) or 'none'})" if filters is not None else ''))

	if args.dry_run:
		print(f"\nImages of each variable over the {'study area' if 'study' in plan.grids else 'ecoregion grid'}:")
		for f, n in recipe.n_images(plan).items():
			print(f"  {f}: {'fixed image' if n is None else n}")
		return

	# Downloads, a grid after the other: directly, in tiles, or by export tasks to Google Drive for the grids with
	# drive: true (or --drive)
	bars = {}

	def progress(file, done, total):
		bar = bars.get(file)
		if bar is None:
			bar = bars[file] = tqdm(total=total, desc=f'  {file}', unit='tile', leave=True)
		bar.update(done - bar.n)
		if done == total:
			bar.close()
			del bars[file]  # the next grid has a file with the same name

	def download(p, out, log=print, wait=not args.no_wait, **kw):
		try:
			return recipe.run(p, out, RECORD, args.overwrite, wait=wait, tile=args.tile, workers=args.workers,
			                  project=project, log=log, progress=progress, **kw)
		except (ValueError, RuntimeError) as e:
			for bar in list(bars.values()):
				bar.close()
			fail(str(e))

	def upgrade(path):
		if os.path.isdir(path) and upgrade_record(path, plan.specs):
			print(f'{rel(path)}: record of the files updated to monitoreo {installed_monitoreo()[0]} (the same '
			      'variables)')

	def stop(signum, frame):  # the app stops a download with SIGTERM: clean up as with Ctrl+C (finally)
		sys.exit(f'Stopped (signal {signum})')

	if hasattr(signal, 'SIGTERM'):
		signal.signal(signal.SIGTERM, stop)
	done, study_tasks = True, None
	if 'study' in plan.grids:
		upgrade(dirs['study'])
		p = recipe.Plan(rec, plan.specs, {'study': plan.grids['study']}, plan.area)
		# with the ecoregion grid too, the export tasks of the study grid run while that grid is downloaded and are
		# waited for after it
		later = 'ecoregion' in plan.grids and not args.no_wait

		def study_log(message):
			if later and message.startswith('The export tasks go on in Earth Engine'):
				message = 'The export tasks go on in Earth Engine: they are waited for after the ecoregion grid.'
			print(message)

		res = download(p, outdir, log=study_log, wait=not args.no_wait and not later, drive_all=args.drive,
		               drive_dir=args.drive_dir)
		if later and res['pending']:
			study_tasks = p
			done = not (res['in_drive'] or res['failed'])
		else:
			done = res['done']
	if 'ecoregion' in plan.grids:
		# the folder keeps the area of its files until all those of the new area are there (start_area)
		eco_dir = dirs['ecoregion']
		pending, ecoregions = start_area(eco_dir, plan, area)
		upgrade(eco_dir)
		check = polygon_guard(eco_dir, pending)

		def eco_log(message):
			check()
			if area and message.startswith('\nEcoregions ('):  # the label of monitoreo for this grid
				message = '\nComparison area (' + message[len('\nEcoregions ('):]
			print(message)

		res = None
		try:
			res = download(recipe.Plan(rec, plan.specs, {'ecoregion': plan.grids['ecoregion']}, plan.area, plan.eco,
			                           ecoregions), outdir, log=eco_log, drive_all=args.drive, drive_dir=args.drive_dir)
		finally:
			check()
			end_area(eco_dir, plan, area, bool(res and res['done']))
		done = done and res['done']
		if res['pending'] and not args.no_wait:  # Ctrl+C stopped waiting for its export tasks: do not wait for others
			study_tasks = None
			done = False
	if study_eco:
		print(f'\nStudy area at the pixel of the ecoregion grid ({study_eco[1]}) -> {rel(dirs_eco)}')
		upgrade(dirs_eco)
		directly = dict(rec, study=dict(rec.get('study') or {}, drive=False))  # never through Google Drive: it is small
		res = download(recipe.Plan(directly, plan.specs, {'study': study_eco[0]}, plan.area), dirs_eco,
		               drive_all=False, prefix=re.sub(r'[^A-Za-z0-9_-]', '_', folder) + '_studyeco_')
		done = done and res['done']
	if study_tasks:
		quiet = re.compile(r'  \S.*: (up to date|export task already started)$')
		res = download(study_tasks, outdir, log=lambda m: quiet.match(m) or print(m), drive_all=args.drive,
		               drive_dir=args.drive_dir)
		done = done and res['done']
	if not done:
		return
	print(f"\nDone: variables in {hft_state.where(folder)} (NoData -9999).")
	if not hft_state.raster_files(folder):
		print(f"{hft_state.where(folder)} has only the ecoregion grid: it is a comparison area for the similarity of "
		      f"the HFTs (ECOSIM_AREA = '{folder}/ecoregion' in config.py), not a set of variables for the HFT model.")
		if args.use:
			print(f'INVARS_FOLDER not changed: {hft_state.where(folder)} has no variables of the study area.')
	elif args.use:
		hft_state.write_option('INVARS_FOLDER', repr(folder))
		print(f"config.py: INVARS_FOLDER = {folder!r}. Next: python3 run_hft.py")
	elif hft_state.invars_folder(cfg) != folder:
		print(f"To use them, set INVARS_FOLDER = {folder!r} in config.py (or run with --use), then python3 run_hft.py")


if __name__ == '__main__':
	main()
