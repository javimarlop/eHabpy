#### eHabitat+ setup: GRASS GIS project and config.py from the input data
#### License: CC BY-SA 3.0

"""
Run once from the segm_optim folder, after copying the input variables
(GeoTIFF) to ../inVars and the vector file of the study area to pa/:

    python3 setup.py

With INVARS_FOLDER in config.py it takes the variables of that subfolder of
../inVars instead (and run_hft.py runs it again when INVARS_FOLDER changes).

It
  1. reads the CRS, resolution and NoData value of every ../inVars/*.tif
     (../inVars/<INVARS_FOLDER>/*.tif),
  2. creates the GRASS project GRASSDB/<study area> with the CRS of the
     variables (if it does not exist yet): the project and the map of the
     study area are both named after the study area file,
  3. imports the variables and the study area into its PERMANENT mapset,
  4. writes palist.csv with every ID of the study area (if it does not exist),
  5. fills GRASSDB, GRASSLOC, ENV_VARS0, STUDY_AREA, COL_ID, RESOLUTION and
     ENV_VARS in config.py; the other options are kept as they are,
  6. records what it imported in run_state.json (see hft_state.py).

The map name of each variable is its file name, preceded by INVARS_FOLDER
for the files of a subfolder (sin_twi/ndvi.tif -> sin_twi_ndvi), so the
variables of several folders can be in the project at the same time.
Running it again only imports the maps that are not in the project yet and
those whose file changed since setup.py imported it (recorded in
run_state.json); --overwrite imports all of them again.
Run 'python3 setup.py --help' for the other options.
"""

import argparse
import glob
import math
import os
import posixpath
import re
import subprocess
import sys

from osgeo import gdal, ogr, osr

import hft_state

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

HERE = os.path.dirname(os.path.abspath(__file__))
INVARS = os.path.normpath(os.path.join(HERE, '..', 'inVars'))
PA_DIR = os.path.join(HERE, 'pa')
CONFIG = os.path.join(HERE, 'config.py')
PALIST = os.path.join(HERE, 'palist.csv')

RASTER_EXT = ('.tif', '.tiff')
VECTOR_EXT = ('.gpkg', '.shp', '.geojson')
SAME_CRS = ['IGNORE_DATA_AXIS_TO_SRS_AXIS_MAPPING=YES', 'CRITERION=EQUIVALENT']


def fail(msg):
	sys.exit('ERROR: ' + msg)


def grass_name(path, folder=''):
	"""GRASS map name from a file name (after FOLDER, the subfolder of ../inVars, if any): letters, digits
	and _, starting with a letter."""
	stem = os.path.splitext(os.path.basename(path))[0]
	name = re.sub(r'\W', '_', f'{folder}_{stem}' if folder else stem, flags=re.ASCII)
	return name if name[:1].isalpha() else 'v_' + name


def fmt_num(v):
	"""Number as written in config.py (-9999 rather than -9999.0)."""
	if v is None:
		return 'None'
	if math.isnan(v):
		return "float('nan')"
	return str(int(v)) if float(v).is_integer() else repr(float(v))


def same_crs(a, b):
	return bool(a.IsSame(b, SAME_CRS)) or a.ExportToProj4() == b.ExportToProj4()


def read_variables(folder):
	"""CRS, resolution and NoData value of every GeoTIFF in ../inVars/FOLDER."""
	indir = os.path.join(INVARS, folder)
	files = sorted(f for f in glob.glob(os.path.join(indir, '*')) if f.lower().endswith(RASTER_EXT) and os.path.isfile(f))
	if not files:
		fail(f'no GeoTIFF files (.tif) found in {indir}')

	variables = {}
	for f in files:
		name = grass_name(f, folder)
		if name in variables:
			fail(f"{variables[name]['file']} and {os.path.basename(f)} give the same map name '{name}': rename one of them")
		ds = gdal.Open(f)
		srs = ds.GetSpatialRef()
		if srs is None:
			fail(f'{f} has no CRS')
		# RESOLUTION and PA_BUFFER are in map units, and in degrees the pixels are neither square nor of equal area
		if srs.IsGeographic():
			fail(f'{f} has a geographic CRS ({srs.GetName()}): reproject the variables to a projected CRS '
			     'in metres (e.g. UTM), for example with gdalwarp -t_srs EPSG:<code>')
		gt = ds.GetGeoTransform()
		if ds.RasterCount > 1:
			print(f'WARNING: {os.path.basename(f)} has {ds.RasterCount} bands, only the first one is used')
		variables[name] = {
			'file': posixpath.join(folder, os.path.basename(f)),  # relative to ../inVars
			'path': f,
			'srs': srs,
			'res': (abs(gt[1]), abs(gt[5])),
			'nodata': ds.GetRasterBand(1).GetNoDataValue(),
		}

	# The model reads the GeoTIFFs directly with a single pixel size (RESOLUTION) for all of them
	first = next(iter(variables.values()))
	res = first['res'][0]
	for name, v in variables.items():
		if not same_crs(v['srs'], first['srs']):
			fail(f"{v['file']} and {first['file']} have different CRS: reproject all variables to the same CRS")
		if not all(math.isclose(r, res, rel_tol=1e-6) for r in v['res']):
			fail(f"pixel size of {v['file']} is {v['res'][0]} x {v['res'][1]} and of {first['file']} {res} x {res}: "
			     'resample all variables to the same square pixel size')
		if v['nodata'] is None:
			print(f"WARNING: {v['file']} has no NoData value: all its pixels are taken as valid data")
	return variables, res


def find_study_area(path):
	"""Vector file of the study area: the one given, else the only vector file in pa/."""
	if path:
		if not os.path.exists(path):
			fail(f'study area file not found: {path}')
		return path
	files = sorted(f for f in glob.glob(os.path.join(PA_DIR, '*')) if f.lower().endswith(VECTOR_EXT))
	if len(files) > 1:
		fail(f'several vector files in {PA_DIR}: {", ".join(map(os.path.basename, files))}\n'
		     'Choose the study area with --study-area')
	return files[0] if files else None


def check_layer(path, layer):
	"""Layer to import from the study area file (None if it only has one)."""
	ds = ogr.Open(path)
	layers = [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]
	if layer and layer not in layers:
		fail(f'layer {layer} not found in {path}; its layers are: {", ".join(layers)}')
	if not layer and len(layers) > 1:
		fail(f'{path} has several layers ({", ".join(layers)}): choose one with --layer')
	return layer


def find_grass():
	"""GRASS GIS executable and installation folder (GISBASE)."""
	for exe in ('grass', 'grass85', 'grass84', 'grass83', 'grass82', 'grass80'):
		try:
			out = subprocess.check_output([exe, '--config', 'path'], text=True,
			                              shell=(os.name == 'nt'), stderr=subprocess.DEVNULL)
			return exe, out.strip().splitlines()[-1]
		except Exception:
			continue
	fail('GRASS GIS not found: activate the environment first (conda activate ehabpy)')


set_option = hft_state.set_option


def set_env_vars(text, entries):
	"""Replace the entries of the ENV_VARS dictionary, keeping its commented-out lines."""
	pattern = re.compile(r'^(ENV_VARS\s*=\s*\{[^\n]*\n)(.*?)(^\})', re.M | re.S)
	m = pattern.search(text)
	if not m:
		return text.rstrip('\n') + '\nENV_VARS = {\n' + ''.join(e + '\n' for e in entries) + '}\n'
	comments = [l + '\n' for l in m.group(2).splitlines() if l.lstrip().startswith('#')]
	body = ''.join(e + '\n' for e in entries) + ''.join(comments)
	return text[:m.start(2)] + body + text[m.end(2):]


def main():
	sys.stdout.reconfigure(line_buffering=True)  # keep our messages in order with those of GRASS
	if not os.path.isfile(CONFIG):
		fail(f'{CONFIG} not found')
	import config

	parser = argparse.ArgumentParser(
		description='Create the GRASS GIS project from ../inVars and pa/ and fill config.py.')
	parser.add_argument('--grassdb', default=getattr(config, 'GRASSDB', '~/grassdata'),
	                    help='GRASS database folder (default: GRASSDB in config.py)')
	parser.add_argument('--study-area', metavar='FILE',
	                    help=f'vector file of the study area (default: the only {"/".join(VECTOR_EXT)} file in pa/); '
	                         'its name is also the name of the GRASS project')
	parser.add_argument('--layer', help='layer of the study area file, if it has several')
	parser.add_argument('--col-id', default=getattr(config, 'COL_ID', 'cat'),
	                    help='integer column with the study area IDs (default: COL_ID in config.py)')
	parser.add_argument('--overwrite', action='store_true',
	                    help='import again the maps that are already in the project')
	args = parser.parse_args()

	# 1. Check the input data before creating anything
	folder = hft_state.invars_folder(hft_state.read_config())
	indir = os.path.join(INVARS, folder)
	if not os.path.isdir(indir) or os.path.commonpath([os.path.realpath(indir), os.path.realpath(INVARS)]) != os.path.realpath(INVARS):
		fail(f'INVARS_FOLDER = {folder!r} in config.py: folder {indir} not found; it must be a subfolder of {INVARS} '
		     "('' for the files directly in it)")
	variables, res = read_variables(folder)
	old_order = list(getattr(config, 'ENV_VARS', {}))
	order = [n for n in old_order if n in variables] + [n for n in variables if n not in old_order]
	ref = variables[order[0]]
	print(f'Variables in {indir}: {", ".join(order)}')
	print(f'Resolution: {fmt_num(res)} | CRS: {ref["srs"].GetName()}')

	sa_file = find_study_area(args.study_area)
	layer = check_layer(sa_file, args.layer) if sa_file else None
	# The GRASS project and the map of the study area take the name of the study area
	study_area = grass_name(sa_file) if sa_file else getattr(config, 'STUDY_AREA', '')

	grassdb = os.path.expanduser(args.grassdb)
	project = os.path.join(grassdb, study_area)
	exists = bool(study_area) and os.path.isdir(os.path.join(project, 'PERMANENT'))
	if not sa_file and not exists:
		fail(f'no vector file of the study area found in {PA_DIR}: copy it there or use --study-area')

	# 2. GRASS project with the CRS of the variables
	exe, gisbase = find_grass()
	if not exists:
		if os.path.exists(project):
			fail(f'{project} exists but is not a GRASS project: move it or rename the study area file')
		os.makedirs(grassdb, exist_ok=True)
		print(f'Creating GRASS project {project}')
		subprocess.run([exe, '-c', ref['path'], '-e', project], check=True, shell=(os.name == 'nt'))

	os.environ['GISBASE'] = gisbase
	sys.path.append(os.path.join(gisbase, 'etc', 'python'))
	import grass.script as gs
	import grass.script.setup as gsetup
	gsetup.init(grassdb, study_area, 'PERMANENT')

	try:
		if exists:
			proj_srs = osr.SpatialReference()
			proj_srs.ImportFromWkt(gs.read_command('g.proj', flags='w'))
			if not same_crs(proj_srs, ref['srs']):
				fail(f'the GRASS project {project} has a different CRS than the variables: '
				     'move it or rename the study area file')
			print(f'Using existing GRASS project {project}')

		# 3. Import the variables and the study area. A map already in the project is imported again
		#    when its file changed since setup.py imported it (size and time recorded in run_state.json)
		imported = hft_state.get_setup(study_area) or {}
		old_sigs = imported.get('variables', {})
		for name in order:
			v = variables[name]
			sig = hft_state.file_sig(v['path'])
			if gs.find_file(name, element='raster', mapset='PERMANENT')['name'] and not args.overwrite:
				if old_sigs.get(v['file'], sig) == sig:
					print(f'  {name}: already in the project')
					continue
				print(f"  {name}: {v['file']} changed since it was imported, importing it again")
			else:
				print(f"  {name}: importing {v['file']}")
			gs.run_command('r.import', input=v['path'], output=name, band=1, overwrite=True, quiet=True)

		if sa_file:
			sa_sigs = hft_state.vector_sigs(sa_file)
			if gs.find_file(study_area, element='vector', mapset='PERMANENT')['name'] and not args.overwrite:
				if imported.get('study_area') in (None, sa_sigs):
					print(f'  {study_area}: already in the project')
				else:
					print(f'  {study_area}: {os.path.basename(sa_file)} changed since it was imported, importing it again')
					gs.run_command('v.import', input=sa_file, layer=layer, output=study_area, overwrite=True, quiet=True)
			else:
				print(f'  {study_area}: importing {os.path.basename(sa_file)}')
				gs.run_command('v.import', input=sa_file, layer=layer, output=study_area, overwrite=True, quiet=True)
		else:
			if not gs.find_file(study_area, element='vector', mapset='PERMANENT')['name']:
				fail(f'no vector file of the study area found in {PA_DIR}: copy it there or use --study-area')
			print(f'  {study_area}: no vector file in pa/, using the study area already in the project')

		# 4. Study area IDs
		columns = gs.vector_columns(study_area)
		if args.col_id not in columns:
			fail(f'column {args.col_id} not found in {study_area}; its columns are: {", ".join(columns)}\n'
			     'Choose the ID column with --col-id')
		if 'INT' not in columns[args.col_id]['type'].upper():
			fail(f'column {args.col_id} of {study_area} is {columns[args.col_id]["type"]}, but the IDs must be integers '
			     '(use --col-id cat for one ID per polygon)')
		values = gs.vector_db_select(study_area, columns=args.col_id)['values'].values()
		ids = sorted({int(v[0]) for v in values if v[0] != ''})
	finally:
		gsetup.finish()

	if not os.path.isfile(PALIST) or not open(PALIST).read().split():
		with open(PALIST, 'w') as f:
			f.write(''.join(f'{i}\n' for i in ids))
		print(f'palist.csv written with the {len(ids)} IDs of {study_area} ({args.col_id})')
	else:
		unknown = [i for i in open(PALIST).read().split() if i not in {str(x) for x in ids}]
		if unknown:
			print(f'WARNING: palist.csv has IDs that are not in {study_area} ({args.col_id}): {", ".join(unknown)}\n'
			      '         Delete palist.csv and run setup.py again to list all IDs of the study area')
		else:
			print('palist.csv kept as it is')

	# 5. config.py
	entries = [f"    {'' if i == 0 else ','}{n!r}: {{'file': {variables[n]['file']!r}, 'nodata': {fmt_num(variables[n]['nodata'])}}}"
	           for i, n in enumerate(order)]
	with open(CONFIG) as f:
		text = f.read()
	for name, value in [('GRASSDB', repr(args.grassdb)), ('GRASSLOC', repr(study_area)),
	                    ('ENV_VARS0', f"[{','.join(order)!r}]"), ('STUDY_AREA', repr(study_area)),
	                    ('COL_ID', repr(args.col_id)), ('RESOLUTION', fmt_num(res))]:
		text = set_option(text, name, value)
		print(f'config.py: {name} = {value}')
	text = set_env_vars(text, entries)
	with open(CONFIG, 'w') as f:
		f.write(text)
	print('config.py: ENV_VARS = {\n' + '\n'.join(entries) + '\n}')

	# 6. What was imported, so that the other steps know when the input data change. The files of the
	#    other folders of ../inVars imported before are kept, to know whether their maps are up to date
	others = {f: s for f, s in old_sigs.items() if posixpath.dirname(f) != folder}
	hft_state.set_setup(study_area, {
		'invars_folder': folder,
		'grassdb': args.grassdb,
		'study_area_file': os.path.abspath(sa_file) if sa_file else None,
		'layer': layer,
		'col_id': args.col_id,
		'variables': dict(others, **{variables[n]['file']: hft_state.file_sig(variables[n]['path']) for n in order}),
		'study_area': hft_state.vector_sigs(sa_file) if sa_file else None,
	})

	print('\nSetup done. Check the other options in config.py and run: python3 run_hft.py')


if __name__ == '__main__':
	main()
