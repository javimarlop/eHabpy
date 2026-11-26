#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / GRASS GIS 8 (2024)
#### License: CC BY-SA 3.0

from multiprocessing import cpu_count, Pool
import multiprocessing
import subprocess
import numpy as np
from scipy import ndimage
import math
import re
import os
import sys
import shutil

import hft_state
import pca_axes

# Importar configuración global de variables y opciones
try:
    from config import ENV_VARS0, RESOLUTION, COL_ID, GRASSDB, GRASSLOC, STUDY_AREA
except ImportError:
    print("ERROR: No se encuentra config.py. Asegúrate de crearlo en el mismo directorio.")
    sys.exit(1)

# OPCIONES DE CONFIGURACIÓN
CLIP_TO_PA = getattr(sys.modules['config'], 'CLIP_TO_PA', False)
FORCE_RESTART = getattr(sys.modules['config'], 'FORCE_RESTART', False) or os.environ.get('EHAB_FORCE_RESTART', '') not in ('', '0')
PA_BUFFER = getattr(sys.modules['config'], 'PA_BUFFER', 0)
SMALL_SEGM_RADIUS = getattr(sys.modules['config'], 'SMALL_SEGM_RADIUS', 5)
MIN_SEGM_SIZE = int(getattr(sys.modules['config'], 'MIN_SEGM_SIZE', 0) or 0)  # 0 = sqrt of the pixels with data
# 'frequent' (the rule up to eHabpy 1.5, also without the option) or 'similar', read as hft_state records it
try:
	SMALL_SEGM_MERGE = hft_state.merge_rule(vars(sys.modules['config']))
except ValueError as e:
	sys.exit(f'ERROR: {e}')
# Thresholds of i.segment, from the lowest to the highest (k = 1..n, their position, names their files), and the layout of segm_id
try:
	SEGM_THRESHOLDS = hft_state.parse_thresholds(getattr(sys.modules['config'], 'SEGM_THRESHOLDS', hft_state.DEFAULT_THRESHOLDS))
except ValueError as e:
	sys.exit(f'ERROR: {e}')
K_DIGITS, CAT_DIGITS = hft_state.segm_id_layout(SEGM_THRESHOLDS)
SEGM_PCA = getattr(sys.modules['config'], 'SEGM_PCA', True)
try:  # the axes of the PCA: a number of axes (SEGM_PCA_AXES); else 0 = the first 3, a % = the first ones that explain at
	# least that % of the variance (SEGM_PCA_VARIANCE)
	SEGM_PCA_AXES = pca_axes.axes_option(vars(sys.modules['config']))
	SEGM_PCA_VARIANCE = pca_axes.variance_option(vars(sys.modules['config']))
except ValueError as e:
	sys.exit(f'ERROR: {e}')
CLIP_PERCENTILE = float(getattr(sys.modules['config'], 'CLIP_PERCENTILE', 0) or 0)
if not 0 <= CLIP_PERCENTILE < 50:
	sys.exit(f'ERROR: CLIP_PERCENTILE in config.py must be from 0 to less than 50, not {CLIP_PERCENTILE}')

SEGM_VARS = [v.strip() for v in ','.join(ENV_VARS0).split(',') if v.strip()]
FIRST_VAR = SEGM_VARS[0]  # its valid pixels define the analysed area

GRASSDBASE = os.environ.get('GRASSDBASE', os.path.expanduser(GRASSDB))
MYLOC = os.environ.get('GRASSLOC', GRASSLOC)
NPROC = int(os.environ.get('EHAB_NPROC', max(1, cpu_count() - 1)))


def find_gisbase():
	gb = os.environ.get('GISBASE')
	if gb and os.path.isdir(os.path.join(gb, 'etc')):
		return gb
	for exe in ('grass', 'grass85', 'grass84', 'grass83', 'grass82', 'grass80'):
		try:
			# shell=True on Windows so that grass84.bat & co. are found
			out = subprocess.check_output([exe, '--config', 'path'], text=True,
			                              shell=(os.name == 'nt'), stderr=subprocess.DEVNULL)
			return out.strip().splitlines()[-1]
		except Exception:
			continue
	return None


def add_grass_to_path():
	gisbase = find_gisbase()
	if gisbase:
		os.environ['GISBASE'] = gisbase
		py = os.path.join(gisbase, 'etc', 'python')
		if py not in sys.path:
			sys.path.append(py)
	return gisbase


def fix_grass_env(gisbase):
	"""GRASS 8.5 on Windows, started from Python with setup.init(): init() changes GISBASE to
	<GISBASE>\\lib\\grass85, which does not exist, and the modules then do not find etc (without its color tables
	r.neighbors, r.colors... fail), so GISBASE goes back; and the modules of vectors need HOME, which Windows does not
	define (grass.bat sets it, init() does not). Elsewhere it does nothing."""
	if gisbase and not os.path.isdir(os.path.join(os.environ.get('GISBASE', ''), 'etc')):
		os.environ['GISBASE'] = gisbase
	if os.name == 'nt' and not os.environ.get('HOME'):
		os.environ['HOME'] = os.path.expanduser('~')


def init_grass(gisdbase, location, mapset):
	gisbase = add_grass_to_path()
	import grass.script as grass
	import grass.script.setup as gsetup
	try:
		gsetup.init(gisdbase, location, mapset)
	except TypeError:
		gsetup.init(gisbase, gisdbase, location, mapset)
	fix_grass_env(gisbase)
	return grass, gsetup


def clean_previous_results():
	"""Borra archivos de control y csv/shp/tiffs (results/ no: moranvar_plots.R mueve los resultados
	anteriores del escenario a results/<STUDY_AREA>/<SCENARIO>_previous cuando escribe los nuevos)."""
	control_files = ['csv/segm_done.csv', 'ongoing.csv', 'done.csv']
	for f in control_files:
		if os.path.exists(f):
			try:
				os.remove(f)
			except Exception as e:
				print(f"No se pudo eliminar {f}: {e}")
	hft_state.reset('segm')

	directories = ['csv', 'shp', 'tiffs']
	for folder in directories:
		if os.path.exists(folder):
			for filename in os.listdir(folder):
				file_path = os.path.join(folder, filename)
				try:
					if os.path.isfile(file_path) or os.path.islink(file_path):
						os.unlink(file_path)
					elif os.path.isdir(file_path):
						shutil.rmtree(file_path)
				except Exception as e:
					print(f"No se pudo eliminar {file_path}: {e}")
		else:
			os.makedirs(folder, exist_ok=True)
			
	print("✓ Limpieza completada exitosamente.\n")


csvname1 = 'csv/segm_done.csv'
csvong1 = 'ongoing.csv'
csvong2 = 'done.csv'
# csv, shp y tiffs los comparten todas las áreas de estudio: este archivo dice de cuál son
csvarea = 'csv/study_area.csv'


def segm_files(pa, k):
	"""Files of the segments of area PA and threshold K read by the statistics step (ehab_segm_core.py): the segm_id
	of every segment, the table cat -> segm_id and the segment raster (labels = cat, 0 outside the segments)."""
	return (f'csv/park_{pa}_{k}.csv', f'csv/segm_{pa}_{k}_ids.csv', f'tiffs/segm_{pa}_{k}.tif')


def _layer_means(data, layers, n):
	"""Mean of every layer (pixels x layers, NaN = no data) in every segment of DATA (0 = none): n x layers."""
	flat = data.ravel()
	out = np.full((n, layers.shape[1]), np.nan)
	for j in range(layers.shape[1]):
		v = layers[:, j]
		ok = (flat > 0) & np.isfinite(v)
		cnt = np.bincount(flat[ok], minlength=n)
		with np.errstate(invalid='ignore', divide='ignore'):
			out[:, j] = np.bincount(flat[ok], weights=v[ok], minlength=n) / cnt
	return out


def merge_small_segments(grass, pa3, minarea, layers=None):
	"""Une cada segmento de pa3 con menos de minarea celdas a uno de los mayores que él a
	SMALL_SEGM_RADIUS píxeles o menos (distancia entre centros de celda): con SMALL_SEGM_MERGE =
	'similar', al más parecido (distancia euclídea entre las medias de las capas LAYERS de la
	segmentación, escaladas a su rango como hace i.segment), y con 'frequent', al más frecuente
	en ese radio (el criterio anterior, que no mira si se parecen). En pasadas repetidas hasta
	que no cambia nada: los grupos de fragmentos cercanos acaban en un mismo segmento. Los que
	siguen siendo pequeños (sin ningún segmento mayor tan cerca, p. ej. fragmentos aislados entre
	píxeles sin datos) se unen después al segmento de minarea celdas o más más cercano, a
	cualquier distancia.
	Se hace en memoria y no con un mapa por segmento: con datos enmascarados hay miles de
	fragmentos aislados y r.patch agotaba los archivos abiertos."""
	import grass.script.array as garray

	data = np.asarray(garray.array(mapname=pa3, null=0, dtype=np.int32))
	ids, counts = np.unique(data[data > 0], return_counts=True)
	if not (counts < minarea).any():
		print(f"{pa3}: no segment < {minarea} cells")
		return
	similar = SMALL_SEGM_MERGE == 'similar' and layers
	if similar:
		# the layers of the segmentation, each scaled to its range in the region (as i.segment does)
		vals = []
		for name in layers:
			a = np.asarray(garray.array(mapname=name, null=np.nan), dtype=np.float32).ravel()
			lo, hi = np.nanmin(a), np.nanmax(a)
			vals.append((a - lo) / (hi - lo) if hi > lo else np.zeros_like(a))
		vals = np.stack(vals, axis=1)
	margin = int(SMALL_SEGM_RADIUS) + 1  # cells around the segment that can be within the radius
	npass = 0
	while True:
		ids, counts = np.unique(data[data > 0], return_counts=True)
		small = ids[counts < minarea]
		if npass == 0:
			nsmall = len(small)
		size = np.zeros(data.max() + 1, dtype=np.int64)
		size[ids] = counts
		lut = np.arange(data.max() + 1, dtype=np.int32)
		boxes = ndimage.find_objects(data)
		means = _layer_means(data, vals, data.max() + 1) if similar and len(small) else None
		for s in small:
			rows, cols = boxes[s - 1]
			sub = data[max(rows.start - margin, 0):rows.stop + margin, max(cols.start - margin, 0):cols.stop + margin]
			inside = sub == s
			near = sub[(ndimage.distance_transform_edt(~inside) <= SMALL_SEGM_RADIUS) & ~inside]
			near = near[near > 0]
			# Only larger segments (the larger ID if equal), so two segments never swap their IDs
			near = near[(size[near] > size[s]) | ((size[near] == size[s]) & (near > s))]
			if near.size:
				ids_near, n = np.unique(near, return_counts=True)
				if means is not None and np.isfinite(means[s]).all():
					d = ((means[ids_near] - means[s]) ** 2).sum(axis=1)
					d[~np.isfinite(d)] = np.inf  # a candidate without data in some layer: the last choice
					lut[s] = ids_near[np.argmin(d)] if np.isfinite(d).any() else ids_near[np.argmax(n)]
				else:
					lut[s] = ids_near[np.argmax(n)]
		if (lut == np.arange(len(lut))).all():
			break
		npass += 1
		# s -> t -> u in the same pass: s ends up in u
		while not (lut[lut] == lut).all():
			lut = lut[lut]
		data = lut[data]
	print(f"{pa3}: {nsmall} segments < {minarea} cells merged into the "
	      f"{'most similar' if similar else 'most frequent'} larger one within {SMALL_SEGM_RADIUS} pixels in {npass} "
	      f"passes; {len(small)} left without a larger segment so close")

	# Those left join the nearest segment with minarea cells or more, at any distance
	changed = npass > 0
	big = np.isin(data, ids[counts >= minarea])
	if len(small) and big.any():
		dist, (nrow, ncol) = ndimage.distance_transform_edt(~big, return_indices=True)
		lut = np.arange(data.max() + 1, dtype=np.int32)
		farthest = 0
		for s in small:
			box = boxes[s - 1]
			inside = data[box] == s
			d = dist[box][inside]
			k = np.argmin(d)
			lut[s] = data[nrow[box][inside][k], ncol[box][inside][k]]
			farthest = max(farthest, d[k])
		print(f"{pa3}: {len(small)} small segments joined the nearest segment >= {minarea} cells "
		      f"(the farthest at {farthest:.1f} pixels)")
		data = lut[data]
		changed = True
	if not changed:
		return
	out = garray.array(dtype=np.int32)
	out[...] = data
	out.write(mapname=pa3, null=0, overwrite=True)


def fsegm(pa):
	"""Segmentation of area PA at the thresholds of SEGM_THRESHOLDS (0.1-0.9 by default; or up to the fixed one): for
	every threshold k exported (k = 1..n, its position in the list), shp/park_segm_<pa>_<k>.shp and the files of
	segm_files(): one raster of all the segments, which the statistics step reads at once (no map or file per
	segment). Every area is 1000 units of progress: half the i.segment of the thresholds, half the export."""
	pa_list_done = np.genfromtxt(csvname1, dtype=str) if os.path.exists(csvname1) else np.array([])
	if pa in pa_list_done:
		hft_state.progress('DONE', 1000)  # an area done before (every area is 1000 units: see below)
	if pa not in pa_list_done:
		current = multiprocessing.current_process()
		mn = current._identity[0]

		# PERMANENT always exists in a GRASS project; each worker then works in its own mapset
		grass, gsetup = init_grass(GRASSDBASE, MYLOC, 'PERMANENT')

		mapset2 = 'm' + str(mn)
		shutil.rmtree(os.path.join(GRASSDBASE, MYLOC, mapset2), ignore_errors=True)
		col = COL_ID
		grass.run_command('g.mapset', mapset=mapset2, flags='c')

		gisbase = os.environ['GISBASE']
		try:
			gsetup.init(GRASSDBASE, MYLOC, mapset2)
		except TypeError:
			gsetup.init(gisbase, GRASSDBASE, MYLOC, mapset2)
		fix_grass_env(gisbase)

		source = STUDY_AREA
		reps = list(SEGM_THRESHOLDS)
		# Fixed threshold: each threshold starts from the segments of the previous one, so the chain is run up to
		# it (same segments as with all the thresholds), but only the chosen one is converted, exported and summarised
		only = hft_state.fixed_level()
		if only:
			reps = reps[:only]

		# Remove leftovers of a previous interrupted run of this area: its segment rasters and their lists
		for k in range(1, hft_state.MAX_THRESHOLDS + 1):
			for leftover in segm_files(pa, k) + ('tiffs/segm_' + str(pa) + '_' + str(k) + '.tif.aux.xml',):
				if os.path.exists(leftover):
					os.remove(leftover)

		pa44x = 'pax_' + str(pa)
		pa0 = 'v0_' + pa
		opt1 = col + '=' + pa
		grass.run_command('v.extract', input=source, output=pa0, where=opt1, overwrite=True)

		# Región: bounding box del área más el buffer opcional, alineada con la malla de las variables
		# (align: mismos bordes de píxel y resolución que FIRST_VAR, así los píxeles son los de ../inVars).
		# grow va en celdas: PA_BUFFER (unidades del mapa) se pasa a celdas redondeando hacia arriba
		region_kwargs = {'vector': pa0, 'align': FIRST_VAR}
		if PA_BUFFER > 0:
			region_kwargs['grow'] = int(math.ceil(PA_BUFFER / RESOLUTION))

		grass.run_command('g.region', **region_kwargs)

		same = pa + 'v2_= const'
		rndmap = 'rndseed=rand(1,10000000000000000000000000)'
		rndname = 'tiffs/rndseed_' + str(pa) + '.tif'
		grass.run_command('r.mapcalc', expression=f'const = if(isnull({FIRST_VAR}),null(),1)', overwrite=True)
		grass.run_command('r.mapcalc', expression=same, overwrite=True)
		grass.run_command('r.mapcalc', seed=10, expression=rndmap, overwrite=True)
		# -c: without a colour table, which GTiff only takes for Byte or UInt16 bands (GDAL printed an ERROR 6)
		grass.run_command('r.out.gdal', input='rndseed', output=rndname, flags='c', overwrite=True)

		a = grass.read_command('r.stats', input='const', flags='nc', separator='\n').splitlines()
		if len(a) == 0: a = [1, 625]
		minarea = MIN_SEGM_SIZE if MIN_SEGM_SIZE > 0 else int(np.sqrt(int(a[1])))

		# Input of i.segment, which scales each map to its range: axes of the PCA (PC1-PC3, the first SEGM_PCA_AXES, or
		# those that explain SEGM_PCA_VARIANCE % of the variance: pca_axes.py), or copies of the variables on the working
		# region, so that their range is also that of the region and not of the whole map
		if SEGM_PCA:
			out = grass.read_command('i.pca', flags='n', input=ENV_VARS0, output=pa44x, overwrite=True)
			print('Eigen values, (vectors), and [percent importance]:\n' + out.rstrip(), flush=True)
			axes = pca_axes.parse(out)
			# (exceptions, not sys.exit: this runs in a worker of the Pool, and pool.map raises them in the main process)
			if re.search(r'\bnan\b', out, re.I):
				# i.pca divides each variable by its standard deviation: nan when one has the same value in every pixel,
				# and the axes would be empty maps
				raise RuntimeError(f'the PCA of the variables of area {pa} could not be computed (i.pca gives nan): a '
				                   'variable has the same value in every pixel of the area (e.g. all 0)? Take it out of the '
				                   'set, or segment without PCA (SEGM_PCA = False)')
			if len(axes) == len(SEGM_VARS):
				k = pca_axes.choose([p for _, _, p in axes], SEGM_PCA_VARIANCE, SEGM_PCA_AXES)
				pca_axes.write(os.path.join('csv', f'pca_axes_{pa}.csv'), axes, k)
				print(f'Segmentation of {pa}: {pca_axes.summary(axes, k, SEGM_PCA_VARIANCE, SEGM_PCA_AXES)}', flush=True)
			elif SEGM_PCA_VARIANCE and not SEGM_PCA_AXES:
				raise RuntimeError(f'the % of the variance of the axes of the PCA of area {pa} could not be read from i.pca '
				                   f'({len(axes)} axes of {len(SEGM_VARS)}), needed by SEGM_PCA_VARIANCE = '
				                   f'{SEGM_PCA_VARIANCE:g}: set it to 0 (the first 3 axes), or give a number of axes '
				                   '(SEGM_PCA_AXES)')
			else:  # an i.pca that prints its table otherwise: the first 3 axes as before (or SEGM_PCA_AXES), without their %
				k = min(SEGM_PCA_AXES or 3, len(SEGM_VARS))
				print(f'WARNING: the % of the variance of the axes of the PCA of {pa} could not be read from i.pca: the '
				      f'segmentation uses the first {k} axes', flush=True)
			segm_input = [f"{pa44x}.{i}" for i in range(1, k + 1)]
		else:
			segm_input = []
			for i, var in enumerate(SEGM_VARS, 1):
				grass.run_command('r.mapcalc', expression=f'{pa44x}_{i} = {var}', overwrite=True)
				segm_input.append(f"{pa44x}_{i}")
		# CLIP_PERCENTILE: each layer clipped to its percentiles in the working region, so that its range (with
		# which i.segment scales it) is not set by a few extreme pixels; NULL cells stay NULL
		if CLIP_PERCENTILE > 0:
			clipped = []
			for i, layer in enumerate(segm_input, 1):
				q = grass.read_command('r.quantile', input=layer, percentiles=f'{CLIP_PERCENTILE},{100 - CLIP_PERCENTILE}', quiet=True)
				lo, hi = [float(line.split(':')[2]) for line in q.split()]
				out = f'{pa44x}_clip{i}'
				grass.run_command('r.mapcalc', expression=f'{out} = max({lo!r}, min({hi!r}, "{layer}"))', overwrite=True)
				clipped.append(out)
			segm_input = clipped
		grass.run_command('i.group', group='segm', input=','.join(segm_input))

		j = 0
		first_export = True  # the first threshold exported has the most segments (they only grow with the threshold)
		for thr in reps:
			pa2 = pa + 'v2_' + str(j)
			pa2s = pa + 'v2_' + str(j - 1)
			pa3 = pa + 'v3'
			pa4 = 'paa_' + pa

			grass.run_command('g.region', **region_kwargs)
			j += 1

			# the first threshold of the list starts from the random seeds with the minimum size; every next one, from
			# the segments of the previous one (so the segments only grow, and the list itself changes them)
			if j == 1:
				grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', minsize=minarea, similarity='euclidean', memory='10000', iterations='20', seeds='rndseed', overwrite=True)
			else:
				grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', similarity='euclidean', memory='10000', iterations='20', seeds=pa2s, overwrite=True)

			hft_state.progress('DONE', f'{500 / len(reps):.4f}')  # i.segment of the thresholds: half of the area

			if only and j < only:
				continue

			if CLIP_TO_PA:
				grass.run_command('r.mask', vector=source, where=opt1)

			opt2 = pa3 + '=' + pa2
			grass.run_command('r.mapcalc', expression=opt2, overwrite=True)

			if CLIP_TO_PA:
				try:
					grass.run_command('r.mask', flags='r')
				except Exception:
					pass

			merge_small_segments(grass, pa3, minarea, segm_input)

			grass.run_command('r.to.vect', input=pa3, output=pa4, type='area', flags='v', overwrite=True)
			# segm_id = area ID + threshold index k + segment cat, both zero-padded: 1 + 5 digits with the nine thresholds
			# of 0.1-0.9 (as up to eHabpy 1.5), else 2 + 7 (hft_state.segm_id_layout); fixed-width suffix so IDs never
			# collide between thresholds or areas, and are the same on reruns
			maxcat = int(grass.read_command('v.db.select', map=pa4, columns='max(cat)', flags='c').strip() or 0)
			# (the length of the area ID is checked before the segmentation: hft_state.segm_id_error)
			if maxcat >= 10 ** CAT_DIGITS:
				raise ValueError(f"Area {pa}: {maxcat} segments at threshold {thr}, segm_id allows at most {10 ** CAT_DIGITS - 1}"
				                 + (" (with SEGM_THRESHOLDS different from the nine of 0.1-0.9: 9999999)" if CAT_DIGITS == 5 else ''))
			# VARCHAR(254): the width of a text field of the shapefile (with VARCHAR, 1000, v.out.ogr warned that it cut it)
			grass.run_command('v.db.addcolumn', map=pa4, columns='cat_pa VARCHAR(254), segm_id numeric')
			grass.run_command('v.db.update', map=pa4, column='cat_pa', value=pa)
			grass.run_command('v.db.update', map=pa4, column='segm_id',
			                  query_column=f"cat_pa || '{j:0{K_DIGITS}d}' || printf('%0{CAT_DIGITS}d', cat)")

			name = 'shp/park_segm_' + str(pa) + '_' + str(j)
			grass.run_command('v.out.ogr', input=pa4, output=name + '.shp', output_layer=os.path.basename(name), format='ESRI_Shapefile', type='area', overwrite=True)

			# cat (the label of the segment in pa3: r.to.vect -v) and segm_id of every segment
			pairs = sorted({tuple(r.split()) for r in grass.read_command(
				'v.db.select', map=pa4, columns='cat,segm_id', separator='space', flags='c').splitlines() if r.strip()},
				key=lambda p: int(p[0]))
			spa_list = np.unique([sid for _, sid in pairs])  # sorted as text, as before
			parkcsv, ids_csv, segm_tif = segm_files(pa, j)

			# One raster of all the segments of the threshold (the statistics step reads it, ehab_segm_core.py):
			# the labels of pa3 (= cat of the vector), 0 outside them, on the working region (the grid of the
			# variables), masked to the pixels with data in FIRST_VAR with CLIP_TO_PA, as the masks of every
			# segment were (tiffs/pa_<segm_id>.tif, one per segment, made with v.extract and v.to.rast up to 1.5)
			segm_map = 'sgr_' + pa
			if CLIP_TO_PA:
				grass.run_command('r.mask', raster=FIRST_VAR)
			grass.run_command('r.mapcalc', expression=f'{segm_map} = {pa3}', overwrite=True)
			if CLIP_TO_PA:
				try:
					grass.run_command('r.mask', flags='r')
				except Exception:
					pass
			grass.run_command('r.null', map=segm_map, null=0)
			# -c: without a colour table, which GTiff only takes for Byte or UInt16 bands (GDAL printed an ERROR 6)
			grass.run_command('r.out.gdal', input=segm_map, output=segm_tif, type='Int32', createopt='COMPRESS=LZW',
			                  flags='c', overwrite=True)

			# The lists read by the statistics step, each written in one go: segm_id of every segment
			# (csv/park_<id>_<k>.csv), cat -> segm_id (csv/segm_<id>_<k>_ids.csv) and k once per segment
			# (csv/ecoregs<k>.csv, shared by the areas)
			with open(ids_csv, 'w') as wb:
				wb.write('cat segm_id\n' + ''.join(f'{cat} {sid}\n' for cat, sid in pairs))
			if len(spa_list):
				with open(parkcsv, 'w') as wb:
					wb.write(''.join(spa + '\n' for spa in spa_list))
			# the areas run in parallel and append to the same csv/ecoregs<k>.csv: a single write() call (O_APPEND)
			# so that the lines of two areas are never mixed
			fd = os.open('csv/ecoregs' + str(j) + '.csv', os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
			try:
				os.write(fd, f'{j}\n'.encode() * len(spa_list))
			finally:
				os.close(fd)
			print(f"Area {pa}, threshold {thr}: {len(spa_list)} segments -> {segm_tif}")
			# the classification into HFTs compares every pair of segments of the threshold with the most segments: with
			# more than hft_state.LOWMEM_SEGMENTS, with little memory (ward_lowmem.py) but a time that grows with segments²
			if first_export and len(spa_list) > hft_state.LOWMEM_SEGMENTS:
				print(f"Area {pa}: {len(spa_list)} segments at threshold {thr}: the classification into HFTs "
				      f"(moranvar_plots.R) will classify the thresholds with more than {hft_state.LOWMEM_SEGMENTS} segments "
				      f"with little memory (ward_lowmem.py; it needs about {hft_state.class_memory_gb(len(spa_list)):.1f} GB "
				      "at most), in a time that grows with segments² (about a minute with 55,000 segments, four with "
				      "166,000).")
			first_export = False

			# the export of the segments (merge, vectors and segment raster) is the other half of the area,
			# shared by the thresholds that are exported
			hft_state.progress('DONE', f'{500 if only else 500 / len(reps):.4f}')

		shutil.rmtree(os.path.join(GRASSDBASE, MYLOC, mapset2), ignore_errors=True)

		with open(csvname1, 'a') as wb:
			wb.write(str(pa) + '\n')


if __name__ == '__main__':
	if not os.path.isdir(os.path.join(GRASSDBASE, MYLOC, 'PERMANENT')):
		print(f"ERROR: GRASS project not found: {os.path.join(GRASSDBASE, MYLOC, 'PERMANENT')}")
		print("GRASSDB in config.py must be the GRASS database folder that contains your project")
		print("(e.g. '~/grassdata'), not the GRASS installation folder given by 'grass --config path'.")
		sys.exit(1)

	hft_state.check_ready('segm')  # the input data have not changed since setup.py imported them; the options are valid
	# and the segm_id of every area of palist.csv fits (hft_state.options_error), before anything is deleted or run

	# Resume only an interrupted run with the same options, palist.csv and input data (run_state.json);
	# otherwise csv, shp and tiffs (shared by all the study areas) are started from scratch
	rec = hft_state.record('segm')
	if FORCE_RESTART:
		print("\n⚠️ FORCE_RESTART activado: Limpiando resultados anteriores para empezar desde cero...")
		clean_previous_results()
	elif rec is None or rec['chain'] != hft_state.chain('segm'):
		why = 'no record of a previous segmentation' if rec is None else '; '.join(hft_state.reasons('segm'))
		print(f"\n⚠️ Starting the segmentation from scratch ({why})...")
		clean_previous_results()
	elif rec['done']:
		print("The segmentation is up to date with config.py, palist.csv and the input data: nothing to do.\n"
		      "Set FORCE_RESTART = True in config.py to run it again anyway.")
		sys.exit(0)
	else:
		print("\nResuming the interrupted segmentation: the IDs in csv/segm_done.csv are skipped")
	hft_state.start('segm')

	for folder in ['csv', 'shp', 'tiffs']:
		os.makedirs(folder, exist_ok=True)
	with open(csvarea, 'w') as wb:
		wb.write(STUDY_AREA + '\n')
	# The value of every threshold k (csv/segm_thresholds.csv), for the labels of the later steps
	hft_state.write_run_thresholds(SEGM_THRESHOLDS, 'csv/segm_thresholds.csv')

	for f in [csvname1, csvong1, csvong2]:
		if not os.path.isfile(f):
			open(f, 'a').close()

	print("Extracting list of PAs")
	pa_list0 = np.genfromtxt('palist.csv', dtype=str)
	pa_list = np.unique(pa_list0)

	# Without the grass package here, a GRASS error in a worker cannot be passed back and the pool hangs
	add_grass_to_path()
	hft_state.progress('UNITS', len(pa_list) * 1000)  # 1000 units for every area
	pool = Pool(NPROC)
	pool.map(fsegm, pa_list)
	pool.close()
	pool.join()
	hft_state.done('segm')
