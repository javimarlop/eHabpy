#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / GRASS GIS 8 (2024)
#### License: CC BY-SA 3.0

from multiprocessing import cpu_count, Pool
import multiprocessing
import subprocess
import numpy as np
import os
import sys
import shutil

# Importar configuración global de variables y opciones
try:
    from config import ENV_VARS0, RESOLUTION, COL_ID, GRASSDB, GRASSLOC, STUDY_AREA
except ImportError:
    print("ERROR: No se encuentra config.py. Asegúrate de crearlo en el mismo directorio.")
    sys.exit(1)

# OPCIONES DE CONFIGURACIÓN
CLIP_TO_PA = getattr(sys.modules['config'], 'CLIP_TO_PA', False)
FORCE_RESTART = getattr(sys.modules['config'], 'FORCE_RESTART', False)
PA_BUFFER = getattr(sys.modules['config'], 'PA_BUFFER', 0)

GRASSDBASE = os.environ.get('GRASSDBASE', os.path.expanduser(GRASSDB))
MYLOC = os.environ.get('GRASSLOC', GRASSLOC)
NPROC = int(os.environ.get('EHAB_NPROC', max(1, cpu_count() - 1)))


def find_gisbase():
	gb = os.environ.get('GISBASE')
	if gb:
		return gb
	for exe in ('grass', 'grass84', 'grass83', 'grass82', 'grass80'):
		try:
			return subprocess.check_output([exe, '--config', 'path'], text=True).strip()
		except Exception:
			continue
	return None


def init_grass(gisdbase, location, mapset):
	gisbase = find_gisbase()
	if gisbase:
		os.environ['GISBASE'] = gisbase
		py = os.path.join(gisbase, 'etc', 'python')
		if py not in sys.path:
			sys.path.append(py)
	import grass.script as grass
	import grass.script.setup as gsetup
	try:
		gsetup.init(gisdbase, location, mapset)
	except TypeError:
		gsetup.init(gisbase, gisdbase, location, mapset)
	return grass, gsetup


def clean_previous_results():
	"""Borra archivos de control y limpia el contenido de las carpetas de salida."""
	print("\n⚠️ FORCE_RESTART activado: Limpiando resultados anteriores para empezar desde cero...")
	
	control_files = ['csv/segm_done.csv', 'ongoing.csv', 'done.csv']
	for f in control_files:
		if os.path.exists(f):
			try:
				os.remove(f)
			except Exception as e:
				print(f"No se pudo eliminar {f}: {e}")

	directories = ['csv', 'shp', 'tiffs', 'results']
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


def fsegm(pa):
	pa_list_done = np.genfromtxt(csvname1, dtype=str) if os.path.exists(csvname1) else np.array([])
	if pa not in pa_list_done:
		current = multiprocessing.current_process()
		mn = current._identity[0]

		mapset = 'm'
		grass, gsetup = init_grass(GRASSDBASE, MYLOC, mapset)

		mapset2 = 'm' + str(mn)
		os.system('rm -rf ' + os.path.join(GRASSDBASE, MYLOC, mapset2))
		col = COL_ID
		grass.run_command('g.mapset', mapset=mapset2, project=MYLOC, dbase=GRASSDBASE, flags='c')

		try:
			gsetup.init(GRASSDBASE, MYLOC, mapset2)
		except TypeError:
			gsetup.init(os.environ['GISBASE'], GRASSDBASE, MYLOC, mapset2)

		source = STUDY_AREA
		reps = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

		pa44x = 'pax_' + str(pa)
		pa0 = 'v0_' + pa
		opt1 = col + '=' + pa
		grass.run_command('v.extract', input=source, output=pa0, where=opt1, overwrite=True)

		# Configuración de región con buffer opcional alrededor del Bounding Box
		region_kwargs = {'vector': pa0, 'res': RESOLUTION}
		if PA_BUFFER > 0:
			region_kwargs['grow'] = PA_BUFFER
		
		grass.run_command('g.region', **region_kwargs)

		same = pa + 'v2_= const'
		rndmap = 'rndseed=rand(1,10000000000000000000000000)'
		rndname = 'tiffs/rndseed_' + str(pa) + '.tif'
		grass.run_command('r.mapcalc', expression='const = if(precip>=0,1,null())', overwrite=True) # change to first variable used
		grass.run_command('r.mapcalc', expression=same, overwrite=True)
		grass.run_command('r.mapcalc', seed=10, expression=rndmap, overwrite=True)
		grass.run_command('r.out.gdal', input='rndseed', output=rndname, overwrite=True)

		a = grass.read_command('r.stats', input='const', flags='nc', separator='\n').splitlines()
		if len(a) == 0: a = [1, 625]
		minarea = int(np.sqrt(int(a[1])))

		grass.run_command('i.pca', flags='n', input=ENV_VARS0, output=pa44x, overwrite=True)
		pcas = f"{pa44x}.1,{pa44x}.2,{pa44x}.3"
		grass.run_command('i.group', group='segm', input=pcas)

		j = 0
		for thr in reps:
			pa2 = pa + 'v2_' + str(j)
			pa2s = pa + 'v2_' + str(j - 1)
			pa3 = pa + 'v3'
			pa4 = 'paa_' + pa
			aleat = np.random.randint(1, 1001)

			grass.run_command('g.region', **region_kwargs)
			j += 1

			if thr == 0.1:
				grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', minsize=minarea, similarity='euclidean', memory='10000', iterations='20', seeds='rndseed', overwrite=True)
			else:
				grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', similarity='euclidean', memory='10000', iterations='20', seeds=pa2s, overwrite=True)

			if CLIP_TO_PA:
				grass.run_command('r.mask', vector=source, where=opt1)

			opt2 = pa3 + '=' + pa2
			grass.run_command('r.mapcalc', expression=opt2, overwrite=True)

			if CLIP_TO_PA:
				try:
					grass.run_command('r.mask', flags='r')
				except Exception:
					pass

			b = grass.read_command('r.stats', input=pa3, flags='nc', separator='\n').splitlines()
			clean = None
			c = pa3
			for g in np.arange(1, len(b), 2):
				if int(b[g]) < minarea:
					c2 = 'old' + str(b[g - 1])
					c22 = c2 + 'b10km'
					c3 = 'new' + str(b[g - 1])
					oper1 = c2 + '=' + 'if(' + pa3 + '==' + str(b[g - 1]) + ',1,null())'
					grass.run_command('r.mapcalc', expression=oper1, overwrite=True)
					grass.run_command('r.buffer', input=c2, output=c22, distances=3, units='kilometers', overwrite=True)
					grass.run_command('r.mask', raster=c22, maskcats='2')
					buff = grass.read_command('r.stats', input=pa3, flags='nc', sort='desc', separator='\n').splitlines()
					try:
						grass.run_command('r.mask', flags='r')
					except Exception:
						pass
					if len(buff) > 0:
						clean = 'T'
						oper1 = c3 + '=' + 'if(' + c2 + '==1,' + str(buff[0]) + ',null())'
						c = c3 + ',' + c
						grass.run_command('r.mapcalc', expression=oper1, overwrite=True)

			if clean == 'T':
				grass.run_command('r.patch', input=c, output=pa3, overwrite=True)

			grass.run_command('r.to.vect', input=pa3, output=pa4, type='area', flags='v', overwrite=True)
			grass.run_command('v.db.addcolumn', map=pa4, columns='cat_pa VARCHAR, aleat VARCHAR, segm_id numeric')
			grass.run_command('v.db.update', map=pa4, column='cat_pa', value=pa)
			grass.run_command('v.db.update', map=pa4, column='aleat', value=aleat)
			grass.run_command('v.db.update', map=pa4, column='segm_id', query_column='cat_pa || cat || aleat')

			name = 'shp/park_segm_' + str(pa) + '_' + str(j)
			grass.run_command('v.out.ogr', input=pa4, output=name + '.shp', output_layer=os.path.basename(name), format='ESRI_Shapefile', type='area', overwrite=True)

			spa_list0 = grass.read_command('v.db.select', map=pa4, column='segm_id').splitlines()
			spa_list = np.unique(spa_list0)
			sn = len(spa_list) - 1
			econ = 'csv/ecoregs' + str(j) + '.csv'
			
			if not os.path.exists(econ):
				open(econ, 'a').close()

			for spx in range(0, sn):
				spa = spa_list[spx]
				if spa == 'segm_id':
					continue
				spa2 = 'svv' + spa
				spa4 = 'spa_' + spa
				spa5 = 'tiffs/pa_' + spa + '.tif'
				spa0 = 'sv0_' + spa
				sopt1 = 'segm_id = ' + spa
				print(spx)
				print("Extracting PA:" + spa)

				grass.run_command('v.extract', input=pa4, output=spa0, where=sopt1, overwrite=True)
				grass.message("setting up the working region")
				grass.run_command('g.region', vector=spa0, res=RESOLUTION)
				grass.run_command('v.to.rast', input=spa0, output=spa0, use='val')
				soptt = spa4 + '=' + spa0

				if CLIP_TO_PA:
					grass.run_command('r.mask', raster='precip')
				
				grass.run_command('r.mapcalc', expression=soptt, overwrite=True)

				if CLIP_TO_PA:
					try:
						grass.run_command('r.mask', flags='r')
					except Exception:
						pass

				grass.run_command('r.null', map=spa4, null=0)
				econame = 'csv/park_' + str(pa) + '_' + str(j) + '.csv'
				eco = str(j)

				grass.run_command('r.out.gdal', input=spa4, output=spa5, overwrite=True)

				with open(econame, 'a') as wb:
					wb.write(spa + '\n')

				with open(econ, 'a') as wb:
					wb.write(eco + '\n')

		os.system('rm -rf ' + os.path.join(GRASSDBASE, MYLOC, mapset2))

		with open(csvname1, 'a') as wb:
			wb.write(str(pa) + '\n')


if __name__ == '__main__':
	if FORCE_RESTART:
		clean_previous_results()

	for f in [csvname1, csvong1, csvong2]:
		if not os.path.isfile(f):
			open(f, 'a').close()

	print("Extracting list of PAs")
	pa_list0 = np.genfromtxt('palist.csv', dtype=str)
	pa_list = np.unique(pa_list0)

	pool = Pool(NPROC)
	pool.map(fsegm, pa_list)
	pool.close()
	pool.join()
