#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / GRASS GIS 8 (2024)
#### License: CC BY-SA 3.0
#### Control files: 'csv/segm_done.csv'; ongoing.csv; done.csv
#### Inputs variables: 9 input variables in GRASS GIS; palist.csv;
#### Outputs: 9 segmentation shapefiles for each PA in shp folder; raster files based on segments for all PAs in tiffs folder; segmentation (ecoregs) and park segments csv files in csv folder;
#### NOtes: configure the GRASS GIS environment/database below (or via the
####        GISBASE / GRASSDBASE / GRASSLOC environment variables) and the number
####        of processors (Pool(n) near the bottom of the file).

from multiprocessing import cpu_count, Pool, Lock
import multiprocessing
import subprocess
from datetime import datetime
import numpy as np
import os
import sys
import csv
import gc

# ----------------------------------------------------------------------------
# GRASS GIS 8 configuration.
#
# Set these to match your installation. They can also be provided through
# environment variables so the script does not need editing:
#   GISBASE    -> GRASS installation dir (``grass --config path``)
#   GRASSDBASE -> GRASS database (GISDBASE) directory
#   GRASSLOC   -> location/project name (Mollweide global location)
# ----------------------------------------------------------------------------
GRASSDBASE = os.environ.get('GRASSDBASE', os.path.expanduser('/Users/javier/grassdata')) #'~/grassdata/ehabgrassdb'))
MYLOC = os.environ.get('GRASSLOC', 'global_MW')
NPROC = int(os.environ.get('EHAB_NPROC', max(1, cpu_count() - 1))) # '12' # was Pool(2); 9 in production


def find_gisbase():
	"""Locate the GRASS installation (GISBASE)."""
	gb = os.environ.get('GISBASE')
	if gb:
		return gb
	for exe in ('grass', 'grass84', 'grass83', 'grass82', 'grass80'):
		try:
			return subprocess.check_output([exe, '--config', 'path'], text=True).strip()
		except Exception:
			continue
	return None


def init_grass(gisdbase, location, mapset, create_mapset=False):
	"""Start a GRASS GIS 8 session and return the grass.script module.

	Falls back to the GRASS 7 ``init`` signature if needed.
	"""
	gisbase = find_gisbase()
	if gisbase:
		os.environ['GISBASE'] = gisbase
		py = os.path.join(gisbase, 'etc', 'python')
		if py not in sys.path:
			sys.path.append(py)
	import grass.script as grass
	import grass.script.setup as gsetup
	try:
		# GRASS 8 signature: init(gisdbase, location, mapset)
		gsetup.init(gisdbase, location, mapset)
	except TypeError:
		# GRASS 7 signature: init(gisbase, gisdbase, location, mapset)
		gsetup.init(gisbase, gisdbase, location, mapset)
	return grass, gsetup


print("Extracting list of PAs")
pa_list0 = np.genfromtxt('palist.csv', dtype=str)
pa_list = np.unique(pa_list0)
print(pa_list)

csvname1 = 'csv/segm_done.csv'
csvong1 = 'ongoing.csv'
csvong2 = 'done.csv'
if os.path.isfile(csvname1) == False:
	os.system('touch ' + str(csvname1))
if os.path.isfile(csvong1) == False:
	os.system('touch ' + str(csvong1))
if os.path.isfile(csvong2) == False:
	os.system('touch ' + str(csvong2))


def fsegm(pa):

	pa_list_done = np.genfromtxt(csvname1, dtype=str)
	if pa not in pa_list_done:
		current = multiprocessing.current_process()
		mn = current._identity[0]
		print('running:', mn)

		mapset = 'm'
		grass, gsetup = init_grass(GRASSDBASE, MYLOC, mapset)

		mapset2 = 'm' + str(mn)  # ehabitat'
		os.system('rm -rf ' + os.path.join(GRASSDBASE, MYLOC, mapset2))
		col = 'wdpaid'
		grass.run_command('g.mapset', mapset=mapset2, project=MYLOC, dbase=GRASSDBASE, flags='c')
		os.system('rm csv/*_' + str(pa) + '_*')
		os.system('rm shp/*_' + str(pa) + '_*')

		# Re-init the session on the freshly created worker mapset.
		try:
			gsetup.init(GRASSDBASE, MYLOC, mapset2)
		except TypeError:
			gsetup.init(os.environ['GISBASE'], GRASSDBASE, MYLOC, mapset2)
		print(pa, mapset2, grass.gisenv())
		ong = str(pa) + str(mapset2) + str(grass.gisenv())
		grass.run_command('g.mapsets', mapset='rasterized_parks', operation='add') # ehabplus_cs,javier
		source = 'wdpa_snapshot_mollweide' # 'cspas'  # 'wdpa_aug14_100km2_moll'

		wb = open(csvong1, 'a')
		wb.write(ong)
		wb.write('\n')
		wb.close()

		reps = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

		print(pa)
		pa44 = 'pa_' + str(pa)
		pa44x = 'pax_' + str(pa)
		pa0 = 'v0_' + pa
		opt1 = col + '=' + pa
		grass.run_command('v.extract', input=source, output=pa0, where=opt1, overwrite=True)  # check inital region from which to copy from!
		pa2 = pa + 'v2_'
		pa3 = pa + 'v3'
		pa4 = 'paa_' + pa
		pa5 = pa4 + '.txt'
		same = pa2 + '= const'
		rndmap = 'rndseed=rand(1,10000000000000000000000000)'
		rndname = 'tiffs/rndseed_' + str(pa) + '.tif'
		grass.run_command('g.region', vector=pa0, res=1000)
		grass.run_command('r.mapcalc', expression='const = if(gcmask>=0,1,null())', overwrite=True)
		grass.run_command('r.mapcalc', expression=same, overwrite=True)
		grass.run_command('r.mapcalc', seed=10, expression=rndmap, overwrite=True)
		grass.run_command('r.out.gdal', input='rndseed', output=rndname, overwrite=True)
		a = grass.read_command('r.stats', input='const', flags='nc', separator='\n').splitlines()
		if len(a) == 0: a = [1, 625]
		minarea = int(np.sqrt(int(a[1])))  # /2 #10
		minaream = minarea  # *1000
		grass.run_command('i.pca', flags='n', input='pre,epr,slope,tree,herb,ndwi,ndvi,ndvi_range,bio', output=pa44x, overwrite=True)  # dem
		pca1 = pa44x + '.1'
		pca2 = pa44x + '.2'
		pca3 = pa44x + '.3'
		pcas = pca1 + ',' + pca2 + ',' + pca3
		grass.run_command('i.group', group='segm', input=pcas)
		os.system('cat ' + os.path.join(GRASSDBASE, MYLOC, mapset2, 'group/segm/REF'))
		j = 0
		for thr in reps:
			pa2 = pa + 'v2_' + str(j)
			pa2s = pa + 'v2_' + str(j - 1)
			aleat = np.random.randint(1, 1001)
			grass.run_command('g.region', vector=pa0, res=1000)
			j = j + 1
			if thr == 0.1:
				grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', minsize=minarea, similarity='euclidean', memory='10000', iterations='20', seeds='rndseed', overwrite=True)  # ,seed=pa2i minsize=minarea,
			else:
				grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', similarity='euclidean', memory='10000', iterations='20', seeds=pa2s, overwrite=True)  # minsize=minarea
			grass.run_command('r.mask', vector=source, where=opt1)
			opt2 = pa3 + '=' + pa2
			grass.run_command('r.mapcalc', expression=opt2, overwrite=True)  # usar const como mapa para crear plantilla de PA con unos y ceros
			#grass.run_command('r.mask', flags='r')  # drop the mask (was: g.rename MASK,masc)
			try:
			    grass.run_command('r.mask', flags='r')
			except Exception:
			    pass  
			print('minarea: ', minarea)

			b = grass.read_command('r.stats', input=pa3, flags='nc', separator='\n').splitlines()
			print(b)
			clean = None
			c = pa3
			for g in np.arange(1, len(b), 2):
				if int(b[g]) < minarea:  # /10: # lower the threshold if omitting min area!
					print('Cleaning small segments I...')
					print('cleaning cat ' + str(b[g - 1]))
					c2 = 'old' + str(b[g - 1])
					c22 = c2 + 'b10km'
					c3 = 'new' + str(b[g - 1])
					oper1 = c2 + '=' + 'if(' + pa3 + '==' + str(b[g - 1]) + ',1,null())'
					grass.run_command('r.mapcalc', expression=oper1, overwrite=True)
					grass.run_command('r.buffer', input=c2, output=c22, distances=3, units='kilometers', overwrite=True)
					grass.run_command('r.mask', raster=c22, maskcats='2')
					buff = grass.read_command('r.stats', input=pa3, flags='nc', sort='desc', separator='\n').splitlines()
					#grass.run_command('r.mask', flags='r')
					try:
					    grass.run_command('r.mask', flags='r')
					except Exception:
					    pass  
					if len(buff) > 0:
						clean = 'T'
						print('New: ' + str(buff[0]))
						oper1 = c3 + '=' + 'if(' + c2 + '==1,' + str(buff[0]) + ',null())'
						c = c3 + ',' + c
						grass.run_command('r.mapcalc', expression=oper1, overwrite=True)
			if clean == 'T':
				print(c)
				grass.run_command('r.patch', input=c, output=pa3, overwrite=True)
				bv = grass.read_command('r.stats', input=pa3, flags='nc', separator='\n').splitlines()
				print(bv)

			b = grass.read_command('r.stats', input=pa3, flags='nc', separator='\n').splitlines()
			print(b)
			clean = None
			c = pa3
			for g in np.arange(1, len(b), 2):
				if int(b[g]) < minarea:  # /10: # lower the threshold if omitting min area!
					print('Cleaning small segments II...')
					print('cleaning cat ' + str(b[g - 1]))
					c2 = 'old' + str(b[g - 1])
					c22 = c2 + 'b10km'
					c3 = 'new' + str(b[g - 1])
					oper1 = c2 + '=' + 'if(' + pa3 + '==' + str(b[g - 1]) + ',1,null())'
					grass.run_command('r.mapcalc', expression=oper1, overwrite=True)
					grass.run_command('r.buffer', input=c2, output=c22, distances=10, units='kilometers', overwrite=True)
					grass.run_command('r.mask', raster=c22, maskcats='2')
					buff = grass.read_command('r.stats', input=pa3, flags='nc', sort='desc', separator='\n').splitlines()
					#grass.run_command('r.mask', flags='r')
					try:
					    grass.run_command('r.mask', flags='r')
					except Exception:
					    pass  
					if len(buff) > 0:
						clean = 'T'
						print('New: ' + str(buff[0]))
						oper1 = c3 + '=' + 'if(' + c2 + '==1,' + str(buff[0]) + ',null())'
						c = c3 + ',' + c
						grass.run_command('r.mapcalc', expression=oper1, overwrite=True)
			if clean == 'T':
				print(c)
				grass.run_command('r.patch', input=c, output=pa3, overwrite=True)
				bv = grass.read_command('r.stats', input=pa3, flags='nc', separator='\n').splitlines()
				print(bv)

			b = grass.read_command('r.stats', input=pa3, flags='nc', sort='desc', separator='\n').splitlines()
			print(b)
			for g in np.arange(1, len(b), 2):
				if int(b[g]) < minarea:  # /10: # lower the threshold if omitting min area!
					print('Cleaning small segments III...')
					print('cleaning cat ' + str(b[g - 1]))
					oper1 = pa3 + '=' + 'if(' + pa3 + '==' + str(b[g - 1]) + ',' + str(b[0]) + ',' + pa3 + ')'
					grass.run_command('r.mapcalc', expression=oper1, overwrite=True)
					bv = grass.read_command('r.stats', input=pa3, flags='nc', separator='\n').splitlines()
					print(bv)

			grass.run_command('r.to.vect', input=pa3, output=pa4, type='area', flags='v', overwrite=True)
			grass.run_command('v.db.addcolumn', map=pa4, columns='wdpaid_pa VARCHAR')
			grass.run_command('v.db.update', map=pa4, column='wdpaid_pa', value=pa)
			grass.run_command('v.db.addcolumn', map=pa4, columns='aleat VARCHAR')
			grass.run_command('v.db.update', map=pa4, column='aleat', value=aleat)
			pa44 = pa4
			pa442 = pa44 + '_diss'
			grass.run_command('v.db.addcolumn', map=pa44, columns='segm_id numeric')  # VARCHAR')
			grass.run_command('v.db.update', map=pa44, column='segm_id', query_column='wdpaid_pa || cat || aleat')
			name = 'shp/park_segm_' + str(pa) + '_' + str(j)
			if os.path.isfile(name + '.shp') == False:
				grass.run_command('v.out.ogr', input=pa44, output=name + '.shp', output_layer=os.path.basename(name), format='ESRI_Shapefile', type='area')
			else:
				grass.run_command('v.out.ogr', flags='a', input=pa44, output=name + '.shp', output_layer=os.path.basename(name), format='ESRI_Shapefile', type='area')

			grass.message("Done1")
			spa_list0 = grass.read_command('v.db.select', map=pa44, column='segm_id').splitlines()
			spa_list = np.unique(spa_list0)
			print(spa_list)
			# save it as a csv excluding last item!

			grass.message("omitting previous masks")
			#grass.run_command('r.mask', flags='r')
			try:
			    grass.run_command('r.mask', flags='r')
			except Exception:
			    pass  
			sn = len(spa_list) - 1  # there is also a segm_id element!
			for spx in range(0, sn):  # 0
				spa = spa_list[spx]
				spa2 = 'svv' + spa
				spa4 = 'spa_' + spa
				spa5 = 'tiffs/pa_' + spa + '.tif'
				spa0 = 'sv0_' + spa
				sopt1 = 'segm_id = ' + spa
				print(spx)
				print("Extracting PA:" + spa)
				grass.run_command('v.extract', input=pa44, output=spa0, where=sopt1, overwrite=True)
				# try to crop PAs shapefile with coastal line or input vars
				grass.message("setting up the working region")
				grass.run_command('g.region', vector=spa0, res=1000)
				grass.run_command('v.to.rast', input=spa0, output=spa0, use='val')  # use='cat',labelcol='segm_id')
				soptt = spa4 + '=' + spa0
				grass.run_command('r.mask', raster='pre')  # new to crop parks to where we have indicators information
				grass.run_command('r.mapcalc', expression=soptt, overwrite=True)  # opt3
				#grass.run_command('r.mask', flags='r')
				try:
				    grass.run_command('r.mask', flags='r')
				except Exception:
				    pass  
				grass.run_command('r.null', map=spa4, null=0)
				econame = 'csv/park_' + str(pa) + '_' + str(j) + '.csv'
				eco = str(j)
				econ = 'csv/ecoregs' + str(j) + '.csv'
				grass.run_command('r.out.gdal', input=spa4, output=spa5, overwrite=True)
				wb = open(econame, 'a')
				wb.write(spa)
				wb.write('\n')
				wb.close()
				wb = open(econ, 'a')
				wb.write(eco)
				wb.write('\n')
				wb.close()
		grass.message("Deleting tmp layers")
		os.system('rm -rf ' + os.path.join(GRASSDBASE, MYLOC, mapset2))

		wb = open(csvong2, 'a')
		wb.write(ong)
		wb.write('\n')
		wb.close()

		wb = open(csvname1, 'a')
		var = str(pa)
		wb.write(var)
		wb.write('\n')
		wb.close()


if __name__ == '__main__':
	pool = Pool(NPROC)  # 9
	pool.map(fsegm, pa_list)
	pool.close()
	pool.join()
