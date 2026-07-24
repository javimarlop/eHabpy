#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / GRASS GIS 8 (2024)
#### License: CC BY-SA 3.0
#### Reqired parameters: class shapefile name without extension ('python subpas_loop_optim.py shapefilename')
#### Control files: pas_segm_tiff_done.csv
#### Inputs variables: ecoregs_moll and pre in GRASS GIS
#### Outputs: single HFT raster files to be used by ehab_optim.py; single ecoregion csv files that contain the HFTs; ecoregs.csv (list of ecoregions processed)
#### NOtes: configure the GRASS GIS environment/database below (or via the
####        GISBASE / GRASSDBASE / GRASSLOC environment variables).

import os
import sys
import subprocess
import numpy as np
from tqdm import tqdm

# ----------------------------------------------------------------------------
# GRASS GIS 8 configuration (see segmentation_pca_par.py for details).
# ----------------------------------------------------------------------------
gisdbase = os.environ.get('GRASSDBASE', os.path.expanduser('~/grassdata/ehabgrassdb'))
location = os.environ.get('GRASSLOC', 'global_MW')
mapset = os.environ.get('GRASSMAPSET', 'm')  # "rasterized_parks"


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


gisbase = find_gisbase()
if gisbase:
	os.environ['GISBASE'] = gisbase
	_py = os.path.join(gisbase, 'etc', 'python')
	if _py not in sys.path:
		sys.path.append(_py)

import grass.script as grass
import grass.script.setup as gsetup

try:
	gsetup.init(gisdbase, location, mapset)
except TypeError:
	gsetup.init(gisbase, gisdbase, location, mapset)

print(grass.gisenv())


def drop_mask():
	"""Remove the current raster mask, ignoring the 'no mask present' error."""
	try:
		grass.run_command('r.mask', flags='r')
	except Exception:
		pass


optres = np.genfromtxt('../segm_optim/results/overall_optim_thresholds.csv', delimiter=" ")
nrpas = optres.shape[0]
for el in tqdm(range(0, nrpas)):
	pa_id = int(optres[el][0])
	ot = int(optres[el][1] * 10)
	source = 'park_segm_' + str(pa_id) + '_' + str(ot) + '_class'
	print(source)
	grass.run_command('v.in.ogr', flags='oe', input='../segm_optim/results/', layer=source, output=source, overwrite=True)
	grass.message("Extracting list of HFTs")
	pa_list0 = grass.read_command('v.db.select', map=source, column='hclst_m').splitlines()
	pa_list = np.unique(pa_list0)
	print(pa_list)
	# save it as a csv excluding last item!

	grass.message("Deleting tmp layers")
	grass.run_command('g.remove', type='raster', pattern='*v3', flags='f')
	grass.run_command('g.remove', type='raster', pattern='*v2', flags='f')
	grass.run_command('g.remove', type='raster', pattern='v0_*', flags='f')
	grass.run_command('g.remove', type='raster', pattern='v0_*', flags='f')
	grass.run_command('g.remove', type='raster', pattern='vv*', flags='f')

	grass.message("omitting previous masks")
	drop_mask()

	csvname1 = 'pas_segm_tiff_done.csv'
	if os.path.isfile(csvname1) == False:
		wb = open(csvname1, 'a')
		wb.write('None')
		wb.write('\n')
		wb.close()

	pa_list_done = np.genfromtxt(csvname1, dtype=str)
	n = len(pa_list) - 1  # there is also a segm_id element!
	for px in tqdm(range(0, n)):  # 0
		pa = pa_list[px]
		paf = str(pa_id) + '_' + pa
		pa2 = 'vv' + pa
		pa4 = 'pa_' + pa
		pa5 = 'pa_' + paf + '.tif'
		pa0 = 'v0_' + pa
		opt1 = 'hclst_m = ' + pa
		if pa not in pa_list_done:
			print(px)
			print("Extracting PA:" + pa)
			grass.run_command('v.extract', input=source, output=pa0, where=opt1, overwrite=True)
			# try to crop PAs shapefile with coastal line or input vars
			grass.message("setting up the working region")
			grass.run_command('g.region', vector=pa0, res=1000)
			grass.run_command('v.to.rast', input=pa0, output=pa0, use='val', value=5)  # use='cat',labelcol='segm_id')
			optt = pa4 + '=' + pa0
			grass.run_command('r.mask', raster='pre')  # new to crop parks to where we have indicators information
			grass.run_command('r.mapcalc', expression=optt, overwrite=True)  # opt3
			drop_mask()  # new
			grass.run_command('r.null', map=pa4, null=0)
			eco_list = grass.read_command('r.stats', input='ecoregs_moll', sort='desc').splitlines()
			print(eco_list)
			eco = eco_list[0]
			if eco == '*' or eco == '-9999' or eco == '-9998':
				if len(eco_list) > 1: eco = eco_list[1]
			if eco == '*' or eco == '-9999' or eco == '-9998':
				grass.run_command('g.region', res=10)
				eco_list = grass.read_command('r.stats', input='ecoregs_moll', sort='desc').splitlines()
				eco = eco_list[0]
				grass.run_command('g.region', res=1000)
			if eco == '*' or eco == '-9999' or eco == '-9998':
				c22 = pa0 + 'b50km'
				grass.run_command('g.region', flags='d')
				grass.run_command('r.buffer', input=pa0, output=c22, distances=50, units='kilometers', overwrite=True)
				grass.run_command('g.region', zoom=c22, res=1000)
				grass.run_command('r.mask', raster=c22, maskcats='2')
				eco_list = grass.read_command('r.stats', input='ecoregs_moll', sort='desc').splitlines()
				eco = eco_list[0]
				print(eco_list)
				if eco == '*' or eco == '-9999' or eco == '-9998':
					eco = 'noterreco'
					if len(eco_list) > 1: eco = eco_list[1]
				drop_mask()  # new
			print('eco: ' + eco)
			econame = str(eco) + '.csv'
			grass.run_command('g.region', res=1000)
			grass.run_command('r.out.gdal', input=pa4, output=pa5, overwrite=True)
			grass.message("Deleting tmp layers")
			grass.run_command('g.remove', type='raster', pattern='*v3', flags='f')
			grass.run_command('g.remove', type='raster', pattern='*v2', flags='f')
			grass.run_command('g.remove', type='raster', pattern='v0_*', flags='f')
			grass.run_command('g.remove', type='raster', pattern='v0_*', flags='f')
			grass.run_command('g.remove', type='raster', pattern='vv*', flags='f')
			grass.run_command('g.remove', type='raster', pattern='pa_*', flags='f')
			grass.run_command('g.remove', type='raster', pattern='*b50km', flags='f')
			wb = open(econame, 'a')
			wb.write(paf)
			wb.write('\n')
			wb.close()
			wb = open('ecoregs.csv', 'a')
			wb.write(eco)
			wb.write('\n')
			wb.close()
			wb = open(csvname1, 'a')
			var = str(paf)
			wb.write(var)
			wb.write('\n')
			wb.close()
			grass.message("Done")
			print("Done PA:" + pa)

# try to paralellize it?
