#### Author: Javier Martinez-Lopez (UTF-8)
#### Modernized for Python 3 / GRASS GIS 8 (2024)
#### License: CC BY-SA 3.0
#### Builds the per-ecoregion buffered raster masks (eco_<id>.tif) used by ehab_optim.py.
#### NOtes: configure the GRASS GIS environment/database below (or via the
####        GISBASE / GRASSDBASE / GRASSLOC environment variables).

import os
import sys
import subprocess
import numpy as np
from tqdm import tqdm

gisdbase = os.environ.get('GRASSDBASE', os.path.expanduser('~/grassdata/ehabgrassdb'))
location = os.environ.get('GRASSLOC', 'global_MW')
mapset = os.environ.get('GRASSMAPSET', 'ehabitat')  # "rasterized_parks"


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


source = 'ecoregs_moll'
grass.message("Extracting list of ecoregs")
list0 = grass.read_command('v.db.select', map=source, column='eco_id').splitlines()
ecolist = np.unique(list0)
# save it as a csv excluding last item!

grass.message("omitting previous masks")
drop_mask()

n = len(ecolist) - 2
print(ecolist)
for px in tqdm(range(0, n)):  # 0,n

	eco = abs(int(ecolist[px]))
	eco2 = int(ecolist[px])
	eco0 = 'v0_' + str(eco)
	opt1 = 'eco_id = ' + str(eco2)
	print(px)
	print("Extracting ECO:" + str(eco))
	grass.run_command('v.extract', input=source, output=eco0, where=opt1, overwrite=True)
	pa2 = 'vv' + str(eco)
	pa3 = str(eco) + 'v3'
	pa4 = 'eco_' + str(eco)
	pa5 = 'eco_' + str(eco) + '.tif'
	pa6 = str(eco) + 'v4'

	grass.message("setting up the working region")
	grass.run_command('g.region', vector=eco0, res=1000)
	grass.run_command('v.to.rast', input=eco0, output=eco0, use='cat', label_column='eco_id', overwrite=True)
	opt3 = pa2 + '= @' + eco0
	opt4 = pa2 + '= round(' + pa2 + ')'
	grass.run_command('r.mapcalc', expression=opt3, overwrite=True)
	grass.run_command('r.mapcalc', expression=opt4, overwrite=True)
	grass.run_command('r.mask', vector=eco0, where=opt1)
	opt2 = pa4 + '=' + pa2
	grass.run_command('r.mapcalc', expression=opt2, overwrite=True)
	drop_mask()
	grass.run_command('g.region', flags='d')
	grass.run_command('r.mask', raster='pre')
	grass.run_command('r.buffer', input=pa4, output=pa3, distances=250, units='kilometers', overwrite=True)
	grass.run_command('g.region', zoom=pa3)
	optk = pa6 + '= if(' + pa3 + '!=0,1,0)'
	grass.run_command('r.mapcalc', expression=optk, overwrite=True)
	grass.run_command('r.null', map=pa6, null=0)
	grass.run_command('r.out.gdal', input=pa6, output=pa5, overwrite=True)
	drop_mask()
	grass.message("Deleting tmp layers")
	grass.run_command('g.remove', type='raster', pattern='*v3', flags='f')
	grass.run_command('g.remove', type='raster', pattern='*v2', flags='f')
	grass.run_command('g.remove', type='raster', pattern='v0_*', flags='f')
	grass.run_command('g.remove', type='vector', pattern='v0_*', flags='f')
	grass.run_command('g.remove', type='raster', pattern='vv*', flags='f')
	grass.message("Done")
	print("Done ECO:" + str(eco))
print("FINISHED")
