import os
import sys
import numpy as np
from tqdm import *

gisbase = os.environ['GISBASE'] = "/usr/local/grass-7.0.svn"
gisdbase = os.path.join("/local1/majavie/hanksgrass7")
location = "global_MW"
mapset   = "rasterized_parks"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as grass
import grass.script.setup as gsetup
 
gsetup.init(gisbase,
            gisdbase, location, mapset)
 
print grass.gisenv()

source = 'ecoregs_moll'
grass. message ("Extracting list of ecoregs")
pa_list0 = grass. read_command ('v.db.select', map=source,column='eco_id'). splitlines ()
pa_list = np.unique(pa_list0)
# save it as a csv excluding last item!

grass. message("omitting previous masks")
grass.run_command('g.remove', rast='MASK')

n = len(pa_list)-2
#print pa_list[1]
for px in tqdm(range(0,n)):

#for pa in pa_list:
 pa = pa_list[px]
 pa0 = 'v0_'+pa
 opt1 = 'wdpa_id = '+pa
 print px
 print "Extracting PA:"+pa
 grass.run_command('v.extract', input=source, out=pa0, where = opt1,overwrite=True)
 pa2 = 'vv'+pa
 #pa3 = pa+'v3'
 pa4 = 'pa_'+pa
 pa5 = 'pa_'+pa+'.tif'
 
 grass. message ("setting up the working region")
 grass.run_command('g.region',vect=pa0,res=1000)
 grass.run_command('v.to.rast',input=pa0,out=pa0,use='cat',labelcol='wdpa_id')
 opt3 = pa2+'= @'+pa0
 opt4 = pa2+'= round('+pa2+')'
 grass.run_command('r.mapcalc',expression=opt3,overwrite=True)
 grass.run_command('r.mapcalc',expression=opt4,overwrite=True)
 grass.run_command('r.mask', vector=pa0, where=opt1)
 opt2 = pa4+'='+pa2
 grass.run_command('r.mapcalc',expression=opt2,overwrite=True)
 grass.run_command('g.remove', rast='MASK')
 grass.run_command('r.null',map=pa4,null=0)
 grass.run_command('r.out.gdal',input=pa4,out=pa5)
 grass. message ("Deleting tmp layers") 
 grass.run_command('g.mremove',rast='*v3',flags='f') 
 grass.run_command('g.mremove',rast='*v2',flags='f') 
 grass.run_command('g.mremove',rast='v0_*',flags='f') 
 grass.run_command('g.mremove',vect='v0_*',flags='f') 
 grass.run_command('g.mremove',rast='vv*',flags='f') 
 grass. message ("Done")
 print "Done PA:"+pa

# try to paralellize it?