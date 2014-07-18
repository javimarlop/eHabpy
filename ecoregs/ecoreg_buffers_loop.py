import os
import sys
import numpy as np
from tqdm import *

#gisbase = os.environ['GISBASE'] = "/usr/local/grass-7.0.svn"
gisbase = os.environ['GISBASE'] = "/home/majavie/grass_last/new/grass7_trunk/dist.x86_64-unknown-linux-gnu"
gisdbase = os.path.join("/local1/majavie/hanksgrass7")
location = "global_MW"
mapset   = "ehabitat"#"rasterized_parks"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as grass
import grass.script.setup as gsetup
 
gsetup.init(gisbase,
            gisdbase, location, mapset)
 
print grass.gisenv()

source = 'ecoregs_moll'
grass. message ("Extracting list of ecoregs")
list0 = grass. read_command ('v.db.select', map=source,column='eco_id'). splitlines ()
list = np.unique(list0)
# save it as a csv excluding last item!

grass. message("omitting previous masks")
grass.run_command('g.remove', rast='MASK')

n = len(list)-2
print list
for px in tqdm(range(2,n)):

#for pa in pa_list:
 eco = abs(int(list[px]))
 eco0 = 'v0_'+str(eco)
 opt1 = 'eco_id = '+str(eco)
 print px
 print "Extracting ECO:"+str(eco)
 grass.run_command('v.extract', input=source, out=eco0, where = opt1,overwrite=True)
 pa2 = 'vv'+str(eco)
 pa3 = str(eco)+'v3'
 pa4 = 'eco_'+str(eco)
 pa5 = 'eco_'+str(eco)+'.tif'
 pa6 = str(eco)+'v4'
 
 grass. message ("setting up the working region")
 grass.run_command('g.region',vect=eco0,res=1000)
 grass.run_command('v.to.rast',input=eco0,out=eco0,use='cat',labelcol='eco_id',overwrite=True)
 opt3 = pa2+'= @'+eco0
 opt4 = pa2+'= round('+pa2+')'
 grass.run_command('r.mapcalc',expression=opt3,overwrite=True)
 grass.run_command('r.mapcalc',expression=opt4,overwrite=True)
 grass.run_command('r.mask', vector=eco0, where=opt1)
 opt2 = pa4+'='+pa2
 grass.run_command('r.mapcalc',expression=opt2,overwrite=True)
 grass.run_command('g.remove', rast='MASK')
 grass.run_command('g.region', flags='d')
 grass.run_command('r.mask', rast='pre')
 grass.run_command('r.buffer',input=pa4,output=pa3,distances=250,units='kilometers',overwrite=True)
 grass.run_command('g.region',zoom=pa3)
 optk = pa6+'= if('+pa3+'!=0,1,0)'
 grass.run_command('r.mapcalc',expression=optk,overwrite=True)
 grass.run_command('r.null',map=pa6,null=0)
 grass.run_command('r.out.gdal',input=pa6,out=pa5,overwrite=True)
 grass.run_command('g.remove', rast='MASK')
 grass. message ("Deleting tmp layers") 
 grass.run_command('g.mremove',rast='*v3',flags='f') 
 grass.run_command('g.mremove',rast='*v2',flags='f') 
 grass.run_command('g.mremove',rast='v0_*',flags='f') 
 grass.run_command('g.mremove',vect='v0_*',flags='f') 
 grass.run_command('g.mremove',rast='vv*',flags='f')
 grass. message ("Done")
 print "Done ECO:"+str(eco)
print "FINISHED"
