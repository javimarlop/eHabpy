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

#source = 'wdpa_snapshot_new_mollweide@javier'
source='parks_80601'
grass. message ("Extracting list of PAs")
pa_list0 = grass. read_command ('v.db.select', map=source,column='wdpa_id'). splitlines ()
pa_list = np.unique(pa_list0)

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
 pa2 = pa+'v2'
 pa3 = pa+'v3'
 pa4 = 'pa_'+pa
  
 grass. message ("setting up the working region")
 grass.run_command('g.region',vect=pa0,res=1000)
 grass. message ("cropping the segments")
 grass.run_command('r.mask', vector=source, where=opt1)
 grass. message ("segmenting the park")
 grass.run_command('i.segment', group='segm', output=pa2, threshold='0.7', method='region_growing', similarity='euclidean', minsize='25', memory='10000', iterations='20',overwrite=True)
 #opt2 = pa3+'='+pa2
 #grass.run_command('r.mapcalc',expression=opt2,overwrite=True)
 grass.run_command('g.remove', rast='MASK')
 grass. message ("Number of cells per segment")
 grass.run_command('r.stats',input=pa2,flags='nc')
 grass. message ("converting to vector")
 grass.run_command('r.to.vect', input=pa2,out=pa4,type ='area',flags='v',overwrite=True)
 grass. message ("adding labels to segments")
 grass.run_command('v.db.addcolumn', map=pa4,col='wdpaid VARCHAR')
 grass.run_command('v.db.addcolumn', map=pa4,col='wdpa_id VARCHAR')
 grass.run_command('v.db.update', map=pa4,col='wdpaid',value=pa)
 grass.run_command('v.db.update', map=pa4,col='wdpa_id',qcol='wdpaid || cat')
 grass. message ("Exporting shapefile")
 if px == 1:
  grass.run_command('v.out.ogr',input=pa4,ola='parks_segmented',dsn='.') 
 else:
  grass.run_command('v.out.ogr',flags='a',input=pa4,ola='parks_segmented',dsn='.') 
 grass. message ("Deleting tmp layers") 
 grass.run_command('g.mremove',rast='*v3',flags='f') 
 grass.run_command('g.mremove',rast='*v2',flags='f') 
 grass.run_command('g.mremove',rast='v0_*',flags='f') 
 grass.run_command('g.mremove',rast='pa_*',flags='f') 
 grass.run_command('g.mremove',vect='v0_*',flags='f') 
 grass.run_command('g.mremove',vect='pa_*',flags='f') 
 grass. message ("Done")
 print "Done PA:"+pa

# try to paralellize it?