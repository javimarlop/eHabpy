# Code by Javier Martinez-Lopez
import os
import sys
import csv
import numpy as np
from tqdm import *
from datetime import datetime

#gisbase = os.environ['GISBASE'] = "/usr/local/grass-7.0.svn"
gisbase = os.environ['GISBASE'] = "/home/majavie/grass7_source/g71/grass7_trunk/dist.x86_64-unknown-linux-gnu"
gisdbase = os.path.join("/local1/majavie/hanksgrass7")
location = "global_MW"
mapset   = "rasterized_parks"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as grass
import grass.script.setup as gsetup
 
gsetup.init(gisbase,
            gisdbase, location, mapset)
 
print grass.gisenv()

source='ecoregs_moll'
#grass. message ("Extracting list of eco")
#pa_list0 = grass. read_command ('v.db.select', map=source,column='eco_id'). splitlines ()
#pa_list = np.unique(pa_list0)
pa_list = '80601'


csvname1 = 'eco_segm_done.csv'
if os.path.isfile(csvname1) == False:
 wb = open(csvname1,'a')
 wb.write('None')
 wb.write('\n')
 wb.close()

pa_list_done = np.genfromtxt(csvname1,dtype='string')
#print pa_list_done

grass. message("omitting previous masks")
grass.run_command('g.remove', rast='MASK')

n = len(pa_list)-1#2
#print pa_list[1]
for px in tqdm(range(0,n)):

#for pa in pa_list:
 pa = pa_list[px]
 if pa not in pa_list_done:
  pa0 = 'v0_'+pa
  opt1 = 'eco_id = '+pa
  print px
  print "Extracting ECO:"+pa
  grass.run_command('v.extract', input=source, out=pa0, where = opt1,overwrite=True)
  pa2 = pa+'v2'
  pa3 = pa+'v3'
  pa4 = 'eco_'+pa
  pa5 = pa4+'.txt'
  same = pa2+'= const'
  
  grass. message ("setting up the working region")
  grass.run_command('g.region',vect=pa0,res=1000)
  #grass.run_command('r.mask', vector=source, where=opt1)
  grass.run_command('r.mapcalc',expression='const = if(gcmask>=0,1,null())',overwrite=True)
  grass.run_command('r.mapcalc',expression=same,overwrite=True)
  a = grass.read_command('r.stats',input='const',flags='nc',separator='\n').splitlines()
  if len(a)==0: a = [1, 625]
  #grass.run_command('g.remove', rast='MASK')
  print a
  minarea = np.sqrt(int(a[1]))#/1000
  print minarea
  grass. message ("segmenting the ecoregion")
  grass.run_command('i.segment', group='segm', output=pa2, threshold='0.7', method='region_growing', similarity='euclidean', minsize=minarea, memory='100000', iterations='20',overwrite=True)
  grass. message ("cropping the segments")
  grass.run_command('r.mask', vector=source, where=opt1)
  opt2 = pa3+'='+pa2
  grass.run_command('r.mapcalc',expression=opt2,overwrite=True) # usar const como mapa para crear plantilla de PA con unos y ceros
  grass.run_command('g.remove', rast='MASK')
  grass. message ("Number of cells per segment")
  grass.run_command('r.stats',input=pa3,out=pa5,overwrite=True) # flags='nc'
  grass. message ("converting to vector")
  grass.run_command('r.to.vect', input=pa3,out=pa4,type ='area',flags='v',overwrite=True)
  grass. message ("adding labels to segments")
  grass.run_command('v.db.addcolumn', map=pa4,col='ecoid VARCHAR')
  grass.run_command('v.db.addcolumn', map=pa4,col='eco_id VARCHAR')
  grass.run_command('v.db.update', map=pa4,col='ecoid',value=pa)
  grass.run_command('v.db.update', map=pa4,col='eco_id',qcol='ecoid || cat')
  #grass. message ("Checking shapefile")
  grass. message ("Exporting shapefile")
  if px == 0:
   grass.run_command('v.out.ogr',input=pa4,ola='eco_segmented',dsn='.') 
  else:
   grass.run_command('v.out.ogr',flags='a',input=pa4,ola='eco_segmented',dsn='.') 
  grass. message ("Deleting tmp layers") 
  grass.run_command('g.mremove',rast='*v3',flags='f') 
  grass.run_command('g.mremove',rast='*v2',flags='f') 
  grass.run_command('g.mremove',rast='v0_*',flags='f') 
  grass.run_command('g.mremove',rast='pa_*',flags='f') 
  grass.run_command('g.mremove',vect='v0_*',flags='f') 
  grass.run_command('g.mremove',vect='pa_*',flags='f') 
  grass. message ("Done")
  print "Done ECO:"+pa
  wb = open(csvname1,'a')
  var = str(pa)
  wb.write(var)
  wb.write('\n')
  wb.close() 
# try to paralellize it?

print str(datetime.now())
print 'END'
