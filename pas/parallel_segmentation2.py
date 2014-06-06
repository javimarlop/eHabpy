# Code by Javier Martinez-Lopez
from multiprocessing import cpu_count,Pool,Lock
import subprocess
import os
import sys
import csv
import numpy as np
from datetime import datetime

GISBASE = os.environ['GISBASE'] = '/home/majavie/grass_last/grass-7.1.svn'
GRASSDBASE = os.path.join("/local1/majavie/hanksgrass7")
MYLOC = "global_MW"
mapset   = "rasterized_parks"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as grass
import grass.script.setup as gsetup
 
gsetup.init(GISBASE,
            GRASSDBASE, MYLOC, mapset)

print grass.gisenv()

source = 'wdpa_snapshot_new_mollweide@javier'
grass. message ("Extracting list of PAs")
pa_list0 = grass. read_command ('v.db.select', map=source,column='wdpa_id'). splitlines ()
pa_list2 = np.unique(pa_list0)
n = len(pa_list2)-1#2
#pa_list = pa_list2[0:n]
pa_list = '257','101922','2017','11','68175','643','555542456' # testing set

csvname1 = 'pas_segm_done.csv'
if os.path.isfile(csvname1) == False:
 wb = open(csvname1,'a')
 wb.write('None')
 wb.write('\n')
 wb.close()

pa_list_done = np.genfromtxt(csvname1,dtype='string')


def function(pa):
 #print pa_list
 mymapset = 'm'+str(pa)
 grass.run_command('g.mapset',mapset=mymapset,flags='c')
 gge = grass.gisenv()
 spn0= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)
 spn= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/SEARCH_PATH'
 wb = open(spn,'a')
 wb.write('PERMANENT')
 wb.write('\n') 
 wb.write('rasterized_parks')
 wb.write('\n')
 wb.write(str(mymapset))
 wb.write('\n')
 wb.close()
 # .gislock remains?

def function2(pa):
 mymapset = 'm'+str(pa)
 grass.run_command('g.mapset',mapset=mymapset)
 pa0 = 'v0_'+pa
 grass.run_command('g.region',vect=pa0)
 grass.run_command('g.remove', rast='MASK') 
 if pa not in pa_list_done:
  opt1 = 'wdpa_id = '+pa
  pa2 = pa+'v2'
  pa3 = pa+'v3'
  pa4 = 'paa_'+pa
  pa44 = 'pa_'+pa
  pa5 = pa4+'.txt'
  same = pa2+'= const'
  grass.run_command('r.mapcalc',expression='const = if(gcmask>=0,1,null())',overwrite=True)
  grass.run_command('r.mapcalc',expression=same,overwrite=True)
  a = grass.read_command('r.stats',input='const',flags='nc',separator='\n').splitlines()
  if len(a)==0: a = [1, 625]
  minarea = np.sqrt(int(a[1]))#/1000
  grass.run_command('i.segment', group='segm', output=pa2, threshold='0.7', method='region_growing', similarity='euclidean', minsize=minarea, memory='100000', iterations='20',overwrite=True)
  grass.run_command('r.mask', vector=source, where=opt1)
  opt2 = pa3+'='+pa2
  grass.run_command('r.mapcalc',expression=opt2,overwrite=True) # usar const como mapa para crear plantilla de PA con unos y ceros
  grass.run_command('g.remove', rast='MASK')
  grass.run_command('r.stats',input=pa3,out=pa5,overwrite=True) # flags='nc'
  grass.run_command('r.to.vect', input=pa3,out=pa4,type ='area',flags='v',overwrite=True)
  grass.run_command('v.db.addcolumn', map=pa4,col='wdpaid VARCHAR')
  grass.run_command('v.db.update', map=pa4,col='wdpaid',value=pa)
  grass.run_command('v.clean', input=pa4,out=pa44,tool='rmarea',thres=minarea,overwrite=True)
  grass.run_command('v.db.addcolumn', map=pa44,col='wdpa_id VARCHAR')
  grass.run_command('v.db.update', map=pa44,col='wdpa_id',qcol='wdpaid || cat')
  grass.run_command('v.out.ogr',flags='a',input=pa44,ola='parks_segmented2',dsn='.') 
  print "Done PA:"+pa
  wb = open(csvname1,'a')
  var = str(pa)
  wb.write(var)
  wb.write('\n')
  wb.close()

pool = Pool()
pool.map(function,pa_list)
pool.close()
pool.join()

pool = Pool()
pool.map(function2,pa_list)
pool.close()
pool.join()

print 'END'
