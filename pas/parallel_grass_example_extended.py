# Code by Javier Martinez-Lopez
from multiprocessing import cpu_count,Pool,Lock
import subprocess
import os
import sys
import csv
import numpy as np
from datetime import datetime

GISBASE = os.environ['GISBASE'] = "/usr/local/grass-7.1.svn"
GRASSDBASE = "/home/javier/Desktop/grassdb"
MYLOC = "newloc_moll"
mapset = 'javier'
col= 'WDPA_ID'

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)

source = 'pas_africa'
#pa_list = 'Spain','France'
grass. message ("Extracting list of PAs")
pa_list0 = grass. read_command ('v.db.select', map=source,column=col). splitlines ()
pa_list2 = np.unique(pa_list0)
pa_list = pa_list2[0:5]
print pa_list

csvname1 = 'pas_segm_done.csv'
if os.path.isfile(csvname1) == False:
 wb = open(csvname1,'a')
 wb.write('None')
 wb.write('\n')
 wb.close()

pa_list_done = np.genfromtxt(csvname1,dtype='string')

def function(pa):
 if pa not in pa_list_done:
  mymapset = 'm'+str(pa)
  spn00= str(GRASSDBASE)+'/'+str(MYLOC)+'/PERMANENT/DEFAULT_WIND'
  spn0x= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)
  spn01= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/WIND'
  spn= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/SEARCH_PATH'
  comm0 = 'mkdir '+spn0x
  comm = 'cp '+spn00+' '+spn01
  os.system(comm0)
  os.system(comm)
  wb = open(spn,'a')
  wb.write('PERMANENT')
  wb.write('\n') 
  wb.write(str(mapset)) # javier
  wb.write('\n')
  wb.write(str(mymapset))
  wb.write('\n')
  wb.close()

  gsetup.init(GISBASE, GRASSDBASE, MYLOC, mymapset)
  pa0 = 'v0_'+pa
# opt1 = 'NAME=\''+pa+'\''
  opt1 = col+'='+pa
  grass.run_command('v.extract', input=source, out=pa0, where=opt1,overwrite=True)
  grass.run_command('g.region',vect=pa0)
  pa5 = pa0+'.tif'
  grass.run_command('r.mask', vector=source,where=opt1)
  #grass.run_command('v.to.rast',input=pa0,typ='area',out='MASK',use='val',val='1',overwrite=True)
  opt2 = pa0+'= const'
  grass.run_command('r.mapcalc',expression=opt2,overwrite=True)
  grass.run_command('r.out.gdal',input=pa0,out=pa5)

  comm1 = 'rm -rf '+spn0x
  os.system(comm1)

pool = Pool()
pool.map(function,pa_list)
pool.close()
pool.join()

print 'END'
