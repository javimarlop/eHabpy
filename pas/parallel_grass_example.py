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
MYLOC = "ecad5_grassdata_ll"
sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

source = 'country_boundaries'
pa_list = 'Spain','France'

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
  wb.write('user1')
  wb.write('\n')
  wb.write(str(mymapset))
  wb.write('\n')
  wb.close()

  gsetup.init(GISBASE, GRASSDBASE, MYLOC, mymapset)
  pa0 = 'v0_'+pa
  opt1 = 'NAME=\''+pa+'\''
  grass.run_command('v.extract', input=source, out=pa0, where=opt1,overwrite=True)
  grass.run_command('g.region',vect=pa0)
  pa5 = pa+'.tif'
  grass.run_command('r.mask', vector=pa0)
  opt2 = pa+'= const'
  grass.run_command('r.mapcalc',expression=opt2,overwrite=True)
  grass.run_command('r.out.gdal',input=pa,out=pa5)

  print "Done:"+pa
  comm1 = 'rm -rf '+spn0x
  os.system(comm1)

pool = Pool()
pool.map(function,pa_list)
pool.close()
pool.join()

print 'END'
