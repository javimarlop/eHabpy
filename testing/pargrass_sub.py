# Code by Javier Martinez-Lopez
from multiprocessing import cpu_count,Pool,Lock
import subprocess
import os
import sys
import csv
import numpy as np
from datetime import datetime

GISBASE = os.environ['GISBASE'] = "/usr/local/grass-7.1.svn"
GRASSDBASE = os.path.join("/home/javier/Desktop/grassdb")
MYLOC = "nc_spm_08_reduced"
mapset = "user1"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as grass
import grass.script.setup as gsetup
 
gsetup.init(GISBASE,
            GRASSDBASE, MYLOC, mapset)

print grass.gisenv()

def function(elem):
 print elem
 mymapset = 'm'+str(elem)
 komm = grass.run_command('g.mapset',mapset=mymapset,flags='c')
 p = subprocess.check_call(komm,shell=T,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 #out, err = p.communicate()
 if p.returncode == 0:
  gge = grass.gisenv()
  spn0= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)
  spn= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/SEARCH_PATH'
  wb = open(spn,'a')
  wb.write('PERMANENT')
  wb.write('\n') 
  wb.write('user1')
  wb.write('\n')
  wb.write(str(mymapset))
  wb.write('\n')
  wb.close()
  pa0 = 's'+str(elem)
  comm2 = 'cat = '+str(elem)
  grass.run_command('g.region',rast='elevation')
  grass.run_command('g.region',res=elem)
  varx = grass.read_command ('g.region',flags='g'). splitlines ()
  wb = open('results.txt','a')
  var = str(elem)+' '+str(gge)+' '+str(varx)
  wb.write(var)
  wb.write('\n')
  wb.close()
  elem=None
  mymapset=None

elems = '100','200','300','400'

pool = Pool()
pool.map(function,elems)
pool.close()
pool.join()


print 'END'
