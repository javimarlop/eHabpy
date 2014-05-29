# Code by Javier Martinez-Lopez
from multiprocessing import cpu_count,Pool,Lock
import os
import sys
import csv

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
 os.environ['GRASS_REGION'] = grass.region_env(res=elem)
 varx = grass.read_command ('g.region',flags='g'). splitlines ()
 wb = open('results.txt','a')
 var = str(elem)+' '+str(varx)
 wb.write(var)
 wb.write('\n')
 wb.close()
 os.environ.pop('GRASS_REGION')

elems = '100','200','300','400'

pool = Pool()
pool.map(function,elems)
pool.close()
pool.join()

print 'END'
