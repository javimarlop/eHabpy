# Code by Javier Martinez-Lopez
from multiprocessing import cpu_count,Pool,Lock
import os
import sys

GISBASE = os.environ['GISBASE'] = '/home/majavie/grass_last/grass-7.1.svn'
GRASSDBASE = os.path.join("/home/majavie/grassdb")
MYLOC = "nc_spm_08_reduced"
mapset = "user1"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as grass
import grass.script.setup as gsetup
 
gsetup.init(GISBASE,
            GRASSDBASE, MYLOC, mapset)

print grass.gisenv()

def function(elem):
 outf='res_'+str(elem)+'.txt'
 os.environ['GRASS_REGION'] = grass.region_env(res=elem)
 varx = grass.read_command ('g.region',flags='g'). splitlines ()
 grass.run_command('r.univar',map='elevation',flags='g',out=outf,overwrite=True)
 os.environ.pop('GRASS_REGION')

#elems = '100','200','300','400'
elems = '100'

function(elems)
#pool = Pool()
#pool.map(function,elems)
#pool.close()
#pool.join()


print 'END'
