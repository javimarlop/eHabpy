# Code by Javier Martinez-Lopez
from multiprocessing import cpu_count,Pool,Lock
import os
import sys
import csv
import numpy as np
from datetime import datetime

GISBASE = os.environ['GISBASE'] = "/usr/lib/grass70"
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
 #l = Lock()
 #l.acquire()
 #rndnr = int(np.random.random_sample(1,)*1000)
 mymapset = 'm'+str(elem)
 #l.release()
 #print mymapset
 grass.run_command('g.mapset',mapset=mymapset,flags='c')
 gge = grass.gisenv()
 spn0= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)
 spn= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/SEARCH_PATH'
 #cmm3 = 'mkdir '+spn0
 #print cmm3
 #os.system(cmm3)
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
 #grass.run_command('v.extract', input='segm', out=pa0, where = comm2,overwrite=True)
 grass.run_command('r.mask', vector='segm', where=comm2)
 opt2 = pa0+' = elevation'
 grass.run_command('r.mapcalc',expression=opt2,overwrite=True)
 grass.run_command('g.region',zoom=pa0)
 varx = grass.read_command ('g.region',flags='g'). splitlines ()
 #varx = grass.read_command ('v.db.select', map='segm',column='newarea2',where=comm2). splitlines ()
 wb = open('results.txt','a')
 var = str(elem)+' '+str(gge)+' '+str(varx)
 wb.write(var)
 wb.write('\n')
 wb.close()
 grass.run_command('g.remove', rast='MASK')
 #comm = 'rm /home/javier/Desktop/grassdb/nc_spm_08_reduced/'+str(mymapset)+'/.gislock'
 ##os.system(comm)
 #rndnr=None
 mymapset=None

elems = '1','2','3','4'
#elems = '1'#,'2','3'

pool = Pool()
pool.map(function,elems)
pool.close()
pool.join()


print 'END'
