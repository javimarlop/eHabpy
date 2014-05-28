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
 print grass.gisenv()
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
 comm2 = 'cat = '+str(elem)
 varx = grass.read_command ('v.db.select', map='segm',column='newarea2',where=comm2). splitlines ()
 wb = open('results.csv','a')
 var = str(elem)+' '+str(varx)
 wb.write(var)
 wb.write('\n')
 wb.close()
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
