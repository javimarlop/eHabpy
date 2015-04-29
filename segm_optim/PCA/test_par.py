# Code by Javier Martinez-Lopez (utf-8)
from multiprocessing import cpu_count,Pool,Lock
import multiprocessing
import subprocess
from datetime import datetime
import numpy as np
from tqdm import *
import os
import sys
import csv
import gc

#GISBASE = os.environ['GISBASE'] = "/home/majavie/hierba_hanks/grass-7.1.svn"#/home/majavie/grass_last/new/grass7_trunk/dist.x86_64-unknown-linux-gnu"
#GRASSDBASE = "/home/majavie/hanksgrass7"
#MYLOC = "global_MW"
#global mapset
#mapset = 'm'#ehabitat'
#col= 'wdpaid'

#sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
#import grass.script as grass
#import grass.script.setup as gsetup

#gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)
#source = 'wdpa_aug14_100km2_moll'
#grass. message ("Extracting list of PAs")
#pa_list0 = grass. read_command ('v.db.select', map=source,column=col). splitlines ()
pa_list0 = np.genfromtxt('palist_100km2.csv',dtype='string')
pa_list2 = np.unique(pa_list0)
n = len(pa_list2)
pa_list = pa_list2[0:50] # testing 5 first!
#pa_list_tmp = ['6317','555523378','555545976','555545971','555546231','71213','4840','151315','257','101922','2017','11','68175','643','555542456','4328', '124389','2006','900883','93294','198301','61612','555577555','555556047','555538721','19297','555580364','555540045','6317','2579','372151','979','2013','983','55557661','555555854','132','19984','35002','10908','1111192','984','9632','1083','801','555523378','555545976','555545971','555546231','71213','4840','151315','257','101922','2017','11','68175','643','555542456','4328', '124389','2006','900883','93294','198301','61612','555577555','555556047','555538721','19297','555580364','555540045'] # ]# '95786'
#pa_list = np.unique(pa_list_tmp)
print pa_list

def fsegm(pa):

 ### mp.current_process().cnt += 1
 current = multiprocessing.current_process()
 mn = current._identity[0]
 print 'running:', mn

 GISBASE = os.environ['GISBASE'] = "/home/majavie/hierba_hanks/grass-7.1.svn"
 GRASSDBASE = "/home/majavie/hanksgrass7"
 MYLOC = "global_MW"
 mapset = 'm'
 sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
 import grass.script as grass
 import grass.script.setup as gsetup
 gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)

 mapset2 = 'm'+str(mn)#ehabitat'
 os.system ('rm -rf /home/majavie/hanksgrass7/global_MW/'+mapset2)
 grass.run_command('g.mapset',mapset=mapset2,location='global_MW',gisdbase='/home/majavie/hanksgrass7',flags='c')

 gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset2)
 print mapset2, grass.gisenv()
 print pa
 grass.run_command('g.mapsets',mapset='ehabitat,rasterized_parks',operation='add')
 grass. message ("Deleting tmp layers")
 os.system ('rm -rf /home/majavie/hanksgrass7/global_MW/'+mapset2)
 gc.collect()

pool = Pool(3)#9
pool.map(fsegm,pa_list)
pool.close()
pool.join()

#print str(datetime.now())
#print 'END'
