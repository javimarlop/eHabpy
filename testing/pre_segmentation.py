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

source = 'wdpa_snapshot_new_mollweide@javier'
grass. message ("Extracting list of PAs")
pa_list0 = grass. read_command ('v.db.select', map=source,column='wdpa_id'). splitlines ()
#pa_list = np.unique(pa_list0)
pa_list = '257','101922','2017','11','68175','643','555542456' # testing set


csvname1 = 'pas_pre_segm_done.csv'
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
  opt1 = 'wdpa_id = '+pa
  print px
  print "Extracting PA:"+pa
  grass.run_command('v.extract', input=source, out=pa0, where = opt1,overwrite=True)
  print "Done PA:"+pa
  wb = open(csvname1,'a')
  var = str(pa)
  wb.write(var)
  wb.write('\n')
  wb.close() 
# try to paralellize it?

print str(datetime.now())
print 'END'
