# Code by Javier Martinez-Lopez (utf-8)
from multiprocessing import cpu_count,Pool,Lock
from datetime import datetime
import numpy as np
import os
import sys
import csv

GISBASE = os.environ['GISBASE'] = "/home/majavie/grass7_source/g71/grass7_trunk/dist.x86_64-unknown-linux-gnu"
GRASSDBASE = "/local1/majavie/hanksgrass7"
MYLOC = "global_MW"
mapset = 'rasterized_parks'
col= 'wdpa_id'

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)

source = 'wdpa_snapshot_new_mollweide@javier'
grass. message ("Extracting list of PAs")
pa_list0 = grass. read_command ('v.db.select', map=source,column=col). splitlines ()
pa_list2 = np.unique(pa_list0)
n = len(pa_list2)
pa_list = pa_list2[0:n-2] # testing 5 first!
#pa_list = '257','101922','2017','11','68175','643','555542456'
print pa_list

csvname1 = 'pas_segm_done_pca.csv'
if os.path.isfile(csvname1) == False:
 wb = open(csvname1,'a')
 wb.write('None')
 wb.write('\n')
 wb.close()

pa_list_done = np.genfromtxt(csvname1,dtype='string')

def segmentation(pa):
 if pa not in pa_list_done:
  mymapset = 'm'+str(pa)
  pa44 = 'pa_'+str(pa)
  spn00= str(GRASSDBASE)+'/'+str(MYLOC)+'/PERMANENT/DEFAULT_WIND'
  spn0x= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)
  spn0xgr= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/group/segm'
  #spn0xgr2= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mapset)+'/group/segm/REF'
  spn0xgr22= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/group/segm/REF'
  spn01= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/WIND'
  spn= str(GRASSDBASE)+'/'+str(MYLOC)+'/'+str(mymapset)+'/SEARCH_PATH'
  #comm0 = 'mkdir '+spn0x
  comm0gr = 'mkdir -p '+spn0xgr
  comm = 'cp '+spn00+' '+spn01
  #commgr = 'cp '+spn0xgr2+' '+spn0xgr22
  #os.system(comm0)
  os.system(comm0gr)
  os.system(comm)
  #os.system(commgr)
  wb = open(spn,'a')
  wb.write('PERMANENT')
  wb.write('\n')
  wb.write('ehabitat')
  wb.write('\n') 
  wb.write(str(mapset)) # javier
  wb.write('\n')
  wb.write(str(mymapset))
  wb.write('\n')
  wb.close()
#  wb = open(spn0xgr22,'a')
#  wb.write(pca1)
#  wb.close()
  gsetup.init(GISBASE, GRASSDBASE, MYLOC, mymapset)
  pa0 = 'v0_'+pa
  opt1 = col+'='+pa
  grass.run_command('v.extract', input=source, out=pa0, where = opt1,overwrite=True) # check inital region from which to copy from!
  pa2 = pa+'v2'
  pa3 = pa+'v3'
  pa4 = 'paa_pca_'+pa
  pa5 = pa4+'.txt'
  same = pa2+'= const'
  #grass. message ("setting up the working region")
  grass.run_command('g.region',vect=pa0,res=1000)
  grass.run_command('r.mapcalc',expression='const = if(gcmask>=0,1,null())',overwrite=True)
  grass.run_command('r.mapcalc',expression=same,overwrite=True)
  a = grass.read_command('r.stats',input='const',flags='nc',separator='\n').splitlines()
  if len(a)==0: a = [1, 625]
  #print a
  minarea = np.sqrt(int(a[1]))#/1000
  #print minarea
  grass.run_command('i.pca', flags='n', input='pre,epr,slope,tree,herb,ndwi,ndvimax2,ndvimin,bio', output=pa44, overwrite=True) # dem
  pca1 = pa44+'.1'
  pca2 = pa44+'.2'
  pca3 = pa44+'.3'
  pcas = pca1+','+pca2+','+pca3
  grass.run_command('i.group',gr='segm',input=pcas)
  #grass. message ("segmenting the park")
  grass.run_command('i.segment', group='segm', output=pa2, threshold='0.8', method='region_growing', similarity='euclidean', memory='100000', minsize=minarea, iterations='20',overwrite=True) # 
  #grass. message ("cropping the segments")
  grass.run_command('r.mask', vector=source, where=opt1)
  opt2 = pa3+'='+pa2
  grass.run_command('r.mapcalc',expression=opt2,overwrite=True) # usar const como mapa para crear plantilla de PA con unos y ceros
  grass.run_command('g.remove', rast='MASK')
  #grass. message ("Number of cells per segment")
  #grass.run_command('r.stats',input=pa3,out=pa5,overwrite=True) # flags='nc'
  #grass. message ("converting to vector")
  grass.run_command('r.to.vect', input=pa3,out=pa4,type ='area',flags='v',overwrite=True)
  #grass. message ("adding labels to segments")
  grass.run_command('v.db.addcolumn', map=pa4,col='wdpaid VARCHAR')
  grass.run_command('v.db.update', map=pa4,col='wdpaid',value=pa)
  #grass. message ("Checking shapefile")
  grass.run_command('v.clean', input=pa4,out=pa44,tool='rmarea',thres=minarea,overwrite=True)
  grass.run_command('v.db.addcolumn', map=pa44,col='wdpa_id VARCHAR')
  grass.run_command('v.db.update', map=pa44,col='wdpa_id',qcol='wdpaid || cat')
  #grass. message ("Exporting shapefile")
  if os.path.isfile('parks_segmented_pca_par.shp') == False:
   grass.run_command('v.out.ogr',input=pa44,ola='parks_segmented_pca_par',type='area',dsn='.') 
  else:
   grass.run_command('v.out.ogr',flags='a',input=pa44,ola='parks_segmented_pca_par',type='area',dsn='.') 
  comm1 = 'rm -rf '+spn0x
  os.system(comm1)
  wb = open(csvname1,'a')
  var = str(pa)
  wb.write(var)
  wb.write('\n')
  wb.close() 

pool = Pool(10)
pool.map(segmentation,pa_list)
pool.close()
pool.join()

print str(datetime.now())
print 'END'
