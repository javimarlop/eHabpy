# Code by Javier Martinez-Lopez (utf-8)
#from multiprocessing import cpu_count,Pool,Lock
from datetime import datetime
import numpy as np
from tqdm import *
import os
import sys
import csv

GISBASE = os.environ['GISBASE'] = "/home/majavie/grass_last/new/grass7_trunk/dist.x86_64-unknown-linux-gnu"
GRASSDBASE = "/local1/majavie/hanksgrass7"
MYLOC = "global_MW"
mapset = 'ehabitat'
col= 'WDPA_ID'#'wdpaid'

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)
# add rnd number to cat!!!
source = 'wdpa_snapshot_new_mollweide'#'wdpa_aug14_100km2_moll'
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

n2 = len(pa_list)
#print pa_list[1]
for px in tqdm(range(0,n2)):
 aleat = np.random.random_integers(100)
#for pa in pa_list:
 pa = pa_list[px]
 if pa not in pa_list_done:
  if os.path.isfile('/local1/majavie/hanksgrass7/global_MW/ehabitat/group/segm/REF') == True:
   os.system ('rm /local1/majavie/hanksgrass7/global_MW/ehabitat/group/segm/REF')
  print pa
  pa44 = 'pa_'+str(pa)
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
  minaream = minarea*1000
  #print minarea
  grass.run_command('i.pca', flags='n', input='pre,eprsqrt,slope,tree,herb,ndwi,ndvimax2,ndvimin,bio', output=pa44, overwrite=True) # dem
  pca1 = pa44+'.1'
  pca2 = pa44+'.2'
  pca3 = pa44+'.3'
  pcas = pca1+','+pca2+','+pca3
  grass.run_command('i.group',gr='segm',input=pcas)
  os.system ('cat /local1/majavie/hanksgrass7/global_MW/ehabitat/group/segm/REF')
  #grass. message ("segmenting the park")
  grass.run_command('i.segment', group='segm', output=pa2, threshold='0.8', method='region_growing', similarity='euclidean', memory='100000', minsize=minarea, iterations='20',overwrite=True) # 
  #grass. message ("cropping the segments")
  grass.run_command('r.mask', vector=source, where=opt1)
  opt2 = pa3+'='+pa2
  grass.run_command('r.mapcalc',expression=opt2,overwrite=True) # usar const como mapa para crear plantilla de PA con unos y ceros
  grass.run_command('g.remove', rast='MASK')

  b = grass.read_command('r.stats',input=pa3,flags='nc',separator='\n').splitlines()
  pa3b = pa3
  c = '-1'
  clean = None
  for g in np.arange(1,len(b),2):
   if b[g] < minarea:
    clean = 'T'
    c = c + ' || '+str(b[g-1])
  if clean == 'T':
   pa3b = pa3+'b'
   optt = pa3b+' = if('+old+'=='+c+', nmode('+old+'[-1,0],'+old+'[-1,1],'+old+'[-1,-1],'+old+'[0,1],'+old+'[0,-1],'+old+'[1,0],'+old+'[1,1],'+old+'[1,-1]),'+old+')'
   grass.run_command('r.mapcalc',expression=optt,overwrite=True)
  #grass. message ("Number of cells per segment")
  #grass.run_command('r.stats',input=pa3,out=pa5,overwrite=True) # flags='nc'
  #grass. message ("converting to vector")
  grass.run_command('r.to.vect', input=pa3b,out=pa4,type ='area',flags='v',overwrite=True)
  #grass. message ("adding labels to segments")
  grass.run_command('v.db.addcolumn', map=pa4,col='wdpaid_pa VARCHAR')
  grass.run_command('v.db.update', map=pa4,col='wdpaid_pa',value=pa)
  grass.run_command('v.db.addcolumn', map=pa4,col='aleat VARCHAR')
  grass.run_command('v.db.update', map=pa4,col='aleat',value=aleat)
  #grass. message ("Checking shapefile")
  grass.run_command('v.clean', input=pa4,out=pa44,tool='rmarea',thres=minaream,overwrite=True)
  grass.run_command('v.db.addcolumn', map=pa44,col='segm_id VARCHAR')
  grass.run_command('v.db.update', map=pa44,col='segm_id',qcol='wdpaid_pa || cat || aleat')
  #grass. message ("Exporting shapefile")
  if os.path.isfile('parks_segmented_pca.shp') == False:
   grass.run_command('v.out.ogr',input=pa44,ola='parks_segmented_pca',type='area',dsn='.') 
  else:
   grass.run_command('v.out.ogr',flags='a',input=pa44,ola='parks_segmented_pca',type='area',dsn='.')
  grass. message ("Deleting tmp layers") 
  grass.run_command('g.mremove',typ='rast',patt='*v3',flags='f') 
  grass.run_command('g.mremove',typ='rast',patt='*v2',flags='f') 
  grass.run_command('g.mremove',typ='rast',patt='v0_*',flags='f') 
  grass.run_command('g.mremove',typ='rast',patt='pa_*',flags='f') 
  grass.run_command('g.mremove',typ='vect',patt='v0_*',flags='f') 
  grass.run_command('g.mremove',typ='vect',patt='pa_*',flags='f') 
  grass.run_command('g.mremove',typ='vect',patt='paa_*',flags='f') 
  grass.run_command('g.mremove',typ='vect',patt='paa_pca*',flags='f')
  os.system ('rm /local1/majavie/hanksgrass7/global_MW/ehabitat/group/segm/REF')
  grass. message ("Done")
  print "Done PA:"+pa 
  wb = open(csvname1,'a')
  var = str(pa)
  wb.write(var)
  wb.write('\n')
  wb.close() 

print str(datetime.now())
print 'END'
