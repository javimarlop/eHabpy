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

GISBASE = os.environ['GISBASE'] = "/home/majavie/grass_last/new/grass7_trunk/dist.x86_64-unknown-linux-gnu"
GRASSDBASE = "/home/majavie/hanksgrass7"
MYLOC = "global_MW"
mapset = 'm4'#ehabitat'
col= 'wdpaid'

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)
source = 'wdpa_aug14_100km2_moll'
grass. message ("Extracting list of PAs")
pa_list0 = grass. read_command ('v.db.select', map=source,column=col). splitlines ()
pa_list2 = np.unique(pa_list0)
n = len(pa_list2)
#pa_list = pa_list2[0:n-2] # testing 5 first!
pa_list_tmp = ['6317','555523378','555545976','555545971','555546231','71213','4840','151315','257','101922','2017','11','68175','643','555542456','4328', '124389','2006','900883','93294','198301','61612','555577555','555556047','555538721','19297','555580364','555540045','6317','2579','372151','979','2013','983','55557661','555555854','132','19984','35002','10908','1111192','984','9632','1083','801','555523378','555545976','555545971','555546231','71213','4840','151315','257','101922','2017','11','68175','643','555542456','4328', '124389','2006','900883','93294','198301','61612','555577555','555556047','555538721','19297','555580364','555540045'] # ]# '95786'
pa_list = np.unique(pa_list_tmp)
print pa_list

csvname1 = 'csv/segm_done.csv'
if os.path.isfile(csvname1) == False:
 os.system('touch '+str(csvname1))
 #wb = open(csvname1,'a')
 ##wb.write('None')
 #wb.write('\n')
 #wb.close()

def fsegm(pa):

 pa_list_done = np.genfromtxt(csvname1,dtype='string')
 current = multiprocessing.current_process()
 mn = current._identity[0]
 print 'running:', mn

 if pa not in pa_list_done:

  GISBASE = os.environ['GISBASE'] = "/home/majavie/grass_last/new/grass7_trunk/dist.x86_64-unknown-linux-gnu"
  GRASSDBASE = "/home/majavie/hanksgrass7"
  MYLOC = "global_MW"
  mapset = 'm'+str(mn)#ehabitat'
  col= 'wdpaid'

  sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
  import grass.script as grass
  import grass.script.setup as gsetup

  gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)
  source = 'wdpa_aug14_100km2_moll'

  grass.run_command('g.mremove',typ='vect',patt='*',flags='f') 
  grass.run_command('g.mremove',typ='rast',patt='*',flags='f')

  reps = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

  if os.path.isfile('/home/majavie/hanksgrass7/global_MW/'+mapset+'/group/segm/REF') == True:
   os.system ('rm /home/majavie/hanksgrass7/global_MW/'+mapset+'/group/segm/REF')
  print pa
  pa44 = 'pa_'+str(pa)
  pa44x = 'pax_'+str(pa)
  pa0 = 'v0_'+pa
  opt1 = col+'='+pa
  grass.run_command('v.extract', input=source, out=pa0, where = opt1,overwrite=True) # check inital region from which to copy from!
  pa2 = pa+'v2_'
  pa3 = pa+'v3'
  pa4 = 'paa_'+pa
  pa5 = pa4+'.txt'
  same = pa2+'= const'
  #grass. message ("setting up the working region")
  grass.run_command('g.region',vect=pa0,res=1000)
  grass.run_command('r.mapcalc',expression='const = if(gcmask>=0,1,null())',overwrite=True)
  grass.run_command('r.mapcalc',expression=same,overwrite=True)
  a = grass.read_command('r.stats',input='const',flags='nc',separator='\n').splitlines()
  if len(a)==0: a = [1, 625]
  #print a
  minarea = np.sqrt(int(a[1]))#/2 #10
  minaream = minarea#*1000
  grass.run_command('i.pca', flags='n', input='pre,eprsqrt,slope,tree,herb,ndwi,ndvimax2,ndvimin,bio', output=pa44x, overwrite=True) # dem
  pca1 = pa44x+'.1'
  pca2 = pa44x+'.2'
  pca3 = pa44x+'.3'
  pcas = pca1+','+pca2+','+pca3
  grass.run_command('i.group',gr='segm',input=pcas)
  os.system ('cat /home/majavie/hanksgrass7/global_MW/'+mapset+'/group/segm/REF')
  j = 0
  for thr in reps:
  	  pa2 = pa+'v2_'+str(j)
  	  pa2s = pa+'v2_'+str(j-1)
	  aleat = np.random.random_integers(1000)
  	  grass.run_command('g.region',vect=pa0,res=1000)
	  j= j + 1
	  if thr==0.1:
	  	grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', minsize=minarea, similarity='euclidean', memory='100000', iterations='20',overwrite=True) #  ,seed=pa2i minsize=minarea,
	  else:
	  	grass.run_command('i.segment', group='segm', output=pa2, threshold=thr, method='region_growing', similarity='euclidean', memory='100000', iterations='20',seed=pa2s,overwrite=True) # minsize=minarea
	  #grass. message ("cropping the segments")
	  grass.run_command('r.mask', vector=source, where=opt1)
	  opt2 = pa3+'='+pa2
	  grass.run_command('r.mapcalc',expression=opt2,overwrite=True) # usar const como mapa para crear plantilla de PA con unos y ceros
	  grass.run_command('g.remove', rast='MASK')
	  print minarea

	  b = grass.read_command('r.stats',input=pa3,flags='nc',separator='\n').splitlines()
	  print b
	  clean = None
	  c = pa3
	  for g in np.arange(1,len(b),2):
	   if np.int(b[g]) < minarea:#/10: # lower the threshold if omitting min area!
	    print 'Cleaning small segments I...'
	    print 'cleaning cat '+ str(b[g-1])
	    c2 = 'old'+ str(b[g-1])
	    c22 = c2+'b10km'
	    c3 = 'new'+ str(b[g-1])
	    oper1 = c2+'='+'if('+pa3+'=='+str(b[g-1])+',1,null())'
	    grass.run_command('r.mapcalc',expression=oper1,overwrite=True)
	    grass.run_command('r.buffer',input=c2,output=c22,distances=3,units='kilometers',overwrite=True)
	    grass.run_command('r.mask', raster=c22,maskc=2)
	    buff = grass.read_command('r.stats',input=pa3,flags='nc',sort='desc',separator='\n').splitlines()
	    grass.run_command('g.remove', rast='MASK')
	    if len(buff) > 0:
	     clean = 'T'
	     print 'New: '+str(buff[0])
	     oper1 = c3+'='+'if('+c2+'==1,'+str(buff[0])+',null())'
	     c = c3 + ',' + c
	     grass.run_command('r.mapcalc',expression=oper1,overwrite=True)
	  if clean=='T':
	   print c
	   grass.run_command('r.patch',input=c,out=pa3,overwrite=True)
	   bv = grass.read_command('r.stats',input=pa3,flags='nc',separator='\n').splitlines()
	   print bv

	  b = grass.read_command('r.stats',input=pa3,flags='nc',separator='\n').splitlines()
	  print b
	  clean = None
	  c = pa3
	  for g in np.arange(1,len(b),2):
	   if np.int(b[g]) < minarea:#/10: # lower the threshold if omitting min area!
	    print 'Cleaning small segments II...'
	    print 'cleaning cat '+ str(b[g-1])
	    c2 = 'old'+ str(b[g-1])
	    c22 = c2+'b10km'
	    c3 = 'new'+ str(b[g-1])
	    oper1 = c2+'='+'if('+pa3+'=='+str(b[g-1])+',1,null())'
	    grass.run_command('r.mapcalc',expression=oper1,overwrite=True)
	    grass.run_command('r.buffer',input=c2,output=c22,distances=10,units='kilometers',overwrite=True)
	    grass.run_command('r.mask', raster=c22,maskc=2)
	    buff = grass.read_command('r.stats',input=pa3,flags='nc',sort='desc',separator='\n').splitlines()
	    grass.run_command('g.remove', rast='MASK')
	    if len(buff) > 0:
	     clean = 'T'
	     print 'New: '+str(buff[0])
	     oper1 = c3+'='+'if('+c2+'==1,'+str(buff[0])+',null())'
	     c = c3 + ',' + c
	     grass.run_command('r.mapcalc',expression=oper1,overwrite=True)
	  if clean=='T':
	   print c
	   grass.run_command('r.patch',input=c,out=pa3,overwrite=True)
	   bv = grass.read_command('r.stats',input=pa3,flags='nc',separator='\n').splitlines()
	   print bv

	  b = grass.read_command('r.stats',input=pa3,flags='nc',sort='desc',separator='\n').splitlines()
	  print b
	  for g in np.arange(1,len(b),2):
	   if np.int(b[g]) < minarea:#/10: # lower the threshold if omitting min area!
	    print 'Cleaning small segments III...'
	    print 'cleaning cat '+ str(b[g-1])
	    oper1 = pa3+'='+'if('+pa3+'=='+str(b[g-1])+','+str(b[0])+','+pa3+')'
	    grass.run_command('r.mapcalc',expression=oper1,overwrite=True)
	    bv = grass.read_command('r.stats',input=pa3,flags='nc',separator='\n').splitlines()
	    print bv

	  #grass. message ("Number of cells per segment")
	  #grass.run_command('r.stats',input=pa3,out=pa5,overwrite=True) # flags='nc'
	  #grass. message ("converting to vector")
	  grass.run_command('r.to.vect', input=pa3,out=pa4,type ='area',flags='v',overwrite=True)
	  #grass. message ("adding labels to segments")
	  grass.run_command('v.db.addcolumn', map=pa4,col='wdpaid_pa VARCHAR')
	  grass.run_command('v.db.update', map=pa4,col='wdpaid_pa',value=pa)
	  grass.run_command('v.db.addcolumn', map=pa4,col='aleat VARCHAR')
	  grass.run_command('v.db.update', map=pa4,col='aleat',value=aleat)
	  #grass. message ("Checking shapefile")
	  pa44 = pa4
	  pa442 = pa44+'_diss'
	  #grass.run_command('v.clean', input=pa4,out=pa44,tool='rmarea',thres=minaream,overwrite=True)
	  grass.run_command('v.db.addcolumn', map=pa44,col='segm_id numeric')#VARCHAR')
	  grass.run_command('v.db.update', map=pa44,col='segm_id',qcol='wdpaid_pa || cat || aleat')
	  #grass. message ("Exporting shapefile")
	  name = 'shp/park_segm_'+str(pa)+'_'+str(j)
	  if os.path.isfile(name+'.shp') == False:
	   grass.run_command('v.out.ogr',input=pa44,ola=name,type='area',dsn='.') 
	  else:
	   grass.run_command('v.out.ogr',flags='a',input=pa44,ola=name,type='area',dsn='.')

	  grass. message ("Done1")
	  #grass.run_command('v.in.ogr',flags='oe',dsn='.',lay=name,out=name,overwrite=True)
    	  spa_list0 = grass. read_command ('v.db.select', map=pa44,column='segm_id'). splitlines ()
	  spa_list = np.unique(spa_list0)
	  print spa_list
	# save it as a csv excluding last item!

	  grass. message("omitting previous masks")
	  grass.run_command('g.remove', rast='MASK')

	#pa_list_done = np.genfromtxt(csvname1,dtype='string')
	  sn = len(spa_list)-1 # there is also a segm_id element!
	  for spx in range(0,sn): # 0
	   spa = spa_list[spx]
	   spa2 = 'svv'+spa
	 #pa3 = pa+'v3'
	   spa4 = 'spa_'+spa
	   spa5 = 'tiffs/pa_'+spa+'.tif'
	   spa0 = 'sv0_'+spa
	   sopt1 = 'segm_id = '+spa
	 #if pa not in pa_list_done:
	   print spx
	   print "Extracting PA:"+spa
	   grass.run_command('v.extract', input=pa44, out=spa0, where = sopt1,overwrite=True)
	  # try to crop PAs shapefile with coastal line or input vars
	   grass. message ("setting up the working region")
	   grass.run_command('g.region',vect=spa0,res=1000)
	   grass.run_command('v.to.rast',input=spa0,out=spa0,use='val')#use='cat',labelcol='segm_id')
	   soptt = spa4+'='+spa0
	   grass.run_command('r.mask', rast='pre') # new to crop parks to where we have indicators information
	   grass.run_command('r.mapcalc',expression=soptt,overwrite=True) # opt3
	   grass.run_command('g.remove', rast='MASK') # new
	   grass.run_command('r.null',map=spa4,null=0)
	   econame = 'csv/park_'+str(pa)+'_'+str(j)+'.csv'
	   eco = str(j)
	   econ = 'csv/ecoregs'+str(j)+'.csv'
	   grass.run_command('r.out.gdal',input=spa4,out=spa5,overwrite=True)
#	   grass. message ("Deleting tmp layers") 
	   wb = open(econame,'a')
	   wb.write(spa)
	   wb.write('\n')
	   wb.close() 
	   wb = open(econ,'a')
	   wb.write(eco)
	   wb.write('\n')
	   wb.close() 
 grass. message ("Deleting tmp layers")
# grass.run_command('g.mremove',typ='rast',patt='sv3',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='sv2',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='sv0_*',flags='f') 
# grass.run_command('g.mremove',typ='vect',patt='sv0_*',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='svv*',flags='f')
# grass.run_command('g.mremove',typ='rast',patt='spa_*',flags='f')
# grass.run_command('g.mremove',typ='rast',patt='*b50km',flags='f')
# grass.run_command('g.mremove',typ='rast',patt='old*',flags='f')
# grass.run_command('g.mremove',typ='rast',patt='new*',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='*b10km',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='*v3',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='*v2',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='v0_*',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='pa*',flags='f') 
# grass.run_command('g.mremove',typ='vect',patt='v0_*',flags='f') 
# grass.run_command('g.mremove',typ='vect',patt='pa*',flags='f') 
# grass.run_command('g.mremove',typ='vect',patt='paa_*',flags='f') 
# grass.run_command('g.mremove',typ='vect',patt='paa_pca*',flags='f')
# grass.run_command('g.mremove',typ='rast',patt='*v3',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='*v2',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='v0_*',flags='f') 
# grass.run_command('g.mremove',typ='vect',patt='v0_*',flags='f') 
# grass.run_command('g.mremove',typ='rast',patt='vv*',flags='f')

 grass.run_command('g.mremove',typ='vect',patt='*',flags='f') 
 grass.run_command('g.mremove',typ='rast',patt='*',flags='f')

 wb = open(csvname1,'a')
 var = str(pa)
 wb.write(var)
 wb.write('\n')
 wb.close() 

	  #grass.run_command('v.to.rast', input=pa44,out=pa44,use='cat',overwrite=True)
	  #os.system('ogr2ogr parks_segmented1_dissolved.shp parks_segmented1.shp -dialect sqlite -sql "SELECT ST_Union(geometry), segm_id FROM parks_segmented1 GROUP BY segm_id"')
	  #os.system('rm parks_segmented1.*')
	  #grass.run_command('r.univar',input='pre',zones=pa44,flags='t',separator='space').splitlines()
	  #r.univar -t map=pre zones=paa_pca_19297 sep=space|cut -f10 -d' '
#print "Done PA1:"+pa 

pool = Pool(9)
pool.map(fsegm,pa_list)
pool.close()
pool.join()

print str(datetime.now())
print 'END'
