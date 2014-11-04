import os
import sys
import numpy as np
from tqdm import *

#gisbase = os.environ['GISBASE'] = "/usr/local/grass-7.0.svn"
gisbase = os.environ['GISBASE'] = "/home/majavie/grass_last/new/grass7_trunk/dist.x86_64-unknown-linux-gnu"
gisdbase = os.path.join("/local1/majavie/hanksgrass7")
location = "global_MW"
mapset   = "ehabitat"#"rasterized_parks"

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as grass
import grass.script.setup as gsetup
 
gsetup.init(gisbase,
            gisdbase, location, mapset)
 
print grass.gisenv()

source = 'parks_segmented_filter@ehabitat'
#grass.run_command('v.in.ogr',flags='oe',dsn='.',lay=source,out=source,overwrite=True)
grass. message ("Extracting list of PAs")
pa_list0 = grass. read_command ('v.db.select', map=source,column='segm_id'). splitlines ()
pa_list = np.unique(pa_list0)
print pa_list
# save it as a csv excluding last item!

grass. message ("Deleting tmp layers") 
grass.run_command('g.mremove',typ='rast',patt='*v3',flags='f') 
grass.run_command('g.mremove',typ='rast',patt='*v2',flags='f') 
grass.run_command('g.mremove',typ='rast',patt='v0_*',flags='f') 
grass.run_command('g.mremove',typ='vect',patt='v0_*',flags='f') 
grass.run_command('g.mremove',typ='rast',patt='vv*',flags='f')

grass. message("omitting previous masks")
grass.run_command('g.remove', rast='MASK')

csvname1 = 'pas_segm_tiff_done.csv'
if os.path.isfile(csvname1) == False:
 wb = open(csvname1,'a')
 wb.write('None')
 wb.write('\n')
 wb.close()

pa_list_done = np.genfromtxt(csvname1,dtype='string')
n = len(pa_list)-2 # there is also a segm_id element!
for px in tqdm(range(0,n)): # 0
 pa = pa_list[px]
 pa2 = 'vv'+pa
 #pa3 = pa+'v3'
 pa4 = 'pa_'+pa
 pa5 = 'pa_'+pa+'.tif'
 pa0 = 'v0_'+pa
 opt1 = 'segm_id = '+pa
 if pa not in pa_list_done:
  print px
  print "Extracting PA:"+pa
  grass.run_command('v.extract', input=source, out=pa0, where = opt1,overwrite=True)
  # try to crop PAs shapefile with coastal line or input vars
  grass. message ("setting up the working region")
  grass.run_command('g.region',vect=pa0,res=1000)
  grass.run_command('v.to.rast',input=pa0,out=pa0,use='val')#use='cat',labelcol='segm_id')
  #opt3 = pa2+'= @'+pa0
  #opt4 = pa2+'= round('+pa2+')'
  #grass.run_command('r.mapcalc',expression=opt3,overwrite=True)
  #grass.run_command('r.mapcalc',expression=opt4,overwrite=True)
  #grass.run_command('r.mask', vector=pa0, where=opt1)
  #opt2 = pa4+'='+pa2
  #opt22 = pa4+'='+pa4
  optt = pa4+'='+pa0
  #grass.run_command('r.mapcalc',expression=opt2,overwrite=True)
  #grass.run_command('g.remove', rast='MASK')
  grass.run_command('r.mask', rast='pre') # new to crop parks to where we have indicators information
  grass.run_command('r.mapcalc',expression=optt,overwrite=True) # opt3
  grass.run_command('g.remove', rast='MASK') # new
  grass.run_command('r.null',map=pa4,null=0)
  eco_list = grass.read_command ('r.stats', input='ecoregs_moll',sort='desc'). splitlines ()
  print eco_list
  eco = eco_list[0]
  if len(eco_list)>1 and eco == '*': eco = eco_list[1]
  if len(eco_list)>1 and eco == '*':
   grass.run_command('g.region',res=10)
   eco_list = grass.read_command ('r.stats', input='ecoregs_moll',sort='desc'). splitlines ()
   eco = eco_list[0]
   if len(eco_list)>1 and eco == '*': eco = eco_list[1]
  print eco
  econame = str(eco)+'.csv'
  grass.run_command('g.region',res=1000)
  grass.run_command('r.out.gdal',input=pa4,out=pa5,overwrite=True)
  grass. message ("Deleting tmp layers") 
  grass.run_command('g.mremove',typ='rast',patt='*v3',flags='f') 
  grass.run_command('g.mremove',typ='rast',patt='*v2',flags='f') 
  grass.run_command('g.mremove',typ='rast',patt='v0_*',flags='f') 
  grass.run_command('g.mremove',typ='vect',patt='v0_*',flags='f') 
  grass.run_command('g.mremove',typ='rast',patt='vv*',flags='f')
  wb = open(econame,'a')
  wb.write(pa)
  wb.write('\n')
  wb.close() 
  wb = open('ecoregs.csv','a')
  wb.write(eco)
  wb.write('\n')
  wb.close() 
  wb = open(csvname1,'a')
  var = str(pa)
  wb.write(var)
  wb.write('\n')
  wb.close() 
  grass. message ("Done")
  print "Done PA:"+pa

# try to paralellize it?
