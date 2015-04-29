import numpy as np
from tqdm import *
import os
import sys
import csv

GISBASE = os.environ['GISBASE'] = "/home/majavie/hierba_hanks/grass-7.1.svn"#grass_last/new/grass7_trunk/dist.x86_64-unknown-linux-gnu"
GRASSDBASE = "/home/majavie/hanksgrass7"
MYLOC = "global_MW"
mapset = 'm1'#ehabitat'
col= 'wdpaid'

sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)

#a = grass.read_command('r.stats',input='MASK',flags='c',separator='\n').splitlines()
#print a
rmk = 'MASK = if( MASK > 0 , MASK , 0)'
print rmk
#grass.run_command('r.mapcalc',expression=rmk,overwrite=True)
grass.run_command('r.mask',flags='r')
#b = grass.read_command('r.stats',input='MASK',flags='c',separator='\n').splitlines()
#print b
