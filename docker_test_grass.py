## python grass gis test docker
import os
import sys

GISBASE = os.environ['GISBASE'] = "/usr/lib/grass70"
sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as grass
import grass.script.setup as gsetup

GRASSDBASE = "/home/majavie/hanksgrass7/"#"/opt/grassdb"
MYLOC = "global_MW"
mapset = 'm365152'
gsetup.init(GISBASE, GRASSDBASE, MYLOC, mapset)

print grass.gisenv()
grass.run_command('g.region',flags='p')
#grass.run_command('g.list',typ='vect')
grass.run_command('r.mapcalc',expr='docker = const',overwrite=True)
