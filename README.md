eHabitat+
============

<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10612.svg)](https://doi.org/10.5281/zenodo.10612)

**eHabitat+** GRASS GIS 7 scripts and Python 2.7 library for automatic delineation of habitats within protected areas (PA) and calculation of maps of probabilities to find areas presenting similar ecological characteristics to those found in PA within the corresponding ecoregion. A habitat similarity index (HSI) is computed based on the ratio between the extent of similar areas around the PA and the PA extent, as well as some  landscape metrics and indices to characterize similar areas to PA. Processed results are being updated and can be accessed through the [DOPA Explorer](https://dopa.jrc.ec.europa.eu/en).

## OS setup

- Install Ubuntu 14.04: http://releases.ubuntu.com/trusty/ 
- Add the ubuntugis stable repository: https://launchpad.net/~ubuntugis/+archive/ubuntu/ppa
- Add "deb https://cloud.r-project.org/bin/linux/ubuntu trusty-cran35/" to /etc/apt/sources.list
- Install gdal command line utilities
- GRASS GIS 7 requirements:
	- https://askubuntu.com/questions/474767/installing-grass-gis-7-0-on-ubuntu-14-04
	- https://grasswiki.osgeo.org/wiki/Compile_and_Install_Ubuntu#Current_stable_Ubuntu_version 
- Install GRASS GIS 7.0.6 (compiled from source): https://grass.osgeo.org/grass70/source/
	- Use the conf_grass7eHabplus.sh config file:

```
sh conf_grass7eHabplus.sh # edit target folder
make -j2 # 2 is the number of available processors
sudo make install
```

- Install python-pysal, scipy, scikit-learn, numpy, gdal, python-fiona
- Install R and required libraries:
	- vegan
	- ggplot2
	- rgdal
	- ade4
	- reshape2
	- RColorBrewer

## Running

### Segmentation

```
python segmentation_pca_par.py

python # opens a python environment
from getmevar import *
run_batch_all()
exit()

python moranvar.py
```

### Similarity

```
python subpas_loop_segm_optim.py # pas folder

python # opens a python environment
from ehab_optim import *
run_batch()
exit()
```

