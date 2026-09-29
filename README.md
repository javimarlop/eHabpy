eHabitat+
============

<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5643271.svg)](https://doi.org/10.5281/zenodo.5643271)

[**eHabitat+**](https://www.sciencedirect.com/science/article/pii/S157495412300119X) GRASS GIS scripts and Python library for automatic delineation of habitats within protected areas (PA) and calculation of maps of probabilities to find areas presenting similar ecological characteristics to those found in PA within the corresponding ecoregion. A habitat similarity index (HSI) is computed based on the ratio between the extent of similar areas around the PA and the PA extent, as well as some  landscape metrics and indices to characterize similar areas to PA. Processed results are being updated and can be accessed through the [DOPA Explorer](https://dopa.jrc.ec.europa.eu/en).

> **Modernized (2024):** the code now runs on **Python 3.10–3.12**, **GRASS GIS 8**,
> **GDAL 3**, and a current **R** stack (`sf` instead of the retired `rgdal`,
> `libpysal`/`esda` instead of PySAL 1.x). See [`MODERNIZATION.md`](MODERNIZATION.md)
> for the full list of changes. The old Python 2.7 / GRASS 7 sources remain in
> the git history.

## OS setup

The recommended way to get the full geospatial stack (GDAL, GRASS, R and all
Python packages) is **conda / mamba**:

You need to install [Anaconda](https://www.anaconda.com/download) in your system.

```
conda env create -f environment.yml
conda env create -f environment_win.yml # Windows users
# conda activate ehabpy
```

Alternatively, install the pieces yourself:

- **GRASS GIS 8** (8.3+) — from your distribution's packages, from
  https://grass.osgeo.org/download/ , or `conda install -c conda-forge grass`.
- **GDAL 3** command-line utilities and Python bindings (`osgeo`).
- **Python 3.10–3.12** packages (see `requirements.txt`):
  `pip install -r requirements.txt`
  (numpy, scipy, joblib, tqdm, geopandas, libpysal, esda; GDAL via conda/system).
- **R 4.x** with: `sf`, `terra`, `vegan`, `ade4`, `ggplot2`, `reshape2`,
  `RColorBrewer`.

Before running the GRASS steps, point the scripts at your GRASS database, e.g.:

First, open the terminal and change to your `segm_optim` folder using `cd`.

```
conda activate ehabpy

export GISBASE=$(grass --config path)
export PYTHONPATH=$GISBASE/etc/python:$PYTHONPATH

# export GRASSDBASE=/path/to/grassdata     # your GISDBASE
# export GRASSLOC=global_MW                # the Mollweide location
```

(The obsolete `conf_grass7eHabplus.sh`, which compiled GRASS 7.0.6 from source
on Ubuntu 14.04, is no longer needed and is kept only for reference.)

## Running

1. You need to create the following folders within the `segm_optim` folder:

- csv
- shp
- tiffs

2. Point the script to the GRASS GIS database and location.

3. Create a `palist.csv` file with the list of IDs that you will process.

4. Edit the list of input variables in `segmentation_pca_par.py`, if necessary.


### Segmentation (_segm_optim_ folder)

```
ulimit -n 8192

python3 segmentation_pca_par.py # it uses parallel processing

export PYTHONNOUSERSITE=1

python3 # opens a python environment
from getmeanvar import * # alternative: from getmedianvar import *
run_batch_all()
exit()

python3 moranvar.py    # calls Rscript moranvar_plots.R at the end
```

### Similarity (to be tested)

```
python3 subpas_loop_segm_optim.py # move to 'pas' folder

python3 # opens a python environment
from ehab_optim import * # alternative: from ehab_optim_median import *
run_batch() # it runs using all available processors in parallel
exit()

Rscript thdi.R
```

