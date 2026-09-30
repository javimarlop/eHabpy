eHabitat+
============

<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5643271.svg)](https://doi.org/10.5281/zenodo.5643271)

[**eHabitat+**](https://www.sciencedirect.com/science/article/pii/S157495412300119X) GRASS GIS scripts and Python library for automatic delineation of habitats within protected areas (PA) and calculation of maps of probabilities to find areas presenting similar ecological characteristics to those found in PA within the corresponding ecoregion. A habitat similarity index (HSI) is computed based on the ratio between the extent of similar areas around the PA and the PA extent, as well as some  landscape metrics and indices to characterize similar areas to PA.

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
```

## Running

1. You need to create the following folders within the `segm_optim` folder:
	- csv
	- shp
	- tiffs
	- pa
2. Locate the vector file with the study area under the `pa` folder.

3. Locate all input variables under the `../inVars` folder.

4. Import all input variables plus the vector file of the study area into a GRASS GIS location (PERMANENT mapset).

5. Point the script to the GRASS GIS database and location.

6. Create a `palist.csv` file with the list of IDs that you will process.

7. Edit the `config.py` file with the required and optional parameters.


### HFTs (_segm_optim_ folder)

```
ulimit -n 8192

python3 segmentation_pca_par.py # it uses parallel processing

export PYTHONNOUSERSITE=1

python3 # opens a python environment
from get_mean_median_var import * 
run_batch_all()
exit()

python3 moranvar.py    # calls Rscript moranvar_plots.R at the end
```

### Similarity (to be tested under the new environment)

```
python3 subpas_loop_segm_optim.py # move to 'pas' folder

python3 # opens a python environment
from ehab_optim import * # alternative: from ehab_optim_median import *
run_batch() # it runs using all available processors in parallel
exit()

Rscript thdi.R
```

