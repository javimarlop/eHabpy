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
Python packages) is **conda / mamba**.

> **Windows users:** GRASS GIS is not available from conda-forge for Windows,
> so the scripts cannot find the `grass` Python package in a native Windows
> conda environment (`ModuleNotFoundError: No module named 'grass'`).
> **Install WSL2 first** (see [Windows: WSL2 setup](#windows-wsl2-setup) below)
> and then follow the Linux instructions from inside the Ubuntu terminal.

### Linux / macOS

You need to install [Anaconda](https://www.anaconda.com/download) (or
[Miniforge](https://github.com/conda-forge/miniforge)) in your system.

```
conda env create -f environment.yml
```

### Windows: WSL2 setup

WSL2 runs a real Linux (Ubuntu) inside Windows 10 (version 2004+) / 11, so the
Linux environment (including GRASS GIS) works unchanged.

1. Open **PowerShell as Administrator**, install WSL2 with Ubuntu and restart
   the computer when asked:

   ```
   wsl --install -d Ubuntu
   ```

2. Open **Ubuntu** from the Start menu and create your Linux user and password.
   All the following commands are run in this Ubuntu terminal, not in
   PowerShell.

3. Install Miniforge (conda) **inside Ubuntu** (a Windows Anaconda install is
   not visible from WSL), accept the defaults, and close and reopen the
   terminal:

   ```
   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
   bash Miniforge3-Linux-x86_64.sh
   ```

4. Get the code into your Linux home folder, either cloning the repository or
   copying an existing Windows folder (Windows drives are under `/mnt/c`,
   `/mnt/d`, ...):

   ```
   sudo apt update && sudo apt install -y git
   git clone https://github.com/javimarlop/eHabpy.git ~/eHabpy
   # or: cp -r /mnt/c/Users/<your_windows_user>/Documents/eHabpy ~/eHabpy
   cd ~/eHabpy
   conda env create -f environment.yml
   ```

5. Keep the code, input data and GRASS database inside the Linux file system
   (e.g. `~/eHabpy`, `~/grassdata`): working on `/mnt/c/...` is much slower.
   Set `GRASSDB` in `config.py` to a Linux path, e.g. `'~/grassdata'`.
   You can browse these files from the Windows Explorer at
   `\\wsl$\Ubuntu\home\<your_linux_user>` (e.g. to open results in QGIS).
   On Windows 11 the GRASS GIS GUI also runs from Ubuntu with `grass`.

From here on, follow the same steps as Linux users.

### Manual install

Alternatively, install the pieces yourself:

- **GRASS GIS 8** (8.3+) — from your distribution's packages, from
  https://grass.osgeo.org/download/ , or `conda install -c conda-forge grass`.
- **GDAL 3** command-line utilities and Python bindings (`osgeo`).
- **Python 3.10–3.12** packages (see `requirements.txt`):
  `pip install -r requirements.txt`
  (numpy, scipy, joblib, tqdm, geopandas, libpysal, esda; GDAL via conda/system).
- **R 4.x** with: `sf`, `terra`, `vegan`, `ade4`, `ggplot2`, `reshape2`,
  `RColorBrewer`.

### Before running

Before running the GRASS steps, point the scripts at your GRASS database, e.g.:

First, open the terminal (the Ubuntu terminal on Windows) and change to your
`segm_optim` folder using `cd`.

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


### Habitat Functional Types (HFTs) (_segm_optim_ folder)

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

### HFTs similarity and landscape metrics (to be tested under the new environment)

```
python3 subpas_loop_segm_optim.py # move to 'pas' folder

python3 # opens a python environment
from ehab_optim import * # alternative: from ehab_optim_median import *
run_batch() # it runs using all available processors in parallel
exit()

Rscript thdi.R
```

