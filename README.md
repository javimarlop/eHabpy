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
   `/mnt/d`, ...). Do **not** use `sudo` for `git clone`, `cp` or `conda`:
   files created with `sudo` belong to `root` and you cannot edit them later.

   ```
   sudo apt update && sudo apt install -y git
   git clone https://github.com/javimarlop/eHabpy.git ~/eHabpy
   # or: cp -r /mnt/c/Users/<your_windows_user>/Documents/eHabpy ~/eHabpy
   cd ~/eHabpy
   conda env create -f environment.yml
   ```

5. Keep the code, input data and GRASS database inside the Linux file system
   (e.g. `~/eHabpy`, `~/grassdata`): working on `/mnt/c/...` is much slower.
   The GRASS database is just an empty folder that will hold your GRASS
   projects; it does not exist until you create it:

   ```
   mkdir -p ~/grassdata
   ```

   This is the folder `GRASSDB` in `config.py` points to (`'~/grassdata'`).
   Do not confuse it with the GRASS *installation* folder printed by
   `grass --config path` (e.g. `~/miniforge3/envs/ehabpy/lib/grass85`), which
   the scripts find by themselves.

6. To open these files from Windows (e.g. to copy input data or to open the
   results in QGIS), run `explorer.exe .` in the Ubuntu terminal from inside
   `~/eHabpy`, or browse to `\\wsl$\Ubuntu\home\<your_linux_user>`.
   Do not use the copy under `\\wsl$\Ubuntu\mnt\wslg\distro\...`: it is a
   read-only view of the same files used internally by WSL, and there is no
   need to remount it read-write or change its owner with `chown`.
   On Windows 11 the GRASS GIS GUI also runs from Ubuntu with `grass`.

**"Permission denied" or "Read-only file system" when creating or editing
files:** check with `pwd` that you are working in
`/home/<your_linux_user>/eHabpy` and not under `/mnt/wslg/distro/...`. If
`ls -l ~` shows `root` as the owner of `eHabpy` (something was run with
`sudo`), take ownership back with (the same for `~/grassdata` if needed):

```
sudo chown -R $USER:$USER ~/eHabpy
```

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

Open the terminal (the Ubuntu terminal on Windows), change to your
`segm_optim` folder using `cd` and activate the environment:

```
conda activate ehabpy
```

The scripts find the GRASS GIS *installation* by themselves. Only if they fail
with `No module named 'grass'`, also run (this is not your GRASS database):

```
export GISBASE=$(grass --config path)
export PYTHONPATH=$GISBASE/etc/python:$PYTHONPATH
```

## Running

1. Copy all input variables (GeoTIFF, all with the same CRS and pixel size)
   to an `inVars` folder next to `segm_optim` (`../inVars` seen from
   `segm_optim`). The map name of each variable is its file name
   (`precip.tif` → `precip`).

2. Copy the vector file of the study area (`.gpkg`, `.shp` or `.geojson`) to
   the `segm_optim/pa` folder. Its file name is also the name of the GRASS
   project and of the study area map (`pa/guajares.gpkg` → `guajares`), so
   rename the file if you want another name.

3. From the `segm_optim` folder run the setup script:

   ```
   python3 setup.py
   ```

   It creates the GRASS GIS project (called *location* before GRASS 8.4)
   named after the study area in the GRASS database folder `GRASSDB` of
   `config.py` (by default `~/grassdata`), with the CRS of the variables,
   imports the variables and the study area into its `PERMANENT` mapset,
   writes `palist.csv` with all the IDs of the study area and fills the
   required options of `config.py` (`GRASSDB`, `GRASSLOC`, `ENV_VARS0`,
   `STUDY_AREA`, `COL_ID`, `RESOLUTION` and `ENV_VARS`) from the data. Use
   another GRASS database folder with, e.g.:

   ```
   python3 setup.py --grassdb ~/grassdata
   ```

   `python3 setup.py --help` lists the other options: `--study-area` when
   there are several vector files in `pa`, `--col-id` for the column with the
   IDs of the study areas (default `cat`, one ID per polygon) and
   `--overwrite` to import the maps again after replacing a file. Running it
   again only imports what is not in the project yet and keeps an existing
   `palist.csv`. If you were given a ready-made GRASS project, copy its whole
   folder into your GRASS database (e.g. `~/grassdata/guajares`): it is used
   when its name is that of the study area file.

4. Edit `palist.csv` to process only some of the IDs, and check the other
   options in `config.py` (`PA_BUFFER`, `SIMIL`, `STAT_AGG`, ...). The output
   folders of `segm_optim` (`csv`, `shp`, `tiffs`, `results`) are created
   automatically by `segmentation_pca_par.py`.

5. The results of each study area (HFT maps and shapefiles, plots,
   `overall_optim_thresholds.csv` and `processing_metadata.txt`) are saved in
   their own folder, `results/<study area>` (e.g. `results/guajares`), and are
   kept when another study area is processed. The intermediate folders `csv`,
   `shp` and `tiffs` only hold the last study area processed.


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

