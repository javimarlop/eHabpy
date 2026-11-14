<p align="center"><img src="app/ehabitat_app/ui/icon.png" width="160" alt="eHabitat+ logo"></p>

eHabitat+
============

<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.595844.svg)](https://doi.org/10.5281/zenodo.595844)

[**eHabitat+**](https://www.sciencedirect.com/science/article/pii/S157495412300119X) GRASS GIS scripts and Python library for automatic delineation of habitats within protected areas (PA) and calculation of maps of probabilities to find areas presenting similar ecological characteristics to those found in PA within the corresponding ecoregion. A habitat similarity index (HSI) is computed based on the ratio between the extent of similar areas around the PA and the PA extent, as well as some  landscape metrics and indices to characterize similar areas to PA.

**Documentation: [wiki](https://github.com/javimarlop/eHabpy/wiki)** (installation, the desktop app for macOS,
Windows and Ubuntu, the tutorial of the Habitat Functional Types, their similarity in the ecoregion with the landscape
metrics and the THDI, the outputs and how the model works). The diagram of the whole workflow is on the
[home page of the wiki](https://github.com/javimarlop/eHabpy/wiki#workflow).

**What's new: [changelog](#changelog)**, what each release brings, at the end of this page.

Cite eHabitat+ as https://doi.org/10.5281/zenodo.595844 (the DOI of all its versions, which leads to the last one;
the DOI of each release is on its page of Zenodo)

---

## Changelog

What each release brings, newest first. A release is a tag of GitHub with its
[release](https://github.com/javimarlop/eHabpy/releases), archived in Zenodo. Each entry names the parts of eHabitat+
it reaches (the **Desktop app**, its **Installers**, the **HFT model**: `setup.py`, `run_hft.py` and their steps, the
**Similarity**: `similarity.py` and `thdi.R`, the **Earth Engine variables**: `get_vars_gee.py`, and the
**Documentation**), and **(fix)** marks an error that was fixed. The desktop app updates itself from `master` every
time it opens, so it already has what is under **Unreleased** before its release.

### Unreleased

- **Desktop app** (fix): The date of the version shows also in a copy that was updated before the app kept it: the app reads it from GitHub when it is missing.
- **Documentation**: Every release in Zenodo has the same title, "eHabitat+ model: mapping habitat functional types and similar areas", its author with ORCID and the licence of eHabitat+, CC BY-SA 3.0 (`.zenodo.json`), instead of the name of the release, the authors of GitHub and CC BY 4.0.
- **Documentation**: The README and the wiki cite eHabitat+ with the DOI of all its versions, https://doi.org/10.5281/zenodo.595844, which leads to the last one, instead of the DOI of version 1.3.

### [1.5](https://github.com/javimarlop/eHabpy/releases/tag/1.5) · 2026-11-14

- **Desktop app**: The version, at the bottom of the side bar, shows the date its commit was published next to the commit and the branch (e.g. `def0624 · master · 2026-11-14`).
- **Documentation**: The changelog of eHabitat+ is the section Changelog at the end of the README, which every change updates, instead of `changelog.txt`.
- **Desktop app**: In Download from MonitorEO, the button **GEE Community Catalog…** adds to the recipe any of the 1,068 images and collections of images of the [GEE Community Catalog](https://gee-community-catalog.org), found by theme or by words of their title, ID, provider or tags. **Read** shows the bands of any ID of Earth Engine, whether it is a fixed image or a collection with dates, and the scale, offset and unit of the band when its record gives them; a collection can have its own period, and the unit of the values goes into the recipe.
- **Desktop app** (fix): Editing a variable of another collection or image of Earth Engine keeps its own season, pixel filters and conditions, which were lost.
- **Desktop app, Earth Engine variables**: A land cover mask leaves out of every variable the pixels of some land covers (water, trees, grassland, wetlands, cropland, shrubland, built-up, bare ground, snow and ice) after a land cover map and a year of your choice: ESA WorldCover (2020, 2021), Dynamic World (every year since 2016), CORINE Land Cover (Europe, 1990 to 2018) or MODIS MCD12Q1 (2001 to 2025). In the recipe it is `where: {landcover: CORINE, year: 2018, exclude: [water, built]}`; the classes of ESA WorldCover 2021 of before (`exclude_landcover`) still work.
- **HFT model**: `hft_means_<id>_<k>.csv` and `hft_stats_<id>_<k>.csv` have the pixel size of the variables (column `resolution`): HFTs made with other pixel sizes are not directly comparable.
- **HFT model, Desktop app**: Tables of the HFTs to compare study areas, in the folder of the scenario: `hft_means_<id>_<k>.csv`, the mean of the pixels of every HFT in each variable, and `hft_stats_<id>_<k>.csv`, their n, mean, sd, quartiles, minimum and maximum, in the original units, so that the tables of several study areas can be joined for an NMDS or a tree of all their HFTs (example in the wiki). The app shows `hft_means` in Results.
- **Desktop app**: Every map has the basemap Sentinel-2 cloudless 2023 (EOX).
- **HFT model, Similarity**: Every plot comes with a `.csv` of the same name with the data it draws.
- **Similarity** (fix): Global ecoregion files (e.g. all of RESOLVE) are read only near the grid, as they could not be reprojected far from it; ecoregions without an ID field get the same IDs in the maps and the tables, and `--id` keeps the other IDs of a folder of results.
- **Similarity, Earth Engine variables**: The HFTs can be compared with a comparison area of your own, inside or outside the study area (`ECOSIM_AREA`, `--area`): a folder of variables with the same names, downloaded with its own grid, and with the study area at its pixel when the area does not contain it.
- **Desktop app**: New page Similarity: the HFTs of a scenario and cut, the pixel size, the options, the run with its log, and the results (THDI, landscape metrics, the map of the similarity and of the similar areas of every HFT, plots).
- **Similarity**: `thdi.R` takes the results of `similarity.py` as arguments (`--hri`, `--hfts`, `--ecoregions`, `--marine`, `--out`), with the same indices, ecoregions in any CRS and NA values.
- **Similarity** (fix): Compared with `ehab_optim.py`, the pixels outside the ecoregion have no data instead of a similarity, the HFT and the ecoregion pixels are on the same grid, no random noise is added to the values of the HFTs (a variable constant in an HFT is left out of the distance, and the run says which), and the values of an HFT are no longer inherited from the previous one.
- **Similarity**: The similarity can use a coarser pixel (`ECOSIM_RESOLUTION`, a multiple of that of the ecoregion variables, aggregated by their mean), and the mean and covariance of every HFT come from the pixels of the similarity grid it covers, so that the HFT and the ecoregion are measured with the same pixel and product.
- **Similarity**: New `segm_optim/similarity.py`, after `run_hft.py`, which replaces `pas/subpas_loop_segm_optim.py` and `ehab_optim.py`: it compares every HFT of a cut with every pixel of its ecoregion (the same variables at the ecoregion scale), writes the similarity map of every HFT, the map of its similar areas and `hri_results.csv` with the landscape metrics of eHabitat+, and runs `thdi.R` and `similarity_plots.R`. New options `ECOSIM_*` in `config.py`.
- **HFT model** (fix): On Windows, `setup.py` and `segmentation_pca_par.py` set `GISBASE` and `HOME` as GRASS GIS needs them: its modules missed their colour tables and the vector modules failed.
- **Installers**: The installer for Windows works on a new Windows 11 (the Visual C++ runtime, short paths for the packages, terra from CRAN and OpenMP for GRASS GIS 8.5). GRASS GIS 8.5 of conda-forge for Windows cannot open databases yet, so the HFT model runs in WSL2 with the installer for Ubuntu, where the app opens the web browser and the folders of Windows.
- **Desktop app**: The study areas to process (the IDs of `palist.csv`) are chosen in HFT mapping, with their names and areas, and the results of each one are shown with its outline highlighted on the map.
- **Installers**: Installer for Ubuntu LTS and other Linux distributions on x86_64 (`app/installer/linux/install.sh`), without sudo, with entries in the applications menu and on the Desktop; without a toolkit for its window, the app opens in the web browser.
- **Desktop app**: The app updates itself from GitHub every time it opens, keeping the files you changed (`config.py`, `gee_vars.yaml`, ...), and updates its environment when `environment.yml` changes; a button checks for updates and a switch turns them off (or `EHAB_NO_UPDATE=1`).
- **Documentation**: The documentation moved from the README to the [wiki](https://github.com/javimarlop/eHabpy/wiki), with a page for the desktop app.
- **Earth Engine variables**: Without an Earth Engine project, `get_vars_gee.py` uses the shared account of MonitorEO, which only downloads directly (not through Google Drive).
- **HFT model**: `moranvar_plots.R` writes the colour of every HFT of each cut to `hft_colours_<id>_<k>[_<cut>].csv`, for the map of the desktop app.
- **Desktop app, Installers**: New desktop app (`app/`): the input variables at the plot and ecoregion scales with their checks, the download of variables from MonitorEO, every option of `config.py` with its default, the runs of `setup.py` and `run_hft.py` with their log, and the plots, tables and an interactive map of the HFTs of each cut. Its installers for macOS and Windows create its environment with micromamba, without administrator rights, with a shortcut on the Desktop.
- **HFT model** (fix): R no longer prints "was not checked for suspicious field assignments" when terra loads (`environment.yml` adds r-codetools), and the radar plots no longer print "Ignoring unknown labels: size".
- **HFT model** (fix): `moranvar_plots.R` no longer fails with "Write error" when it writes the shapefile of the HFTs with variable names that only differ at the end.
- **Earth Engine variables**: `get_vars_gee.py` uses the Python package of [MonitorEO](https://github.com/monitoreo-obsnev/monitoreo) instead of its own copy, so it computes the variables as the QGIS plugin, the R package and the MCP server of MonitorEO do; every run installs its last version from GitHub (`--no-update` skips it).

### [1.4](https://github.com/javimarlop/eHabpy/releases/tag/1.4) · 2026-10-27

- **Earth Engine variables**: New `segm_optim/get_vars_gee.py` and its recipe `gee_vars.yaml`: it downloads the input variables from Google Earth Engine with the catalog of MonitorEO, for the study area with a fine pixel (`inVars/<folder>`) and for the RESOLVE ecoregions it touches with a coarser one (`inVars/<folder>/ecoregion`): variables of the catalog, other collections or fixed images (e.g. a DEM with `terrain: slope`), with their period, season, pixel filters (clouds, water, snow, aerosols) and one or several metrics. Only the files whose recipe or grid changed are downloaded again, and large grids can go through Google Drive (`drive: true`).
- **HFT model**: New `CLIP_PERCENTILE` option in `config.py` (default 0): every layer of the segmentation and every variable of the classification is clipped to its percentiles, so that a few extreme pixels do not set how much each one weighs; the statistics, the threshold selection and the plots keep the original values.
- **HFT model**: New plot `cor_<id>_vars.png`: the correlation between every pair of variables in the pixels of the segments, to see which variables repeat the same gradient.
- **HFT model**: New `INVARS_FOLDER` option in `config.py`: a subfolder of `inVars` with another set of variables, to try it without moving files; `setup.py` imports it with the name of the folder in the names of the maps, and an empty `SCENARIO` takes its name, so that every set of variables has its own results.
- **HFT model**: `moranvar_plots.R` no longer draws the NMDS of the segments, whose r² of the variables was easily misread as their importance: the eta² of the heatmap gives it, and the PCA biplot shows how the HFTs are distributed.
- **HFT model**: New tree of the HFTs of the cuts with `SIMIL = 3` (`hft_tree_<id>_<k>.png` and `.csv`): which HFTs of a cut belong to the same HFT in the broader cuts and how each one splits in the finer ones.
- **HFT model**: `moranvar_plots.R` suggests a cut, the one whose number of HFTs is closest to the bend of the R² curve of the selected threshold (beyond it, more HFTs add little explained variance); its files carry `_optim` after the name of the cut and its figures say "suggested".
- **HFT model**: Hierarchical colours of the HFTs with `SIMIL = 3` ("tree colors" in HCL): every HFT of the broadest cut gets a part of the colour circle and the HFTs inside it share it out, distinguishable with up to about 10 HFTs, also with colour blindness; the legends follow the order of the dendrogram.
- **HFT model**: With `SIMIL = 3` the tree is also cut at 10 and 90 % of its height (`p10` and `p90`, besides `p25`, `p50` and `p75`).
- **HFT model**: New `CLASS_METHOD` option in `config.py`: `ward_area` (default), Ward's method weighted by the pixels of each segment, which gave HFTs with a higher R² and silhouette than the original method in the Guajares data; `ward_area_kmeans`, followed by k-means; or `ward`, the original.
- **HFT model**: New plot `r2_<id>_hfts.png`: the variance of the pixel values explained by the HFTs when the tree of every threshold is cut into 2 to 30 or more HFTs; the plot of the thresholds shows the score of the method and the threshold selected.
- **HFT model**: New `THRESHOLD_METHOD`, `THRESHOLD_ALPHA` and `MORAN_NEIGHBOURS` options in `config.py`: `moranvar.py` chooses the segmentation threshold with the F-function of Johnson et al. (2015) (`f`, the default), the Global Score of Espindola et al. (2006) (`gs`, `gs_fixed`), the original crossing of the curves (`crossing`, `crossing_first`) or a fixed threshold, and writes the score of every threshold to `csv/<id>_threshold.csv`.
- **HFT model** (fix): `ehab_segm_core.py` no longer adds a random noise to the pixel values before computing the statistics of the segments, which raised their means and their variance.
- **HFT model**: New `SEGM_PCA` option in `config.py` (default `True`): with `False` the variables themselves are segmented instead of the first three axes of their PCA.
- **HFT model**: Scenarios: the `SCENARIO` option of `config.py` writes the results to `results/<STUDY_AREA>/<SCENARIO>`; `run_hft.py` shows it, lets you change it (`--scenario`) and warns when the scenario has results made with other options.
- **Documentation**: Tutorial of the HFT module, the interpretation of all the outputs and how the model works.
- **HFT model**: `moranvar_plots.R` moves the previous results of the study area to `results/<STUDY_AREA>_previous`, so that the files of runs with different options are not mixed; `moranvar.py` stops with an error when `moranvar_plots.R` fails.
- **HFT model**: `setup.py` stops when the variables have a geographic CRS, and imports again the files that changed since it imported them.
- **HFT model** (fix): The statistics of the segments are computed in double precision without rounding, the same on every run, and the pixel offsets are rounded instead of truncated.
- **HFT model** (fix): The working region and the segments are aligned with the grid of the first variable, so that the masks of the segments fall exactly on its pixels.
- **HFT model** (fix): `PA_BUFFER` is in map units, as documented: it grew the region by that number of pixels (5 km instead of 500 m with 10 m pixels).
- **HFT model**: `segmentation_pca_par.py` resumes an interrupted run only when its inputs did not change, and does nothing when it is up to date; `FORCE_RESTART` is now `False` and only forces a new segmentation.
- **HFT model** (fix): `get_mean_median_var.py` and `moranvar.py` can run again without duplicating the rows of their tables.
- **HFT model**: New `run_hft.py`: it runs only the steps that are out of date with `config.py`, `palist.csv`, the input files and the code, and those after them (`--status` says what would run and why, `--from` repeats from a step); a script run by hand says what to run first.
- **HFT model**: New plots of the HFTs: boxplots of the values of all the pixels of every HFT in each variable, heatmaps of their profiles (the median of every HFT in each variable, with the variables ordered by the share of their variance explained by the HFTs, eta²) and PCA biplots with the HFTs and the variables.
- **HFT model**: The HFT maps are drawn with ggplot2 at 300 dpi, with the legend outside the map, coordinates, the outline of the study area in black and a title with the study area, threshold, cut and number of HFTs; the radar plots have the same title, every HFT in numeric order with the colours of the map, and the names of the variables shortened and wrapped instead of cut.
- **HFT model** (fix): `ehab_segm_core.py` no longer leaves out of the HFTs the segments that start in the first row or column of the variables (19 % of the area in the Guajares data).
- **HFT model**: The segments smaller than `minarea` are merged in memory, much faster and without "Too many open files", with the segments within `SMALL_SEGM_RADIUS` pixels (new option, default 5) in repeated passes, and then with the nearest large segment, so that none is left.
- **HFT model**: The results of every study area are saved in `results/<STUDY_AREA>`, and `FORCE_RESTART` only removes those of the study area processed.
- **HFT model**: The valid pixels are those with data in the first variable of `ENV_VARS0` (before, those of `precip` of 0 or more), and the pixel size comes from `RESOLUTION` in `config.py` (before, 10 m).
- **HFT model**: New `segm_optim/setup.py`: it creates the GRASS GIS project from the GeoTIFFs of `inVars` and the study area of `pa/`, writes `palist.csv` and fills the required options of `config.py` from the data.
- **HFT model**: New `SIMIL` option in `config.py`: `SIMIL = 3` cuts the tree of the HFTs at 25, 50 and 75 % of its height, each cut in its folder and with its name in the files, and `processing_metadata.txt` records the options of the run.
- **HFT model** (fix): The optimal threshold is interpolated on the values of the thresholds and is always one of those tested, and `overall_optim_thresholds.csv` keeps one row per study area.
- **HFT model** (fix): `segmentation_pca_par.py` starts from the PERMANENT mapset, creates its output folders and stops with an error instead of hanging when a command of GRASS GIS fails.
- **HFT model** (fix): Resuming with `FORCE_RESTART = False` first removes the leftovers of an interrupted study area.
- **HFT model** (fix): The IDs of the segments no longer collide between thresholds or study areas, and are the same on every run.
- **eHabitat+**: eHabitat+ runs with Python 3.10 to 3.12, GRASS GIS 8, GDAL 3, NumPy up to 2.2 and R 4 with sf instead of rgdal (see `MODERNIZATION.md`).

### [1.3](https://github.com/javimarlop/eHabpy/releases/tag/1.3) · 2021-11-03

- **Similarity**: New habitat diversity and irreplaceability modules and indices.

### [1.2](https://github.com/javimarlop/eHabpy/releases/tag/1.2) · 2021-11-02

- **Similarity**: New script for computing the THDI, the habitat diversity and irreplaceability indices.
- **Similarity**: New script that finds the landscape patches similar to the reference HFTs with their median values.
- **HFT model**: New script that finds the optimal segmentation threshold and the final HFTs with the median values.
- **HFT model**: New options in `moranvar_plots.R`.
- **HFT model**: A fixed seed map for each protected area in the first segmentation step, so that the results are the same with different versions of GRASS GIS.
- **HFT model**: The sub-areas of the protected areas are processed by themselves.
- **HFT model**: The dissolve operation is done in Python.
- **Similarity**: New column names and a new column in the results of the similarity (`hriresults`).
- **Similarity** (fix): The similarity indices are right when a mask is used.
- **Documentation**: Every script is documented, and the README has installation instructions for each operating system and a tutorial.

### [1.1](https://github.com/javimarlop/eHabpy/releases/tag/1.1) · 2016-05-23

- **eHabitat+**: Enhanced segmentation and similarity maps for DOPA, the Digital Observatory for Protected Areas.

### [1.0](https://github.com/javimarlop/eHabpy/releases/tag/1.0) · 2014-06-20

- **eHabitat+**: First release.
