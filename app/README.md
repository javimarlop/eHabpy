# eHabitat+ desktop app

A window (pywebview) over the scripts of `segm_optim`: input variables at the plot and ecoregion scales, download from
MonitorEO (the RESOLVE ecoregions that the study area touches or a comparison area of your own, with the study area at
the pixel of that grid when the area does not contain it; `monitoreo` is kept at its last version of GitHub when the
app opens and its catalog is read again without restarting), options of the HFT model, run and results with an
interactive map of the HFTs, and the similarity of the HFTs in their ecoregion or in a comparison area
(`similarity.py`: pixels of each HFT at the chosen pixel size, maps of the similarity and of the similar areas,
landscape metrics and THDI).

- Install: `installer/macos/install.command` (macOS), `bash installer/linux/install.sh` (Ubuntu and other Linux on
  x86_64, and Windows in WSL2) or `installer\windows\install.bat` (native Windows, where the HFT model does not run
  yet: GRASS GIS 8.5.0 of conda-forge for Windows cannot open its databases).
- Open: the shortcut **eHabitat+** on the Desktop, or `python app/ehabitat.py` in the environment of `environment.yml`.

Documentation: [Desktop app](https://github.com/javimarlop/eHabpy/wiki/Desktop-app) in the wiki.
