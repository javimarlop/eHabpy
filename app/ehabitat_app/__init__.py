"""eHabitat+ app: a window (pywebview) over the scripts of eHabpy.

    paths     folders of the copy of eHabpy the app works on
    envsetup  activation of the conda environment when the app is not opened from an activated one
    jobs      background jobs (scripts of eHabpy and other long tasks) with their log, progress and cancellation
    inputs    input variables (GeoTIFFs of the study area and ecoregion scales) and the study area
    model     parameters of config.py, palist.csv and the steps of the HFT model (run_hft.py)
    gee       MonitorEO: catalog (read again when monitoreo is updated), Earth Engine connection and downloads of the
              study area and its ecoregions or a comparison area (get_vars_gee.py)
    results   results of the HFT model: plots, tables and the HFT map
    ecosim    similarity of the HFTs in their ecoregion or a comparison area (similarity.py): options, plan, runs and
              their maps
    server    local HTTP server of the interface and the result files (and of the API without pywebview)
    api       functions called by the interface
    main      opens the window
"""

__version__ = '0.1.0'
