"""Python core of MonitorEO (OBSNEV): the variables of the MonitorEO app of Google Earth Engine
(https://ee-ugr-default-project.projects.earthengine.app/view/monitoreo-obsnev), computed with the Earth Engine
Python API and downloaded as GeoTIFF on a given grid.

    catalog    catalog of variables of MonitorEO (catalog.json)
    sensors    collections of each sensor: quality masks, scale and offset, formula of the variable
    composite  temporal aggregation: one metric over the period, or by year, month or 16 days and then between years
    variable   Earth Engine image of a variable described by a dictionary (a 'spec')
    areas      study area, ecoregions (RESOLVE Ecoregions 2017)
    grid       grids aligned to the pixel size and download of an image in tiles
    drive      export of an image to Google Drive (batch task) and import of the exported files onto a grid

It needs the Earth Engine API (earthengine-api) and an Earth Engine project: run ee.Initialize(project=...) first.
"""
