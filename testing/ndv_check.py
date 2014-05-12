from osgeo import gdal
import numpy as np

# Read the input raster into a Numpy array
infile = "results/2607_eco.tif"
data0   = gdal.Open(infile)
data = data0.GetRasterBand(1)
#arr    = data.ReadAsArray()

# Do some processing....

# Save out to a GeoTiff

# First of all, gather some information from the original file
#[cols,rows] = arr.shape
trans       = data0.GetGeoTransform()
print trans
proj        = data0.GetProjection()
print proj
ndv = data.GetNoDataValue()
print ndv
