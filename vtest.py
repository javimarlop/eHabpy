from osgeo import ogr
import os

# Get the input Layer
inShapefile = "pas/ecoregs_moll_80601.shp"
inDriver = ogr.GetDriverByName("ESRI Shapefile")
inDataSource = inDriver.Open(inShapefile, 0)
inLayer = inDataSource.GetLayer()

# Create the output Layer
outShapefile = "parks_centroids.shp"
outDriver = ogr.GetDriverByName("ESRI Shapefile")

# Remove output shapefile if it already exists
if os.path.exists(outShapefile):
    outDriver.DeleteDataSource(outShapefile)

# Create the output shapefile
outDataSource = outDriver.CreateDataSource(outShapefile)
outLayer = outDataSource.CreateLayer("parks_centroids", geom_type=ogr.wkbPoint)

# Add input Layer Fields to the output Layer
inLayerDefn = inLayer.GetLayerDefn()
for i in range(0, inLayerDefn.GetFieldCount()):
    fieldDefn = inLayerDefn.GetFieldDefn(i)
    outLayer.CreateField(fieldDefn)

# Get the output Layer's Feature Definition
outLayerDefn = outLayer.GetLayerDefn()

# Add features to the ouput Layer
for i in range(0, inLayer.GetFeatureCount()):
    # Get the input Feature
    inFeature = inLayer.GetFeature(i)
    # Create output Feature
    outFeature = ogr.Feature(outLayerDefn)
    # Add field values from input Layer
    for i in range(0, outLayerDefn.GetFieldCount()):
        outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
    # Set geometry as centroid
    geom = inFeature.GetGeometryRef()
    centroid = geom.Centroid()
    outFeature.SetGeometry(centroid)
    # Add new feature to output Layer
    outLayer.CreateFeature(outFeature)

# Close DataSources
inDataSource.Destroy()
outDataSource.Destroy()
