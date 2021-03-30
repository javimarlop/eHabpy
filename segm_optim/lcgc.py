from gdalconst import GA_ReadOnly
from datetime import datetime
from osgeo import ogr,gdal
import numpy as np
import os, sys, csv
from tqdm import *

def lcgc(bhcat='',cluster=''):

 pa_list = np.genfromtxt('results/overall_optim_thresholds2.csv',dtype='int',skip_header=0)#,usecols=(1))#np.genfromtxt('wdpaids_marine.csv',dtype='int')
 #pa_list = pa_listt[,0]
 #res_list = pa_listt[,1]#np.genfromtxt('results/overall_optim_thresholds2.csv',dtype='int',skip_header=0,usecols=(2))
 n = len(pa_list)
 #pa_list = pa_list2[0:n] #n-2
 #print pa_list
 #n2 = len(pa_list)
 print n
 field_name_target = 'hclst_m'
 idd = roads = mval = None

 csvname = 'results/gclc_aver_in.csv' # +str(bhcat)+
 if os.path.isfile(csvname) == False:
  wb = open(csvname,'a')
  wb.write('wdpaid cluster gclc aver')
  wb.write('\n')
  wb.close()
 #f=open('','w')
 #f.write('wdpaid, roads, tot, perc\n')
 csvname1 = 'csv/gclc_'+str(bhcat)+'_in_done_100km2.csv'
 if os.path.isfile(csvname1) == False:
  wb = open(csvname1,'a')
  wb.write('None')
  wb.write('\n')
  wb.close()

 gtif_file = '/home/majavie/DOPA/consensus_full_class_LC/consensus_full_class_'+str(bhcat)+'_moll.tif' # 1.tif
 gtif_dataset = gdal.Open(gtif_file, GA_ReadOnly)

 transform = gtif_dataset.GetGeoTransform()
 xOrigin = transform[0]
 yOrigin = transform[3]
 pixelWidth = transform[1]
 pixelHeight = transform[5]
 banddataraster = gtif_dataset.GetRasterBand(1)
 #segmlist = range(1,9)

 pa_list_done = np.genfromtxt(csvname1,dtype='int')

 for px in range(0,n):
    idd2 = pa_list[px,0]
    res = pa_list[px,1]
    #print idd2,res
    if idd2 not in pa_list_done:
	    #clusterlistcsv = 'csv/park_'+idd2+'_'+res+'.csv'
	    #segmlist = np.genfromtxt(segmlistcsv,dtype='int')
	    #n2 = len(segmlist)
	    #for px in range(0,n2):
	    #idd = segmlist[px]

	    inShapefile = 'results/park_segm_'+str(idd2)+'_'+str(res)+'_class.shp'
	    inDriver = ogr.GetDriverByName("ESRI Shapefile")
	    inDataSource = inDriver.Open(inShapefile, 0)
	    inLayer = inDataSource.GetLayer()
	    
	    #inLayer.GetFieldAsString(field_name_target)
 	    #for cluster in segmlist: # indent1
	    # Get the input Layer
	    print idd2,res,'cluster:',cluster
	    opt = field_name_target+" = '"+str(cluster)+"'"
	    #print opt
	    inLayer.SetAttributeFilter(opt)
	    if inLayer.GetFeatureCount() > 0: # indent2
		    ola = 'pa_tmpfs_'+str(bhcat)+'_'+str(idd2)+'_'+str(cluster)+'.shp'
		    # Create the output LayerS
		    outShapefile = os.path.join( os.path.split( inShapefile )[0], ola )
		    #print outShapefile
		    outDriver = ogr.GetDriverByName("ESRI Shapefile")

		    # Remove output shapefile if it already exists
		    if os.path.exists(outShapefile):
		        outDriver.DeleteDataSource(outShapefile)

		    # Create the output shapefile
		    outDataSource = outDriver.CreateDataSource(outShapefile)
		    out_lyr_name = os.path.splitext( os.path.split( outShapefile )[1] )[0]
		    outLayer = outDataSource.CreateLayer( out_lyr_name, geom_type=ogr.wkbMultiPolygon )

		    # Add input Layer Fields to the output Layer if it is the one we want
		    inLayerDefn = inLayer.GetLayerDefn()
		    for i in range(0, inLayerDefn.GetFieldCount()):
		        fieldDefn = inLayerDefn.GetFieldDefn(i)
		        fieldName = fieldDefn.GetName()
		        if fieldName not in field_name_target:
		            continue
		        outLayer.CreateField(fieldDefn)

		    # Get the output Layer's Feature Definition
		    outLayerDefn = outLayer.GetLayerDefn()

		    # Add features to the ouput Layer
		    for inFeature in inLayer:
		        # Create output Feature
		        outFeature = ogr.Feature(outLayerDefn)

		        # Add field values from input Layer
		        for i in range(0, outLayerDefn.GetFieldCount()):
		            fieldDefn = outLayerDefn.GetFieldDefn(i)
		            fieldName = fieldDefn.GetName()
		            if fieldName not in field_name_target:
		                continue

		            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(),
		                inFeature.GetField(i))

		        # Set geometry as centroid
		        geom = inFeature.GetGeometryRef()
		        outFeature.SetGeometry(geom.Clone())
		        # Add new feature to output Layer
		        outLayer.CreateFeature(outFeature)

		    # Close DataSources
		    inDataSource.Destroy()
		    outDataSource.Destroy()

		    pixel_size = 1000
		    NoData_value = 255

		    # Filename of input OGR file
		    vector_fn = ola#'sierranevada2.shp'

		    # Filename of the raster Tiff that will be created
		    #raster_fn = 'sierranevada2.tif'

		    # Open the data source and read in the extent
		    source_ds = ogr.Open(outShapefile)#vector_fn)
		    source_layer = source_ds.GetLayer()
		    x_min, x_max, y_min, y_max = source_layer.GetExtent()
		    print str(x_min)+" "+str(x_max)+" "+str(y_min)+" "+str(y_max)
		    tx = abs(x_min)-abs(x_max)
		    ty = abs(y_min)-abs(y_max)
		    ggo = 1
		    #if tx<100 or ty<100: ggo=0
		    #print 'tx:',tx
		    #print 'ty:',ty
		    #print 'ggo:',ggo
	#	    if tx<1000: tx=1000
	#	    if ty<1000: ty=1000
		    if ggo ==1:
			    # Create the destination data source
			    raster_fn = 'tmpf/pa_'+str(bhcat)+'_'+str(idd2)+'_'+str(cluster)+'.tif'
			    x_res = int((x_max - x_min) / pixel_size)
			    y_res = int((y_max - y_min) / pixel_size)
			    print 'xres:',x_res
			    print 'yres:',y_res
	#		    if x_res<1000: x_res=1000
	#		    if y_res<1000: y_res=1000
			    target_ds = gdal.GetDriverByName('MEM').Create('', x_res, y_res, 1, gdal.GDT_Int32)
			    #target_ds2 = gdal.GetDriverByName('GTiff').Create(raster_fn, x_res, y_res, 1, gdal.GDT_Int32)
			    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
			    #target_ds2.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
			    band = target_ds.GetRasterBand(1)
			    #band2 = target_ds.GetRasterBand(1)
			    band.SetNoDataValue(NoData_value)
			    #band2.SetNoDataValue(NoData_value)

			    # Rasterize
			    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])# options = [""]
			    #gdal.RasterizeLayer(target_ds2, [1], source_layer, burn_values=[1])#)

			    array = band.ReadAsArray().astype(np.int32)


			    xoff = int((x_min - xOrigin)/pixelWidth)
			    yoff = int((yOrigin - y_max)/pixelWidth)
			    xcount = int((x_max - x_min)/pixelWidth)#+1
			    ycount = int((y_max - y_min)/pixelWidth)#+1
			    tcount = xcount + ycount
			    print 'tcount:',tcount

			    if xoff < 0:
			        xoff = 0
			    if yoff< 0:
			        yoff =0
		            sumx = xoff + xcount
		            sumy = yoff + ycount
		            go = 1
		            if sumx > gtif_dataset.RasterXSize: go = 0
		            if sumy > gtif_dataset.RasterYSize: go = 0
		            if tcount == 0: go = 0
		            print 'go:',go
		            if go == 1:
				    dataraster = banddataraster.ReadAsArray(xoff, yoff, xcount, ycount)#.astype(np.float)
				    print 'unique:',np.unique(array)
				    mask = np.where(array == 1,dataraster,(0)) 
				    #roads = np.sum(mask)
				    mval = np.mean(dataraster)
				    #tot = np.sum(array)
				    #perc = float(roads/tot)*100
				    print 'pa','cl','cat','mn'
				    var = str(idd2)+" "+str(cluster)+" "+str(bhcat)+" "+str(mval)
				    print var
				    dele = 'rm tmpf/pa_'+str(bhcat)+'*' # '_'+str(idd)+
				    #os.system(dele)
				    #f.write('{0}, {1}, {2}, {3}\n'.format(idd, roads, tot, perc))
				    wb = open(csvname,'a')
				    wb.write(var)
				    wb.write('\n')
				    wb.close()
				    target_ds = dataraster = array = outShapefile = outFeature = outDataSource = vector_fn = outLayer = source_ds = opt = ola = source_layer = inDataSource = None
    # = None
    dele = 'rm results/*tmpfs*' # '_'+str(idd)+
    os.system(dele)
    wb = open(csvname1,'a')
    varx = str(idd)
    wb.write(varx)
    wb.write('\n')
    wb.close() 


def run_batch():
 segmlist = range(1,9)
 lclist = range(1,13)
 for bhcat0 in lclist:
  for clusterx in segmlist:
   lcgc(bhcat=bhcat0,cluster=clusterx)
 print 'DONE'

 
