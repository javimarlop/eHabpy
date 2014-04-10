from __future__ import division
from time import clock
t0 = clock() 
import mkl
import numpy as np
import scipy
from scipy.linalg import cholesky, solve_triangular
from sklearn.externals.joblib import Parallel, delayed
from multiprocessing import cpu_count
import csv
import os
import sys
from osgeo import ogr,gdal
from tqdm import *
mkl.set_num_threads(1)
#except ImportError:
#pass

#adfGeoTransform[0] /* top left x */
#adfGeoTransform[1] /* w-e pixel resolution */
#adfGeoTransform[2] /* rotation, 0 if image is "north up" */
#adfGeoTransform[3] /* top left y */
#adfGeoTransform[4] /* rotation, 0 if image is "north up" */
#adfGeoTransform[5] /* n-s pixel resolution */ 

# Thanks to Sturla Molden for the Mahalanobis functions

####
def _schedule(n, nproc):
	 """ guided scheduler """
	 start = 0
	 size = (n - start) // nproc
	 while size > 100:
		 yield slice(start,start+size)
		 start += size
		 size = (n - start) // nproc
	 yield slice(start,n+1)
	 return

def _mahalanobis_distances(m, L, X):
	 cX = X - m[np.newaxis,:]
	 tmp = solve_triangular(L, cX.T, lower=True).T
	 tmp **= 2
	 #return np.sqrt(tmp.sum(axis=1))
	 return tmp.sum(axis=1)

def mahalanobis_distances(m, S, X, parallel=True):
	 L = cholesky(S, lower=True)
	 n = X.shape[0]
	 if parallel:
		 nproc = cpu_count()
		 res = (Parallel(n_jobs=-1)
				(delayed(_mahalanobis_distances)
				  (m, L, X[s,:])
					for s in _schedule(n,nproc)))
		 return np.hstack(res)
	 else:
		 return _mahalanobis_distances(m, L, X)


# scipy.spatial.distance.mahalanobis for comparison

from scipy.spatial import distance

def _mahalanobis_distances_scipy(m, SI, X):
	 n = X.shape[0]
	 mahal = np.zeros(n)
	 for i in xrange(X.shape[0]):
		 x = X[i,:]
		 mahal[i] = distance.mahalanobis(x,m,SI)
	 return mahal

def mahalanobis_distances_scipy(m, S, X, parallel=True):
	 SI = np.linalg.inv(S)
	 n = X.shape[0]
	 if parallel:
		 nproc = cpu_count()
		 res = (Parallel(n_jobs=-1)
				(delayed(_mahalanobis_distances_scipy)
				 (m, SI, X[s,:])
				   for s in _schedule(n,nproc)))
		 return np.hstack(res)
	 else:
		 return _mahalanobis_distances_scipy(m, SI, X)
###

print "scripts loaded"

#from osgeo import ogr,gdal

indir = 'inVars'

herbf = 'herb.tif'
treef = 'tree.tif'
ndvif = 'ndvi.tif'
ndwif = 'ndwi.tif'
slopef = 'slope.tif'
demf = 'dem.tif'
biof = 'bio.tif'
eprf = 'epr.tif'
pref = 'pre.tif'

demf_globalfile=indir+'/'+demf
src_ds_dem_global = gdal.Open(demf_globalfile)
dem_global = src_ds_dem_global.GetRasterBand(1)
gt_dem_global = src_ds_dem_global.GetGeoTransform()

print 'dem'

biof_globalfile=indir+'/'+biof
src_ds_bio_global = gdal.Open(biof_globalfile)
bio_global = src_ds_bio_global.GetRasterBand(1)
gt_bio_global = src_ds_bio_global.GetGeoTransform()

print 'bio'

pref_globalfile=indir+'/'+pref
src_ds_pre_global = gdal.Open(pref_globalfile)
pre_global = src_ds_pre_global.GetRasterBand(1)
gt_pre_global = src_ds_pre_global.GetGeoTransform()

print 'pre'

eprf_globalfile=indir+'/'+eprf
src_ds_epr_global = gdal.Open(eprf_globalfile)
epr_global = src_ds_epr_global.GetRasterBand(1)
gt_epr_global = src_ds_epr_global.GetGeoTransform()

print 'epr'

herbf_globalfile=indir+'/'+herbf
src_ds_herb_global = gdal.Open(herbf_globalfile)
herb_global = src_ds_herb_global.GetRasterBand(1)
gt_herb_global = src_ds_herb_global.GetGeoTransform()

print 'herb'

ndvif_globalfile=indir+'/'+ndvif
src_ds_ndvi_global = gdal.Open(ndvif_globalfile)
ndvi_global = src_ds_ndvi_global.GetRasterBand(1)
gt_ndvi_global = src_ds_ndvi_global.GetGeoTransform()

print 'ndvi'

ndwif_globalfile=indir+'/'+ndwif
src_ds_ndwi_global = gdal.Open(ndwif_globalfile)
ndwi_global = src_ds_ndwi_global.GetRasterBand(1)
gt_ndwi_global = src_ds_ndwi_global.GetGeoTransform()

print 'ndwi'

slopef_globalfile=indir+'/'+slopef
src_ds_slope_global = gdal.Open(slopef_globalfile)
slope_global = src_ds_slope_global.GetRasterBand(1)
gt_slope_global = src_ds_slope_global.GetGeoTransform()

print 'slope'

treef_globalfile=indir+'/'+treef
src_ds_tree_global = gdal.Open(treef_globalfile)
tree_global = src_ds_tree_global.GetRasterBand(1)
gt_tree_global = src_ds_tree_global.GetGeoTransform()

print 'tree'

print "Global variables imported"

# LOOP ECOREGIONS

ecofile='ecoregs/ecoreg_80601f0.tif'
src_ds_eco = gdal.Open(ecofile)
eco = src_ds_eco.GetRasterBand(1)
eco_mask0 = eco.ReadAsArray(0,0,eco.XSize,eco.YSize).astype(np.int32)
eco_mask = eco_mask0.flatten()
gt_eco = src_ds_eco.GetGeoTransform()

print 'eco mask'

xoff = int((gt_eco[0]-gt_dem_global[0])/1000)
yoff = int((gt_dem_global[3]-gt_eco[3])/1000)
dem_eco_bb0 = dem_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
dem_eco_bb = dem_eco_bb0.flatten()
dem_eco = np.where(eco_mask == 1, (dem_eco_bb),(None))

print 'eco dem'

xoff = int((gt_eco[0]-gt_epr_global[0])/1000)
yoff = int((gt_epr_global[3]-gt_eco[3])/1000)
epr_eco_bb0 = epr_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
epr_eco_bb = epr_eco_bb0.flatten()
epr_eco = np.where(eco_mask == 1, (epr_eco_bb),(0))

print 'eco epr'

xoff = int((gt_eco[0]-gt_slope_global[0])/1000)
yoff = int((gt_slope_global[3]-gt_eco[3])/1000)
slope_eco_bb0 = slope_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
slope_eco_bb = slope_eco_bb0.flatten()
slope_eco = np.where(eco_mask == 1, (slope_eco_bb),(0))

print 'eco slope'

xoff = int((gt_eco[0]-gt_ndvi_global[0])/1000)
yoff = int((gt_ndvi_global[3]-gt_eco[3])/1000)
ndvi_eco_bb0 = ndvi_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
ndvi_eco_bb = ndvi_eco_bb0.flatten()
ndvi_eco = np.where(eco_mask == 1, (ndvi_eco_bb),(0))

print 'eco ndvi'

xoff = int((gt_eco[0]-gt_ndwi_global[0])/1000)
yoff = int((gt_ndwi_global[3]-gt_eco[3])/1000)
ndwi_eco_bb0 = ndwi_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
ndwi_eco_bb = ndwi_eco_bb0.flatten()
ndwi_eco = np.where(eco_mask == 1, (ndwi_eco_bb),(0))

print 'eco ndwi'

xoff = int((gt_eco[0]-gt_pre_global[0])/1000)
yoff = int((gt_pre_global[3]-gt_eco[3])/1000)
pre_eco_bb0 = pre_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
pre_eco_bb = pre_eco_bb0.flatten()
pre_eco = np.where(eco_mask == 1, (pre_eco_bb),(0))

print 'eco pre'

xoff = int((gt_eco[0]-gt_bio_global[0])/1000)
yoff = int((gt_bio_global[3]-gt_eco[3])/1000)
bio_eco_bb0 = bio_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
bio_eco_bb = bio_eco_bb0.flatten()
bio_eco = np.where(eco_mask == 1, (bio_eco_bb),(0))

print 'eco bio'

xoff = int((gt_eco[0]-gt_tree_global[0])/1000)
yoff = int((gt_tree_global[3]-gt_eco[3])/1000)
tree_eco_bb0 = tree_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
tree_eco_bb = tree_eco_bb0.flatten()
tree_eco = np.where(eco_mask == 1, (tree_eco_bb),(0))

print 'eco tree'

xoff = int((gt_eco[0]-gt_herb_global[0])/1000)
yoff = int((gt_herb_global[3]-gt_eco[3])/1000)
herb_eco_bb0 = herb_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
herb_eco_bb = herb_eco_bb0.flatten()
herb_eco = np.where(eco_mask == 1, (herb_eco_bb),(0))

print 'eco herb'

ind_eco = np.column_stack((dem_eco,bio_eco,pre_eco,epr_eco,herb_eco,ndvi_eco,ndwi_eco,slope_eco,tree_eco))

print ind_eco.shape

print 'ecovars stacked'

pa_list0 = np.genfromtxt('parks.csv',dtype='int') # crear este archivo en subpas!
pa_list = np.unique(pa_list0)

n = len(pa_list)
#print pa_list[1]
for px in range(0,n):

#for pa in pa_list:
 pa = pa_list[px]
 print pa
 outfile = 'results/'+str(pa)+'_eco.tif'
 pa4 = 'pas/pa_'+str(pa)+'.tif'
 
 num_bands=src_ds_eco.RasterCount
 driver = gdal.GetDriverByName("GTiff")
 dst_ds = driver.Create( outfile,src_ds_eco.RasterXSize,src_ds_eco.RasterYSize,num_bands,gdal.GDT_Float32)
 dst_ds.SetGeoTransform( src_ds_eco.GetGeoTransform())
 dst_ds.SetProjection( src_ds_eco.GetProjectionRef())
 
 pafile=pa4
 src_ds_pa = gdal.Open(pafile)
 par = src_ds_pa.GetRasterBand(1)
 pa_mask0 = par.ReadAsArray(0,0,par.XSize,par.YSize).astype(np.int32)
 pa_mask = pa_mask0.flatten()
 ind = pa_mask > 0
 gt_pa = src_ds_pa.GetGeoTransform()

 xoff = int((gt_pa[0]-gt_dem_global[0])/1000)
 yoff = int((gt_dem_global[3]-gt_pa[3])/1000)
 dem_pa_bb0 = dem_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 dem_pa_bb = dem_pa_bb0.flatten()
 dem_pa = dem_pa_bb[ind]
 
 
 print 'pa dem'
 print dem_pa.min()
 print dem_pa.max()
 
 xoff = int((gt_pa[0]-gt_tree_global[0])/1000)
 yoff = int((gt_tree_global[3]-gt_pa[3])/1000)
 tree_pa_bb0 = tree_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 tree_pa_bb = tree_pa_bb0.flatten()
 tree_pa = tree_pa_bb[ind]
 
 print 'pa tree'
 print tree_pa.min()
 print tree_pa.max()
  
 xoff = int((gt_pa[0]-gt_epr_global[0])/1000)
 yoff = int((gt_epr_global[3]-gt_pa[3])/1000)
 epr_pa_bb0 = epr_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 epr_pa_bb = epr_pa_bb0.flatten()
 epr_pa = epr_pa_bb[ind]
 
 print 'pa epr'
 print epr_pa.min()
 print epr_pa.max()
 
 xoff = int((gt_pa[0]-gt_pre_global[0])/1000)
 yoff = int((gt_pre_global[3]-gt_pa[3])/1000)
 pre_pa_bb0 = pre_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 pre_pa_bb = pre_pa_bb0.flatten()
 pre_pa = pre_pa_bb[ind]
 
 print 'pa pre'
 print pre_pa.min()
 print pre_pa.max()
 
 xoff = int((gt_pa[0]-gt_bio_global[0])/1000)
 yoff = int((gt_bio_global[3]-gt_pa[3])/1000)
 bio_pa_bb0 = bio_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 bio_pa_bb = bio_pa_bb0.flatten()
 bio_pa = bio_pa_bb[ind]
 
 print 'pa bio'
 print bio_pa.min()
 print bio_pa.max()
 
 xoff = int((gt_pa[0]-gt_slope_global[0])/1000)
 yoff = int((gt_slope_global[3]-gt_pa[3])/1000)
 slope_pa_bb0 = slope_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 slope_pa_bb = slope_pa_bb0.flatten()
 slope_pa = slope_pa_bb[ind]
 
 print 'pa slope'
 print slope_pa.min()
 print slope_pa.max()
 
 xoff = int((gt_pa[0]-gt_ndwi_global[0])/1000)
 yoff = int((gt_ndwi_global[3]-gt_pa[3])/1000)
 ndwi_pa_bb0 = ndwi_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 ndwi_pa_bb = ndwi_pa_bb0.flatten()
 ndwi_pa = ndwi_pa_bb[ind]
 
 print 'pa ndwi'
 print ndwi_pa.min()
 print ndwi_pa.max()
 
 xoff = int((gt_pa[0]-gt_ndvi_global[0])/1000)
 yoff = int((gt_ndvi_global[3]-gt_pa[3])/1000)
 ndvi_pa_bb0 = ndvi_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 ndvi_pa_bb = ndvi_pa_bb0.flatten()
 ndvi_pa = ndvi_pa_bb[ind]
 
 print 'pa ndvi'
 print ndvi_pa.min()
 print ndvi_pa.max()
 
 xoff = int((gt_pa[0]-gt_herb_global[0])/1000)
 yoff = int((gt_herb_global[3]-gt_pa[3])/1000)
 herb_pa_bb0 = herb_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
 herb_pa_bb = herb_pa_bb0.flatten()
 herb_pa = herb_pa_bb[ind]
 
 print 'pa herb'
 print herb_pa.min()
 print herb_pa.max()
 
 ind_pa = np.column_stack((dem_pa,bio_pa,pre_pa,epr_pa,herb_pa,ndvi_pa,ndwi_pa,slope_pa,tree_pa))

 print ind_pa.shape
 
 print "PA masked"

 Ymean = np.mean(ind_pa,axis=0)
 print "Ymean ok"
 Ycov = np.cov(ind_pa,rowvar=False)
 print "Ycov ok"

 mh = mahalanobis_distances(Ymean, Ycov, ind_eco, parallel=False)
# mh = mahalanobis_distances(Ymean, Ycov, ind_eco, parallel=True)
# mh = mahalanobis_distances_scipy(Ymean, Ycov, ind_eco, parallel=True)
# mh = mahalanobis_distances_scipy(Ymean, Ycov, ind_eco, parallel=False)
 print "mh ok"

 from scipy.stats import chisqprob
 pmh = chisqprob(mh,9).reshape((eco.YSize,eco.XSize))
 pmhh = np.where(pmh <= 1e-10,None, pmh)
 print "pmh ok" # quitar valores muy bajos!

 dst_ds.GetRasterBand(1).WriteArray(pmhh)

# calculate single HRI 0.5 value
 pmh2 = pmhh.flatten()
 hr1 = np.where(pmh2 >= 0.5, 1,0)
 hr2 = sum(hr1) #- ind_pa.shape[0] # PROBLEM WITH PA_in pixels
 print hr2
 hr3 = float(hr2/ind_pa.shape[0])
 print hr3
 wb = open('results/hri_results.csv','a')
 var = str(pa)+' '+str(hr3)
 wb.write(var)
 wb.write('\n')
 wb.close()

#before the loop! file_writer.writerow(['wdpa_id', 'ap', 'wn', 'time1', 'time2']) 
 print "results exported"
 
t1 = clock()
print("Time spent: %f min" % ((t1-t0)/60,))
print "DONE"