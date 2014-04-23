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
eco_list0 = np.genfromtxt('pas/ecoregs.csv',dtype='int') # crear este archivo en subpas!
eco_list = np.unique(eco_list0)

m = len(eco_list)
#print pa_list[1]
for pm in tqdm(range(3,m)): # 3,m # 0 without the negative ecoregs!
 eco = eco_list[pm]
 print eco
 ecofile='ecoregs/eco_'+str(eco)+'.tif'
 ecoparksf = 'pas/'+str(eco)+'.csv'
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
 dem_eco = np.where(eco_mask == 1, (dem_eco_bb),(0))
 maskdem = np.isnan(dem_eco)
 dem_eco[maskdem] = np.interp(np.flatnonzero(maskdem), np.flatnonzero(~maskdem), dem_eco[~maskdem])

 print 'eco dem'

 xoff = int((gt_eco[0]-gt_epr_global[0])/1000)
 yoff = int((gt_epr_global[3]-gt_eco[3])/1000)
 epr_eco_bb0 = epr_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 epr_eco_bb = epr_eco_bb0.flatten()
 epr_eco = np.where(eco_mask == 1, (epr_eco_bb),(0))
 maskepr = np.isnan(epr_eco)
 epr_eco[maskepr] = np.interp(np.flatnonzero(maskepr), np.flatnonzero(~maskepr), epr_eco[~maskepr])

 print 'eco epr'

 xoff = int((gt_eco[0]-gt_slope_global[0])/1000)
 yoff = int((gt_slope_global[3]-gt_eco[3])/1000)
 slope_eco_bb0 = slope_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 slope_eco_bb = slope_eco_bb0.flatten()
 slope_eco = np.where(eco_mask == 1, (slope_eco_bb),(0))
 maskslope = np.isnan(slope_eco)
 slope_eco[maskslope] = np.interp(np.flatnonzero(maskslope), np.flatnonzero(~maskslope), slope_eco[~maskslope])

 print 'eco slope'

 xoff = int((gt_eco[0]-gt_ndvi_global[0])/1000)
 yoff = int((gt_ndvi_global[3]-gt_eco[3])/1000)
 ndvi_eco_bb0 = ndvi_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 ndvi_eco_bb = ndvi_eco_bb0.flatten()
 ndvi_eco = np.where(eco_mask == 1, (ndvi_eco_bb),(0))
 maskndvi = np.isnan(ndvi_eco)
 ndvi_eco[maskndvi] = np.interp(np.flatnonzero(maskndvi), np.flatnonzero(~maskndvi), ndvi_eco[~maskndvi])

 print 'eco ndvi'

 xoff = int((gt_eco[0]-gt_ndwi_global[0])/1000)
 yoff = int((gt_ndwi_global[3]-gt_eco[3])/1000)
 ndwi_eco_bb0 = ndwi_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 ndwi_eco_bb = ndwi_eco_bb0.flatten()
 ndwi_eco = np.where(eco_mask == 1, (ndwi_eco_bb),(0))
 maskndwi = np.isnan(ndwi_eco)
 ndwi_eco[maskndwi] = np.interp(np.flatnonzero(maskndwi), np.flatnonzero(~maskndwi), ndwi_eco[~maskndwi])

 print 'eco ndwi'

 xoff = int((gt_eco[0]-gt_pre_global[0])/1000)
 yoff = int((gt_pre_global[3]-gt_eco[3])/1000)
 pre_eco_bb0 = pre_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 pre_eco_bb = pre_eco_bb0.flatten()
 pre_eco = np.where(eco_mask == 1, (pre_eco_bb),(0))
 maskpre = np.isnan(pre_eco)
 pre_eco[maskpre] = np.interp(np.flatnonzero(maskpre), np.flatnonzero(~maskpre), pre_eco[~maskpre])

 print 'eco pre'

 xoff = int((gt_eco[0]-gt_bio_global[0])/1000)
 yoff = int((gt_bio_global[3]-gt_eco[3])/1000)
 bio_eco_bb0 = bio_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 bio_eco_bb = bio_eco_bb0.flatten()
 bio_eco = np.where(eco_mask == 1, (bio_eco_bb),(0))
 maskbio = np.isnan(bio_eco)
 bio_eco[maskbio] = np.interp(np.flatnonzero(maskbio), np.flatnonzero(~maskbio), bio_eco[~maskbio])

 print 'eco bio'

 xoff = int((gt_eco[0]-gt_tree_global[0])/1000)
 yoff = int((gt_tree_global[3]-gt_eco[3])/1000)
 tree_eco_bb0 = tree_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 tree_eco_bb = tree_eco_bb0.flatten()
 tree_eco = np.where(eco_mask == 1, (tree_eco_bb),(0))
 masktree = np.isnan(tree_eco)
 tree_eco[masktree] = np.interp(np.flatnonzero(masktree), np.flatnonzero(~masktree), tree_eco[~masktree])

 print 'eco tree'

 xoff = int((gt_eco[0]-gt_herb_global[0])/1000)
 yoff = int((gt_herb_global[3]-gt_eco[3])/1000)
 herb_eco_bb0 = herb_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
 herb_eco_bb = herb_eco_bb0.flatten()
 herb_eco = np.where(eco_mask == 1, (herb_eco_bb),(0))
 maskherb = np.isnan(herb_eco)
 herb_eco[maskherb] = np.interp(np.flatnonzero(maskherb), np.flatnonzero(~maskherb), herb_eco[~maskherb])

 print 'eco herb'

 ind_eco = np.column_stack((dem_eco,bio_eco,pre_eco,epr_eco,herb_eco,ndvi_eco,ndwi_eco,slope_eco,tree_eco))

 print ind_eco.shape

 print 'ecovars stacked'

 pa_list0 = np.genfromtxt(ecoparksf,dtype='int') # crear este archivo en subpas!
 pa_list = np.unique(pa_list0)

 n = len(pa_list)
#print pa_list[1]
 for px in range(0,n): # 0,n

#for pa in pa_list:
  pa = pa_list[px]
  print pa
  outfile = 'results/'+str(pa)+'_eco.tif'
  pa4 = 'pas/pa_'+str(pa)+'.tif'
  
  import os.path
  done = os.path.isfile(outfile)
  
  if done == False:
  
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
	  mask2dem = np.isnan(dem_pa)
	  dem_pa[mask2dem] = np.interp(np.flatnonzero(mask2dem), np.flatnonzero(~mask2dem), dem_pa[~mask2dem])
	 
	  print 'pa dem'
	  print dem_pa.min()
	  print dem_pa.max()
	 
	  xoff = int((gt_pa[0]-gt_tree_global[0])/1000)
	  yoff = int((gt_tree_global[3]-gt_pa[3])/1000)
	  tree_pa_bb0 = tree_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  tree_pa_bb = tree_pa_bb0.flatten()
	  tree_pa = tree_pa_bb[ind]
	  mask2tree = np.isnan(tree_pa)
	  tree_pa[mask2tree] = np.interp(np.flatnonzero(mask2tree), np.flatnonzero(~mask2tree), tree_pa[~mask2tree])
	 
	  print 'pa tree'
	  print tree_pa.min()
	  print tree_pa.max()
	  
	  xoff = int((gt_pa[0]-gt_epr_global[0])/1000)
	  yoff = int((gt_epr_global[3]-gt_pa[3])/1000)
	  epr_pa_bb0 = epr_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  epr_pa_bb = epr_pa_bb0.flatten()
	  epr_pa = epr_pa_bb[ind]
	  mask2epr = np.isnan(epr_pa)
	  epr_pa[mask2epr] = np.interp(np.flatnonzero(mask2epr), np.flatnonzero(~mask2epr), epr_pa[~mask2epr])
	 
	  print 'pa epr'
	  print epr_pa.min()
	  print epr_pa.max()
	 
	  xoff = int((gt_pa[0]-gt_pre_global[0])/1000)
	  yoff = int((gt_pre_global[3]-gt_pa[3])/1000)
	  pre_pa_bb0 = pre_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  pre_pa_bb = pre_pa_bb0.flatten()
	  pre_pa = pre_pa_bb[ind]
	  mask2pre = np.isnan(pre_pa)
	  pre_pa[mask2pre] = np.interp(np.flatnonzero(mask2pre), np.flatnonzero(~mask2pre), pre_pa[~mask2pre])
	 
	  print 'pa pre'
	  print pre_pa.min()
	  print pre_pa.max()
	 
	  xoff = int((gt_pa[0]-gt_bio_global[0])/1000)
	  yoff = int((gt_bio_global[3]-gt_pa[3])/1000)
	  bio_pa_bb0 = bio_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  bio_pa_bb = bio_pa_bb0.flatten()
	  bio_pa = bio_pa_bb[ind]
	  mask2bio = np.isnan(bio_pa)
	  bio_pa[mask2bio] = np.interp(np.flatnonzero(mask2bio), np.flatnonzero(~mask2bio), bio_pa[~mask2bio])

	  print 'pa bio'
	  print bio_pa.min()
	  print bio_pa.max()
	 
	  xoff = int((gt_pa[0]-gt_slope_global[0])/1000)
	  yoff = int((gt_slope_global[3]-gt_pa[3])/1000)
	  slope_pa_bb0 = slope_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  slope_pa_bb = slope_pa_bb0.flatten()
	  slope_pa = slope_pa_bb[ind]
	  mask2slope = np.isnan(slope_pa)
	  slope_pa[mask2slope] = np.interp(np.flatnonzero(mask2slope), np.flatnonzero(~mask2slope), slope_pa[~mask2slope])

	  print 'pa slope'
	  print slope_pa.min()
	  print slope_pa.max()
	 
	  xoff = int((gt_pa[0]-gt_ndwi_global[0])/1000)
	  yoff = int((gt_ndwi_global[3]-gt_pa[3])/1000)
	  ndwi_pa_bb0 = ndwi_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  ndwi_pa_bb = ndwi_pa_bb0.flatten()
	  ndwi_pa = ndwi_pa_bb[ind]
	  mask2ndwi = np.isnan(ndwi_pa)
	  ndwi_pa[mask2ndwi] = np.interp(np.flatnonzero(mask2ndwi), np.flatnonzero(~mask2ndwi), ndwi_pa[~mask2ndwi])

	  print 'pa ndwi'
	  print ndwi_pa.min()
	  print ndwi_pa.max()
	 
	  xoff = int((gt_pa[0]-gt_ndvi_global[0])/1000)
	  yoff = int((gt_ndvi_global[3]-gt_pa[3])/1000)
	  ndvi_pa_bb0 = ndvi_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  ndvi_pa_bb = ndvi_pa_bb0.flatten()
	  ndvi_pa = ndvi_pa_bb[ind]
	  mask2ndvi = np.isnan(ndvi_pa)
	  ndvi_pa[mask2ndvi] = np.interp(np.flatnonzero(mask2ndvi), np.flatnonzero(~mask2ndvi), ndvi_pa[~mask2ndvi])

	  print 'pa ndvi'
	  print ndvi_pa.min()
	  print ndvi_pa.max()
	 
	  xoff = int((gt_pa[0]-gt_herb_global[0])/1000)
	  yoff = int((gt_herb_global[3]-gt_pa[3])/1000)
	  herb_pa_bb0 = herb_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
	  herb_pa_bb = herb_pa_bb0.flatten()
	  herb_pa = herb_pa_bb[ind]
	  mask2herb = np.isnan(herb_pa)
	  herb_pa[mask2herb] = np.interp(np.flatnonzero(mask2herb), np.flatnonzero(~mask2herb), herb_pa[~mask2herb])

	  print 'pa herb'
	  print herb_pa.min()
	  print herb_pa.max()
	 
	  #tot = tree_pa.max() + herb_pa.min()
	  #print tot
	  #if  tot == 100: tree_pa = np.random.random_sample(len(tree_pa),) + tree_pa
	  
	  tree_pa = np.random.random_sample(len(tree_pa),) + tree_pa
	  herb_pa = np.random.random_sample(len(herb_pa),) + herb_pa
	  
	  ind_pa = np.column_stack((dem_pa,bio_pa,pre_pa,epr_pa,herb_pa,ndvi_pa,ndwi_pa,slope_pa,tree_pa))

	  print ind_pa.shape
	 
	  print "PA masked"

	  Ymean = np.mean(ind_pa,axis=0)
	  #print Ymean
	  print "Ymean ok"
	  Ycov = np.cov(ind_pa,rowvar=False)
	  #print 'Ycov'
	  #print Ycov
	  print "Ycov ok"

	#mh = mahalanobis_distances(Ymean, Ycov, ind_eco, parallel=False)
	  #mh = mahalanobis_distances(Ymean, Ycov, ind_eco, parallel=True)
	  mh2 = mahalanobis_distances_scipy(Ymean, Ycov, ind_eco, parallel=True)
	  mh = mh2*mh2
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
