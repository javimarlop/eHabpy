# Code by Javier Martinez-Lopez
from __future__ import division
#from time import clock
#t0 = clock() 
import mkl
from datetime import datetime
import numpy as np
import scipy.ndimage as nd
import os.path
import scipy
from scipy.linalg import cholesky, solve_triangular
from sklearn.externals.joblib import Parallel, delayed
from multiprocessing import cpu_count,Pool,Lock
import csv
import os
import sys
from osgeo import ogr,gdal
from tqdm import *
mkl.set_num_threads(1)
#except ImportError:
#pass

s = nd.generate_binary_structure(2,2)

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

def process_pa(pa) #ecor
 pa = pa_list[px]
 print pa
 outfile = 'results/'+str(ecor)+'_'+str(pa)+'.tif'
 #print outfile
 pa4 = 'pas/pa_'+str(pa)+'.tif'
 #print pa4
	  
 dropcols = np.arange(9,dtype=int)
	  
 #import os.path
 done = os.path.isfile(outfile)
 avail2 = os.path.isfile(pa4)
 if done == False and avail2 == True:
	  
	 pafile=pa4
	 #print pafile
	 src_ds_pa = gdal.Open(pafile)
	 par = src_ds_pa.GetRasterBand(1)
	 pa_mask0 = par.ReadAsArray(0,0,par.XSize,par.YSize).astype(np.int32)
	 pa_mask = pa_mask0.flatten()
	 #print pa_mask.max()
	 #print pa_mask.min()
	 ind = pa_mask == int(pa)#> 0
	 sum_pa_mask = sum(pa_mask[ind])/int(pa)
	 print sum_pa_mask
	 sum_pa_mask_inv = len(pa_mask[pa_mask == 0])
	 print sum_pa_mask_inv
	 print len(pa_mask)
	 ratiogeom = 10000
	 if sum_pa_mask > 0: ratiogeom = sum_pa_mask_inv/sum_pa_mask
	 print ratiogeom
	 gt_pa = src_ds_pa.GetGeoTransform()
	 
	 xoff = int((gt_pa[0]-gt_pre_global[0])/1000)
	 yoff = int((gt_pre_global[3]-gt_pa[3])/1000)
	 #print xoff
	 #print yoff
	 
	 if xoff>0 and yoff>0:

		 num_bands=src_ds_eco.RasterCount
		 driver = gdal.GetDriverByName("GTiff")
		 dst_options = ['COMPRESS=LZW']
		 dst_ds = driver.Create( outfile,src_ds_eco.RasterXSize,src_ds_eco.RasterYSize,num_bands,gdal.GDT_Float32,dst_options)
		 dst_ds.SetGeoTransform( src_ds_eco.GetGeoTransform())
		 dst_ds.SetProjection( src_ds_eco.GetProjectionRef())
	 
#		 dem_pa_bb0 = dem_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
#		 dem_pa_bb = dem_pa_bb0.flatten()
#		 dem_pa0 = dem_pa_bb[ind]
#		 dem_pa = np.where(dem_pa0 == 65535.0, (float('NaN')),(dem_pa0))
#		 mask2dem = np.isnan(dem_pa)
#		 if mask2dem.all() == True:
#				 dropcols[0] = -1
#		 else:
#				 dem_pa[mask2dem] = np.interp(np.flatnonzero(mask2dem), np.flatnonzero(~mask2dem), dem_pa[~mask2dem])
#				 dem_pa = np.random.random_sample(len(dem_pa),)/1000 + dem_pa
#				 print 'pa dem'
#				 print dem_pa.min()
#				 print dem_pa.max()
			 
		 xoff = int((gt_pa[0]-gt_tree_global[0])/1000)
		 yoff = int((gt_tree_global[3]-gt_pa[3])/1000)
		 tree_pa_bb0 = tree_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 tree_pa_bb = tree_pa_bb0.flatten()
		 tree_pa0 = tree_pa_bb[ind]
		 tree_pa = np.where(tree_pa0 == 255.0, (float('NaN')),(tree_pa0))
		 mask2tree = np.isnan(tree_pa)
		 if mask2tree.all() == True:
				 dropcols[8] = -8
		 else:
				 tree_pa[mask2tree] = np.interp(np.flatnonzero(mask2tree), np.flatnonzero(~mask2tree), tree_pa[~mask2tree])
				 tree_pa = np.random.random_sample(len(tree_pa),)/1000 + tree_pa
				 print 'pa tree'
				 print tree_pa.min()
				 print tree_pa.max()
			
		 xoff = int((gt_pa[0]-gt_epr_global[0])/1000)
		 yoff = int((gt_epr_global[3]-gt_pa[3])/1000)
		 epr_pa_bb0 = epr_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 epr_pa_bb = epr_pa_bb0.flatten()
		 epr_pa0 = epr_pa_bb[ind]
		 epr_pa = np.where(epr_pa0 == 65535.0, (float('NaN')),(epr_pa0))
		 mask2epr = np.isnan(epr_pa)
		 if mask2epr.all() == True:
				 dropcols[3] = -3
		 else:
				 epr_pa[mask2epr] = np.interp(np.flatnonzero(mask2epr), np.flatnonzero(~mask2epr), epr_pa[~mask2epr])
				 epr_pa = np.random.random_sample(len(epr_pa),)/1000 + epr_pa
				 print 'pa epr'
				 print epr_pa.min()
				 print epr_pa.max()
			 
		 xoff = int((gt_pa[0]-gt_pre_global[0])/1000)
		 yoff = int((gt_pre_global[3]-gt_pa[3])/1000)
		 pre_pa_bb0 = pre_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 pre_pa_bb = pre_pa_bb0.flatten()
		 pre_pa0 = pre_pa_bb[ind]
		 pre_pa = np.where(pre_pa0 == 65535.0, (float('NaN')),(pre_pa0))
		 mask2pre = np.isnan(pre_pa)
		 if mask2pre.all() == True:
				 dropcols[2] = -2
		 else:
				 pre_pa[mask2pre] = np.interp(np.flatnonzero(mask2pre), np.flatnonzero(~mask2pre), pre_pa[~mask2pre])
				 pre_pa = np.random.random_sample(len(pre_pa),)/1000 + pre_pa
				 print 'pa pre'
				 print pre_pa.min()
				 print pre_pa.max()
			 
		 xoff = int((gt_pa[0]-gt_bio_global[0])/1000)
		 yoff = int((gt_bio_global[3]-gt_pa[3])/1000)
		 bio_pa_bb0 = bio_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 bio_pa_bb = bio_pa_bb0.flatten()
		 bio_pa0 = bio_pa_bb[ind]
		 bio_pa = np.where(bio_pa0 == 65535.0, (float('NaN')),(bio_pa0))
		 mask2bio = np.isnan(bio_pa)
		 if mask2bio.all() == True:
				 dropcols[1] = -1
		 else:
				 bio_pa[mask2bio] = np.interp(np.flatnonzero(mask2bio), np.flatnonzero(~mask2bio), bio_pa[~mask2bio])
				 bio_pa = np.random.random_sample(len(bio_pa),)/1000 + bio_pa
				 print 'pa bio'
				 print bio_pa.min()
				 print bio_pa.max()
			 
		 xoff = int((gt_pa[0]-gt_slope_global[0])/1000)
		 yoff = int((gt_slope_global[3]-gt_pa[3])/1000)
		 slope_pa_bb0 = slope_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 slope_pa_bb = slope_pa_bb0.flatten()
		 slope_pa0 = slope_pa_bb[ind]
		 slope_pa = np.where(slope_pa0 == 65535.0, (float('NaN')),(slope_pa0))
		 mask2slope = np.isnan(slope_pa)
		 if mask2slope.all() == True:
				 dropcols[7] = -7
		 else:
				 slope_pa[mask2slope] = np.interp(np.flatnonzero(mask2slope), np.flatnonzero(~mask2slope), slope_pa[~mask2slope])
				 slope_pa = np.random.random_sample(len(slope_pa),)/1000 + slope_pa
				 print 'pa slope'
				 print slope_pa.min()
				 print slope_pa.max()
			 
		 xoff = int((gt_pa[0]-gt_ndwi_global[0])/1000)
		 yoff = int((gt_ndwi_global[3]-gt_pa[3])/1000)
		 ndwi_pa_bb0 = ndwi_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 ndwi_pa_bb = ndwi_pa_bb0.flatten()
		 ndwi_pa0 = ndwi_pa_bb[ind]
		 ndwi_pa = np.where(ndwi_pa0 == 255.0, (float('NaN')),(ndwi_pa0))
		 mask2ndwi = np.isnan(ndwi_pa)
		 if mask2ndwi.all() == True:
				 dropcols[6] = -6
		 else:
				 ndwi_pa[mask2ndwi] = np.interp(np.flatnonzero(mask2ndwi), np.flatnonzero(~mask2ndwi), ndwi_pa[~mask2ndwi])
				 ndwi_pa = np.random.random_sample(len(ndwi_pa),)/1000 + ndwi_pa
				 print 'pa ndwi'
				 print ndwi_pa.min()
				 print ndwi_pa.max()
			 
		 xoff = int((gt_pa[0]-gt_ndvimax_global[0])/1000)
		 yoff = int((gt_ndvimax_global[3]-gt_pa[3])/1000)
		 ndvimax_pa_bb0 = ndvimax_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 ndvimax_pa_bb = ndvimax_pa_bb0.flatten()
		 ndvimax_pa0 = ndvimax_pa_bb[ind]
		 ndvimax_pa = np.where(ndvimax_pa0 == 65535.0, (float('NaN')),(ndvimax_pa0))
		 mask2ndvimax = np.isnan(ndvimax_pa)
		 if mask2ndvimax.all() == True:
				 dropcols[5] = -5
		 else:
				 ndvimax_pa[mask2ndvimax] = np.interp(np.flatnonzero(mask2ndvimax), np.flatnonzero(~mask2ndvimax), ndvimax_pa[~mask2ndvimax])
				 ndvimax_pa = np.random.random_sample(len(ndvimax_pa),)/1000 + ndvimax_pa
				 print 'pa ndvimax'
				 print ndvimax_pa.min()
				 print ndvimax_pa.max()

		 xoff = int((gt_pa[0]-gt_ndvimin_global[0])/1000)
		 yoff = int((gt_ndvimin_global[3]-gt_pa[3])/1000)
		 ndvimin_pa_bb0 = ndvimin_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 ndvimin_pa_bb = ndvimin_pa_bb0.flatten()
		 ndvimin_pa0 = ndvimin_pa_bb[ind]
		 ndvimin_pa = np.where(ndvimin_pa0 == 65535.0, (float('NaN')),(ndvimin_pa0))
		 mask2ndvimin = np.isnan(ndvimin_pa)
		 if mask2ndvimin.all() == True:
				 dropcols[5] = -5
		 else:
				 ndvimin_pa[mask2ndvimin] = np.interp(np.flatnonzero(mask2ndvimin), np.flatnonzero(~mask2ndvimin), ndvimin_pa[~mask2ndvimin])
				 ndvimin_pa = np.random.random_sample(len(ndvimin_pa),)/1000 + ndvimin_pa
				 print 'pa ndvimin'
				 print ndvimin_pa.min()
				 print ndvimin_pa.max()
			 
		 xoff = int((gt_pa[0]-gt_herb_global[0])/1000)
		 yoff = int((gt_herb_global[3]-gt_pa[3])/1000)
		 herb_pa_bb0 = herb_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
		 herb_pa_bb = herb_pa_bb0.flatten()
		 herb_pa0 = herb_pa_bb[ind]
		 herb_pa = np.where(herb_pa0 == 255.0, (float('NaN')),(herb_pa0))
		 mask2herb = np.isnan(herb_pa)
		 if mask2herb.all() == True:
				 dropcols[4] = -4
		 else:
				 herb_pa[mask2herb] = np.interp(np.flatnonzero(mask2herb), np.flatnonzero(~mask2herb), herb_pa[~mask2herb])
				 herb_pa = np.random.random_sample(len(herb_pa),)/1000 + herb_pa
				 print 'pa herb'
				 print herb_pa.min()
				 print herb_pa.max()
			 
		 cols = dropcols[dropcols>=0]
		 #tot = tree_pa.max() + herb_pa.min()
		 #print tot
		 #if  tot == 100: tree_pa = np.random.random_sample(len(tree_pa),) + tree_pa
			  
		 ind_pa0 = np.column_stack((bio_pa,pre_pa,epr_pa,herb_pa,ndvimax_pa,ndvimin_pa,ndwi_pa,slope_pa,tree_pa))
			  
		 ind_pa = ind_pa0[:,cols]
		 ind_eco = ind_eco0[:,cols]

		 print ind_pa.shape
			 
		 print "PA masked"

		 Ymean = np.mean(ind_pa,axis=0)
			  #print Ymean
		 print "Ymean ok"
		 Ycov = np.cov(ind_pa,rowvar=False)
			  #print 'Ycov'
			  #print Ycov
		 print "Ycov ok"

		 mh2 = mahalanobis_distances(Ymean, Ycov, ind_eco, parallel=False)
			  #mh = mahalanobis_distances(Ymean, Ycov, ind_eco, parallel=True)
		 #mh2 = mahalanobis_distances_scipy(Ymean, Ycov, ind_eco, parallel=True)
		 mh = mh2*mh2
			# mh = mahalanobis_distances_scipy(Ymean, Ycov, ind_eco, parallel=False)
		 print "mh ok"

		 from scipy.stats import chisqprob
		 pmh = chisqprob(mh,9).reshape((eco.YSize,eco.XSize))
		 pmhh = np.where(pmh <= 1e-10,None, pmh)
		 print "pmh ok" # quitar valores muy bajos!

		 dst_ds.GetRasterBand(1).WriteArray(pmhh)
		 dst_ds = None
		 
		 #pmh2 = pmhh.flatten()
		 hr11 = np.where(pmhh >= 0.5, 1,0)
		 hr1 = hr11.flatten()
		 hr1sum = sum(hr1)
		 print hr1sum
		 hr1insumaver = hr1insum = 0
		 hr1sumaver = hr1sum
		 hr1averpa = hr3aver = num_featuresaver = hr1medianpa = None
		 labeled_array, num_features = nd.label(hr11, structure=s)
		 
	 
		 # read back the tif file and clip park values
		 #hrif = outfile
		 src_ds_sim = gdal.Open(outfile)
		 sim = src_ds_sim.GetRasterBand(1)
		 gt_sim = src_ds_sim.GetGeoTransform()
		 
		 #print gt_pa[0]
		 #print gt_sim[0]
		 
		 #print gt_sim[3]
		 #print gt_pa[3]
		 
		 #print gt_pa[1]
		 #print gt_pa[5]
		 
		 #print gt_sim[1]
		 #print gt_sim[5]
		 
		 xoff = int((gt_pa[0]-gt_sim[0])/1000)
		 #print xoff
		 yoff = int((gt_sim[3]-gt_pa[3])/1000)
		 #print yoff
		 
		 xextentpa = xoff + par.XSize
		 yextentpa = yoff + par.YSize
		 xless = sim.XSize - xextentpa
		 yless = sim.YSize - yextentpa

		 xsize = par.XSize
		 ysize = par.YSize
		 if xoff>0 and yoff>0 and ratiogeom < 100: # also check if results are not empty?
			 if xless < 0: xsize = xsize + xless
			 if yless < 0: ysize = ysize + yless
			 hri_pa_bb0 = sim.ReadAsArray(xoff,yoff,xsize,ysize).astype(np.float32)
			 hri_pa_bb = hri_pa_bb0.flatten()
			 indd = hri_pa_bb > 0
			 hri_pa0 = hri_pa_bb[indd]
			 #print hri_pa0.max()
			 #print hri_pa0.min()
			 hr1averpa = np.mean(hri_pa0[~np.isnan(hri_pa0)])
			 hr1medianpa = np.median(hri_pa0[~np.isnan(hri_pa0)])
			 #hr1p25pa = np.percentile(hri_pa0[~np.isnan(hri_pa0)],25)
			 print 'mean similarity in the park is '+str(hr1averpa)
			 hr1insum = sum(np.where(hri_pa0 >= 0.5, 1,0)) # use hr1averpa as threshold instead!		 
			 hr1inaver = np.where(hri_pa0 >= hr1averpa, 1,0)
			 hr1insumaver = sum(hr1inaver)
			 print hr1insum
			 hr1averr = np.where(pmhh >= hr1averpa, 1,0)
			 hr1aver = hr1averr.flatten()
			 labeled_arrayaver, num_featuresaver = nd.label(hr1averr, structure=s)
			 hr1sumaver = sum(hr1aver)
			 hr2aver = hr1sumaver - hr1insumaver
			 hr3aver = float(hr2aver/ind_pa.shape[0])
                                   aggregation = float(hr2aver/num_featuresaver)
		 
			# calculate single HRI 0.5 value
		 
		 #print len(pmh2)
		 #pmh2in = pmh2[ind]
		 #print len(pmh2in)

		 #hr1in = hr1[ind]#np.where(pmh2in >= 0.5, 1,0)
		 #hr1insum = sum(hr1in)
		 #hr2 = sum(hr1) #- hr1insum#- ind_pa.shape[0] # PROBLEM WITH PA_in pixels
		 hr2 = hr1sum - hr1insum
		 #print sum(hr1)
		 #print hr1insum
		 print hr2
		 hr3 = float(hr2/ind_pa.shape[0])
		 
		 print hr3
		 wb = open(csvname,'a')
		 var = str(ecor)+' '+str(pa)+' '+str(hr3)+' '+str(hr1averpa)+' '+str(hr3aver)+' '+str(hr1medianpa)+' '+str(num_features)+' '+str(num_featuresaver)+' '+str(aggregation) # exclude PA! #+' '+str(hr1p25pa)#
		 wb.write(var)
		 wb.write('\n')
		 wb.close()

			#before the loop! file_writer.writerow(['wdpa_id', 'ap', 'wn', 'time1', 'time2']) 
		 print "results exported"
	  wb = open('results/ecoregs_done.csv','a')
	  var = str(ecor)
	  wb.write(var)
	  wb.write('\n')
	  wb.close() 

###
print "scripts loaded"

#from osgeo import ogr,gdal

indir = 'inVars'

herbf = 'herb.tif'
treef = 'tree.tif'
ndvimaxf = 'ndvimax.tif'
ndviminf = 'ndvimin.tif'
ndwif = 'ndwi.tif'
slopef = 'slope.tif'
#demf = 'dem.tif'
biof = 'bio.tif'
eprf = 'epr.tif'
pref = 'pre.tif'

#demf_globalfile=indir+'/'+demf
#src_ds_dem_global = gdal.Open(demf_globalfile)
#dem_global = src_ds_dem_global.GetRasterBand(1)
#gt_dem_global = src_ds_dem_global.GetGeoTransform()

#print 'dem'

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

ndvimaxf_globalfile=indir+'/'+ndvimaxf
src_ds_ndvimax_global = gdal.Open(ndvimaxf_globalfile)
ndvimax_global = src_ds_ndvimax_global.GetRasterBand(1)
gt_ndvimax_global = src_ds_ndvimax_global.GetGeoTransform()

print 'ndvimax'

ndviminf_globalfile=indir+'/'+ndviminf
src_ds_ndvimin_global = gdal.Open(ndviminf_globalfile)
ndvimin_global = src_ds_ndvimin_global.GetRasterBand(1)
gt_ndvimin_global = src_ds_ndvimin_global.GetGeoTransform()

print 'ndvimin'

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

#eco_list_done = 'None'

csvname1 = 'results/ecoregs_done.csv'
if os.path.isfile(csvname1) == False:
 wb = open(csvname1,'a')
 wb.write('None')
 wb.write('\n')
 wb.close()

csvname = 'results/hri_results.csv'
if os.path.isfile(csvname) == False:
 wb = open(csvname,'a')
 wb.write('ecoregion wdpaid hri50 averpasim hriaver medianpasim nrpatches nrpatchesaver aggregation')
 wb.write('\n')
 wb.close()

# LOOP ECOREGIONS
eco_list0 = np.genfromtxt('pas/ecoregs.csv',dtype='int') # crear este archivo en subpas!
eco_list = np.unique(eco_list0)
print eco_list
eco_list_done = np.genfromtxt('results/ecoregs_done.csv',dtype='int')
print eco_list_done

m = len(eco_list)
#print pa_list[1]
for pm in tqdm(range(0,m)): # 3,m # 0 without the negative ecoregs!
 ecor = eco_list[pm]
 print ecor
 if ecor not in eco_list_done:
  ecofile='ecoregs/eco_'+str(ecor)+'.tif'
  avail = os.path.isfile(ecofile)
  if avail == True:
	  ecoparksf = 'pas/'+str(ecor)+'.csv'
	  src_ds_eco = gdal.Open(ecofile)
	  eco = src_ds_eco.GetRasterBand(1)
	  eco_mask0 = eco.ReadAsArray(0,0,eco.XSize,eco.YSize).astype(np.int32)
	  eco_mask = eco_mask0.flatten()
	  gt_eco = src_ds_eco.GetGeoTransform()
	  print 'eco mask'
	  
	  #if gt_eco[0]>gt_dem_global[0] and gt_eco[3]>gt_dem_global[3]:

#	  xoff = int((gt_eco[0]-gt_dem_global[0])/1000)
#	  yoff = int((gt_dem_global[3]-gt_eco[3])/1000)
#	  dem_eco_bb0 = dem_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
#	  dem_eco_bb = dem_eco_bb0.flatten()
#	  dem_eco0 = np.where(eco_mask == 1, (dem_eco_bb),(0))
#	  dem_eco = np.where(dem_eco0 == 65535.0, (float('NaN')),(dem_eco0))
#	  maskdem = np.isnan(dem_eco)
#	  dem_eco[maskdem] = np.interp(np.flatnonzero(maskdem), np.flatnonzero(~maskdem), dem_eco[~maskdem])

#	  print 'eco dem'

	  xoff = int((gt_eco[0]-gt_epr_global[0])/1000)
	  yoff = int((gt_epr_global[3]-gt_eco[3])/1000)
	  epr_eco_bb0 = epr_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  epr_eco_bb = epr_eco_bb0.flatten()
	  epr_eco0 = np.where(eco_mask == 1, (epr_eco_bb),(0))
	  epr_eco = np.where(epr_eco0 == 65535.0, (float('NaN')),(epr_eco0))
	  maskepr = np.isnan(epr_eco)
	  epr_eco[maskepr] = np.interp(np.flatnonzero(maskepr), np.flatnonzero(~maskepr), epr_eco[~maskepr])

	  print 'eco epr'

	  xoff = int((gt_eco[0]-gt_slope_global[0])/1000)
	  yoff = int((gt_slope_global[3]-gt_eco[3])/1000)
	  slope_eco_bb0 = slope_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  slope_eco_bb = slope_eco_bb0.flatten()
	  slope_eco0 = np.where(eco_mask == 1, (slope_eco_bb),(0))
	  slope_eco = np.where(slope_eco0 == 65535.0, (float('NaN')),(slope_eco0))
	  maskslope = np.isnan(slope_eco)
	  slope_eco[maskslope] = np.interp(np.flatnonzero(maskslope), np.flatnonzero(~maskslope), slope_eco[~maskslope])

	  print 'eco slope'

	  xoff = int((gt_eco[0]-gt_ndvimax_global[0])/1000)
	  yoff = int((gt_ndvimax_global[3]-gt_eco[3])/1000)
	  ndvimax_eco_bb0 = ndvimax_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  ndvimax_eco_bb = ndvimax_eco_bb0.flatten()
	  ndvimax_eco0 = np.where(eco_mask == 1, (ndvimax_eco_bb),(0))
	  ndvimax_eco = np.where(ndvimax_eco0 == 65535.0, (float('NaN')),(ndvimax_eco0))
	  maskndvimax = np.isnan(ndvimax_eco)
	  ndvimax_eco[maskndvimax] = np.interp(np.flatnonzero(maskndvimax), np.flatnonzero(~maskndvimax), ndvimax_eco[~maskndvimax])

	  print 'eco ndvimax'

	  xoff = int((gt_eco[0]-gt_ndvimin_global[0])/1000)
	  yoff = int((gt_ndvimin_global[3]-gt_eco[3])/1000)
	  ndvimin_eco_bb0 = ndvimin_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  ndvimin_eco_bb = ndvimin_eco_bb0.flatten()
	  ndvimin_eco0 = np.where(eco_mask == 1, (ndvimin_eco_bb),(0))
	  ndvimin_eco = np.where(ndvimin_eco0 == 65535.0, (float('NaN')),(ndvimin_eco0))
	  maskndvimin = np.isnan(ndvimin_eco)
	  ndvimin_eco[maskndvimin] = np.interp(np.flatnonzero(maskndvimin), np.flatnonzero(~maskndvimin), ndvimin_eco[~maskndvimin])

	  print 'eco ndvimin'

	  xoff = int((gt_eco[0]-gt_ndwi_global[0])/1000)
	  yoff = int((gt_ndwi_global[3]-gt_eco[3])/1000)
	  ndwi_eco_bb0 = ndwi_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  ndwi_eco_bb = ndwi_eco_bb0.flatten()
	  ndwi_eco0 = np.where(eco_mask == 1, (ndwi_eco_bb),(0))
	  ndwi_eco = np.where(ndwi_eco0 == 255.0, (float('NaN')),(ndwi_eco0))
	  maskndwi = np.isnan(ndwi_eco)
	  ndwi_eco[maskndwi] = np.interp(np.flatnonzero(maskndwi), np.flatnonzero(~maskndwi), ndwi_eco[~maskndwi])

	  print 'eco ndwi'

	  xoff = int((gt_eco[0]-gt_pre_global[0])/1000)
	  yoff = int((gt_pre_global[3]-gt_eco[3])/1000)
	  pre_eco_bb0 = pre_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  pre_eco_bb = pre_eco_bb0.flatten()
	  pre_eco0 = np.where(eco_mask == 1, (pre_eco_bb),(0))
	  pre_eco = np.where(pre_eco0 == 65535.0, (float('NaN')),(pre_eco0))
	  maskpre = np.isnan(pre_eco)
	  pre_eco[maskpre] = np.interp(np.flatnonzero(maskpre), np.flatnonzero(~maskpre), pre_eco[~maskpre])

	  print 'eco pre'

	  xoff = int((gt_eco[0]-gt_bio_global[0])/1000)
	  yoff = int((gt_bio_global[3]-gt_eco[3])/1000)
	  bio_eco_bb0 = bio_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  bio_eco_bb = bio_eco_bb0.flatten()
	  bio_eco0 = np.where(eco_mask == 1, (bio_eco_bb),(0))
	  bio_eco = np.where(bio_eco0 == 65535.0, (float('NaN')),(bio_eco0))
	  maskbio = np.isnan(bio_eco)
	  bio_eco[maskbio] = np.interp(np.flatnonzero(maskbio), np.flatnonzero(~maskbio), bio_eco[~maskbio])

	  print 'eco bio'

	  xoff = int((gt_eco[0]-gt_tree_global[0])/1000)
	  yoff = int((gt_tree_global[3]-gt_eco[3])/1000)
	  tree_eco_bb0 = tree_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  tree_eco_bb = tree_eco_bb0.flatten()
	  tree_eco0 = np.where(eco_mask == 1, (tree_eco_bb),(0))
	  tree_eco = np.where(tree_eco0 == 255.0, (float('NaN')),(tree_eco0))
	  masktree = np.isnan(tree_eco)
	  tree_eco[masktree] = np.interp(np.flatnonzero(masktree), np.flatnonzero(~masktree), tree_eco[~masktree])

	  print 'eco tree'

	  xoff = int((gt_eco[0]-gt_herb_global[0])/1000)
	  yoff = int((gt_herb_global[3]-gt_eco[3])/1000)
	  herb_eco_bb0 = herb_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
	  herb_eco_bb = herb_eco_bb0.flatten()
	  herb_eco0 = np.where(eco_mask == 1, (herb_eco_bb),(0))
	  herb_eco = np.where(herb_eco0 == 255.0, (float('NaN')),(herb_eco0))
	  maskherb = np.isnan(herb_eco)
	  herb_eco[maskherb] = np.interp(np.flatnonzero(maskherb), np.flatnonzero(~maskherb), herb_eco[~maskherb])

	  print 'eco herb'

	  ind_eco0 = np.column_stack((bio_eco,pre_eco,epr_eco,herb_eco,ndvimax_eco,ndvimin_eco,ndwi_eco,slope_eco,tree_eco))

	  #print ind_eco.shape

	  print 'ecovars stacked'

	  pa_list0 = np.genfromtxt(ecoparksf,dtype='int') # crear este archivo en subpas!
	  pa_list = np.unique(pa_list0)

	  pool = Pool()
	  pool.map(process_pa,pa_list)
	  pool.close()
	  pool.join()



#t1 = clock()
print str(datetime.now())
#print("Time spent: %f min" % ((t1-t0)/60,))
print "END"
