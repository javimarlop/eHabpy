#### Script by Javier Martinez-Lopez (UTF-8)

from __future__ import division
from datetime import datetime
import numpy as np
import scipy.ndimage as nd
import os.path
import scipy
from scipy.linalg import cholesky, solve_triangular
from scipy.spatial import distance
#from scipy.stats import chisqprob
from scipy.stats.distributions import chi2
from sklearn.externals.joblib import Parallel, delayed
from multiprocessing import cpu_count
import csv
import os
import sys
from osgeo import ogr,gdal

# ## Create  directory safely - don't disturb it if it exists already - LB
def safelyMakeDir(d):
	try:
		os.makedirs(d)
		return True
	except OSError:
		if os.path.isdir(d):
			print("Can't create directory: %s - a directory with this name already exists. It will be used for the results of the analysis." % d)	
			return True
		else: 
			print("Can't create directory: %s." % d)	
			return False

def _schedule(n, nproc):
	"""guided scheduler""" 
	start = 0
	size = (n - start) // nproc
	while size > 100:
		yield slice(start, start + size)
		start += size
		size = (n - start) // nproc
	yield slice(start, n + 1)
	return
	
def _mahalanobis_distances(m, L, X):
	cX = X - m[np.newaxis, :]
	tmp = solve_triangular(L, cX.T, lower=True).T
	tmp **= 2
	# return np.sqrt(tmp.sum(axis=1))
	return tmp.sum(axis=1)

def mahalanobis_distances(m, S, X, parallel=True):
	L = cholesky(S, lower=True)
	n = X.shape[0]
	if parallel:
		nproc = cpu_count()
		res = (Parallel(n_jobs= -1) # -1
				(delayed(_mahalanobis_distances)
				  (m, L, X[s, :])
					for s in _schedule(n, nproc)))
		return np.hstack(res)
	else:
		return _mahalanobis_distances(m, L, X)


# scipy.spatial.distance.mahalanobis for comparison

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
		res = (Parallel(n_jobs= -1) #-1
				(delayed(_mahalanobis_distances_scipy)
				 (m, SI, X[s,:])
				   for s in _schedule(n,nproc)))
		return np.hstack(res)
	else:
		return _mahalanobis_distances_scipy(m, SI, X)
###

gmaps = 0
nwpath = ''

def initglobalmaps():
	
	#	SHARED FOLDER PATH OR LOCAL DIRECTORY
	indir = os.path.join(os.path.sep, nwpath, 'inVars')
	print indir
	herbf = 'herb.tif'
	treef = 'tree.tif'
	ndvimaxf = 'ndvimax.tif'
	ndviminf = 'ndvimin.tif'
	ndwif = 'ndwi.tif'
	slopef = 'slope.tif'
	biof = 'bio.tif'
	eprf = 'eprsqrt2.tif'
	pref = 'pre.tif'
	
	biof_globalfile = os.path.join(os.path.sep, indir, biof)
	print biof_globalfile
	global	src_ds_bio_global
	src_ds_bio_global = gdal.Open(biof_globalfile)
	global	bio_global
	bio_global = src_ds_bio_global.GetRasterBand(1)
	global	gt_bio_global
	gt_bio_global = src_ds_bio_global.GetGeoTransform()
	print 'bio'
	
	pref_globalfile = os.path.join(os.path.sep, indir, pref)
	global	src_ds_pre_global
	src_ds_pre_global = gdal.Open(pref_globalfile)
	global	pre_global
	pre_global = src_ds_pre_global.GetRasterBand(1)
	global	gt_pre_global
	gt_pre_global = src_ds_pre_global.GetGeoTransform()
	print 'pre'
	
	eprf_globalfile = os.path.join(os.path.sep, indir, eprf)
	global	src_ds_epr_global
	src_ds_epr_global = gdal.Open(eprf_globalfile)
	global	epr_global
	epr_global = src_ds_epr_global.GetRasterBand(1)
	global	gt_epr_global
	gt_epr_global = src_ds_epr_global.GetGeoTransform()
	print 'epr'
	
	herbf_globalfile = os.path.join(os.path.sep, indir, herbf)
	global	src_ds_herb_global
	src_ds_herb_global = gdal.Open(herbf_globalfile)
	global	herb_global
	herb_global = src_ds_herb_global.GetRasterBand(1)
	global	gt_herb_global
	gt_herb_global = src_ds_herb_global.GetGeoTransform()
	print 'herb'
	
	ndvimaxf_globalfile = os.path.join(os.path.sep, indir, ndvimaxf)
	global	src_ds_ndvimax_global
	src_ds_ndvimax_global = gdal.Open(ndvimaxf_globalfile)
	global	ndvimax_global
	ndvimax_global = src_ds_ndvimax_global.GetRasterBand(1)
	global	gt_ndvimax_global
	gt_ndvimax_global = src_ds_ndvimax_global.GetGeoTransform()
	print 'ndvimax'
	
	ndviminf_globalfile = os.path.join(os.path.sep, indir, ndviminf)
	global	src_ds_ndvimin_global
	src_ds_ndvimin_global = gdal.Open(ndviminf_globalfile)
	global	ndvimin_global
	ndvimin_global = src_ds_ndvimin_global.GetRasterBand(1)
	global	gt_ndvimin_global
	gt_ndvimin_global = src_ds_ndvimin_global.GetGeoTransform()
	print 'ndvimin'
	
	ndwif_globalfile = os.path.join(os.path.sep, indir, ndwif)
	global	src_ds_ndwi_global
	src_ds_ndwi_global =	gdal.Open(ndwif_globalfile)
	global	ndwi_global
	ndwi_global = src_ds_ndwi_global.GetRasterBand(1)
	global	gt_ndwi_global
	gt_ndwi_global = src_ds_ndwi_global.GetGeoTransform()
	print 'ndwi'
	
	slopef_globalfile = os.path.join(os.path.sep, indir, slopef)
	global	src_ds_slope_global
	src_ds_slope_global = gdal.Open(slopef_globalfile)
	global	slope_global
	slope_global = src_ds_slope_global.GetRasterBand(1)
	global	gt_slope_global
	gt_slope_global = src_ds_slope_global.GetGeoTransform()
	print 'slope'
	
	treef_globalfile = os.path.join(os.path.sep, indir, treef)	
	global	src_ds_tree_global
	src_ds_tree_global = gdal.Open(treef_globalfile)
	global	tree_global
	tree_global = src_ds_tree_global.GetRasterBand(1)
	global	gt_tree_global
	gt_tree_global = src_ds_tree_global.GetGeoTransform()
	print 'tree'
	print "Global variables imported"
	global	gmaps
	gmaps = 1

def	ehabitat(ecor,nw,nwpathout):

	global	nwpath
	if nw=='':
		nwpath = os.getcwd()
	else:
		nwpath = nw
		
	if gmaps == 0:
		initglobalmaps()
	if nwpathout=='':
		#outdir = 'results'	#	ToDo: locally	create	folder "results"	if it	does	not	exist!
		outdir = os.path.join(os.path.sep, os.getcwd(), 'results')
		safelyMakeDir(outdir)
	else:
		#outdir = nwpathout+'/results'	#	SHARED	FOLDER	PATH
		outdir = os.path.join(os.path.sep, nwpathout, 'results')
		safelyMakeDir(outdir)
		
	treepamin = treepamax = eprpamin = eprpamax = prepamin = prepamax = biopamin = biopamax = slopepamin = slopepamax = ndwipamin = ndwipamax = ndvimaxpamin = ndvimaxpamax = ndviminpamin = ndviminpamax = hpamin = hpamax = None
	s = nd.generate_binary_structure(2,2)	#	most	restrictive	pattern	for	the	landscape	patches
	#	LOCAL FOLDER
	csvname1 = os.path.join(os.path.sep, outdir, 'ecoregs_done.csv')
	print csvname1
	if os.path.isfile(csvname1) == False:
		wb = open(csvname1,'a')
		wb.write('None')
		wb.write('\n')
		wb.close()
	#	LOCAL FOLDER	
	csvname = os.path.join(os.path.sep, outdir, 'hri_results.csv')
	print csvname
	if os.path.isfile(csvname) == False:
		wb = open(csvname,'a')
		wb.write('ecoregion wdpaid averpasim hr2aver pxpa hr1insumaver hriaver nfeatsaver lpratio lpratio2 numpszok lpmaxsize aggregation treepamin treepamax eprpamin eprpamax prepamin prepamax biopamin biopamax slopepamin slopepamax ndwipamin ndwipamax ndvimaxpamin ndvimaxpamax ndviminpamin ndviminpamax hpamin hpamax treepamean eprpamean prepamean biopamean slopepamean ndwipamean ndvimaxpamean ndviminpamean hpamean')
		wb.write('\n')
		wb.close()
	treepamean = eprpamean = prepamean = biopamean = slopepamean = ndwipamean = ndvimaxpamean = ndviminpamean = hpamean = None
	ef = 'eco_'+str(ecor)+'.tif'
	ecofile = os.path.join(os.path.sep, nwpath, 'ecoregs', ef)
	#ecofile = os.path.join(os.path.sep, nwpath, os.path.sep,'ecoregs', os.path.sep, ef)
	print ecofile
	avail = os.path.isfile(ecofile)
	if avail == True:
		eco_csv = str(ecor)+'.csv'
		print eco_csv
		ecoparksf = os.path.join(os.path.sep, nwpath, 'pas', eco_csv)
		#ecoparksf = os.path.join(os.path.sep, nwpath, os.path.sep, 'pas', os.path.sep, eco_csv)
		print ecoparksf
		#ecoparksf = nwpath+'/pas/'+str(ecor)+'.csv'
		src_ds_eco = gdal.Open(ecofile)
		eco = src_ds_eco.GetRasterBand(1)
		eco_mask0 = eco.ReadAsArray(0,0,eco.XSize,eco.YSize).astype(np.int32)
		eco_mask = eco_mask0.flatten()
		gt_eco = src_ds_eco.GetGeoTransform()
		print 'eco mask'
		xoff = int((gt_eco[0]-gt_epr_global[0])/1000)
		yoff = int((gt_epr_global[3]-gt_eco[3])/1000)
		epr_eco_bb0 = epr_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		epr_eco_bb = epr_eco_bb0.flatten()
		epr_eco0 = np.where(eco_mask == 1,	(epr_eco_bb),(0))
		epr_eco = np.where(epr_eco0 == 65535.0,	(float('NaN')),(epr_eco0))
		maskepr = np.isnan(epr_eco)
		epr_eco[maskepr] = np.interp(np.flatnonzero(maskepr),	np.flatnonzero(~maskepr),	epr_eco[~maskepr])
		print 'eco epr'
		xoff = int((gt_eco[0]-gt_slope_global[0])/1000)
		yoff = int((gt_slope_global[3]-gt_eco[3])/1000)
		slope_eco_bb0 = slope_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		slope_eco_bb = slope_eco_bb0.flatten()
		slope_eco0 = np.where(eco_mask == 1,	(slope_eco_bb),(0))
		slope_eco = np.where(slope_eco0 == 65535.0,	(float('NaN')),(slope_eco0))
		maskslope = np.isnan(slope_eco)
		slope_eco[maskslope] = np.interp(np.flatnonzero(maskslope),	np.flatnonzero(~maskslope),	slope_eco[~maskslope])
		print 'eco slope'
		xoff = int((gt_eco[0]-gt_ndvimax_global[0])/1000)
		yoff = int((gt_ndvimax_global[3]-gt_eco[3])/1000)
		ndvimax_eco_bb0 = ndvimax_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		ndvimax_eco_bb = ndvimax_eco_bb0.flatten()
		ndvimax_eco0 = np.where(eco_mask == 1,	(ndvimax_eco_bb),(0))
		ndvimax_eco = np.where(ndvimax_eco0 == 65535.0,	(float('NaN')),(ndvimax_eco0))
		maskndvimax = np.isnan(ndvimax_eco)
		ndvimax_eco[maskndvimax] = np.interp(np.flatnonzero(maskndvimax),	np.flatnonzero(~maskndvimax),	ndvimax_eco[~maskndvimax])
		print 'eco ndvimax'
		xoff = int((gt_eco[0]-gt_ndvimin_global[0])/1000)
		yoff = int((gt_ndvimin_global[3]-gt_eco[3])/1000)
		ndvimin_eco_bb0 = ndvimin_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		ndvimin_eco_bb = ndvimin_eco_bb0.flatten()
		ndvimin_eco0 = np.where(eco_mask == 1,	(ndvimin_eco_bb),(0))
		ndvimin_eco = np.where(ndvimin_eco0 == 65535.0,	(float('NaN')),(ndvimin_eco0))
		maskndvimin = np.isnan(ndvimin_eco)
		ndvimin_eco[maskndvimin] = np.interp(np.flatnonzero(maskndvimin),	np.flatnonzero(~maskndvimin),	ndvimin_eco[~maskndvimin])
		print 'eco ndvimin'
		xoff = int((gt_eco[0]-gt_ndwi_global[0])/1000)
		yoff = int((gt_ndwi_global[3]-gt_eco[3])/1000)
		ndwi_eco_bb0 = ndwi_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		ndwi_eco_bb = ndwi_eco_bb0.flatten()
		ndwi_eco0 = np.where(eco_mask == 1,	(ndwi_eco_bb),(0))
		ndwi_eco = np.where(ndwi_eco0 == 255.0,	(float('NaN')),(ndwi_eco0))
		maskndwi = np.isnan(ndwi_eco)
		ndwi_eco[maskndwi] = np.interp(np.flatnonzero(maskndwi),	np.flatnonzero(~maskndwi),	ndwi_eco[~maskndwi])
		print 'eco ndwi'
		xoff = int((gt_eco[0]-gt_pre_global[0])/1000)
		yoff = int((gt_pre_global[3]-gt_eco[3])/1000)
		pre_eco_bb0 = pre_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		pre_eco_bb = pre_eco_bb0.flatten()
		pre_eco0 = np.where(eco_mask == 1,	(pre_eco_bb),(0))
		pre_eco = np.where(pre_eco0 == 65535.0,	(float('NaN')),(pre_eco0))
		maskpre = np.isnan(pre_eco)
		pre_eco[maskpre] = np.interp(np.flatnonzero(maskpre),	np.flatnonzero(~maskpre),	pre_eco[~maskpre])
		print 'eco pre'
		xoff = int((gt_eco[0]-gt_bio_global[0])/1000)
		yoff = int((gt_bio_global[3]-gt_eco[3])/1000)
		bio_eco_bb0 = bio_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		bio_eco_bb = bio_eco_bb0.flatten()
		bio_eco0 = np.where(eco_mask == 1,	(bio_eco_bb),(0))
		bio_eco = np.where(bio_eco0 == 65535.0,	(float('NaN')),(bio_eco0))
		maskbio = np.isnan(bio_eco)
		bio_eco[maskbio] = np.interp(np.flatnonzero(maskbio),	np.flatnonzero(~maskbio),	bio_eco[~maskbio])
		print 'eco bio'
		xoff = int((gt_eco[0]-gt_tree_global[0])/1000)
		yoff = int((gt_tree_global[3]-gt_eco[3])/1000)
		tree_eco_bb0 = tree_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		tree_eco_bb = tree_eco_bb0.flatten()
		tree_eco0 = np.where(eco_mask == 1,	(tree_eco_bb),(0))
		tree_eco = np.where(tree_eco0 == 255.0,	(float('NaN')),(tree_eco0))
		masktree = np.isnan(tree_eco)
		tree_eco[masktree] = np.interp(np.flatnonzero(masktree),	np.flatnonzero(~masktree),	tree_eco[~masktree])
		print 'eco tree'
		xoff = int((gt_eco[0]-gt_herb_global[0])/1000)
		yoff = int((gt_herb_global[3]-gt_eco[3])/1000)
		herb_eco_bb0 = herb_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		herb_eco_bb = herb_eco_bb0.flatten()
		herb_eco0 = np.where(eco_mask == 1,	(herb_eco_bb),(0))
		herb_eco = np.where(herb_eco0 == 255.0,	(float('NaN')),(herb_eco0))
		maskherb = np.isnan(herb_eco)
		herb_eco[maskherb] = np.interp(np.flatnonzero(maskherb),	np.flatnonzero(~maskherb),	herb_eco[~maskherb])
		print 'eco herb'
		ind_eco0 = np.column_stack((bio_eco,pre_eco,epr_eco,herb_eco,ndvimax_eco,ndvimin_eco,ndwi_eco,slope_eco,tree_eco))
		print 'ecovars stacked'
		
		print ecoparksf
		pa_list0 = np.genfromtxt(ecoparksf,dtype='string')	# crear este archivo en subpas!
		pa_list = np.unique(pa_list0)
		n = len(pa_list)
		for	px in range(0,n): #	0,n

			pa = pa_list[px]
			print pa

			outfile = os.path.join(os.path.sep, outdir, str(ecor)+'_'+str(pa)+'.tif')
			outfile2 = os.path.join(os.path.sep, outdir, str(ecor)+'_'+str(pa)+'_lp.tif')
			outfile3 = os.path.join(os.path.sep, outdir, str(ecor)+'_'+str(pa)+'_mask.tif')
			#outfile = outdir+'/'+str(ecor)+'_'+str(pa)+'.tif'	#	LOCAL FOLDER
			pa_infile = 'pa_'+str(pa)+'.tif'

			pa4 = os.path.join(os.path.sep, nwpath, 'pas', pa_infile)
			#pa4 = os.path.join(os.path.sep, nwpath, os.path.sep, 'pas', os.path.sep, pa_infile)
			print pa4
			#pa4 = nwpath+'/pas/pa_'+str(pa)+'.tif'

			dropcols = np.arange(9,dtype=int)
			done = os.path.isfile(outfile)
			avail2 = os.path.isfile(pa4)
			if done == False and avail2 == True:
				pafile=pa4
				src_ds_pa = gdal.Open(pafile)
				par = src_ds_pa.GetRasterBand(1)
				pa_mask0 = par.ReadAsArray(0,0,par.XSize,par.YSize).astype(np.int32)
				pa_mask = pa_mask0.flatten()
				ind = pa_mask >	0 #==int(pa)
				go = 1
				sum_pa_mask = sum(pa_mask[ind])#/int(pa)
				if sum_pa_mask < 3: go = 0	#	not	processing	areas	smaller	than	3	pixels
				print sum_pa_mask
				sum_pa_mask_inv = len(pa_mask[pa_mask == 0])
				print sum_pa_mask_inv
				print len(pa_mask)
				ratiogeom = 10000
				if sum_pa_mask > 0: ratiogeom = sum_pa_mask_inv/sum_pa_mask
				#print ratiogeom
				gt_pa = src_ds_pa.GetGeoTransform()
				xoff = int((gt_pa[0]-gt_pre_global[0])/1000)
				yoff = int((gt_pre_global[3]-gt_pa[3])/1000)
				if xoff>0 and yoff>0 and go == 1:
					num_bands=src_ds_eco.RasterCount
					driver = gdal.GetDriverByName("GTiff")
					dst_options = ['COMPRESS=LZW']
					dst_ds = driver.Create(	outfile,src_ds_eco.RasterXSize,src_ds_eco.RasterYSize,num_bands,gdal.GDT_Float32,dst_options)
					dst_ds.SetGeoTransform(	src_ds_eco.GetGeoTransform())
					dst_ds.SetProjection(	src_ds_eco.GetProjectionRef())
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
						tree_pa[mask2tree] = np.interp(np.flatnonzero(mask2tree),	np.flatnonzero(~mask2tree),	tree_pa[~mask2tree])
						tree_pa = np.random.random_sample(len(tree_pa),)/1000 + tree_pa
						print 'pa tree'

						treepamin = round(tree_pa.min(),2)
						treepamax = round(tree_pa.max(),2)
						treepamean = round(np.mean(tree_pa),2)
						print treepamin
						print treepamax
						treediff = abs(tree_pa.min()-tree_pa.max())
						if treediff < 0.001: dropcols[8] = -8

					xoff = int((gt_pa[0]-gt_epr_global[0])/1000)
					yoff = int((gt_epr_global[3]-gt_pa[3])/1000)
					epr_pa_bb0 = epr_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					epr_pa_bb = epr_pa_bb0.flatten()
					epr_pa0 = epr_pa_bb[ind]
					epr_pa = np.where(epr_pa0 == 65535.0,	(float('NaN')),(epr_pa0))
					mask2epr = np.isnan(epr_pa)
					if mask2epr.all() == True:
						dropcols[2] = -2
					else:
						epr_pa[mask2epr] = np.interp(np.flatnonzero(mask2epr),	np.flatnonzero(~mask2epr),	epr_pa[~mask2epr])
						epr_pa = np.random.random_sample(len(epr_pa),)/1000 + epr_pa
						print 'pa epr'

						eprpamin = round(epr_pa.min(),2)
						eprpamax = round(epr_pa.max(),2)
						eprpamean = round(np.mean(epr_pa),2)
						print eprpamin
						print eprpamax
						eprdiff = abs(epr_pa.min()-epr_pa.max())
						if eprdiff < 0.001: dropcols[2] = -2

					xoff = int((gt_pa[0]-gt_pre_global[0])/1000)
					yoff = int((gt_pre_global[3]-gt_pa[3])/1000)
					pre_pa_bb0 = pre_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					pre_pa_bb = pre_pa_bb0.flatten()
					pre_pa0 = pre_pa_bb[ind]
					pre_pa = np.where(pre_pa0 == 65535.0,	(float('NaN')),(pre_pa0))
					mask2pre = np.isnan(pre_pa)
					if mask2pre.all() == True:
						dropcols[1] = -1
					else:
						pre_pa[mask2pre] = np.interp(np.flatnonzero(mask2pre),	np.flatnonzero(~mask2pre),	pre_pa[~mask2pre])
						pre_pa = np.random.random_sample(len(pre_pa),)/1000 + pre_pa
						print 'pa pre'

						prepamin = round(pre_pa.min(),2)
						prepamax = round(pre_pa.max(),2)
						prepamean = round(np.mean(pre_pa),2)
						print prepamin
						print prepamax
						prediff = abs(pre_pa.min()-pre_pa.max())
						if prediff < 0.001: dropcols[1] = -1

					xoff = int((gt_pa[0]-gt_bio_global[0])/1000)
					yoff = int((gt_bio_global[3]-gt_pa[3])/1000)
					bio_pa_bb0 = bio_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					bio_pa_bb = bio_pa_bb0.flatten()
					bio_pa0 = bio_pa_bb[ind]
					bio_pa = np.where(bio_pa0 == 65535.0,	(float('NaN')),(bio_pa0))
					mask2bio = np.isnan(bio_pa)
					if mask2bio.all() == True:
						dropcols[0] = -0
					else:
						bio_pa[mask2bio] = np.interp(np.flatnonzero(mask2bio),	np.flatnonzero(~mask2bio),	bio_pa[~mask2bio])
						bio_pa = np.random.random_sample(len(bio_pa),)/1000 + bio_pa
						print 'pa bio'

						biopamin = round(bio_pa.min(),2)
						biopamax = round(bio_pa.max(),2)
						biopamean = round(np.mean(bio_pa),2)
						print biopamin
						print biopamax
						biodiff = abs(bio_pa.min()-bio_pa.max())
						if biodiff < 0.001: dropcols[0] = -0

					xoff = int((gt_pa[0]-gt_slope_global[0])/1000)
					yoff = int((gt_slope_global[3]-gt_pa[3])/1000)
					slope_pa_bb0 = slope_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					slope_pa_bb = slope_pa_bb0.flatten()
					slope_pa0 = slope_pa_bb[ind]
					slope_pa = np.where(slope_pa0 == 65535.0,	(float('NaN')),(slope_pa0))
					mask2slope = np.isnan(slope_pa)
					if mask2slope.all() == True:
						dropcols[7] = -7
					else:
						slope_pa[mask2slope] = np.interp(np.flatnonzero(mask2slope),	np.flatnonzero(~mask2slope),	slope_pa[~mask2slope])
						slope_pa = np.random.random_sample(len(slope_pa),)/1000 + slope_pa
						print 'pa slope'

						slopepamin = round(slope_pa.min(),2)
						slopepamax = round(slope_pa.max(),2)
						slopepamean = round(np.mean(slope_pa),2)
						print slopepamin
						print slopepamax
						slopediff = abs(slope_pa.min()-slope_pa.max())
						if slopediff < 0.001: dropcols[7] = -7

					xoff = int((gt_pa[0]-gt_ndwi_global[0])/1000)
					yoff = int((gt_ndwi_global[3]-gt_pa[3])/1000)
					ndwi_pa_bb0 = ndwi_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					ndwi_pa_bb = ndwi_pa_bb0.flatten()
					ndwi_pa0 = ndwi_pa_bb[ind]
					ndwi_pa = np.where(ndwi_pa0 == 255.0,	(float('NaN')),(ndwi_pa0))
					mask2ndwi = np.isnan(ndwi_pa)
					if mask2ndwi.all() == True:
						dropcols[6] = -6
					else:
						ndwi_pa[mask2ndwi] = np.interp(np.flatnonzero(mask2ndwi),	np.flatnonzero(~mask2ndwi),	ndwi_pa[~mask2ndwi])
						ndwi_pa = np.random.random_sample(len(ndwi_pa),)/1000 + ndwi_pa
						print 'pa ndwi'

						ndwipamin = round(ndwi_pa.min(),2)
						ndwipamax = round(ndwi_pa.max(),2)
						ndwipamean = round(np.mean(ndwi_pa),2)
						print ndwipamin
						print ndwipamax
						ndwidiff = abs(ndwi_pa.min()-ndwi_pa.max())
						if ndwidiff < 0.001: dropcols[6] = -6

					xoff = int((gt_pa[0]-gt_ndvimax_global[0])/1000)
					yoff = int((gt_ndvimax_global[3]-gt_pa[3])/1000)
					ndvimax_pa_bb0 = ndvimax_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					ndvimax_pa_bb = ndvimax_pa_bb0.flatten()
					ndvimax_pa0 = ndvimax_pa_bb[ind]
					ndvimax_pa = np.where(ndvimax_pa0 == 65535.0,	(float('NaN')),(ndvimax_pa0))
					mask2ndvimax = np.isnan(ndvimax_pa)
					if mask2ndvimax.all() == True:
						dropcols[4] = -4
					else:
						ndvimax_pa[mask2ndvimax] = np.interp(np.flatnonzero(mask2ndvimax),	np.flatnonzero(~mask2ndvimax),	ndvimax_pa[~mask2ndvimax])
						ndvimax_pa = np.random.random_sample(len(ndvimax_pa),)/1000 + ndvimax_pa
						print 'pa ndvimax'

						ndvimaxpamin = round(ndvimax_pa.min(),2)
						ndvimaxpamax = round(ndvimax_pa.max(),2)
						ndvimaxpamean = round(np.mean(ndvimax_pa),2)
						print ndvimaxpamin
						print ndvimaxpamax
						ndvimaxdiff = abs(ndvimax_pa.min()-ndvimax_pa.max())
						if ndvimaxdiff < 0.001: dropcols[4] = -4

					xoff = int((gt_pa[0]-gt_ndvimin_global[0])/1000)
					yoff = int((gt_ndvimin_global[3]-gt_pa[3])/1000)
					ndvimin_pa_bb0 = ndvimin_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					ndvimin_pa_bb = ndvimin_pa_bb0.flatten()
					ndvimin_pa0 = ndvimin_pa_bb[ind]
					ndvimin_pa = np.where(ndvimin_pa0 == 65535.0,	(float('NaN')),(ndvimin_pa0))
					mask2ndvimin = np.isnan(ndvimin_pa)
					if mask2ndvimin.all() == True:
						dropcols[5] = -5
					else:
						ndvimin_pa[mask2ndvimin] = np.interp(np.flatnonzero(mask2ndvimin),	np.flatnonzero(~mask2ndvimin),	ndvimin_pa[~mask2ndvimin])
						ndvimin_pa = np.random.random_sample(len(ndvimin_pa),)/1000 + ndvimin_pa
						print 'pa ndvimin'

						ndviminpamin = round(ndvimin_pa.min(),2)
						ndviminpamax = round(ndvimin_pa.max(),2)
						ndviminpamean = round(np.mean(ndvimin_pa),2)
						print ndviminpamin
						print ndviminpamax
						ndvimindiff = abs(ndvimin_pa.min()-ndvimin_pa.max())
						if ndvimindiff < 0.001: dropcols[5] = -5

					xoff = int((gt_pa[0]-gt_herb_global[0])/1000)
					yoff = int((gt_herb_global[3]-gt_pa[3])/1000)
					herb_pa_bb0 = herb_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					herb_pa_bb = herb_pa_bb0.flatten()
					herb_pa0 = herb_pa_bb[ind]
					herb_pa = np.where(herb_pa0 == 255.0,	(float('NaN')),(herb_pa0))
					mask2herb = np.isnan(herb_pa)
					if mask2herb.all() == True:
						dropcols[3] = -3
					else:
						herb_pa[mask2herb] = np.interp(np.flatnonzero(mask2herb),	np.flatnonzero(~mask2herb),	herb_pa[~mask2herb])
						herb_pa = np.random.random_sample(len(herb_pa),)/1000 + herb_pa
						print 'pa herb'

						hpamin = round(herb_pa.min(),2)
						hpamax = round(herb_pa.max(),2)
						hpamean = round(np.mean(herb_pa),2)
						print hpamin
						print hpamax
						hdiff = abs(herb_pa.min()-herb_pa.max())
						if hdiff < 0.001: dropcols[3] = -3

					cols = dropcols[dropcols>=0]
					ind_pa0 = np.column_stack((bio_pa,pre_pa,epr_pa,herb_pa,ndvimax_pa,ndvimin_pa,ndwi_pa,slope_pa,tree_pa))
					ind_pa = ind_pa0[:,cols]
					ind_eco = ind_eco0[:,cols]
					print ind_pa.shape
					hr1sum = hr1insum = indokpsz = pszok = sumpszok = lpratio2 = numpszok = hr1averpa = hr3aver = hr2aver = pszmax = num_featuresaver = lpratio = hr1medianpa = hr1insumaver = pxpa = aggregation = None
					print "PA masked"
					#print ind_pa
					if ind_pa.shape[0]>4 and ind_pa.shape[1]>1: 
						Ymean = np.mean(ind_pa,axis=0)
						print 'Max. mean value is '+ str(Ymean.max())
						print "Ymean ok"
						Ycov = np.cov(ind_pa,rowvar=False)
						print 'Max. cov value is '+ str(Ycov.max())
						print "Ycov	ok"
						#mh = mahalanobis_distances(Ymean,	Ycov,	ind_eco,	parallel=False)
						#mh2 = mahalanobis_distances(Ymean,	Ycov,	ind_eco,	parallel=True)
						mh2 = mahalanobis_distances_scipy(Ymean,	Ycov,	ind_eco,	parallel=True) # previous working version
						#mh2 = mahalanobis_distances_scipy(Ymean,	Ycov,	ind_eco,	parallel=False)
						maxmh=mh2.max()
						print 'Max. mh value is '+ str(maxmh)
						print 'Max. mh value is nan: '+ str(np.isnan(maxmh))
						mh = mh2*mh2
						print "mh ok"
						pmh = chi2.sf(mh,len(cols)).reshape((eco.YSize,eco.XSize)) # chisqprob
						pmhh = np.where(pmh	<=	0.001,None,	pmh)
						print "pmh ok"	#	quitar	valores	muy	bajos!
						pmhhmax = pmhh.max()
						print 'Max. similarity value is '+ str(pmhhmax)
						dst_ds.GetRasterBand(1).WriteArray(pmhh)
						dst_ds = None
						hr11 = np.where(pmhh>0,1,0) # 0.5
						hr1 = hr11.flatten()
						hr1sum = sum(hr1)
						print 'Number of pixels with similarity higher than 0 is '+str(hr1sum)
						hr1insumaver = hr1insum = 0
						hr1sumaver = hr1sum
						src_ds_sim = gdal.Open(outfile)
						sim = src_ds_sim.GetRasterBand(1)
						gt_sim = src_ds_sim.GetGeoTransform()
						xoff = int((gt_pa[0]-gt_sim[0])/1000)
						yoff = int((gt_sim[3]-gt_pa[3])/1000)
						xextentpa = xoff + par.XSize
						yextentpa = yoff + par.YSize
						xless = sim.XSize - xextentpa
						yless = sim.YSize - yextentpa
						xsize = par.XSize
						ysize = par.YSize
						if xoff>0 and yoff>0 and pmhhmax>0.01 and hr1sum>1 and maxmh!=float('NaN'):#and ratiogeom < 100: #	also	checks	if results	are	not	empty

							# reading the similarity ecoregion without the PA (tmp mask)
							os.system('gdal_merge.py '+str(ecofile)+' '+str(pa4)+' -o '+str(outfile3)+' -ot Int32')
							hri_pa_bb03 = sim.ReadAsArray().astype(np.float32)
							hri_pa_bb3 = hri_pa_bb03.flatten()
							
							src_ds_sim2 = gdal.Open(outfile3)
							sim2 = src_ds_sim2.GetRasterBand(1)
							gt_sim2 = src_ds_sim2.GetGeoTransform()
							hri_pa_bb02 = sim2.ReadAsArray().astype(np.int32)
							#hri_pa_bb2 = hri_pa_bb02.flatten()
							hri_pa_bb02_max = hri_pa_bb02.max()
							print 'PA: '+str(pa)
							print 'PA (= max) value from mask = '+str(hri_pa_bb02_max)
							if hri_pa_bb02.shape == hri_pa_bb03.shape:
							 hri_pa02 = np.where(hri_pa_bb02 == pa,0,hri_pa_bb03) # hri_pa_bb02_max


							 if xless < 0: xsize = xsize + xless
							 if yless < 0: ysize = ysize + yless
							 hri_pa_bb0 = sim.ReadAsArray(xoff,yoff,xsize,ysize).astype(np.float32)
							 hri_pa_bb = hri_pa_bb0.flatten()
							 indd = hri_pa_bb > 0
							 hri_pa0 = hri_pa_bb[indd]
							 print 'Total number of pixels with similarity values in PA: '+str(len(hri_pa0))
							 hr1averpa = round(np.mean(hri_pa0[~np.isnan(hri_pa0)]),2)
							 #print hr1averpa
							 #hr1medianpa = np.median(hri_pa0[~np.isnan(hri_pa0)])
							 print 'mean similarity in the park is '+str(hr1averpa)
							 #hr1insum = sum(np.where(hri_pa0 >= 0.5,	1,0))	#	use	hr1averpa	as	threshold	instead!						
							 hr1inaver = np.where(hri_pa0 >= hr1averpa,	1,0)
							 hr1insumaver = sum(hr1inaver)
							 #print hr1insum
							 ##labeled_arrayin, num_featuresin = nd.label(hr1inaver,	structure=s)
							 hr1averr = np.where(hri_pa02 >= hr1averpa,	1,0) # pmhh
							 hr1aver = hr1averr.flatten()
							 print 'Total number of pixels with similarity values in ECO: '+str(sum(hr1aver))
							 labeled_arrayaver, num_featuresaver = nd.label(hr1averr,	structure=s)
							 print 'Nr of similar patches found: '+str(num_featuresaver)
							 if num_featuresaver > 0:
							  lbls = np.arange(1, num_featuresaver+1)
							  psizes = nd.labeled_comprehension(labeled_arrayaver, labeled_arrayaver, lbls, np.count_nonzero, float, 0) #-1
							  pszmax = psizes.max()#-hr1insumaver
							  dst_ds2 = driver.Create(outfile2,src_ds_eco.RasterXSize,src_ds_eco.RasterYSize,num_bands,gdal.GDT_Int32,dst_options)
							  dst_ds2.SetGeoTransform(src_ds_eco.GetGeoTransform())
							  dst_ds2.SetProjection(src_ds_eco.GetProjectionRef())
							  dst_ds2.GetRasterBand(1).WriteArray(labeled_arrayaver)
							  dst_ds2 = None
							  #num_feats = num_features - num_featuresaver
							  hr1sumaver = sum(hr1aver)
							  hr2aver = hr1sumaver #- hr1insumaver
							  pxpa = ind_pa.shape[0]
							  indokpsz = psizes >= pxpa
							  pszsok = psizes[indokpsz] # NEW
							  sumpszok = sum(pszsok)
							  lpratio=round(float(pszmax/pxpa),2)
							  lpratio2=round(float(sumpszok/pxpa),2)
							  numpszok = len(pszsok)
							  hr3aver = round(float(hr2aver/pxpa),2)
							  aggregation = round(float(hr2aver/num_featuresaver),2)
						#hr2 = hr1sumaver - hr1insumaver
						#print hr2
						#hr3 = float(hr2/ind_pa.shape[0])
						#print hr3
					wb = open(csvname,'a')
					var = str(ecor)+' '+str(pa)+' '+str(hr1averpa)+' '+str(hr2aver)+' '+str(pxpa)+' '+str(hr1insumaver)+' '+str(hr3aver)+' '+str(num_featuresaver)+' '+str(lpratio)+' '+str(lpratio2)+' '+str(numpszok)+' '+str(pszmax)+' '+str(aggregation)+' '+str(treepamin)+' '+str(treepamax)+' '+str(eprpamin)+' '+str(eprpamax)+' '+str(prepamin)+' '+str(prepamax)+' '+str(biopamin)+' '+str(biopamax)+' '+str(slopepamin)+' '+str(slopepamax)+' '+str(ndwipamin)+' '+str(ndwipamax)+' '+str(ndvimaxpamin)+' '+str(ndvimaxpamax)+' '+str(ndviminpamin)+' '+str(ndviminpamax)+' '+str(hpamin)+' '+str(hpamax)+' '+str(treepamean)+' '+str(eprpamean)+' '+str(prepamean)+' '+str(biopamean)+' '+str(slopepamean)+' '+str(ndwipamean)+' '+str(ndvimaxpamean)+' '+str(ndviminpamean)+' '+str(hpamean)#	exclude	PA!	#+' '+str(hr1p25pa)#	'+str(hr3)+'	+' '+str(hr1medianpa)+' '+str(num_features)+' '
					wb.write(var)
					wb.write('\n')
					wb.close()
					print "results exported"
					os.system('rm '+str(outfile3)) 
		wb = open(csvname1,'a')	#	LOCAL	FOLDER
		var = str(ecor)
		wb.write(var)
		wb.write('\n')
		wb.close()	
	print "END ECOREG: " + str(ecor)

def	run_batch():
#	if __name__ == '__main__':
	#from	datetime	import	datetime
	#import	numpy	as	np
	eco_list0 = np.genfromtxt(nwpath + 'pas' + os.path.sep + 'ecoregs.csv',dtype='int')	#	crear	este	archivo	en	subpas!
	eco_list = np.unique(eco_list0)
	#print eco_list
	m = len(eco_list)
	for	pm	in	range(0,m): #	3,m	#	0	without	the	negative	ecoregs!
		ecor = eco_list[pm]
		print ecor
		ehabitat(ecor,'','')
	print str(datetime.now())
	print "BATCH END"

