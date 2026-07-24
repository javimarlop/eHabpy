#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / GDAL 3+ (2024)
#### License: CC BY-SA 3.0
####
#### Shared core for getmeanvar.py and getmedianvar.py.
####
#### The original getmeanvar.py / getmedianvar.py contained nine almost identical
#### copies of the same function (ehabitat1 .. ehabitat9), one per segmentation
#### threshold level (1..9). The nine copies differed ONLY in a numeric suffix used
#### in two CSV file names. getmedianvar.py was in turn a byte-for-byte copy of
#### getmeanvar.py with np.mean replaced by np.median for the per-variable summary.
####
#### This module factors all of that into a single parameterized implementation:
####   * `k`      -> segmentation threshold level (1..9), used in the CSV names
####   * `aggfun` -> np.mean (getmeanvar) or np.median (getmedianvar)
####
#### Control files: 'csv/segm_done.csv', 'csv/park_*_*.csv'
#### Inputs variables: 9 input variables in inVars folder
#### Outputs: park.., park_hri_results and ecoregs_done csv files

from datetime import datetime
import numpy as np
import scipy.ndimage as nd
import os.path
import scipy
from scipy.linalg import cholesky, solve_triangular
from scipy.spatial import distance
from multiprocessing import cpu_count
import csv
import os
import sys
from osgeo import ogr, gdal

gdal.UseExceptions()
ogr.UseExceptions()


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


# ## Mahalanobis functions by Sturla Molden (kept for parity with the original module)
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


def _mahalanobis_distances_scipy(m, SI, X):
	n = X.shape[0]
	mahal = np.zeros(n)
	for i in range(X.shape[0]):
		x = X[i, :]
		mahal[i] = distance.mahalanobis(x, m, SI)
	return mahal


gmaps = 0
nwpath = ''
park = ''


def initglobalmaps():

	#	SHARED FOLDER PATH OR LOCAL DIRECTORY
	indir = os.path.join(nwpath, '../inVars')
	print(indir)
	herbf = 'herb.tif'
	treef = 'tree.tif'
	ndvimaxf = 'ndvimax.tif'
	ndviminf = 'ndvimin.tif'
	ndwif = 'ndwi.tif'
	slopef = 'slope.tif'
	biof = 'bio.tif'
	eprf = 'eprsqrt2.tif'
	pref = 'pre.tif'

	biof_globalfile = os.path.join(indir, biof)
	print(biof_globalfile)
	global	src_ds_bio_global
	src_ds_bio_global = gdal.Open(biof_globalfile)
	global	bio_global
	bio_global = src_ds_bio_global.GetRasterBand(1)
	global	gt_bio_global
	gt_bio_global = src_ds_bio_global.GetGeoTransform()
	print('bio')

	pref_globalfile = os.path.join(indir, pref)
	global	src_ds_pre_global
	src_ds_pre_global = gdal.Open(pref_globalfile)
	global	pre_global
	pre_global = src_ds_pre_global.GetRasterBand(1)
	global	gt_pre_global
	gt_pre_global = src_ds_pre_global.GetGeoTransform()
	print('pre')

	eprf_globalfile = os.path.join(indir, eprf)
	global	src_ds_epr_global
	src_ds_epr_global = gdal.Open(eprf_globalfile)
	global	epr_global
	epr_global = src_ds_epr_global.GetRasterBand(1)
	global	gt_epr_global
	gt_epr_global = src_ds_epr_global.GetGeoTransform()
	print('epr')

	herbf_globalfile = os.path.join(indir, herbf)
	global	src_ds_herb_global
	src_ds_herb_global = gdal.Open(herbf_globalfile)
	global	herb_global
	herb_global = src_ds_herb_global.GetRasterBand(1)
	global	gt_herb_global
	gt_herb_global = src_ds_herb_global.GetGeoTransform()
	print('herb')

	ndvimaxf_globalfile = os.path.join(indir, ndvimaxf)
	global	src_ds_ndvimax_global
	src_ds_ndvimax_global = gdal.Open(ndvimaxf_globalfile)
	global	ndvimax_global
	ndvimax_global = src_ds_ndvimax_global.GetRasterBand(1)
	global	gt_ndvimax_global
	gt_ndvimax_global = src_ds_ndvimax_global.GetGeoTransform()
	print('ndvimax')

	ndviminf_globalfile = os.path.join(indir, ndviminf)
	global	src_ds_ndvimin_global
	src_ds_ndvimin_global = gdal.Open(ndviminf_globalfile)
	global	ndvimin_global
	ndvimin_global = src_ds_ndvimin_global.GetRasterBand(1)
	global	gt_ndvimin_global
	gt_ndvimin_global = src_ds_ndvimin_global.GetGeoTransform()
	print('ndvimin')

	ndwif_globalfile = os.path.join(indir, ndwif)
	global	src_ds_ndwi_global
	src_ds_ndwi_global = gdal.Open(ndwif_globalfile)
	global	ndwi_global
	ndwi_global = src_ds_ndwi_global.GetRasterBand(1)
	global	gt_ndwi_global
	gt_ndwi_global = src_ds_ndwi_global.GetGeoTransform()
	print('ndwi')

	slopef_globalfile = os.path.join(indir, slopef)
	global	src_ds_slope_global
	src_ds_slope_global = gdal.Open(slopef_globalfile)
	global	slope_global
	slope_global = src_ds_slope_global.GetRasterBand(1)
	global	gt_slope_global
	gt_slope_global = src_ds_slope_global.GetGeoTransform()
	print('slope')

	treef_globalfile = os.path.join(indir, treef)
	global	src_ds_tree_global
	src_ds_tree_global = gdal.Open(treef_globalfile)
	global	tree_global
	tree_global = src_ds_tree_global.GetRasterBand(1)
	global	gt_tree_global
	gt_tree_global = src_ds_tree_global.GetGeoTransform()
	print('tree')
	print("Global variables imported")
	global	gmaps
	gmaps = 1


def ehabitat(ecor, nw, nwpathout, k, aggfun=np.mean):
	"""Compute per-HFT summary statistics for segmentation level ``k``.

	``k``      : segmentation threshold level (1..9); used in the CSV names.
	``aggfun`` : np.mean (getmeanvar) or np.median (getmedianvar).
	"""

	global	nwpath
	if nw == '':
		nwpath = os.getcwd()
	else:
		nwpath = nw

	if gmaps == 0:
		initglobalmaps()
	if nwpathout == '':
		outdir = nwpath
	else:
		outdir = nwpath

	treepamin = treepamax = eprpamin = eprpamax = prepamin = prepamax = biopamin = biopamax = slopepamin = slopepamax = ndwipamin = ndwipamax = ndvimaxpamin = ndvimaxpamax = ndviminpamin = ndviminpamax = hpamin = hpamax = treepavar = eprpavar = prepavar = biopavar = slopepavar = ndwipavar = ndvimaxpavar = ndviminpavar = hpavar = treepavar2 = eprpavar2 = prepavar2 = biopavar2 = slopepavar2 = ndwipavar2 = ndvimaxpavar2 = ndviminpavar2 = hpavar2 = None
	#	LOCAL FOLDER
	csvname1 = os.path.join(outdir, 'csv/ecoregs_done' + str(k) + '.csv')
	print(csvname1)
	if os.path.isfile(csvname1) == False:
		wb = open(csvname1, 'a')
		wb.write('None')
		wb.write('\n')
		wb.close()
	#	LOCAL FOLDER
	csvname = os.path.join(outdir, 'csv/park_' + str(park) + '_hri_results' + str(k) + '.csv')
	print(csvname)
	if os.path.isfile(csvname) == False:
		wb = open(csvname, 'a')
		wb.write('ecoregion segm_id treepamean eprpamean prepamean biopamean slopepamean ndwipamean ndvimaxpamean ndviminpamean hpamean treepavar eprpavar prepavar biopavar slopepavar ndwipavar ndvimaxpavar ndviminpavar hpavar sumpamask treepavar2 eprpavar2 prepavar2 biopavar2 slopepavar2 ndwipavar2 ndvimaxpavar2 ndviminpavar2 hpavar2')
		wb.write('\n')
		wb.close()
	treepamean = eprpamean = prepamean = biopamean = slopepamean = ndwipamean = ndvimaxpamean = ndviminpamean = hpamean = None

	eco_csv = 'csv/park_' + str(park) + '_' + str(ecor) + '.csv'
	print(eco_csv)
	ecoparksf = os.path.join(nwpath, eco_csv)

	pa_list0 = np.genfromtxt(ecoparksf, dtype=int)	# crear este archivo en subpas!
	pa_list = np.unique(pa_list0)
	n = len(pa_list)
	for	px in range(0, n):  # 0,n

		pa = pa_list[px]
		print(pa)

		outfile = os.path.join(outdir, 'csv/park' + '_' + str(pa) + '_' + 'hri_results' + str(ecor) + '.csv')
		pa_infile = 'tiffs/pa_' + str(pa) + '.tif'

		pa4 = os.path.join(nwpath, pa_infile)
		print(pa4)

		dropcols = np.arange(9, dtype=int)
		done = os.path.isfile(outfile)
		avail2 = os.path.isfile(pa4)
		if done == False and avail2 == True:
			pafile = pa4
			src_ds_pa = gdal.Open(pafile)
			par = src_ds_pa.GetRasterBand(1)
			pa_mask0 = par.ReadAsArray(0, 0, par.XSize, par.YSize).astype(np.int32)
			pa_mask = pa_mask0.flatten()
			ind = pa_mask > 0  # ==int(pa)
			sum_pa_mask = sum(pa_mask[ind])  # /int(pa)
			print(sum_pa_mask)
			sum_pa_mask_inv = len(pa_mask[pa_mask == 0])
			print(sum_pa_mask_inv)
			print(len(pa_mask))
			ratiogeom = 10000
			if sum_pa_mask > 0: ratiogeom = sum_pa_mask_inv / sum_pa_mask
			gt_pa = src_ds_pa.GetGeoTransform()
			xoff = int((gt_pa[0] - gt_pre_global[0]) / 1000)
			yoff = int((gt_pre_global[3] - gt_pa[3]) / 1000)
			if xoff > 0 and yoff > 0:  # and go == 1:
				xoff = int((gt_pa[0] - gt_tree_global[0]) / 1000)
				yoff = int((gt_tree_global[3] - gt_pa[3]) / 1000)
				tree_pa_bb0 = tree_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				tree_pa_bb = tree_pa_bb0.flatten()
				tree_pa0 = tree_pa_bb[ind]
				tree_pa = np.where(tree_pa0 == 255.0, (float('NaN')), (tree_pa0))
				mask2tree = np.isnan(tree_pa)
				if mask2tree.all() == True:
					dropcols[8] = -8
				else:
					tree_pa[mask2tree] = np.interp(np.flatnonzero(mask2tree), np.flatnonzero(~mask2tree), tree_pa[~mask2tree])
					tree_pa = np.random.random_sample(len(tree_pa),) / 1000 + tree_pa
					print('pa tree')

					treepamin = round(tree_pa.min(), 2)
					treepamax = round(tree_pa.max(), 2)
					treepamean = round(aggfun(tree_pa), 2)
					treepavar = round(np.var(tree_pa), 2)
					treepavar2 = treepavar * sum_pa_mask
					print(treepamin)
					print(treepamax)
					treediff = abs(tree_pa.min() - tree_pa.max())
					if treediff < 0.001: dropcols[8] = -8

				xoff = int((gt_pa[0] - gt_epr_global[0]) / 1000)
				yoff = int((gt_epr_global[3] - gt_pa[3]) / 1000)
				epr_pa_bb0 = epr_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				epr_pa_bb = epr_pa_bb0.flatten()
				epr_pa0 = epr_pa_bb[ind]
				epr_pa = np.where(epr_pa0 == 65535.0, (float('NaN')), (epr_pa0))
				mask2epr = np.isnan(epr_pa)
				if mask2epr.all() == True:
					dropcols[2] = -2
				else:
					epr_pa[mask2epr] = np.interp(np.flatnonzero(mask2epr), np.flatnonzero(~mask2epr), epr_pa[~mask2epr])
					epr_pa = np.random.random_sample(len(epr_pa),) / 1000 + epr_pa
					print('pa epr')

					eprpamin = round(epr_pa.min(), 2)
					eprpamax = round(epr_pa.max(), 2)
					eprpamean = round(aggfun(epr_pa), 2)
					eprpavar = round(np.var(epr_pa), 2)
					eprpavar2 = eprpavar * sum_pa_mask
					print(eprpamin)
					print(eprpamax)
					eprdiff = abs(epr_pa.min() - epr_pa.max())
					if eprdiff < 0.001: dropcols[2] = -2

				xoff = int((gt_pa[0] - gt_pre_global[0]) / 1000)
				yoff = int((gt_pre_global[3] - gt_pa[3]) / 1000)
				pre_pa_bb0 = pre_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				pre_pa_bb = pre_pa_bb0.flatten()
				pre_pa0 = pre_pa_bb[ind]
				pre_pa = np.where(pre_pa0 == 65535.0, (float('NaN')), (pre_pa0))
				mask2pre = np.isnan(pre_pa)
				if mask2pre.all() == True:
					dropcols[1] = -1
				else:
					pre_pa[mask2pre] = np.interp(np.flatnonzero(mask2pre), np.flatnonzero(~mask2pre), pre_pa[~mask2pre])
					pre_pa = np.random.random_sample(len(pre_pa),) / 1000 + pre_pa
					print('pa pre')

					prepamin = round(pre_pa.min(), 2)
					prepamax = round(pre_pa.max(), 2)
					prepamean = round(aggfun(pre_pa), 2)
					prepavar = round(np.var(pre_pa), 2)
					prepavar2 = prepavar * sum_pa_mask
					print(prepamin)
					print(prepamax)
					prediff = abs(pre_pa.min() - pre_pa.max())
					if prediff < 0.001: dropcols[1] = -1

				xoff = int((gt_pa[0] - gt_bio_global[0]) / 1000)
				yoff = int((gt_bio_global[3] - gt_pa[3]) / 1000)
				bio_pa_bb0 = bio_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				bio_pa_bb = bio_pa_bb0.flatten()
				bio_pa0 = bio_pa_bb[ind]
				bio_pa = np.where(bio_pa0 == 65535.0, (float('NaN')), (bio_pa0))
				mask2bio = np.isnan(bio_pa)
				if mask2bio.all() == True:
					dropcols[0] = -0
				else:
					bio_pa[mask2bio] = np.interp(np.flatnonzero(mask2bio), np.flatnonzero(~mask2bio), bio_pa[~mask2bio])
					bio_pa = np.random.random_sample(len(bio_pa),) / 1000 + bio_pa
					print('pa bio')
					biopamin = round(bio_pa.min(), 2)
					biopamax = round(bio_pa.max(), 2)
					biopamean = round(aggfun(bio_pa), 2)
					biopavar = round(np.var(bio_pa), 2)
					biopavar2 = biopavar * sum_pa_mask
					print(biopamin)
					print(biopamax)
					biodiff = abs(bio_pa.min() - bio_pa.max())
					if biodiff < 0.001: dropcols[0] = -0

				xoff = int((gt_pa[0] - gt_slope_global[0]) / 1000)
				yoff = int((gt_slope_global[3] - gt_pa[3]) / 1000)
				slope_pa_bb0 = slope_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				slope_pa_bb = slope_pa_bb0.flatten()
				slope_pa0 = slope_pa_bb[ind]
				slope_pa = np.where(slope_pa0 == 65535.0, (float('NaN')), (slope_pa0))
				mask2slope = np.isnan(slope_pa)
				if mask2slope.all() == True:
					dropcols[7] = -7
				else:
					slope_pa[mask2slope] = np.interp(np.flatnonzero(mask2slope), np.flatnonzero(~mask2slope), slope_pa[~mask2slope])
					slope_pa = np.random.random_sample(len(slope_pa),) / 1000 + slope_pa
					print('pa slope')

					slopepamin = round(slope_pa.min(), 2)
					slopepamax = round(slope_pa.max(), 2)
					slopepamean = round(aggfun(slope_pa), 2)
					slopepavar = round(np.var(slope_pa), 2)
					slopepavar2 = slopepavar * sum_pa_mask
					print(slopepamin)
					print(slopepamax)
					slopediff = abs(slope_pa.min() - slope_pa.max())
					if slopediff < 0.001: dropcols[7] = -7

				xoff = int((gt_pa[0] - gt_ndwi_global[0]) / 1000)
				yoff = int((gt_ndwi_global[3] - gt_pa[3]) / 1000)
				ndwi_pa_bb0 = ndwi_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				ndwi_pa_bb = ndwi_pa_bb0.flatten()
				ndwi_pa0 = ndwi_pa_bb[ind]
				ndwi_pa = np.where(ndwi_pa0 == 255.0, (float('NaN')), (ndwi_pa0))
				mask2ndwi = np.isnan(ndwi_pa)
				if mask2ndwi.all() == True:
					dropcols[6] = -6
				else:
					ndwi_pa[mask2ndwi] = np.interp(np.flatnonzero(mask2ndwi), np.flatnonzero(~mask2ndwi), ndwi_pa[~mask2ndwi])
					ndwi_pa = np.random.random_sample(len(ndwi_pa),) / 1000 + ndwi_pa
					print('pa ndwi')

					ndwipamin = round(ndwi_pa.min(), 2)
					ndwipamax = round(ndwi_pa.max(), 2)
					ndwipamean = round(aggfun(ndwi_pa), 2)
					ndwipavar = round(np.var(ndwi_pa), 2)
					ndwipavar2 = ndwipavar * sum_pa_mask
					print(ndwipamin)
					print(ndwipamax)
					ndwidiff = abs(ndwi_pa.min() - ndwi_pa.max())
					if ndwidiff < 0.001: dropcols[6] = -6

				xoff = int((gt_pa[0] - gt_ndvimax_global[0]) / 1000)
				yoff = int((gt_ndvimax_global[3] - gt_pa[3]) / 1000)
				ndvimax_pa_bb0 = ndvimax_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				ndvimax_pa_bb = ndvimax_pa_bb0.flatten()
				ndvimax_pa0 = ndvimax_pa_bb[ind]
				ndvimax_pa = np.where(ndvimax_pa0 == 65535.0, (float('NaN')), (ndvimax_pa0))
				mask2ndvimax = np.isnan(ndvimax_pa)
				if mask2ndvimax.all() == True:
					dropcols[4] = -4
				else:
					ndvimax_pa[mask2ndvimax] = np.interp(np.flatnonzero(mask2ndvimax), np.flatnonzero(~mask2ndvimax), ndvimax_pa[~mask2ndvimax])
					ndvimax_pa = np.random.random_sample(len(ndvimax_pa),) / 1000 + ndvimax_pa
					print('pa ndvimax')

					ndvimaxpamin = round(ndvimax_pa.min(), 2)
					ndvimaxpamax = round(ndvimax_pa.max(), 2)
					ndvimaxpamean = round(aggfun(ndvimax_pa), 2)
					ndvimaxpavar = round(np.var(ndvimax_pa), 2)
					ndvimaxpavar2 = ndvimaxpavar * sum_pa_mask
					print(ndvimaxpamin)
					print(ndvimaxpamax)
					ndvimaxdiff = abs(ndvimax_pa.min() - ndvimax_pa.max())
					if ndvimaxdiff < 0.001: dropcols[4] = -4

				xoff = int((gt_pa[0] - gt_ndvimin_global[0]) / 1000)
				yoff = int((gt_ndvimin_global[3] - gt_pa[3]) / 1000)
				ndvimin_pa_bb0 = ndvimin_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				ndvimin_pa_bb = ndvimin_pa_bb0.flatten()
				ndvimin_pa0 = ndvimin_pa_bb[ind]
				ndvimin_pa = np.where(ndvimin_pa0 == 65535.0, (float('NaN')), (ndvimin_pa0))
				mask2ndvimin = np.isnan(ndvimin_pa)
				if mask2ndvimin.all() == True:
					dropcols[5] = -5
				else:
					ndvimin_pa[mask2ndvimin] = np.interp(np.flatnonzero(mask2ndvimin), np.flatnonzero(~mask2ndvimin), ndvimin_pa[~mask2ndvimin])
					ndvimin_pa = np.random.random_sample(len(ndvimin_pa),) / 1000 + ndvimin_pa
					print('pa ndvimin')

					ndviminpamin = round(ndvimin_pa.min(), 2)
					ndviminpamax = round(ndvimin_pa.max(), 2)
					ndviminpamean = round(aggfun(ndvimin_pa), 2)
					ndviminpavar = round(np.var(ndvimin_pa), 2)
					ndviminpavar2 = ndviminpavar * sum_pa_mask
					print(ndviminpamin)
					print(ndviminpamax)
					ndvimindiff = abs(ndvimin_pa.min() - ndvimin_pa.max())
					if ndvimindiff < 0.001: dropcols[5] = -5

				xoff = int((gt_pa[0] - gt_herb_global[0]) / 1000)
				yoff = int((gt_herb_global[3] - gt_pa[3]) / 1000)
				herb_pa_bb0 = herb_global.ReadAsArray(xoff, yoff, par.XSize, par.YSize).astype(np.float32)
				herb_pa_bb = herb_pa_bb0.flatten()
				herb_pa0 = herb_pa_bb[ind]
				herb_pa = np.where(herb_pa0 == 255.0, (float('NaN')), (herb_pa0))
				mask2herb = np.isnan(herb_pa)
				if mask2herb.all() == True:
					dropcols[3] = -3
				else:
					herb_pa[mask2herb] = np.interp(np.flatnonzero(mask2herb), np.flatnonzero(~mask2herb), herb_pa[~mask2herb])
					herb_pa = np.random.random_sample(len(herb_pa),) / 1000 + herb_pa
					print('pa herb')

					hpamin = round(herb_pa.min(), 2)
					hpamax = round(herb_pa.max(), 2)
					hpamean = round(aggfun(herb_pa), 2)
					hpavar = round(np.var(herb_pa), 2)
					hpavar2 = hpavar * sum_pa_mask
					print(hpamin)
					print(hpamax)
					hdiff = abs(herb_pa.min() - herb_pa.max())
					if hdiff < 0.001: dropcols[3] = -3

				hr1sum = hr1insum = indokpsz = pszok = sumpszok = lpratio2 = numpszok = hr1averpa = hr3aver = hr2aver = pszmax = num_featuresaver = lpratio = hr1medianpa = hr1insumaver = pxpa = aggregation = None
				print("PA masked")
				wb = open(csvname, 'a')
				var = str(ecor) + ' ' + str(pa) + ' ' + str(treepamean) + ' ' + str(eprpamean) + ' ' + str(prepamean) + ' ' + str(biopamean) + ' ' + str(slopepamean) + ' ' + str(ndwipamean) + ' ' + str(ndvimaxpamean) + ' ' + str(ndviminpamean) + ' ' + str(hpamean) + ' ' + str(treepavar) + ' ' + str(eprpavar) + ' ' + str(prepavar) + ' ' + str(biopavar) + ' ' + str(slopepavar) + ' ' + str(ndwipavar) + ' ' + str(ndvimaxpavar) + ' ' + str(ndviminpavar) + ' ' + str(hpavar) + ' ' + str(sum_pa_mask) + ' ' + str(treepavar2) + ' ' + str(eprpavar2) + ' ' + str(prepavar2) + ' ' + str(biopavar2) + ' ' + str(slopepavar2) + ' ' + str(ndwipavar2) + ' ' + str(ndvimaxpavar2) + ' ' + str(ndviminpavar2) + ' ' + str(hpavar2)
				wb.write(var)
				wb.write('\n')
				wb.close()
				print("results exported")
	wb = open(csvname1, 'a')	#	LOCAL	FOLDER
	var = str(ecor)
	wb.write(var)
	wb.write('\n')
	wb.close()


def run_batch(k, aggfun=np.mean):
	eco_list0 = np.genfromtxt('csv/ecoregs' + str(k) + '.csv', dtype=int)	#	crear	este	archivo	en	subpas!
	eco_list = np.unique(eco_list0)
	m = len(eco_list)
	for	pm	in	range(0, m):  # 3,m	#	0	without	the	negative	ecoregs!
		ecor = eco_list[pm]
		print(ecor)
		ehabitat(ecor, '', '', k, aggfun)
	print(str(datetime.now()))
	print("BATCH END")


def run_batch_all(aggfun=np.mean):
	ecox_list0 = np.genfromtxt('csv/segm_done.csv', dtype=int)	#	crear	este	archivo	en	subpas!
	ecox_list = np.unique(ecox_list0)
	mx = len(ecox_list)
	for	pmx	in	range(0, mx):  # 3,m	#	0	without	the	negative	ecoregs!
		global park
		park = ecox_list[pmx]
		print(park)
		for k in range(1, 10):
			run_batch(k, aggfun)
	print(str(datetime.now()))
	print("BATCH END")


def run_park(paid, aggfun=np.mean):
	global park
	park = paid
	print(park)
	for k in range(1, 10):
		run_batch(k, aggfun)
	print(str(datetime.now()))
	print("BATCH END")
