#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 / PySAL 2.x (libpysal + esda) / GeoPandas (2024)
#### License: CC BY-SA 3.0
#### Control files: 'csv/segm_done.csv'
#### Notes: At the end it calls the 'Rscript moranvar_plots.R' script.
####        PySAL 1.x (monolithic ``pysal``) was split; Moran's I now comes from
####        ``esda`` and spatial weights from ``libpysal``. The per-segment
####        polygon dissolve is done with GeoPandas instead of the old
####        fiona + shapely + itertools.groupby block.
#### Outputs: moran_... and var... csv files in csv folder; All outputs in results folder

import numpy as np
import os
import sys
import csv
import itertools

import libpysal
import esda
import geopandas as gpd

ecox_list0 = np.genfromtxt('csv/segm_done.csv', dtype=int)	#	crear	este	archivo	en	subpas!
ecox_list = np.unique(ecox_list0)  # ['19297'] ['555542456']
mx = len(ecox_list)
for	pmx	in	range(0, mx):
	park = ecox_list[pmx]
	print('park id is:', park)
	csvname = 'csv/' + str(park) + '_movar_segmsum.csv'
	csvnamem = 'csv/' + str(park) + '_moran_segm.csv'
	csvnamev = 'csv/' + str(park) + '_var_segm.csv'
	csvnamem2 = 'csv/' + str(park) + '_moran_mean.csv'
	csvnamev2 = 'csv/' + str(park) + '_var_mean.csv'
	csvname2 = 'csv/' + str(park) + '_movar_thresholds.csv'
	for i in np.arange(2, 11):  # loop by variable
		print('variable is:', i)
		mors = []
		varis = []
		thr = []
		for k in range(1, 10):  # loop by segmentation threshold
			wpn = 0
			layrname = 'park_segm_' + str(park) + '_' + str(k)
			shpname = 'shp/park_segm_' + str(park) + '_' + str(k) + '.shp'
			shpname2 = 'shp/park_segm_' + str(park) + '_' + str(k) + '_diss.shp'
			hriname = 'csv/park_' + str(park) + '_hri_results' + str(k) + '.csv'
			sumareas = np.genfromtxt(hriname, delimiter=' ', skip_header=1, usecols=(20))
			if os.path.isfile(shpname2) == False:
				# Dissolve polygons by segm_id (unary union per group, keep first
				# attributes) -- GeoPandas replacement for the old fiona/itertools code.
				gdf = gpd.read_file(shpname)
				dissolved = gdf.dissolve(by='segm_id', as_index=False, aggfunc='first')
				dissolved.to_file(shpname2)

			w = libpysal.weights.Rook.from_shapefile(shpname2)
			print('w.n is:', w.n)
			wpn = w.n
			if wpn > 1:  # !=
				print('segmentation threshold is:', k)
				sareas = sum(sumareas)  # area of the PA (sum of the segments)
				thr.append(k)
				y30 = np.genfromtxt(hriname, delimiter=' ', skip_header=1, usecols=(i))
				mi = esda.Moran(y30, w)  # , two_tailed=False)
				mm = abs(mi.I)
				print('M.I. is:', mm)
				i2 = i + 19
				y330 = np.genfromtxt(hriname, delimiter=' ', skip_header=1, usecols=(i2))
				wv = sum(y330) / sareas
				print('Sum of the variance is:', wv)
				mors.append(mm)
				varis.append(wv)
				wb = open(csvname2, 'a')
				outxt = str(k) + ' ' + str(w.n) + ' ' + str(sareas)
				wb.write(outxt)
				wb.write('\n')
				wb.close()
		print('list of MIs:', mors)
		print('list of variances:', varis)
		mors = np.asarray(mors, dtype=float)
		varis = np.asarray(varis, dtype=float)
		m3 = (mors - min(mors)) / (max(mors) - min(mors))
		v3 = (max(varis) - varis) / (max(varis) - min(varis))
		tot = (m3 + v3) / 2

		wb = open(csvnamem, 'a')
		for f in np.arange(0, len(m3)):
			wb.write('{},'.format(str(m3[f])))
		wb.write('\n')
		wb.close()

		wb = open(csvnamev, 'a')
		for f in np.arange(0, len(m3)):
			wb.write('{},'.format(str(v3[f])))
		wb.write('\n')
		wb.close()

		wb = open(csvname, 'a')
		for f in np.arange(0, len(tot)):
			wb.write('{},'.format(str(tot[f])))
		wb.write('\n')
		wb.close()
	csvname5 = 'csv/' + str(park) + '_movar_results.csv'
	thrs = np.genfromtxt(csvname2, delimiter=' ', skip_header=0, usecols=(0))  # check delimiter!
	segs = np.unique(thrs)
	h = -1
	f = len(segs)  # +1
	for d in np.arange(0, f):
		h = h + 1
		thr = np.genfromtxt(csvname, delimiter=',', skip_header=0, usecols=(d))
		print('threshold:', thr)
		tm = np.nanmean(thr)
		wb = open(csvname5, 'a')
		var = str(segs[h]) + ' ' + str(tm)
		wb.write(var)
		wb.write('\n')
		wb.close()

		m3d = np.genfromtxt(csvnamem, delimiter=',', skip_header=0, usecols=(d))
		print('morans:', m3d)
		m3m = np.nanmean(m3d)
		wb = open(csvnamem2, 'a')
		varm = str(segs[h]) + ' ' + str(m3m)
		wb.write(varm)
		wb.write('\n')
		wb.close()

		v3d = np.genfromtxt(csvnamev, delimiter=',', skip_header=0, usecols=(d))
		print('variances', v3d)
		v3m = np.nanmean(v3d)
		wb = open(csvnamev2, 'a')
		varv = str(segs[h]) + ' ' + str(v3m)
		wb.write(varv)
		wb.write('\n')
		wb.close()

os.system('Rscript moranvar_plots.R')
print("BATCH END")
