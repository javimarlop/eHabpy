from __future__ import division
import pysal
import numpy as np
import os
import sys
import csv

#os.system('touch optim_thresholds.csv')
ecox_list0 = np.genfromtxt('csv/segm_done.csv',dtype='int')	#	crear	este	archivo	en	subpas!
ecox_list = np.unique(ecox_list0) # ['19297'] ['555542456']
mx = len(ecox_list)
for	pmx	in	range(0,mx):
 park = ecox_list[pmx]
 print 'park id is:',park
 #csvname4 = str(park)+'_movar_segmsum0.csv'
 #csvname40 = str(park)+'_movar_segmsum02.csv'
 csvname = 'csv/'+str(park)+'_movar_segmsum.csv'
 csvnamem = 'csv/'+str(park)+'_moran_segm.csv'
 csvnamev = 'csv/'+str(park)+'_var_segm.csv'
 csvnamem2 = 'csv/'+str(park)+'_moran_mean.csv'
 csvnamev2 = 'csv/'+str(park)+'_var_mean.csv'
 #csvname3 = str(park)+'_movar_segmsumf.csv'
 #csvname30 = str(park)+'_movar_segmsumf0.csv'
 csvname2 = 'csv/'+str(park)+'_movar_thresholds.csv'
 for i in np.arange(2,11): # loop by variable
  print 'variable is:',i
  mors = []
  varis = []
  thr = []
  for k in range(1,10): # loop by segmentation threshold
	wpn = 0
	layrname = 'park_segm_'+str(park)+'_'+str(k)
	shpname = 'shp/park_segm_'+str(park)+'_'+str(k)+'.shp'
	shpname2 = 'shp/park_segm_'+str(park)+'_'+str(k)+'_diss.shp'
	hriname = 'csv/park_'+str(park)+'_hri_results'+str(k)+'.csv'
	sumareas = np.genfromtxt(hriname,delimiter=' ',skip_header=1,usecols=(20))
	if os.path.isfile(shpname2) == False:
	 os.system('ogr2ogr '+shpname2+' '+shpname+' -dialect sqlite -sql "SELECT ST_Union(geometry), segm_id FROM '+layrname+' GROUP BY segm_id"')
	w=pysal.rook_from_shapefile(shpname2)
	print 'w.n is:',w.n
	wpn = w.n
	if wpn >1: #!=
		print 'segmentation threshold is:',k
		sareas = sum(sumareas) # area of the PA (sum of the segments)
		thr.append(k)
		y30 = np.genfromtxt(hriname,delimiter=' ',skip_header=1,usecols=(i))
		mi = pysal.Moran(y30, w)#, two_tailed=False)
		mm = abs(mi.I)
		print 'M.I. is:',mm
		i2=i+19
		y330 = np.genfromtxt(hriname,delimiter=' ',skip_header=1,usecols=(i2))
		wv = sum(y330)/sareas
		print 'Sum of the variance is:',wv
 		mors.append(mm)
 		varis.append(wv)
		wb = open(csvname2,'a')
		outxt = str(k)+' '+str(w.n)+' '+str(sareas)
		wb.write(outxt)
		wb.write('\n')
		wb.close()
  print 'list of MIs:',mors
  print 'list of variances:',varis
  m3 = (mors-min(mors))/(max(mors)-min(mors))
  #v3 = (varis-min(varis))/(max(varis)-min(varis))
  v3 = (max(varis)-varis)/(max(varis)-min(varis))
  #if max(varis) == 0 and min(varis)==0: v3 = 0
  tot = (m3 + v3)/2

  wb = open(csvnamem,'a')
  for f in np.arange(0,len(m3)):
   wb.write('{},'.format(str(m3[f])))
  wb.write('\n')
  wb.close()

  wb = open(csvnamev,'a')
  for f in np.arange(0,len(m3)):
   wb.write('{},'.format(str(v3[f])))
  wb.write('\n')
  wb.close()

  wb = open(csvname,'a')
  #wb.write('{}'.format(str(tot[0:len(m3)])))
  for f in np.arange(0,len(tot)):
   wb.write('{},'.format(str(tot[f])))
  #wb.write('{0},{1},{2},{3},{4},{5},{6},{7},{8}\n'.format(str(tot[0]),str(tot[1]),str(tot[2]),str(tot[3]),str(tot[4]),str(tot[5]),str(tot[6]),str(tot[7]),str(tot[8])))#str(tot))
  wb.write('\n')
  wb.close()
 #os.system('cat '+csvname+' | sed "s/]\+//g" > '+csvname3) # \t
 #os.system('cat '+csvname4+' | sed "s/[\+//g" > '+csvname3) #\t
 #os.system('cat '+csvname30+' | sed "s/ \+//g" > '+csvname3) #\t
 #os.system('rm '+csvname+' '+csvname4)
 csvname5 = 'csv/'+str(park)+'_movar_results.csv'
 thrs = np.genfromtxt(csvname2,delimiter=' ',skip_header=0,usecols=(0)) # check delimiter!
 segs = np.unique(thrs)
 h = -1
 f = len(segs)#+1
 for d in np.arange(0,f):
  h = h +1
  thr = np.genfromtxt(csvname,delimiter=',',skip_header=0,usecols=(d))
  print 'threshold:',thr
  tm = np.nanmean(thr)
  wb = open(csvname5,'a')
  var = str(segs[h])+' '+str(tm)
  wb.write(var)
  wb.write('\n')
  wb.close()

  m3d = np.genfromtxt(csvnamem,delimiter=',',skip_header=0,usecols=(d))
  print 'morans:',m3d
  m3m = np.nanmean(m3d)
  wb = open(csvnamem2,'a')
  varm = str(segs[h])+' '+str(m3m)
  wb.write(varm)
  wb.write('\n')
  wb.close()

  v3d = np.genfromtxt(csvnamev,delimiter=',',skip_header=0,usecols=(d))
  print 'variances',v3d
  v3m = np.nanmean(v3d)
  wb = open(csvnamev2,'a')
  varv = str(segs[h])+' '+str(v3m)
  wb.write(varv)
  wb.write('\n')
  wb.close()

os.system('Rscript moranvar_plots.R')
print "BATCH END"
