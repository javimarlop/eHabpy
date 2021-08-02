library(sf)
library(ggplot2)

hri<-read.table('results/hri_results.csv',sep=' ', header=T)
ots<-read.table('segm_optim/results/overall_optim_thresholds.csv',sep=' ', header=F)
teow<-read_sf('ecoregions/terrwwf_cor_moll_54009_areas.shp')
meow<-read_sf('ecoregions/meow_cor_moll_54009_areas.shp')

# FOR LOOP BY wdpaid
	for(i in 1:1){ # dim(ots)[1]
		paid<-ots$V1[i]
		thr<-ots$V2[i]*10
		source = paste('segm_optim/results/park_segm_',paid,'_',thr,'_class.shp',sep='')
		pa0<-read_sf(source)

		pateow<-st_intersection(teow,pa0)
		nteow<-length(unique(pateow$ECO_ID))

		pameow<-st_intersection(meow,pa0)
		nmeow<-length(unique(pameow$ECO_CODE))
		
		# export as list for each PA

	}
