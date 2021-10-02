library(sf)
library(ggplot2)

hri<-read.table('results/hri_results.csv',sep=' ', header=T) # input
ots<-read.table('segm_optim/results/overall_optim_thresholds.csv',sep=' ', header=F) # control file
teow<-read_sf('ecoregions/terrwwf_cor_moll_54009_areas.shp')
meow<-read_sf('ecoregions/meow_cor_moll_54009_areas.shp')

a<-list()
	for(i in 1:dim(ots)[1]){ # 1 # FOR LOOP BY wdpaid

		paid<-ots$V1[i]
		thr<-ots$V2[i]*10
		source = paste('segm_optim/results/park_segm_',paid,'_',thr,'_class.shp',sep='')
		pa0<-read_sf(source)

		pateow<-st_intersection(teow,pa0)
		nteow<-length(unique(pateow$ECO_ID))

		pameow<-st_intersection(meow,pa0)
		nmeow<-length(unique(pameow$ECO_CODE))
		
		# export as list for each PA
a[i]<-paid

a[[i]][2]<-nteow
a[[i]][3]<-nmeow

neco<-nteow+nmeow

a[[i]][4]<-neco

df<-hri[hri$wdpaid==paid,]

nhft<-dim(df)[1]

a[[i]][5]<-nhft

area<-sum(df$NrPixHFT)

perc<-sapply(df$NrPixHFT, function(x) x/area)

sihftx<-sapply(perc, function(x) x*log(x))

sihft<-(-1)*sum(sihftx)

a[[i]][6]<-sihft

ecohab<-neco*sihft

a[[i]][7]<-ecohab

msim<-median(df$HSR) # TotNrPixAllSimLcpPatEqLargArea (old name)

a[[i]][8]<-msim

thdi<-ecohab/msim

a[[i]][9]<-thdi

	} # end of for loop

#print(a)

res <- data.frame(matrix(unlist(a), nrow=length(a), byrow=TRUE))
names(res)<-c('wdpaid','TerrEco','MarEco','TotEco','NrHFT','SIH','THD','THI','THDI')

write.table(res,'results/thdi_res.csv',sep=' ',row.names=F)
