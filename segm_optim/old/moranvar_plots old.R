# Modernized for R 4.x: rgdal (retired 2023) replaced by sf; ggplot2:::rescale01
# (internal, removed) replaced by a local helper.
#library(pastecs)
library(vegan)
#library(FactoMineR)
library(ade4)
library(ggplot2)
library(reshape2)
library(sf)
library(RColorBrewer)

# ggplot2 no longer exports the internal rescale01(); define it locally.
# Same formula as the original ggplot2:::rescale01. The zero-range guard is an
# intentional robustness improvement: for a constant column ggplot2 returned
# NaN (0/0), which would later crash CreateRadialPlot(); here it maps to 0.
rescale01 <- function(x) {
	rng <- range(x, na.rm = TRUE)
	if (diff(rng) == 0) return(rep(0, length(x)))
	(x - rng[1]) / (rng[2] - rng[1])
}

kols<-brewer.pal(8,'Set1')

is.linear.radar <- function(coord) TRUE

df0<-read.table('csv/segm_done.csv',sep=' ',header=F)
ecox_list <- unique(df0[,1])
mx <- length(ecox_list)
optim<-2 # new option: values 0 (fixed threshold of 0.5), 1 (default optimum), 2 (one step higher number of segments than the optimum threshold)
simil<-2 # new option: values 1 (default value), 2 (correspoding to the number of segments with an upper limit)
upper<-2 # new option: values 1 (default value), 2 (increased upper limit)
maxnclas<-8 # get number of variables -1 from file

for (pmx in 1:mx){
	park = ecox_list[pmx]
	print(park)
	name <- paste('csv/',park,'_movar_results.csv',sep='')
	df0<-read.table(name,sep=' ',header=F)
	df1<-df0[!is.nan(df0$V2),]
	fthr<-0
	if(dim(df1)[1]==0){fthr<-1; print("Using a fixed threshold")}

	namem <- paste('csv/',park,'_moran_mean.csv',sep='')
	df0m<-read.table(namem,sep=' ',header=F)
	namev <- paste('csv/',park,'_var_mean.csv',sep='')
	df0v<-read.table(namev,sep=' ',header=F)

	name2 <- paste('csv/',park,'_movar_thresholds.csv',sep='')
	df02<-read.table(name2,sep=' ',header=F)
	if(upper==1){upplim<-ceiling(log10(df02[1,3]))} # default value
	if(upper==2){upplim<-ceiling(log(df02[1,3],6))} # new value
	print(paste('upper limit is ',upplim,sep=''))

	x1<-df0m$V2
	x2<-df0v$V2

	above<-x1>x2
	# Points always intersect when above=TRUE, then FALSE or reverse
	intersect.points<-which(diff(above)!=0)
	# Find the slopes for each line segment.
	x1.slopes<-x1[intersect.points+1]-x1[intersect.points]
	x2.slopes<-x2[intersect.points+1]-x2[intersect.points]
	# Find the intersection for each segment.
	x.points<-intersect.points + ((x2[intersect.points] - x1[intersect.points]) / (x1.slopes-x2.slopes))
	y.points<-x1[intersect.points] + (x1.slopes*(x.points-intersect.points))

	newtot<-df0$V2#/2
	uppv <- 1# max(c(df0m$V2,df0v$V2,newtot))
	png(paste('results/',park,'_moran_var_mean.png',sep=''))
	plot(df0m$V1,df0m$V2,ylim=c(0,uppv), col=1,typ='o',ylab='Moran I. and Variance',xlab='Similarity threshold (x 0.1)',main=park) # black
	lines(df0v$V1,df0v$V2,col=3,typ='o') # green
	lines(df0$V1,newtot,col=4,typ='o') # blue
	points(x.points,y.points,col='red')
	legend("top", leg=c('M.I','Var','Aver'), col=c(1,3,4), lt = 1)
	dev.off()

	res<-5 # (standard version) 5
	res0<-as.data.frame(cbind(5,NA))
	res2<-res0#[1,]

	if(length(x.points)>0){ #>=
	if(optim==1){res<-round(x.points[1])} # default value #length(x.points)
	if(optim==2){res<-round(x.points[1]) - 1}  #new value #length(x.points)
	res0<-as.data.frame(cbind(x.points,y.points))
	res2<-res0[1,]
	}

	k<-paste('0.',res,sep='')
	res2[1,2]<-as.numeric(k)
	res2[1,1]<-park
	write.table(res0,paste('results/',park,'_optim_thresholds.csv',sep=''),sep=' ',col.names=F,row.names=F)
	write.table(res2,'results/overall_optim_thresholds.csv',append=T,sep=' ',col.names=F,row.names=F)



	namef <- paste('csv/','park_',park,'_hri_results',res,'.csv',sep='')
	print(namef)
	hri<-read.table(namef,sep=' ',header=T)
	skaled <- as.data.frame(lapply(hri[,3:20], rescale01))
	try(dmh<-vegdist(skaled,"euclidean",na.rm=T))#"")
	hclust(dmh,"ward.D2")->hclust_mh #  # ward.D2
	cophval<-0
	try(chcl<-cophenetic(hclust_mh))
	try(cophval<-cor(chcl,dmh))
	try(metaMDS(dmh)->mds_mh)
	if(simil==1){q25<-quantile(hclust_mh$height)[2]} # default value
	if(simil==2){q25<-min(hclust_mh$height) - 0.1} # new value
	print(q25)
	cutree(hclust_mh,h=q25)->hclust_mean #
	print(hclust_mean)
	ncl<-length(unique(hclust_mean))
	print(paste('Number of HFTs is ',ncl,sep=''))

	if(ncl==1){
	as.numeric(factor(hri$segm_id))->hclust_mean
	ncl<-length(unique(hclust_mean))
	}

	#upplim<-8#3
	if(ncl>upplim){
	cutree(hclust_mh,k=upplim)->hclust_mean
	ncl<-length(unique(hclust_mean))
	}

	if(ncl>maxnclas){
	cutree(hclust_mh,k=maxnclas)->hclust_mean
	ncl<-length(unique(hclust_mean))
	}

	if(ncl>1){

	png(paste('results/hclust_',park,'_',res,'_segms_mean.png',sep=''))
	try(plot(hclust_mh,hang=-1,main=park,sub=cophval));try(rect.hclust(hclust_mh,k=ncl))
	dev.off()

	png(paste('results/NMDS_',park,'_',res,'_segms_mean.png',sep=''))
	try(plot(mds_mh))
	try(s.class(mds_mh$points,as.factor(hclust_mean),col=1:length(unique(hclust_mean))))
	dev.off()

	hrin<-cbind(hri[,1:11],hclust_mean)
	names(hrin)[3:11]<-c("Woody","Aridity","Precip","BioTemp","Slope","NDWI","NDVI_MAX","NDVI_MIN","Grassland")
	hri3<-melt(hrin[,3:12],'hclust_mean')
	hri4<-dcast(hri3,hclust_mean ~ variable,mean)

	scaled <- as.data.frame(lapply(hri4[,2:10], rescale01))
	scaled$model <- hri4[,1]#rownames(hri[,3:11])

	scaled2<-cbind(scaled[,10],scaled[,1:9])
	names(scaled2)[1]<-'group'
	source('CreateRadialPlot.R')
	CreateRadialPlot(scaled2,plot.extent.x = 1.5)
	rpn=paste('results/radarplot_',park,'_',res,'_segms_mean.png',sep='')
	ggsave(filename=rpn)

	# rgdal::readOGR -> sf::st_read
	segm_pa<-st_read('shp', layer=paste('park_segm_',park,'_',res,'_diss',sep=''), quiet=TRUE)

	merge(segm_pa,hrin[,c(2,12)],by='segm_id')->segm_pa_class

	scaled0<-scaled2
	kl<-dim(scaled0)[2]+1
	scaled0[,kl]<-rep(park,dim(scaled0)[1])
	names(scaled0)[1]<-'hclust_mean'
	names(scaled0)[kl]<-'wdpaid'
	merge(segm_pa_class,scaled0,by='hclust_mean')->segm_pa_class2

	# rgdal::writeOGR -> sf::st_write. NOTE: the ESRI Shapefile driver truncates
	# field names to 10 characters, so 'hclust_mean' becomes 'hclust_mea'. The
	# downstream scripts (subpas_loop_segm_optim.py, lcgc.py) look up column
	# 'hclst_m'; verify/align this name for your data if needed.
	st_write(segm_pa_class2,
			 file.path('results', paste('park_segm_',park,'_',res,'_class.shp',sep='')),
			 delete_layer=TRUE, quiet=TRUE)

	rpn2=paste('results/map_',park,'_',res,'_segms_hclust.png',sep='')
	png(paste(rpn2))
	plot(st_geometry(segm_pa_class),col=kols[segm_pa_class$hclust_mean],main=park) # col=segm_pa_class$hclust_mean
	legend("bottomright", leg=unique(segm_pa_class$hclust_mean), col=unique(kols[segm_pa_class$hclust_mean]), pch = 19, title = "Legend")
	dev.off()

	# MERGE BACK THE CATEGORIES IN THE DISS SHAPEFILE
}}
