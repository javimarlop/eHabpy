# Load libraries

library(raster)
library(rasterVis)
#library(colorspace)
#library(randomcoloR)
library(gridExtra)
library(sp)

hft1<-raster('../../results/10701_2572_1_lp.tif')
#pa<-read_sf('park_segm_2572_7_class.shp')
pa<-readOGR('park_segm_2572_7_class.shp')

rw0<-raster('af_avoided_potential_removed_soil_mass.tif')
rw<-rw0*0.3254

burundi <- getData('GADM', country='Burundi', level=0)
rwanda <- getData('GADM', country='Rwanda', level=0)

### Palette

#Uniquesrw <- cellStats(rw,stat=unique)
#Uniques.max.rw <- max(Uniquesrw,na.rm=T)
#Uniques.min.rw <- min(Uniquesrw, na.rm=T)
#my.at.rw <- c(Uniques.max.rw,20,10,5,2,1,Uniques.min.rw)

#Uniquesbu <- cellStats(bu,stat=unique)
#Uniques.max.bu <- max(Uniquesbu,na.rm=T)
#Uniques.min.bu <- min(Uniquesbu, na.rm=T)
#my.at.bu <- c(Uniques.max.bu,20,10,5,2,1,Uniques.min.bu)

#my.brks<-c(20,10,5,2)#,2,1)
#myColorkey <- list(at=my.brks, 
#                   labels=list(at=my.brks, labels= c('>20',10,5,'<2',''))) # expression(>= 50)

#rw_plot<-levelplot(rw, par.settings=YlOrRdTheme(),at = my.at.rw, colorkey=myColorkey, margin=FALSE, xlab=NULL, ylab=NULL, main='Rwanda')

#bu_plot<-levelplot(log(bu), par.settings=YlOrRdTheme(),at = my.at.bu, colorkey=myColorkey, margin=FALSE, xlab=NULL, ylab=NULL, main='Burundi') 

png('fig7.png',res=200,width=1440,height=900) 
#tiff('fig7.tif',res=200,width=1440,height=900) 

levelplot(hft1, par.settings=YlOrRdTheme(), margin=FALSE, xlab=NULL, ylab=NULL, main='HFT 1 similar landscape patches') + layer(sp.polygons(pa['hclst_m'])) 

dev.off()

# Export
#png('fig7.png',res=200,width=1680,height=1050) 
#grid.arrange(rw_plot,bu_plot, ncol=2)
#dev.off()

#tiff('fig7.tif',res=200,width=1680,height=1050) 
#grid.arrange(rw_plot,bu_plot, ncol=2)
#dev.off()


