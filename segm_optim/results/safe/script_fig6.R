library(sf)
library(RColorBrewer)

bc = read_sf("ws.shp")
names(bc)[2]<-'Mean flood regulation service surplus/deficit'

tiff('fig6.tif',width=1280,height=1024,res=200)
#png('fig6.png',width=1280,height=1024,res=200)
plot(bc[,2], key.pos = 1, axes = TRUE,cex=2.5)#, graticule = FALSE)
dev.off()

# library(tmap)
# qtm(bc[,2],fill='rps',fill.palette="-Blues")


