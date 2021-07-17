library(rgdal)

kkdu<-readOGR('park_segm_2572_7_class.shp')

library(raster)

hft1<-raster('../../results/10701_2572_1.tif')
hft1_lp<-raster('../../results/10701_2572_1_lp.tif')

plot(hft1)
plot(kkdu, add=True)

library(basemaps)
library(sf)

pa<-read_sf('park_segm_2572_7_class.shp')
library(RColorBrewer)

#options(sf_max.plot=1)

plot(pa['hclst_m'])

#plot(st_geometry(pa), col = sf.colors(4, categorical = TRUE), border = 'grey', axes = TRUE)

# https://cengel.github.io/R-spatial/mapping.html
library(ggplot2)

ggplot(pa) + geom_sf(aes(fill=hclst_m))

# https://github.com/Chrisjb/basemapR

#ggplot(pa) + geom_sf(aes(fill=hclst_m)) + basemap_gglayer(bb) + coord_sf() + scale_fill_identity()

#pa_wgs84 <- st_transform(pa, 4326)
pa_3857 <- st_transform(pa, 3857)

ggplot(pa_3857) + basemap_gglayer(base) + geom_sf(aes(fill=hclst_m)) + coord_sf() + scale_fill_identity()

###

bb<-st_bbox(hft1)

get_maptypes()

set_defaults(map_service = "esri", map_type = "world_physical_map")

bm<-basemap_geotif(bb) # save basemap as geotiff

base<-raster('basemap_kkdu_moll.tif')

plot(base,legend=F)
plot(kkdu,add=T)
plot(hft1_lp,add=T,legend=F)
#plot(hft1,add=T)


