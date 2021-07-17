 library(sf)
 #library(basemaps)
 library(raster)
 library(ggplot2)
 library(ggspatial)
  
 pa<-read_sf('kkdu.shp')#'park_segm_2572_7_class.shp')
 pa_3857 <- st_transform(pa, 3857)
 pa_4326 <- st_transform(pa, 4326)
 
# set_defaults(map_service = "esri", map_type = "world_physical_map")

 hft1<-raster('../../results/10701_2572_1.tif')
 hft1_lp<-raster('../../results/10701_2572_1_lp.tif')
 #base<-st_bbox(hft1)

#ggplot(pa_3857['hclst_m']) + geom_sf(aes(fill=as.factor(hclst_m))) + scale_colour_brewer(palette = 'Set1') + coord_sf() + theme_bw() + labs(title='Kakadu')  #+ scale_fill_identity() #  + scale_colour_brewer(palette = 'Set1') + basemap_gglayer(base)

hft<-ggplot(data = pa) + annotation_map_tile('stamenwatercolor') + geom_sf(aes(fill = as.factor(hclst_m))) + annotation_scale() + scale_colour_brewer(palette = 'Set1') + ggtitle(label = "Kakadu", subtitle = "HFTs") + annotation_north_arrow(location = "br", which_north = "true") + guides(fill=guide_legend("HFT"))

#hft1_lp<-
ggplot(data = pa) + layer_spatial(hft1_lp) + annotation_map_tile('stamenwatercolor') + geom_sf(aes(fill = as.factor(hclst_m))) + annotation_scale() + ggtitle(label = "Kakadu", subtitle = "Similar landscape patches to HFT 1") + annotation_north_arrow(location = "br", which_north = "true") + guides(fill=guide_legend("HFT"))


  
