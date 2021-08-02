 library(sf)
 library(dplyr)
 #library(basemaps)
 library(raster)
 library(ggplot2)
 library(ggspatial)
  
 pa0<-read_sf('park_segm_2572_6_class.shp')#'park_segm_2572_7_class.shp')
 pa <- pa0 %>% group_by(hclst_m) %>% summarize()
 #pa_3857 <- st_transform(pa, 3857)
 #pa_4326 <- st_transform(pa, 4326)
 
# set_defaults(map_service = "esri", map_type = "world_physical_map")

 #hft1<-raster('../../results/10701_2572_1.tif')
 hft1_lp<-raster('../../results/10701_2572_1_lp.tif')
 hft2_lp<-raster('../../results/10701_2572_2_lp.tif')
 hft3_lp<-raster('../../results/10701_2572_3_lp.tif')
 hft4_lp<-raster('../../results/10701_2572_4_lp.tif')
 hft5_lp<-raster('../../results/10701_2572_5_lp.tif')
 hft6_lp<-raster('../../results/10701_2572_6_lp.tif')  
 #base<-st_bbox(hft1)

 hfts_lp<-hft1_lp+hft2_lp+hft3_lp+hft4_lp+hft5_lp+hft6_lp

 #hft1_lp_recl<-reclassify(hft1_lp,c(1,max(getValues(hft1_lp)),1))
 #hft1_lp_final <- clamp(hft1_lp_recl, lower=1, useValues=FALSE)

 reclass_raster<-function(raster){
 hft1_lp_recl<-reclassify(raster,c(1,max(getValues(raster)),1))
 hft1_lp_final <- clamp(hft1_lp_recl, lower=1, useValues=FALSE)
 return(hft1_lp_final)
 }

hfts<-reclass_raster(hfts_lp)

hft1<-reclass_raster(hft1_lp)
hft2<-reclass_raster(hft2_lp)
hft3<-reclass_raster(hft3_lp)
hft4<-reclass_raster(hft4_lp)
hft5<-reclass_raster(hft5_lp)
hft6<-reclass_raster(hft6_lp)

#ggplot(pa_3857['hclst_m']) + geom_sf(aes(fill=as.factor(hclst_m))) + scale_colour_brewer(palette = 'Set1') + coord_sf() + theme_bw() + labs(title='Kakadu')  #+ scale_fill_identity() #  + scale_colour_brewer(palette = 'Set1') + basemap_gglayer(base)

# HFTs
hft<-ggplot(data = pa) + annotation_map_tile('stamenwatercolor', zoom=8) + geom_sf(aes(fill = as.factor(hclst_m))) + annotation_scale() + scale_fill_brewer(palette = 'Set1') + ggtitle(label = "Kakadu", subtitle = "HFTs") + annotation_north_arrow(location = "br", which_north = "true") + guides(fill=guide_legend("HFT"))

ggsave('hft_kakadu.png',hft)

# Similarity using single hft
#ggplot(data = pa) + annotation_map_tile('stamenwatercolor') + geom_sf(fill='transparent') + annotation_scale() + ggtitle(label = "Kakadu", subtitle = "Similar landscape patches to all HFTs") + annotation_north_arrow(location = "br", which_north = "true") + guides(fill=guide_legend("HFT")) + layer_spatial(hft1) + layer_spatial(hft2) + layer_spatial(hft3) + layer_spatial(hft4) + layer_spatial(hft5) + layer_spatial(hft6) + scale_fill_continuous(na.value = "transparent")

# Similarity using sum of hfts
hfts<-ggplot(data = pa) + annotation_map_tile('stamenwatercolor', zoom=8) + geom_sf(fill='transparent') + annotation_scale() + ggtitle(label = "Kakadu", subtitle = "Similar landscape patches to all HFTs") + annotation_north_arrow(location = "br", which_north = "true") + layer_spatial(trim(hfts)) + scale_fill_continuous(na.value = "transparent") + theme(legend.position = "none") # + guides(fill=guide_legend("HFT"))

ggsave('hfts_kakadu.png',hfts) #.tiff

plot_hfts<-function(name=''){
hfts<-ggplot(data = pa) + annotation_map_tile('stamenwatercolor', zoom=8) + geom_sf(fill='transparent') + annotation_scale() + ggtitle(label = "Kakadu", subtitle = "Similar landscape patches") + annotation_north_arrow(location = "br", which_north = "true") + layer_spatial(trim(get(name))) + scale_fill_continuous(na.value = "transparent") + theme(legend.position = "none") # + guides(fill=guide_legend("HFT"))
ggsave(paste(name,'_simil_kakadu.png',sep=''),hfts) #.tiff
return(hfts)
}

phft1<-plot_hfts('hft1')
phft2<-plot_hfts('hft2')
phft3<-plot_hfts('hft3')
phft4<-plot_hfts('hft4')
phft5<-plot_hfts('hft5')
phft6<-plot_hfts('hft6')

png('Kakadu_simil.png',width=1680,height=1050) 
grid.arrange(phft1, phft2, phft3, phft4, phft5,  phft6, ncol=2)
dev.off()

#library(maptiles)

#nc <- st_transform(nc_raw, "EPSG:3857")
## dowload tiles and compose raster (SpatRaster)
#nc_osm <- get_tiles(nc, crop = TRUE)
## display map
#plot_tiles(nc_osm)
