# exclude all rows with NAs
hri<-read.table('hri_results_no_nas.csv',sep=' ',header=T,dec='.')
# Export segments shapefile table to csv (QGIS WIN) as db.csv and import it in R
db<-read.table('db.csv',sep=',',header=T)
# merge with segmentation output shapefile table
names(db)<-c('cat','label','wdpa_id','wdpaid')
duplicated(db$wdpaid)->dupl_index
db0<-db[!dupl_index,]

# normalize values by PA/ECO/variable
hri0<-merge(hri,db0,by='wdpaid')

library(vegan)
vegdist(hri0[,12:38],"euclidean")->vdist_euclid_segms
hclust(vdist_euclid_segms,"ward.D2")->hclust_dist_euclid_segms
library(FactoMineR)
PCA(hri0[,12:38])->PCA_segms
CA(hri0[,12:38])->CA_segms
metaMDS(vdist_euclid_segms)->mds_segms
#rf proximities

