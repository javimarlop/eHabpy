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
rownames(hri0)<-hri0$wdpaid

library(vegan)
hclust(vegdist(hri0[,12:38],"euclidean"),"ward.D2")->hclust_dist_euclid_segms
cutree(hclust_dist_euclid_segms,h=18000)->hclust_dist_euclid_segms_10cl
cutree(hclust_dist_euclid_segms,h=5000)->hclust_dist_euclid_segms_40cl
cutree(hclust_dist_euclid_segms,h=9100)->hclust_dist_euclid_segms_20cl
cutree(hclust_dist_euclid_segms,h=6200)->hclust_dist_euclid_segms_30cl

write.table(as.data.frame(cbind(hri0$wdpaid,hclust_dist_euclid_segms_10cl,hclust_dist_euclid_segms_20cl)),'habitats_10_20cl.csv',row.names=F)

library(FactoMineR)
PCA(hri0[,12:38])->PCA_segms
CA(hri0[,12:38])->CA_segms

plot(CA_segms$row$coord[,1],CA_segms$row$coord[,2],col=hclust_dist_euclid_segms_10cl)
identify(CA_segms$row$coord[,1],CA_segms$row$coord[,2])

metaMDS(vdist_euclid_segms)->mds_segms
#rf proximities
plot(mds_segms$points,col=hclust_dist_euclid_segms_10cl)
identify(mds_segms$points)->examples_cats
#[1]   886  3294  4939  6850  8774 11563 12120 12518 13457 13774
hri0$wdpaid[examples_cats[1]]


table(hclust_dist_euclid_segms_10cl)
#hclust_dist_euclid_segms_10cl
#   1    2    3    4    5    6    7    8    9   10
#1780 2068 2861  597 3417 1760  342 1130  108   41

boxplot(hri0$prepamax ~ hclust_dist_euclid_segms_10cl)

library(vegan)
library(ade4)
png('mds_segms_10cl.png')
plot(mds_segms)
s.class(mds_segms$points,as.factor(hclust_dist_euclid_segms_10cl),col=1:10)
dev.off()

png('mds_segms_20cl.png')
plot(mds_segms)
s.class(mds_segms$points,as.factor(hclust_dist_euclid_segms_20cl),col=1:20)
dev.off()

