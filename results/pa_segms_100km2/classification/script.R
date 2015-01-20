# cat ../hri_results.csv |grep -v None > hri_results.csv
hri<-read.table('hri_results.csv',sep=' ',header=T,dec='.')
# Export segments shapefile table to csv (QGIS WIN) as db.csv and import it in R
#db<-read.table('parks_segmented.csv',sep=',',header=T)
# merge with segmentation output shapefile table
#names(db)<-c('cat','label','wdpa_id','wdpaid')
names(hri)[2]<-c('segm_id')
#duplicated(db$wdpaid)->dupl_index
#db0<-db[!dupl_index,]

# normalize values by PA/ECO/variable
#hri0<-merge(hri,db0,by='wdpaid')
#rownames(hri0)<-hri0$wdpaid

library(vegan)
hclust(vegdist(hri[,14:40],"euclidean"),"ward.D2")->hclust_dist_euclid_segms
cutree(hclust_dist_euclid_segms,h=50000)->hclust_dist_euclid_segms_10cl
cutree(hclust_dist_euclid_segms,h=21000)->hclust_dist_euclid_segms_20cl
cutree(hclust_dist_euclid_segms,h=100000)->hclust_dist_euclid_segms_5cl
cutree(hclust_dist_euclid_segms,h=30000)->hclust_dist_euclid_segms_15cl
#cutree(hclust_dist_euclid_segms,h=7000)->hclust_dist_euclid_segms_25cl
#cutree(hclust_dist_euclid_segms,h=6200)->hclust_dist_euclid_segms_30cl
#cutree(hclust_dist_euclid_segms,h=5500)->hclust_dist_euclid_segms_35cl
#cutree(hclust_dist_euclid_segms,h=5000)->hclust_dist_euclid_segms_40cl

#df<-as.data.frame(cbind(hri0$segm_id,hclust_dist_euclid_segms_5cl,hclust_dist_euclid_segms_10cl,hclust_dist_euclid_segms_15cl,hclust_dist_euclid_segms_20cl,hclust_dist_euclid_segms_25cl))
#names(df)<-c('segment','5cl','10cl','15cl','20cl','25cl')

#write.table(df,'habitats_cl.csv',row.names=F)

library(FactoMineR)
PCA(hri[,14:40])->PCA_segms
CA(hri[,14:40])->CA_segms
decorana(hri[,14:40])->DCA_segms

png('CA_segms.png')
plot(CA_segms)
dev.off()

#plot(CA_segms_no_high_epr_min$row$coord[,1],CA_segms_no_high_epr_min$row$coord[,2],typ='n')
#text(CA_segms_no_high_epr_min$col$coord[,1],CA_segms_no_high_epr_min$col$coord[,2],labels=rownames(CA_segms_no_high_epr_min$col$coord),col=2)

#plot(CA_segms$row$coord[,1],CA_segms$row$coord[,2],col=hclust_dist_euclid_segms_5cl)
#text(CA_segms$col$coord[,1],CA_segms$col$coord[,2],labels=rownames(CA_segms$col$coord),col=2)
#identify(CA_segms_no_high_epr$row$coord[,1],CA_segms_no_high_epr$row$coord[,2])

#plot(CA_segms_no_epr$row$coord[,1],CA_segms_no_epr$row$coord[,2])#,col=hclust_dist_euclid_segms_5cl)
#text(CA_segms_no_epr$col$coord[,1],CA_segms_no_epr$col$coord[,2],labels=rownames(CA_segms_no_epr$col$coord),col=2)

metaMDS(vegdist(hri[,14:40],"euclidean"))->mds_segms
#rf proximities
#plot(mds_segms$points,col=hclust_dist_euclid_segms_10cl)
#identify(mds_segms$points)->examples_cats
##[1]   886  3294  4939  6850  8774 11563 12120 12518 13457 13774
#hri0$wdpaid[examples_cats[1]]

table(hclust_dist_euclid_segms_5cl)
#    1     2     3     4     5
# 7835  4324 10624  1849   458

table(hclust_dist_euclid_segms_10cl)
#   1    2    3    4    5    6    7    8    9   10
#2601 3430  894 5234 4913 5711 1539  310  174  284

table(hclust_dist_euclid_segms_15cl)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
#2601 1204 2226  858 2053 3181   36 3571 5711 1342 1077  462  310  174  284

#png('premax_10cl_boxplots.png')
#boxplot(hri$prepamax ~ hclust_dist_euclid_segms_10cl)
#dev.off()

library(vegan)
library(ade4)

png('mds_segms_5cl.png')
plot(mds_segms)
s.class(mds_segms$points,as.factor(hclust_dist_euclid_segms_5cl),col=1:5)
dev.off()

png('mds_segms_10cl.png')
plot(mds_segms)
s.class(mds_segms$points,as.factor(hclust_dist_euclid_segms_10cl),col=1:10)
dev.off()

png('mds_segms_15cl.png')
plot(mds_segms)
s.class(mds_segms$points,as.factor(hclust_dist_euclid_segms_15cl),col=1:15)
dev.off()

