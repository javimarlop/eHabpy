library(rgdal)
segments<-readOGR(dsn='.',lay='parks_segmented')
segments@data->data
dim(data)
data[!duplicated(data[c(3,5)]),]->data2
dim(data2)
as.data.frame(table(data2$wdpaid_pa))->data3
dim(data3)

max(data3$Freq)
data3[which.max(data3$Freq),]

data3[order(data3$Freq),]->data3_ord
data.frame(table(data3_ord$Freq))->num_segms_freq

library(reshape2)
datax<-read.table('total_areas.csv',sep=',',header=T)
dcast(datax[,c(3,6)], wdpaid_pa ~ ., sum)->data_wdpaid
dcast(datax[,c(5,6)], segm_id ~ ., sum)->data_segm_id
names(data_segm_id)[2]<-'m2'
names(data_wdpaid)[2]<-'m2'

write.table(data_segm_id,'data_segm_id.csv',sep=',',row.names=F)
names(data3_ord)[1]<-'wdpaid_pa'
merge(data3_ord,data_wdpaid,by='wdpaid_pa',all=T)->data3_ord_areas_m2
data3_ord_areas_m2$m2/1000000->data3_ord_areas_m2[,4]
names(data3_ord_areas_m2)[4]<-'km2'

data3_ord_areas_m2[order(data3_ord_areas_m2$km2),]->data3_ord_areas_km2
index150km2<-data3_ord_areas_km2$km2 >= 150
data3_ord_areas_150km2<-data3_ord_areas_km2[index150km2,]

data3_ord_areas_150km2_sqrt<-data3_ord_areas_150km2
sqrt(data3_ord_areas_150km2$km2)->data3_ord_areas_150km2_sqrt[,5]#/10
(data3_ord_areas_150km2_sqrt$Freq*1000/data3_ord_areas_150km2_sqrt$V5)->data3_ord_areas_150km2_sqrt[,6]
hist(data3_ord_areas_150km2_sqrt$V5)
hist(data3_ord_areas_150km2_sqrt$V6)
data3_ord_areas_150km2_sqrt[order(data3_ord_areas_150km2_sqrt$V6),]->data3_ord_areas_150km2_final_sqrt
names(data3_ord_areas_150km2_final_sqrt)[5:6]<-c('sqrt_area','awhd')
tail(data3_ord_areas_150km2_final_sqrt)

write.table(data3_ord_areas_150km2_final_sqrt,'data3_ord_areas_150km2_final_sqrt.csv',sep=',',row.names=F)

index150shp <- segments$wdpaid_pa %in% data3_ord_areas_150km2$wdpaid_pa

segments_filter<-segments[index150shp,]
writeOGR(segments_filter,dsn='.',layer='parks_segmented_filter',driver="ESRI Shapefile")

write.table(num_segms_freq_filter,'num_segms_freq_filter.csv',sep=',',row.names=F)
#write.table(data3_ord,'data3_ord.csv',sep=',',row.names=F)

png('freq_segments_150km2.png')
plot(num_segms_freq$Var1,num_segms_freq$Freq)
dev.off()

data.frame(table(data3_ord_areas_150km2_final_sqrt$Freq))->num_segms_freq_filter

png('freq_segments_150km2_filter.png')
plot(num_segms_freq_filter$Var1,num_segms_freq_filter$Freq)
dev.off()

