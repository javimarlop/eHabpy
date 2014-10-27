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

write.table(num_segms_freq,'num_segms_freq.csv',sep=',',row.names=F)
write.table(data3_ord,'data3_ord.csv',sep=',',row.names=F)

png('freq_segments_150km2.png')
plot(num_segms_freq$Var1,num_segms_freq$Freq)
dev.off()

