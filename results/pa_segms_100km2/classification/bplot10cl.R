bplot10cl<-function(){
for(i in 14:40){
png(paste(names(hri)[i],'_10cl_boxplots.png',sep=''))
boxplot(hri[,i] ~ hclust_dist_euclid_segms_10cl,main=paste(names(hri)[i]))
dev.off()
}}
