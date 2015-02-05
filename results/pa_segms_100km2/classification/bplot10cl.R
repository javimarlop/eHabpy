bplot10cl<-function(){
for(i in 14:40){
png(paste(names(hri)[i],'_10cl_boxplots.png',sep=''))
boxplot(hri[,i] ~ hclust_dist_euclid_segms_10cl,main=paste(names(hri)[i]))
dev.off()
}}

bplot10cl2<-function(){
for(i in 14:40){
png(paste(names(hri)[i],'_10cl_boxplots_mh.png',sep=''))
boxplot(hri[,i] ~ hclust_dist_euclid_segms_10cl2,main=paste(names(hri)[i]))
dev.off()
}}

bplot10cl3<-function(){
for(i in c(3,14:40)){
png(paste(names(hri)[i],'_10cl_boxplots_mh_all.png',sep=''))
boxplot(hri[,i] ~ hclust_dist_euclid_segms_10cl3,main=paste(names(hri)[i]))
dev.off()
}}
