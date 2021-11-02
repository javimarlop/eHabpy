library(FactoMineR)

df<-read.table('thdi_res_names.csv',sep=',',header=T)

pcadf<-PCA(df[,6:8])

png('CA_PA_indices.png')
CA(df[,6:8])
dev.off()

# plot with axes and PAs
hcl<-hclust(dist(df[,6:8]),"ward.D2")

png('hclust_PA_indices.png')
plot(hcl,hang=-1)
dev.off()





