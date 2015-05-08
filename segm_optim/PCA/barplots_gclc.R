require(ggplot2)
#library(RColorBrewer)

#kols<-brewer.pal(8,'Set1')


df0<-read.table('results/gclc_aver_in.csv',sep=' ',header=T)
pas<-unique(df0$wdpaid)
mx<-length(pas)

for (pmx in 1:mx){

	park = pas[pmx]
	print(park)
	dfpa<-df0[df0$wdpaid==park,]

m <- ggplot(dfpa, aes(x=factor(gclc),y=aver,fill=factor(gclc)))

m + geom_bar(stat="identity") + facet_wrap(~cluster)+ labs(x = paste("GCLC for park ",park,sep=''))#+theme_bw()

ggsave(paste('results/gclc_cats_',park,'.png',sep=''))

}
