library(ggplot2)
library(reshape)

hri0<-read.table('hri_results.csv',sep=' ',header=T)

hri1<-hri0[,c(2,3,33:41)]

names(hri1)[2:11]<-c('HFT','Woody','Aridity','Precipitation','Biotemperature','Slope','NDWI','NDVImax','NDVImin','Grassland')

hri<-melt(hri1, id=c('wdpaid','HFT'))

ggplot(hri[hri$wdpaid=='2572',], aes(HFT, value)) + geom_col() + facet_wrap(~ variable, scales='free') + labs(title = "Kakadu")
ggsave('2572_barplots.png')

ggplot(hri[hri$wdpaid=='61612',], aes(HFT, value)) + geom_col() + facet_wrap(~ variable, scales='free') + labs(title = "Canaima")
ggsave('61612_barplots.png')

ggplot(hri[hri$wdpaid=='2017',], aes(HFT, value)) + geom_col() + facet_wrap(~ variable, scales='free') + labs(title = "Virunga")
ggsave('2017_barplots.png')

ggplot(hri[hri$wdpaid=='19297',], aes(HFT, value)) + geom_col() + facet_wrap(~ variable, scales='free') + labs(title = "Udzungwa")
ggsave('19297_barplots.png')

ggplot(hri[hri$wdpaid=='555577555',], aes(HFT, value)) + geom_col() + facet_wrap(~ variable, scales='free') + labs(title = "Okavango")
ggsave('555577555_barplots.png')

ggplot(hri[hri$wdpaid=='555538721',], aes(HFT, value)) + geom_col() + facet_wrap(~ variable, scales='free') + labs(title = "Sierra Nevada")
ggsave('555538721_barplots.png')


