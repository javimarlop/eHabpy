library(ggplot2)
library(reshape2)

df<-read.table('thdi_res.csv',sep=' ',header=T)

###

df2<-melt(df[,-c(2,3,5)],id='wdpaid')

ggplot(df2) + geom_col(aes(factor(wdpaid),value)) + facet_grid(rows=vars(variable), scales='free') + theme_bw() + xlab('Protected area')

ggsave('PA_indices.png')

###

df3<-melt(df[,c(1:3)],id='wdpaid')

ggplot(df3) + geom_col(aes(factor(wdpaid),value,fill=variable)) + scale_fill_brewer(palette='Set2')

ggsave('PA_ecoregs.png')


