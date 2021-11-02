library(ggplot2)
library(reshape2)

df<-read.table('thdi_res.csv',sep=' ',header=T)

###

df2<-melt(df[,-c(2,3,5)],id='wdpaid')

ggplot(df2) + geom_col(aes(factor(wdpaid),value)) + facet_grid(rows=vars(variable), scales='free') + theme_bw() + xlab('Protected area')

ggsave('PA_indices.png')

ggplot(df2) + geom_col(aes(factor(wdpaid),value)) + facet_wrap(vars(variable), scales='free') + theme_bw() + xlab('Protected area') + theme(axis.text.x = element_text(angle = 45, hjust=1)) + theme(axis.title.y = element_blank(), axis.title.x = element_blank()) 

ggsave('PA_indices2.png')

###

df3<-melt(df[,c(1:3)],id='wdpaid')

ggplot(df3) + geom_col(aes(factor(wdpaid),value,fill=variable)) + scale_fill_brewer(palette='Set2')

ggsave('PA_ecoregs.png')

###

dfv2<-read.table('thdi_res_names.csv',sep=',',header=T)

#
ggplot(dfv2, aes(TotEco, THR)) + geom_point() + geom_smooth(method='lm', se = FALSE)
ggsave('totecoregs_thr.png')

summary(lm(dfv2$THR ~ dfv2$TotEco))

#
ggplot(dfv2, aes(TotEco, TotSimEco)) + geom_point()

summary(lm(dfv2$TotSimEco ~ dfv2$TotEco))
summary(lm(dfv2$TotSimEco ~ dfv2$TotEco + I(dfv2$TotEco^2)))

ggplot(dfv2, aes(x=TotEco, y=TotSimEco)) + geom_point()+stat_smooth(se=F, method='lm', formula=y~poly(x,2))
ggsave('totecoregs_totsimeco.png')

