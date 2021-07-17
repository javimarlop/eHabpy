library(ggplot2)

hri<-read.table('hri_results_boxplots.csv',sep=' ',header=T)

ggplot(hri[hri$wdpaid=='2572',], aes(HFT, value)) + geom_boxplot() + facet_wrap(~ variable, scales='free') + labs(title = "Kakadu")
ggsave('2572_boxplots.png')

ggplot(hri[hri$wdpaid=='61612',], aes(HFT, value)) + geom_boxplot() + facet_wrap(~ variable, scales='free') + labs(title = "Canaima")
ggsave('61612_boxplots.png')

ggplot(hri[hri$wdpaid=='2017',], aes(HFT, value)) + geom_boxplot() + facet_wrap(~ variable, scales='free') + labs(title = "Virunga")
ggsave('2017_boxplots.png')

ggplot(hri[hri$wdpaid=='19297',], aes(HFT, value)) + geom_boxplot() + facet_wrap(~ variable, scales='free') + labs(title = "Udzungwa")
ggsave('19297_boxplots.png')

ggplot(hri[hri$wdpaid=='555577555',], aes(HFT, value)) + geom_boxplot() + facet_wrap(~ variable, scales='free') + labs(title = "Okavango")
ggsave('555577555_boxplots.png')

ggplot(hri[hri$wdpaid=='555538721',], aes(HFT, value)) + geom_boxplot() + facet_wrap(~ variable, scales='free') + labs(title = "Sierra Nevada")
ggsave('555538721_boxplots.png')


