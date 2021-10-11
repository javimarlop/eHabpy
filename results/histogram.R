# url: https://ggplot2.tidyverse.org/reference/geom_histogram.html

library(ggplot2)
#library(reshape2)

df<-read.table('hri_results.csv',sep=' ',header=T)

ggplot(df, aes(HSR)) + geom_histogram()
# ggplot(df, aes(HSR)) + geom_histogram() + facet_wrap(vars(wdpaid),scales='free')

ggsave('hsr_hist.png')

##



