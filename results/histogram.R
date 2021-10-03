# url: https://ggplot2.tidyverse.org/reference/geom_histogram.html

library(ggplot2)
#library(reshape2)

df<-read.table('hri_results.csv',sep=' ',header=T)

ggplot(df, aes(HSR)) + geom_histogram()

ggsave('hsr_hist.png')

##



