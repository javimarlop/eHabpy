
ex <- parse(file.path('/Users/javier/Downloads/Claude/eHabpy', 'segm_optim', 'moranvar_plots.R'))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name('<-')) && as.character(e[[2]]) %in% c('std_stats', 'class_tree')) eval(e)
for (stat in c('mean', 'median')) {
	hri <- read.csv(file.path('/Users/javier/Downloads/Claude/eHabpy/feature_comparison/review_tmp', paste0('hri_', stat, '.csv')))
	ct <- class_tree(hri, 'ward_area')
	out <- sapply(2:20, function(k) cutree(ct$tree, k=k))
	write.csv(out, file.path('/Users/javier/Downloads/Claude/eHabpy/feature_comparison/review_tmp', paste0('r_', stat, '.csv')), row.names=FALSE)
	write.csv(ct$z, file.path('/Users/javier/Downloads/Claude/eHabpy/feature_comparison/review_tmp', paste0('rz_', stat, '.csv')), row.names=FALSE)
}
