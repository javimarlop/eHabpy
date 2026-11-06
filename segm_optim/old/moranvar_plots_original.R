# Modernized for R 4.x & Refactored for Dynamic Variables
library(vegan)
library(ade4)
library(ggplot2)
library(reshape2)
library(sf)
library(RColorBrewer)

rescale01 <- function(x) {
	rng <- range(x, na.rm = TRUE)
	if (diff(rng) == 0) return(rep(0, length(x)))
	(x - rng[1]) / (rng[2] - rng[1])
}

kols <- brewer.pal(8, 'Set1')

is.linear.radar <- function(coord) TRUE

df0 <- read.table('csv/segm_done.csv', sep=' ', header=F)
ecox_list <- unique(df0[,1])
mx <- length(ecox_list)

optim <- 2 # 0: fijo, 1: óptimo, 2: un paso superior
simil <- 2 # 1: default, 2: límite inferior
upper <- 2 # 1: default, 2: límite superior incrementado

for (pmx in 1:mx) {
	park = ecox_list[pmx]
	print(paste("Procesando park ID:", park))
	
	name <- paste('csv/', park, '_movar_results.csv', sep='')
	if (!file.exists(name)) next
	
	df0 <- read.table(name, sep=' ', header=F)
	df1 <- df0[!is.nan(df0$V2),]
	fthr <- 0
	if (dim(df1)[1] == 0) { fthr <- 1; print("Using a fixed threshold") }

	namem <- paste('csv/', park, '_moran_mean.csv', sep='')
	df0m <- read.table(namem, sep=' ', header=F)
	namev <- paste('csv/', park, '_var_mean.csv', sep='')
	df0v <- read.table(namev, sep=' ', header=F)

	name2 <- paste('csv/', park, '_movar_thresholds.csv', sep='')
	df02 <- read.table(name2, sep=' ', header=F)
	if (upper == 1) { upplim <- ceiling(log10(df02[1,3])) }
	if (upper == 2) { upplim <- ceiling(log(df02[1,3], 6)) }
	print(paste('Upper limit is ', upplim, sep=''))

	x1 <- df0m$V2
	x2 <- df0v$V2

	above <- x1 > x2
	intersect.points <- which(diff(above) != 0)
	x1.slopes <- x1[intersect.points + 1] - x1[intersect.points]
	x2.slopes <- x2[intersect.points + 1] - x2[intersect.points]
	x.points <- intersect.points + ((x2[intersect.points] - x1[intersect.points]) / (x1.slopes - x2.slopes))
	y.points <- x1[intersect.points] + (x1.slopes * (x.points - intersect.points))

	newtot <- df0$V2
	uppv <- 1
	
	dir.create('results', showWarnings = FALSE)
	png(paste('results/', park, '_moran_var_mean.png', sep=''))
	plot(df0m$V1, df0m$V2, ylim=c(0, uppv), col=1, typ='o', ylab='Moran I. and Variance', xlab='Similarity threshold (x 0.1)', main=park)
	lines(df0v$V1, df0v$V2, col=3, typ='o')
	lines(df0$V1, newtot, col=4, typ='o')
	points(x.points, y.points, col='red')
	legend("top", leg=c('M.I', 'Var', 'Aver'), col=c(1,3,4), lt = 1)
	dev.off()

	res <- 5
	res0 <- as.data.frame(cbind(5, NA))
	res2 <- res0

	if (length(x.points) > 0) {
		if (optim == 1) { res <- round(x.points[1]) }
		if (optim == 2) { res <- round(x.points[1]) - 1 }
		res0 <- as.data.frame(cbind(x.points, y.points))
		res2 <- res0[1,]
	}

	k <- paste('0.', res, sep='')
	res2[1,2] <- as.numeric(k)
	res2[1,1] <- park
	write.table(res0, paste('results/', park, '_optim_thresholds.csv', sep=''), sep=' ', col.names=F, row.names=F)
	write.table(res2, 'results/overall_optim_thresholds.csv', append=T, sep=' ', col.names=F, row.names=F)

	namef <- paste('csv/park_', park, '_hri_results', res, '.csv', sep='')
	print(paste("Reading HRI file:", namef))
	if (!file.exists(namef)) next
	
	hri <- read.table(namef, sep=' ', header=T)

	# --- CÁLCULO DINÁMICO DE VARIABLES ---
	# Estructura de columnas en hri: ecoregion, segm_id, [N medias], [N varianzas], sumpamask, [N var2]
	# Total columnas = 3*N + 3 => N = (ncol - 3) / 3
	num_vars <- (ncol(hri) - 3) / 3
	maxnclas <- num_vars - 1  # Número máximo de clases dinámico

	# Escalar las medias y varianzas de las variables (cols 3 a 2 + 2*N)
	skaled <- as.data.frame(lapply(hri[, 3:(2 + 2 * num_vars)], rescale01))
	
	try(dmh <- vegdist(skaled, "euclidean", na.rm=T))
	hclust(dmh, "ward.D2") -> hclust_mh
	cophval <- 0
	try(chcl <- cophenetic(hclust_mh))
	try(cophval <- cor(chcl, dmh))
	try(metaMDS(dmh) -> mds_mh)
	
	if (simil == 1) { q25 <- quantile(hclust_mh$height)[2] }
	if (simil == 2) { q25 <- min(hclust_mh$height) - 0.1 }
	
	cutree(hclust_mh, h=q25) -> hclust_mean
	ncl <- length(unique(hclust_mean))
	print(paste('Number of HFTs is ', ncl, sep=''))

	if (ncl == 1) {
		as.numeric(factor(hri$segm_id)) -> hclust_mean
		ncl <- length(unique(hclust_mean))
	}

	if (ncl > upplim) {
		cutree(hclust_mh, k=upplim) -> hclust_mean
		ncl <- length(unique(hclust_mean))
	}

	if (ncl > maxnclas) {
		cutree(hclust_mh, k=maxnclas) -> hclust_mean
		ncl <- length(unique(hclust_mean))
	}

	if (ncl > 1) {
		png(paste('results/hclust_', park, '_', res, '_segms_mean.png', sep=''))
		try(plot(hclust_mh, hang=-1, main=park, sub=cophval)); try(rect.hclust(hclust_mh, k=ncl))
		dev.off()

		png(paste('results/NMDS_', park, '_', res, '_segms_mean.png', sep=''))
		try(plot(mds_mh))
		try(s.class(mds_mh$points, as.factor(hclust_mean), col=1:length(unique(hclust_mean))))
		dev.off()

		# Construir tabla hrin con las N variables de media y hclust_mean al final
		hrin <- cbind(hri[, 1:(2 + num_vars)], hclust_mean)
		
		# Extraer y limpiar nombres de variables desde las cabeceras del CSV
		raw_names <- names(hri)[3:(2 + num_vars)]
		clean_names <- sub("pamean$", "", raw_names)
		names(hrin)[3:(2 + num_vars)] <- clean_names

		# Reshape para el Radarplot
		hri3 <- melt(hrin[, 3:(3 + num_vars)], id.vars='hclust_mean')
		hri4 <- dcast(hri3, hclust_mean ~ variable, mean)

		scaled <- as.data.frame(lapply(hri4[, 2:(num_vars + 1), drop=FALSE], rescale01))
		scaled$model <- hri4[, 1]

		scaled2 <- cbind(group = scaled$model, scaled[, 1:num_vars, drop=FALSE])
		
		source('CreateRadialPlot.R')
		CreateRadialPlot(scaled2, plot.extent.x = 1.5)
		rpn <- paste('results/radarplot_', park, '_', res, '_segms_mean.png', sep='')
		ggsave(filename=rpn)

		# Integración en Shapefile
		shp_layer <- paste('park_segm_', park, '_', res, '_diss', sep='')
		if (file.exists(file.path('shp', paste0(shp_layer, '.shp')))) {
			segm_pa <- st_read('shp', layer=shp_layer, quiet=TRUE)

			# Merge segm_id y la columna de clase (posición 2 + num_vars + 1)
			merge(segm_pa, hrin[, c(2, 2 + num_vars + 1)], by='segm_id') -> segm_pa_class

			scaled0 <- scaled2
			kl <- dim(scaled0)[2] + 1
			scaled0[, kl] <- rep(park, dim(scaled0)[1])
			names(scaled0)[1] <- 'hclust_mean'
			names(scaled0)[kl] <- 'wdpaid'
			merge(segm_pa_class, scaled0, by='hclust_mean') -> segm_pa_class2

			st_write(segm_pa_class2,
					 file.path('results', paste('park_segm_', park, '_', res, '_class.shp', sep='')),
					 delete_layer=TRUE, quiet=TRUE)

			rpn2 <- paste('results/map_', park, '_', res, '_segms_hclust.png', sep='')
			png(rpn2)
			plot(st_geometry(segm_pa_class), col=kols[segm_pa_class$hclust_mean], main=park)
			legend("bottomright", leg=unique(segm_pa_class$hclust_mean), col=unique(kols[segm_pa_class$hclust_mean]), pch = 19, title = "Legend")
			dev.off()
		}
	}
}
