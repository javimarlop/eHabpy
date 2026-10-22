# Modernized for R 4.x & Refactored for Dynamic Variables and PA Outline Overlays
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

# Función para buscar y leer dinámicamente cualquier formato vectorial en la carpeta 'pa'
# (con varios, el que tiene el nombre del área de estudio, como en setup.py)
find_and_read_pa_vector <- function(pa_dir = "pa", park_id = NULL, area = NULL) {
	if (!dir.exists(pa_dir)) return(NULL)
	
	files <- list.files(pa_dir, full.names = TRUE)
	if (length(files) == 0) return(NULL)
	
	# Filtro de extensiones vectoriales soportadas por GDAL / sf
	vec_exts <- "\\.(shp|gpkg|geojson|json|kml|gml|sqlite|tab)$"
	valid_files <- files[grep(vec_exts, files, ignore.case = TRUE)]
	
	if (length(valid_files) == 0) return(NULL)
	
	selected_file <- NULL
	if (!is.null(area)) {
		# setup.py cambia los caracteres que no son letras ASCII, números o _ por _ (y añade v_ si empieza por número)
		key <- function(x) gsub('[^A-Za-z0-9]', '', x)
		stems <- key(sub('\\.[^.]*$', '', basename(valid_files)))
		match_file <- valid_files[stems %in% key(c(area, sub('^v_', '', area)))]
		if (length(match_file) > 0) {
			selected_file <- match_file[1]
		}
	}
	if (is.null(selected_file) && !is.null(park_id)) {
		match_file <- valid_files[grep(as.character(park_id), valid_files)]
		if (length(match_file) > 0) {
			selected_file <- match_file[1]
		}
	}
	
	if (is.null(selected_file)) {
		selected_file <- valid_files[1]
	}
	
	pa_sf <- tryCatch({
		st_read(selected_file, quiet = TRUE)
	}, error = function(e) {
		warning(paste("No se pudo cargar la capa vectorial de área de estudio:", selected_file, "-", e$message))
		return(NULL)
	})
	
	return(pa_sf)
}

kols <- brewer.pal(8, 'Set1')

is.linear.radar <- function(coord) TRUE

df0 <- read.table('csv/segm_done.csv', sep=' ', header=F)
ecox_list <- unique(df0[,1])
mx <- length(ecox_list)

# Lee el valor (sin comentario ni comillas) de una opción sencilla (NOMBRE = valor) de config.py
read_config <- function(name) {
	if (!file.exists('config.py')) return(NA)
	hit <- grep(paste0('^\\s*', name, '\\s*='), readLines('config.py', warn=FALSE), value=TRUE)
	if (length(hit) == 0) return(NA)
	gsub("^\\s*['\"]?|['\"]?\\s*$", '', sub('#.*$', '', sub('^[^=]*=', '', hit[length(hit)])))
}

read_config_num <- function(name, default) {
	val <- suppressWarnings(as.numeric(read_config(name)))
	if (is.na(val)) default else val
}

optim <- 2 # 0: fijo, 1: óptimo, 2: un paso superior
simil <- read_config_num('SIMIL', 2) # ver SIMIL en config.py
upper <- 2 # 1: default, 2: límite superior incrementado

if (!simil %in% 1:3) stop("SIMIL en config.py debe ser 1, 2 o 3")

# Resultados en results/<área de estudio>: la de los csv (anotada por segmentation_pca_par.py) o STUDY_AREA de config.py
area <- if (file.exists('csv/study_area.csv')) readLines('csv/study_area.csv', n=1, warn=FALSE) else read_config('STUDY_AREA')
if (is.na(area) || area == '') stop("No se encuentra el área de estudio: falta STUDY_AREA en config.py")
outroot <- file.path('results', area)

# Cortes del árbol: simil = 3 genera una clasificación por corte (25, 50 y 75 % de la altura máxima), cada una en su subcarpeta
if (simil == 3) {
	cuts <- c(p25=0.25, p50=0.50, p75=0.75)
} else {
	cuts <- c(results=NA)
}

# Metadatos de la ejecución: opciones de config.py (sin comentarios) y de este script
cfg <- if (file.exists('config.py')) readLines('config.py', warn=FALSE) else character(0)
cfg <- sub('\\s+#.*$', '', cfg)
cfg <- cfg[!grepl('^\\s*(#|$)', cfg)]
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
	'eHabitat+ processing metadata: options used to produce the files in this folder',
	paste('Date:', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
	paste('Study area:', area),
	paste('IDs processed (csv/segm_done.csv):', paste(ecox_list, collapse=', ')),
	'',
	'## config.py',
	cfg,
	'',
	'## moranvar_plots.R',
	paste('optim =', optim, '# 0: fixed threshold 0.5, 1: first Moran/variance crossing, 2: first crossing - 1'),
	paste('upper =', upper, '# max. HFTs (SIMIL 1 and 2): 1: log10(n. pixels), 2: log6(n. pixels)')
), file.path(outroot, 'processing_metadata.txt'))

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

	# moranvar.py appends to these files on every run: keep one row per threshold (the latest)
	df0 <- df0[!duplicated(df0$V1, fromLast=TRUE), ]
	df0m <- df0m[!duplicated(df0m$V1, fromLast=TRUE), ]
	df0v <- df0v[!duplicated(df0v$V1, fromLast=TRUE), ]
	df0 <- df0[order(df0$V1), ]
	df0m <- df0m[order(df0m$V1), ]
	df0v <- df0v[order(df0v$V1), ]

	name2 <- paste('csv/', park, '_movar_thresholds.csv', sep='')
	df02 <- read.table(name2, sep=' ', header=F)
	if (upper == 1) { upplim <- ceiling(log10(df02[1,3])) }
	if (upper == 2) { upplim <- ceiling(log(df02[1,3], 6)) }
	print(paste('Upper limit is ', upplim, sep=''))

	x1 <- df0m$V2
	x2 <- df0v$V2
	thr <- df0m$V1 # thresholds actually evaluated (some may be missing)

	# Crossing points, interpolated on the threshold values (not on row positions)
	above <- x1 > x2
	intersect.points <- which(diff(above) != 0)
	x1.slopes <- x1[intersect.points + 1] - x1[intersect.points]
	x2.slopes <- x2[intersect.points + 1] - x2[intersect.points]
	t.points <- (x2[intersect.points] - x1[intersect.points]) / (x1.slopes - x2.slopes)
	x.points <- thr[intersect.points] + t.points * (thr[intersect.points + 1] - thr[intersect.points])
	y.points <- x1[intersect.points] + (x1.slopes * t.points)

	newtot <- df0$V2
	uppv <- 1
	
	png(file.path(outroot, paste(park, '_moran_var_mean.png', sep='')))
	plot(df0m$V1, df0m$V2, ylim=c(0, uppv), col=1, typ='o', ylab='Moran I. and Variance', xlab='Similarity threshold (x 0.1)', main=park)
	lines(df0v$V1, df0v$V2, col=3, typ='o')
	lines(df0$V1, newtot, col=4, typ='o')
	points(x.points, y.points, col='red')
	legend("top", leg=c('M.I', 'Var', 'Aver'), col=c(1,3,4), lt = 1)
	dev.off()

	res <- 5
	res0 <- as.data.frame(cbind(5, NA))

	if (length(x.points) > 0) {
		if (optim == 1) { res <- round(x.points[1]) }
		if (optim == 2) { res <- round(x.points[1]) - 1 }
		# Use the closest evaluated threshold (optim = 2 can give 0 when the curves cross near 1)
		res <- thr[which.min(abs(thr - res))]
		res0 <- as.data.frame(cbind(x.points, y.points))
	}

	k <- paste('0.', res, sep='')
	write.table(res0, file.path(outroot, paste(park, '_optim_thresholds.csv', sep='')), sep=' ', col.names=F, row.names=F)

	# One row per park: replace the previous row of this park instead of appending a duplicate
	overall <- file.path(outroot, 'overall_optim_thresholds.csv')
	res2 <- data.frame(V1 = park, V2 = as.numeric(k))
	if (file.exists(overall) && file.size(overall) > 0) {
		prev <- read.table(overall, sep=' ', header=F)
		res2 <- rbind(prev[prev$V1 != park, ], res2)
	}
	write.table(res2, overall, sep=' ', col.names=F, row.names=F)

	namef <- paste('csv/park_', park, '_hri_results', res, '.csv', sep='')
	print(paste("Reading HRI file:", namef))
	if (!file.exists(namef)) { print(paste("WARNING: HRI file not found, park skipped:", namef)); next }
	
	hri <- read.table(namef, sep=' ', header=T)

	num_vars <- (ncol(hri) - 3) / 3
	maxnclas <- num_vars - 1

	skaled <- as.data.frame(lapply(hri[, 3:(2 + 2 * num_vars)], rescale01))
	
	try(dmh <- vegdist(skaled, "euclidean", na.rm=T))
	hclust(dmh, "ward.D2") -> hclust_mh
	cophval <- 0
	try(chcl <- cophenetic(hclust_mh))
	try(cophval <- cor(chcl, dmh))
	try(metaMDS(dmh) -> mds_mh)
	
	for (cn in names(cuts)) {
		outdir <- if (simil == 3) file.path(outroot, cn) else outroot
		tag <- if (simil == 3) paste0('_', cn) else '' # el corte también va en el nombre de cada fichero
		dir.create(outdir, showWarnings = FALSE)

		if (simil == 1) { hcut <- quantile(hclust_mh$height)[2] }
		if (simil == 2) { hcut <- min(hclust_mh$height) - 0.1 }
		if (simil == 3) { hcut <- cuts[[cn]] * max(hclust_mh$height) } # fracción de la altura máxima del árbol

		cutree(hclust_mh, h=hcut) -> hclust_mean
		ncl <- length(unique(hclust_mean))
		print(paste('Number of HFTs is ', ncl, if (simil == 3) paste0(' (', cn, ')'), sep=''))

		if (ncl == 1) {
			as.numeric(factor(hri$segm_id)) -> hclust_mean
			ncl <- length(unique(hclust_mean))
		}

		# simil = 3 no impone límite superior al número de HFTs
		if (simil != 3 && ncl > upplim) {
			cutree(hclust_mh, k=upplim) -> hclust_mean
			ncl <- length(unique(hclust_mean))
		}

		if (simil != 3 && ncl > maxnclas) {
			cutree(hclust_mh, k=maxnclas) -> hclust_mean
			ncl <- length(unique(hclust_mean))
		}

		# Colores del mapa: Set1 hasta 8 HFTs, interpolados si hay más
		kolsn <- if (ncl <= length(kols)) kols else colorRampPalette(kols)(ncl)

		if (ncl > 1) {
			png(file.path(outdir, paste('hclust_', park, '_', res, tag, '_segms_mean.png', sep='')))
			try(plot(hclust_mh, hang=-1, main=park, sub=cophval)); try(rect.hclust(hclust_mh, k=ncl))
			dev.off()

			png(file.path(outdir, paste('NMDS_', park, '_', res, tag, '_segms_mean.png', sep='')))
			try(plot(mds_mh))
			try(s.class(mds_mh$points, as.factor(hclust_mean), col=1:length(unique(hclust_mean))))
			dev.off()

			hrin <- cbind(hri[, 1:(2 + num_vars)], hclust_mean)
			raw_names <- names(hri)[3:(2 + num_vars)]
			clean_names <- sub("pamean$", "", raw_names)
			names(hrin)[3:(2 + num_vars)] <- clean_names

			hri3 <- melt(hrin[, 3:(3 + num_vars)], id.vars='hclust_mean')
			hri4 <- dcast(hri3, hclust_mean ~ variable, mean)

			scaled <- as.data.frame(lapply(hri4[, 2:(num_vars + 1), drop=FALSE], rescale01))
			scaled$model <- hri4[, 1]

			scaled2 <- cbind(group = scaled$model, scaled[, 1:num_vars, drop=FALSE])
			
			source('CreateRadialPlot.R')
			# Mismos colores que el mapa (kolsn indexado por nº de HFT)
			CreateRadialPlot(scaled2, plot.extent.x = 1.5, group.colours = setNames(kolsn[scaled2$group], scaled2$group))
			rpn <- file.path(outdir, paste('radarplot_', park, '_', res, tag, '_segms_mean.png', sep=''))
			ggsave(filename=rpn)

			# Integración en Shapefile y Mapa final
			shp_layer <- paste('park_segm_', park, '_', res, '_diss', sep='')
			if (file.exists(file.path('shp', paste0(shp_layer, '.shp')))) {
				segm_pa <- st_read('shp', layer=shp_layer, quiet=TRUE)

				merge(segm_pa, hrin[, c(2, 2 + num_vars + 1)], by='segm_id') -> segm_pa_class

				scaled0 <- scaled2
				kl <- dim(scaled0)[2] + 1
				scaled0[, kl] <- rep(park, dim(scaled0)[1])
				names(scaled0)[1] <- 'hclust_mean'
				names(scaled0)[kl] <- 'wdpaid'
				merge(segm_pa_class, scaled0, by='hclust_mean') -> segm_pa_class2

				st_write(segm_pa_class2,
						 file.path(outdir, paste('park_segm_', park, '_', res, tag, '_class.shp', sep='')),
						 delete_layer=TRUE, quiet=TRUE)

				# --- DIBUJO DEL MAPA DE HFTs CON EL CONTORNO DE ÁREA DE ESTUDIO (ggplot2) ---
				rpn2 <- file.path(outdir, paste('map_', park, '_', res, tag, '_segms_hclust.png', sep=''))
				hfts <- sort(unique(segm_pa_class$hclust_mean))
				segm_pa_class$HFT <- factor(segm_pa_class$hclust_mean, levels=hfts)

				# Bordes de los segmentos en gris fino para que destaque el contorno negro del área de estudio
				hftmap <- ggplot() +
					geom_sf(data=segm_pa_class, aes(fill=HFT), colour='grey50', linewidth=0.05) +
					scale_fill_manual(values=setNames(kolsn[hfts], hfts)) +
					labs(title=area, subtitle=paste0('ID ', park, ' | threshold ', k, if (simil == 3) paste0(' | cut ', cn), ' | ', ncl, ' HFTs')) +
					theme_bw(base_size=14) +
					theme(axis.title=element_blank())

				# Cargar y dibujar la capa del área de estudio en negro
				pa_outline <- find_and_read_pa_vector("pa", park, area)
				if (!is.null(pa_outline)) {
					# Reproyectar al CRS del mapa si difieren
					if (!is.na(st_crs(segm_pa_class)) && !is.na(st_crs(pa_outline)) && st_crs(segm_pa_class) != st_crs(pa_outline)) {
						pa_outline <- st_transform(pa_outline, st_crs(segm_pa_class))
					}
					hftmap <- hftmap + geom_sf(data=pa_outline, fill=NA, colour='black', linewidth=0.8)
				}
				hftmap <- hftmap + coord_sf(expand=FALSE)

				# 12 pulgadas de ancho a 300 dpi; el alto sigue la proporción del área (más título)
				bb <- st_bbox(segm_pa_class)
				asp <- as.numeric((bb['ymax'] - bb['ymin']) / (bb['xmax'] - bb['xmin']))
				ggsave(rpn2, hftmap, width=12, height=min(max(10 * asp + 1.5, 5), 16), dpi=300)
			}
		}
	}
}
