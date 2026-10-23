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

# Etiquetas del radar: sin las partes (separadas por '_') que comparten al principio y al
# final todos los nombres de variable, y en líneas de 'width' caracteres como máximo
short_labels <- function(x, width = 18) {
	parts <- strsplit(x, '_', fixed = TRUE)
	n <- min(lengths(parts))
	pre <- 0
	suf <- 0
	if (length(x) > 1) {
		while (pre < n - 1 && length(unique(sapply(parts, function(p) p[pre + 1]))) == 1) pre <- pre + 1
		while (suf < n - 1 - pre && length(unique(sapply(parts, function(p) p[length(p) - suf]))) == 1) suf <- suf + 1
	}
	sapply(parts, function(p) paste(strwrap(paste(p[(pre + 1):(length(p) - suf)], collapse = ' '), width = width), collapse = '\n'))
}

# Colores jerárquicos para cortes anidados del mismo árbol (SIMIL = 3): cada HFT del corte con menos
# clases conserva su color (base_cols) y los HFTs de los cortes más finos que contiene son tonos de ese
# color, de oscuro a claro en el orden del dendrograma; un HFT intermedio toma el tono medio de los
# HFTs finos que contiene. cls_list: clasificaciones (una por corte) de la más gruesa a la más fina.
hier_colours <- function(tree, cls_list, base_cols) {
	coarse <- cls_list[[1]]
	fine <- cls_list[[length(cls_list)]]
	ford <- unique(fine[tree$order])
	pos <- numeric(max(fine))
	for (g in unique(coarse)) {
		fg <- ford[ford %in% fine[coarse == g]]
		pos[fg] <- (seq_along(fg) - 0.5) / length(fg)
	}
	mixcol <- function(a, b, w) rgb(t((1 - w) * col2rgb(a) + w * col2rgb(b)), maxColorValue=255)
	lapply(cls_list, function(cl) {
		cols <- character(max(cl))
		for (h in unique(cl)) {
			g <- coarse[cl == h][1]
			ramp <- colorRamp(c(mixcol(base_cols[g], 'black', 0.5), base_cols[g], mixcol(base_cols[g], 'white', 0.65)), space='Lab')
			cols[h] <- rgb(ramp(mean(pos[unique(fine[cl == h])])), maxColorValue=255)
		}
		cols
	})
}

# Archivo (en ../inVars) y NoData de cada variable, de ENV_VARS en config.py
read_env_vars <- function() {
	lines <- readLines('config.py', warn = FALSE)
	lines <- lines[!grepl('^\\s*#', lines)]
	m <- regmatches(lines, regexec("'([^']+)'\\s*:\\s*\\{\\s*'file'\\s*:\\s*'([^']+)'\\s*,\\s*'nodata'\\s*:\\s*([^,}]+)", lines))
	m <- m[lengths(m) == 4]
	data.frame(var = sapply(m, `[`, 2), file = sapply(m, `[`, 3),
	           nodata = suppressWarnings(as.numeric(trimws(sapply(m, `[`, 4)))))
}

# Estadísticos de un boxplot (como geom_boxplot): cuartiles y bigotes hasta 1.5 x IQR
box_stats <- function(x) {
	x <- x[!is.na(x)]
	if (length(x) == 0) return(c(ymin = NA, lower = NA, middle = NA, upper = NA, ymax = NA))
	q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
	iqr <- q[3] - q[1]
	c(ymin = min(x[x >= q[1] - 1.5 * iqr]), lower = q[1], middle = q[2], upper = q[3], ymax = max(x[x <= q[3] + 1.5 * iqr]))
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

# Registro de los pasos (hft_state.py, run_state.json): para si un paso anterior no está al día con
# config.py, palist.csv o los datos, y anota cuándo empieza y termina este
run_state <- function(action) {
	suppressWarnings(system2(Sys.getenv('EHAB_PYTHON', 'python3'), c('hft_state.py', action, 'plots')))
}
st <- run_state('check')
if (st == 3) stop('a step before moranvar_plots.R is out of date (see the message above)', call. = FALSE)
if (st != 0) warning('could not check the steps run before with hft_state.py (is python3 in the PATH?)')
run_state('start')

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
# Los resultados anteriores de esta área se mueven a results/<área>_previous (sustituyendo los que hubiera),
# para no mezclar ficheros de ejecuciones con otras opciones (p. ej. otro umbral en el nombre)
if (length(list.files(outroot, recursive = TRUE)) > 0) {
	prev <- paste0(outroot, '_previous')
	unlink(prev, recursive = TRUE)
	if (file.rename(outroot, prev)) print(paste('Previous results moved to', prev))
}
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
	# metaMDS parte de configuraciones aleatorias: con la misma semilla sale igual en cada ejecución
	mds_mh <- NULL
	set.seed(1)
	try(mds_mh <- metaMDS(dmh, trace=FALSE))

	# Valores de las variables en cada píxel (../inVars), para los boxplots por HFT
	var_names <- sub("pamean$", "", names(hri)[3:(2 + num_vars)])
	ev <- read_env_vars()
	ev <- ev[match(var_names, ev$var), ]
	pixvals <- NULL
	if (anyNA(ev$var) || !all(file.exists(file.path('..', 'inVars', ev$file)))) {
		print("WARNING: variable files of config.py (ENV_VARS) not found in ../inVars, boxplots skipped")
	} else {
		rs <- lapply(file.path('..', 'inVars', ev$file), terra::rast, lyrs = 1)
		pixgrid <- rs[[1]]
		# Todas en la malla de la primera variable
		pixvals <- sapply(rs, function(r) {
			if (!terra::compareGeom(r, pixgrid, stopOnError = FALSE)) r <- terra::resample(r, pixgrid, method = 'near')
			terra::values(r, mat = FALSE)
		})
		for (i in seq_len(nrow(ev))) {
			if (!is.na(ev$nodata[i])) pixvals[which(pixvals[, i] == ev$nodata[i]), i] <- NA
		}
		pixarea <- terra::values(terra::cellSize(pixgrid, unit = 'ha', mask = FALSE), mat = FALSE)
	}

	# SIMIL = 3: colores jerárquicos entre los cortes (ver hier_colours), solo para cortes con más de un HFT
	hier_cols <- list()
	if (simil == 3) {
		cls_list <- lapply(cuts, function(f) cutree(hclust_mh, h=f * max(hclust_mh$height)))
		nk <- sapply(cls_list, function(cl) length(unique(cl)))
		cls_list <- cls_list[nk > 1][order(nk[nk > 1])]
		if (length(cls_list) > 0) {
			nb <- length(unique(cls_list[[1]]))
			hier_cols <- hier_colours(hclust_mh, cls_list, if (nb <= length(kols)) kols else colorRampPalette(kols)(nb))
		}
	}

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

		# Colores del mapa: Set1 hasta 8 HFTs, interpolados si hay más; con SIMIL = 3, jerárquicos
		kolsn <- if (ncl <= length(kols)) kols else colorRampPalette(kols)(ncl)
		if (!is.null(hier_cols[[cn]])) kolsn <- hier_cols[[cn]]
		# Orden de los HFTs en leyendas y ejes: con colores jerárquicos el del dendrograma (los tonos de un
		# mismo color quedan juntos), si no el numérico
		hft_order <- if (!is.null(hier_cols[[cn]])) unique(hclust_mean[hclust_mh$order]) else sort(unique(hclust_mean))

		if (ncl > 1) {
			# Título de las figuras: área de estudio, umbral, corte del árbol y nº de HFTs
			plot_subtitle <- paste0('ID ', park, ' | threshold ', k, if (simil == 3) paste0(' | cut ', cn), ' | ', ncl, ' HFTs')

			png(file.path(outdir, paste('hclust_', park, '_', res, tag, '_segms_mean.png', sep='')))
			try(plot(hclust_mh, hang=-1, main=park, sub=cophval)); try(rect.hclust(hclust_mh, k=ncl))
			dev.off()

			# --- NMDS DE LOS SEGMENTOS (las mismas distancias que la clasificación) ---
			# Cada segmento coloreado por su HFT y con tamaño según su superficie, unido al centroide del
			# HFT; elipses (68 %) con 10 HFTs o menos; variables ajustadas como vectores con envfit
			if (!is.null(mds_mh)) {
				nm <- data.frame(x=mds_mh$points[, 1], y=mds_mh$points[, 2], HFT=hclust_mean)
				nm_hfts <- hft_order[hft_order %in% nm$HFT]
				nm$HFT <- factor(nm$HFT, levels=nm_hfts)
				cell_ha <- if (!is.null(pixvals)) median(pixarea, na.rm=TRUE) else read_config_num('RESOLUTION', 1)^2 / 1e4
				nm$ha <- hri$sumpamask * cell_ha
				nm_ha <- tapply(nm$ha, nm$HFT, sum)
				nm_labels <- paste0(nm_hfts, ': ', format(round(nm_ha), big.mark=',', trim=TRUE), ' ha (', round(100 * nm_ha / sum(nm_ha), 1), ' %)')
				nm_cols <- setNames(kolsn[nm_hfts], nm_hfts)
				cent <- aggregate(cbind(cx=x, cy=y) ~ HFT, nm, mean)
				nm <- merge(nm, cent, by='HFT')
				ext <- max(abs(c(nm$x, nm$y)))

				# Variables (medias de los segmentos) ajustadas a la ordenación: dirección y r2; las no
				# significativas (p >= 0.05, 999 permutaciones) en gris claro
				fit <- NULL
				set.seed(1)
				try(fit <- envfit(mds_mh, hri[, 3:(2 + num_vars), drop=FALSE], permutations=999, na.rm=TRUE))
				if (!is.null(fit)) {
					fv <- data.frame(var=paste0(short_labels(var_names), '\nr2 = ', sprintf('%.2f', fit$vectors$r)),
					                 x=fit$vectors$arrows[, 1] * sqrt(fit$vectors$r), y=fit$vectors$arrows[, 2] * sqrt(fit$vectors$r),
					                 signif=fit$vectors$pvals < 0.05)
					fv[, c('x', 'y')] <- fv[, c('x', 'y')] * 0.9 * ext / max(sqrt(fv$x^2 + fv$y^2), 1e-9)
					flen <- pmax(sqrt(fv$x^2 + fv$y^2), 1e-9)
					fv$lx <- fv$x / flen * 1.15 * ext
					fv$ly <- fv$y / flen * 1.15 * ext
				}

				nmdsp <- ggplot(nm, aes(x=x, y=y)) +
					geom_hline(yintercept=0, colour='grey85') + geom_vline(xintercept=0, colour='grey85') +
					geom_segment(aes(x=cx, y=cy, xend=x, yend=y, colour=HFT), alpha=0.35, linewidth=0.3) +
					(if (length(nm_hfts) <= 10) stat_ellipse(data=nm[nm$HFT %in% names(which(table(nm$HFT) >= 4)), ],
					                                         aes(colour=HFT), type='norm', level=0.68, linewidth=0.6)) +
					geom_point(aes(colour=HFT, size=ha), alpha=0.8)
				if (!is.null(fit)) {
					nmdsp <- nmdsp +
						geom_segment(data=fv, aes(x=x, y=y, xend=lx, yend=ly), colour='grey60', linewidth=0.3, linetype='dotted') +
						geom_segment(data=fv, aes(x=0, y=0, xend=x, yend=y), colour=ifelse(fv$signif, 'grey20', 'grey70'),
						             linewidth=0.5, arrow=arrow(length=unit(0.2, 'cm'))) +
						ggrepel::geom_label_repel(data=fv, aes(x=lx, y=ly, label=var), colour=ifelse(fv$signif, 'black', 'grey55'),
						                          size=3, lineheight=0.9, seed=1, box.padding=0.25, min.segment.length=Inf,
						                          max.overlaps=Inf, label.padding=unit(0.15, 'lines'), fill=alpha('white', 0.9))
				}
				nmdsp <- nmdsp +
					geom_label(data=cent, aes(x=cx, y=cy, label=HFT, fill=HFT), fontface='bold', size=3.5, show.legend=FALSE) +
					scale_colour_manual(values=nm_cols, labels=nm_labels, name='HFT (area)') +
					scale_fill_manual(values=nm_cols) +
					scale_size_area(max_size=7, name='Segment area (ha)') +
					guides(colour=guide_legend(override.aes=list(size=4, alpha=1, linetype=0))) +
					coord_equal() +
					labs(title=area, subtitle=paste0(plot_subtitle, sprintf(' | NMDS of the segments, stress = %.3f', mds_mh$stress),
					                                 if (!is.null(fit)) ' | vectors: variables fitted with envfit (grey: p >= 0.05)'),
					     x='NMDS1', y='NMDS2') +
					theme_bw(base_size=13)
				ggsave(file.path(outdir, paste('NMDS_', park, '_', res, tag, '_segms_mean.png', sep='')), nmdsp, width=12, height=9, dpi=300)
			}

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
			# Etiquetas cortas y en varias líneas, con margen a los lados para que no se corten
			radar <- CreateRadialPlot(scaled2, axis.labels = short_labels(clean_names), axis.label.size = 3.5, axis.label.offset = 1.25,
			                          plot.extent.x.sf = 1.9, plot.extent.y.sf = 1.35,
			                          group.colours = setNames(kolsn[scaled2$group], scaled2$group), group.order = hft_order)
			radar <- radar + labs(title=area, subtitle=plot_subtitle) + theme(plot.title=element_text(size=16), plot.subtitle=element_text(size=14))
			rpn <- file.path(outdir, paste('radarplot_', park, '_', res, tag, '_segms_mean.png', sep=''))
			ggsave(filename=rpn, plot=radar, width=10, height=7.5, dpi=300)

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
				hfts <- hft_order[hft_order %in% segm_pa_class$hclust_mean]
				segm_pa_class$HFT <- factor(segm_pa_class$hclust_mean, levels=hfts)

				# Bordes de los segmentos en gris fino para que destaque el contorno negro del área de estudio
				hftmap <- ggplot() +
					geom_sf(data=segm_pa_class, aes(fill=HFT), colour='grey50', linewidth=0.05) +
					scale_fill_manual(values=setNames(kolsn[hfts], hfts)) +
					labs(title=area, subtitle=plot_subtitle) +
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

				# --- BOXPLOTS DE CADA VARIABLE POR HFT CON LOS VALORES DE LOS PÍXELES ---
				if (!is.null(pixvals)) {
					cls <- terra::values(terra::rasterize(terra::vect(segm_pa_class), pixgrid, field='hclust_mean'), mat=FALSE)
					inhft <- !is.na(cls)
					hftcls <- factor(cls[inhft], levels=hfts)

					# Superficie de cada HFT (ha y % de la superficie clasificada) en la leyenda
					ha <- tapply(pixarea[inhft], hftcls, sum)
					ha[is.na(ha)] <- 0
					hft_labels <- paste0(hfts, ': ', format(round(ha), big.mark=',', trim=TRUE), ' ha (', round(100 * ha / sum(ha), 1), ' %)')

					var_labels <- short_labels(var_names)
					bstats <- do.call(rbind, lapply(seq_along(var_names), function(i) {
						s <- sapply(split(pixvals[inhft, i], hftcls), box_stats)
						data.frame(var=var_labels[i], HFT=hfts, t(s))
					}))
					bstats$var <- factor(bstats$var, levels=var_labels)
					bstats$HFT <- factor(bstats$HFT, levels=hfts)

					# Un panel por variable con su propia escala (valores originales)
					ncolf <- ceiling(sqrt(length(var_names)))
					nrowf <- ceiling(length(var_names) / ncolf)
					boxp <- ggplot(bstats, aes(x=HFT, ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax, fill=HFT)) +
						geom_boxplot(stat='identity', linewidth=0.3, width=0.7) +
						facet_wrap(~ var, ncol=ncolf, scales='free_y') +
						scale_fill_manual(values=setNames(kolsn[hfts], hfts), labels=hft_labels, name='HFT (area)') +
						labs(title=area, subtitle=paste0(plot_subtitle, ' | pixel values, whiskers to 1.5 x IQR'), x='HFT', y=NULL) +
						theme_bw(base_size=13)
					rpn3 <- file.path(outdir, paste('boxplot_', park, '_', res, tag, '_hft_pixels.png', sep=''))
					ggsave(rpn3, boxp, width=4.5 * ncolf + 3.5, height=3.5 * nrowf + 1.5, dpi=300)

					# --- MAPA DE CALOR DE LOS PERFILES DE LOS HFTs ---
					# Mediana de cada HFT como z respecto a la media y la desviación típica de los píxeles del
					# área; variables ordenadas por eta2 (parte de su varianza que explican los HFTs) y HFTs
					# en el orden del dendrograma
					pmean <- colMeans(pixvals[inhft, , drop=FALSE], na.rm=TRUE)
					psd <- apply(pixvals[inhft, , drop=FALSE], 2, sd, na.rm=TRUE)
					eta2 <- sapply(seq_along(var_names), function(i) {
						x <- pixvals[inhft, i]
						ok <- !is.na(x)
						m <- tapply(x[ok], hftcls[ok], mean)
						n <- tapply(x[ok], hftcls[ok], length)
						sum(n * (m - mean(x[ok]))^2, na.rm=TRUE) / sum((x[ok] - mean(x[ok]))^2)
					})
					vi <- match(as.character(bstats$var), var_labels)
					vord <- order(-eta2)
					col_labels <- paste0(var_labels, '\neta2 = ', sprintf('%.2f', eta2))
					hord <- unique(hclust_mean[hclust_mh$order])
					hord <- hord[hord %in% hfts]
					heat <- data.frame(x=match(vi, vord),
					                   row=factor(hft_labels[match(bstats$HFT, hfts)], levels=rev(hft_labels[match(hord, hfts)])),
					                   median=bstats$middle, z=(bstats$middle - pmean[vi]) / psd[vi])
					zlim <- max(abs(heat$z), 0.5, na.rm=TRUE)
					heat$txtcol <- ifelse(abs(heat$z) > 0.6 * zlim, 'white', 'black')
					# Columna estrecha con el color de cada HFT en el mapa
					swatch <- data.frame(row=factor(hft_labels, levels=levels(heat$row)))
					heatp <- ggplot(heat, aes(x=x, y=row)) +
						geom_tile(aes(fill=z), colour='white', linewidth=0.6) +
						geom_text(aes(label=signif(median, 3), colour=txtcol), size=3.4) +
						geom_tile(data=swatch, aes(x=0.25, y=row), fill=kolsn[hfts], width=0.25, height=0.8) +
						scale_colour_identity() +
						scale_fill_gradient2(low='#2166AC', mid='white', high='#B2182B', limits=c(-zlim, zlim), name='z of the\nmedian') +
						scale_x_continuous(breaks=seq_along(vord), labels=col_labels[vord], position='top', expand=c(0, 0.05)) +
						labs(title=area, subtitle=paste0(plot_subtitle, ' | median of the pixels of each HFT, coloured as z = (median - mean) / sd of the study area'),
						     x=NULL, y=NULL) +
						theme_minimal(base_size=13) +
						theme(panel.grid=element_blank(), axis.text.x.top=element_text(size=10))
					rpn4 <- file.path(outdir, paste('heatmap_', park, '_', res, tag, '_hft_profiles.png', sep=''))
					ggsave(rpn4, heatp, width=4.5 + 1.7 * length(var_names), height=2.5 + 0.4 * length(hfts), dpi=300)

					# --- BIPLOT DE LA PCA: CENTROIDES DE LOS HFTs Y VARIABLES COMO VECTORES ---
					# PCA de los píxeles con las variables estandarizadas; los vectores son las correlaciones
					# de cada variable con los ejes (escaladas al tamaño de la nube de HFTs)
					okpix <- inhft & complete.cases(pixvals)
					zpix <- scale(pixvals[okpix, , drop=FALSE])
					eg <- eigen(cor(zpix), symmetric=TRUE)
					ve <- 100 * eg$values / sum(eg$values)
					axes2 <- if (length(var_names) >= 3) 2:3 else 2
					scores <- zpix %*% eg$vectors[, c(1, axes2), drop=FALSE]
					loads <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)))
					hfpix <- factor(cls[okpix], levels=hfts)
					panels <- sprintf('PC1 (%.1f %%) - PC%d (%.1f %%)', ve[1], axes2, ve[axes2])
					cents <- ells <- arrs <- NULL
					for (a in seq_along(axes2)) {
						sc <- scores[, c(1, a + 1)]
						for (h in hfts) {
							s <- sc[hfpix == h, , drop=FALSE]
							if (nrow(s) < 3) next
							mu <- colMeans(s)
							# Elipse que contiene la mitad de los píxeles del HFT (aproximación normal)
							e <- eigen(cov(s), symmetric=TRUE)
							ang <- seq(0, 2 * pi, length.out=100)
							xy <- t(mu + e$vectors %*% diag(sqrt(pmax(e$values, 0) * qchisq(0.5, 2))) %*% rbind(cos(ang), sin(ang)))
							cents <- rbind(cents, data.frame(panel=panels[a], HFT=h, x=mu[1], y=mu[2], ha=ha[as.character(h)], colour=kolsn[h]))
							ells <- rbind(ells, data.frame(panel=panels[a], HFT=h, x=xy[, 1], y=xy[, 2], colour=kolsn[h]))
						}
						arrs <- rbind(arrs, data.frame(panel=panels[a], var=var_labels, x=loads[, 1], y=loads[, axes2[a]]))
					}
					# Con muchos HFTs las elipses se convierten en una maraña: solo se dibujan con 10 o menos
					show_ell <- length(hfts) <= 10
					ext <- if (show_ell) max(abs(c(ells$x, ells$y))) else 1.3 * max(abs(c(cents$x, cents$y)))
					arrscale <- 0.9 * ext / max(sqrt(arrs$x^2 + arrs$y^2))
					arrs$x <- arrs$x * arrscale
					arrs$y <- arrs$y * arrscale
					# Etiquetas de las variables en un anillo por fuera de las elipses, en la dirección de su
					# vector, para que no tapen los centroides (que quedan cerca del centro)
					rlab <- 1.1 * ext
					alen <- pmax(sqrt(arrs$x^2 + arrs$y^2), 1e-9)
					arrs$lx <- arrs$x / alen * rlab
					arrs$ly <- arrs$y / alen * rlab
					for (d in c('cents', 'ells', 'arrs')) assign(d, transform(get(d), panel=factor(panel, levels=panels)))
					biplot <- ggplot() +
						geom_hline(yintercept=0, colour='grey80') + geom_vline(xintercept=0, colour='grey80') +
						# Elipses tenues: los HFTs se solapan mucho en el espacio de los píxeles
						(if (show_ell) geom_path(data=ells, aes(x=x, y=y, group=HFT, colour=colour), linewidth=0.4, alpha=0.45)) +
						geom_segment(data=arrs, aes(x=x, y=y, xend=lx, yend=ly), colour='grey60', linewidth=0.3, linetype='dotted') +
						geom_segment(data=arrs, aes(x=0, y=0, xend=x, yend=y), colour='grey25', linewidth=0.5,
						             arrow=arrow(length=unit(0.2, 'cm'))) +
						geom_point(data=cents, aes(x=x, y=y, fill=colour, size=ha), shape=21, colour='black') +
						geom_text(data=cents, aes(x=x, y=y, label=HFT), size=3, fontface='bold') +
						ggrepel::geom_label_repel(data=arrs, aes(x=lx, y=ly, label=var), size=3, lineheight=0.9, seed=1,
						                          box.padding=0.25, min.segment.length=Inf, max.overlaps=Inf,
						                          label.padding=unit(0.15, 'lines'), fill=alpha('white', 0.9)) +
						scale_colour_identity() + scale_fill_identity() +
						scale_size(range=c(5, 13), name='HFT area (ha)') +
						facet_wrap(~ panel, ncol=2) + coord_equal() +
						labs(title=area, subtitle=paste0(plot_subtitle, ' | PCA of the pixels (standardized variables): HFT centroids',
						                                 if (show_ell) ', ellipses with half of their pixels', ', variables as vectors'),
						     x='PC1', y=paste(paste0('PC', axes2), collapse=' / ')) +
						theme_bw(base_size=13)
					rpn5 <- file.path(outdir, paste('pca_biplot_', park, '_', res, tag, '_hft.png', sep=''))
					ggsave(rpn5, biplot, width=7 * length(axes2) + 2.5, height=8, dpi=300)
				}
			}
		}
	}
}

run_state('done')
