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

# Colores jerárquicos para cortes anidados del mismo árbol (SIMIL = 3), como los "tree colors" de Tennekes
# y de Jonge (2014) en el espacio HCL: los HFTs del corte con menos HFTs se reparten el círculo de tonos, y
# cada HFT reparte su parte entre los HFTs del corte siguiente que contiene (en proporción a los HFTs del
# corte más fino que hay en cada uno y en el orden del dendrograma). El tono de cada HFT es el centro de su
# parte, así que una familia tiene tonos próximos en todos los cortes y distintos de los de las otras; entre
# hermanos alterna además la luminosidad. Un HFT que no cambia de un corte al siguiente conserva su color.
# cls_list: clasificaciones (una por corte), de la más general a la más detallada. frac: parte del círculo
# que usan los HFTs del corte más general (el resto, huecos entre familias); frac_in: lo mismo dentro de
# cada HFT. lum y chroma (HCL) elegidos con el validador de paletas: con hasta unos 10 HFTs los colores
# vecinos en la leyenda se distinguen también con daltonismo
hier_colours <- function(tree, cls_list, frac=0.75, frac_in=1, lum=c(47, 69), chroma=80) {
	fine <- cls_list[[length(cls_list)]]
	out <- list()
	prev <- NULL
	for (lv in seq_along(cls_list)) {
		cl <- cls_list[[lv]]
		ord <- unique(cl[tree$order])
		parent <- if (lv == 1) rep(1, length(ord)) else sapply(ord, function(h) cls_list[[lv - 1]][cl == h][1])
		cur <- list(lo=numeric(max(cl)), hi=numeric(max(cl)), col=character(max(cl)))
		for (p in unique(parent)) {
			kids <- ord[parent == p]
			a <- if (lv == 1) 0 else prev$lo[p]
			b <- if (lv == 1) 360 else prev$hi[p]
			if (lv > 1 && length(kids) == 1) {
				cur$lo[kids] <- a
				cur$hi[kids] <- b
				cur$col[kids] <- prev$col[p]
				next
			}
			nfine <- sapply(kids, function(h) length(unique(fine[cl == h])))
			edges <- a + (b - a) * c(0, cumsum(nfine)) / sum(nfine)
			for (j in seq_along(kids)) {
				pad <- (1 - (if (lv == 1) frac else frac_in)) / 2 * (edges[j + 1] - edges[j])
				cur$lo[kids[j]] <- edges[j] + pad
				cur$hi[kids[j]] <- edges[j + 1] - pad
				L <- if (lv == 1) 60 else lum[(j - 1) %% length(lum) + 1]
				cur$col[kids[j]] <- hcl((cur$lo[kids[j]] + cur$hi[kids[j]]) / 2 + 15, chroma, L, fixup=TRUE)
			}
		}
		out[[lv]] <- cur$col
		prev <- cur
	}
	setNames(out, names(cls_list))
}

# Árbol de los HFTs de los cortes (SIMIL = 3), del corte más general (arriba) al más detallado (abajo): cada HFT es
# un círculo con su número y su color, de tamaño según su superficie, unido a los HFTs del corte siguiente con los
# que comparte segmentos (grosor según la superficie compartida). Así se ve en qué HFT de los cortes más generales
# está cada HFT y en qué HFTs de los más detallados se divide. Los HFTs del corte más detallado van en el orden
# del dendrograma y los de los otros cortes en el centro de los suyos. Con 'ward_area_kmeans' los cortes pueden no
# estar anidados: un HFT puede venir de varios del corte de arriba, y la opacidad de cada unión es la parte de la
# superficie del HFT de abajo que viene del de arriba.
# cls_list: HFT de cada segmento en cada corte, del más general al más detallado; cols_list: colores de los HFTs
# de cada corte (por nº de HFT); ha: superficie de cada segmento; labels: rótulo de cada corte; hl: corte
# resaltado (el sugerido). Devuelve la figura, su tamaño y una tabla con cada combinación de HFTs de los cortes
# (con cortes anidados, una por HFT del corte más detallado) y su superficie
hft_tree <- function(tree, cls_list, cols_list, ha, labels, hl=NA) {
	nl <- length(cls_list)
	cuts <- names(cls_list)
	yl <- nl + 1 - seq_len(nl)
	# Posición horizontal de los HFTs, de abajo arriba: la media de las de los HFTs del corte de abajo, cada uno
	# con la parte de su superficie que viene de este HFT (con cortes anidados, toda o nada)
	ord <- unique(cls_list[[nl]][tree$order])
	xs <- vector('list', nl)
	xs[[nl]] <- setNames(seq_along(ord), ord)
	edges <- NULL
	for (lv in rev(seq_len(nl - 1))) {
		ov <- tapply(ha, list(cls_list[[lv]], cls_list[[lv + 1]]), sum, default=0)
		share <- sweep(ov, 2, colSums(ov), '/')
		xs[[lv]] <- setNames(as.vector(share %*% xs[[lv + 1]][colnames(ov)]) / rowSums(share), rownames(ov))
		e <- which(ov > 0, arr.ind=TRUE)
		edges <- rbind(edges, data.frame(x=xs[[lv]][e[, 1]], y=yl[lv], xend=xs[[lv + 1]][colnames(ov)[e[, 2]]], yend=yl[lv + 1],
		                                 ha=ov[e], share=share[e]))
	}
	nodes <- do.call(rbind, lapply(seq_len(nl), function(lv) {
		a <- tapply(ha, cls_list[[lv]], sum)
		data.frame(HFT=names(a), x=xs[[lv]][names(a)], y=yl[lv], ha=as.vector(a), fill=cols_list[[lv]][as.integer(names(a))])
	}))
	# Número en blanco sobre los colores oscuros
	nodes$txt <- ifelse(colSums(col2rgb(nodes$fill) / 255 * c(0.299, 0.587, 0.114)) < 0.55, 'white', 'black')
	nested <- is.null(edges) || all(abs(edges$share - 1) < 1e-9)

	p <- ggplot() +
		(if (!is.na(hl) && hl %in% cuts) annotate('rect', xmin=-Inf, xmax=Inf, ymin=yl[match(hl, cuts)] - 0.45,
		                                          ymax=yl[match(hl, cuts)] + 0.45, fill='#eb6834', alpha=0.1)) +
		geom_hline(yintercept=yl, colour='grey90') +
		(if (!is.null(edges)) geom_segment(data=edges, aes(x=x, y=y, xend=xend, yend=yend, linewidth=ha, alpha=share),
		                                   colour='grey40', lineend='round')) +
		geom_point(data=nodes, aes(x=x, y=y, size=ha, fill=fill), shape=21, colour='grey20', stroke=0.4) +
		geom_text(data=nodes, aes(x=x, y=y, label=HFT, colour=txt), size=3.2, fontface='bold') +
		scale_fill_identity() + scale_colour_identity() +
		scale_size(range=c(6, 16), name='HFT area (ha)') +
		scale_linewidth(range=c(0.4, 4), guide='none') +
		scale_alpha(range=c(0.15, 1), limits=c(0, 1), breaks=c(0.25, 0.5, 0.75, 1), labels=c('25 %', '50 %', '75 %', '100 %'),
		            name='Part of the lower HFT\nfrom the upper one', guide=if (nested) 'none' else 'legend') +
		scale_x_continuous(expand=expansion(add=0.8)) +
		scale_y_continuous(breaks=yl, labels=labels, limits=c(0.5, nl + 0.5), expand=c(0, 0)) +
		labs(x=NULL, y=NULL) +
		theme_minimal(base_size=13) +
		theme(panel.grid=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=12, hjust=0))

	# Combinaciones de HFTs de los cortes, en el orden del árbol
	comb <- as.data.frame(cls_list, check.names=FALSE)
	comb$segments <- 1
	comb$ha <- ha
	paths <- aggregate(cbind(segments, ha) ~ ., comb, sum)
	paths <- paths[do.call(order, lapply(seq_len(nl), function(lv) xs[[lv]][as.character(paths[[cuts[lv]]])])), ]
	paths$pct <- round(100 * paths$ha / sum(paths$ha), 2)
	paths$ha <- round(paths$ha, 2)
	list(plot=p, paths=paths, width=max(9, 0.5 * length(ord) + 4), height=1.3 * nl + 2.5)
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

# Distancias entre los segmentos del método original ('ward'): estadístico (media o mediana) y varianza
# de cada variable, reescalados a 0-1 entre los segmentos
segm_dist <- function(hri) {
	nv <- (ncol(hri) - 3) / 3
	vegdist(as.data.frame(lapply(hri[, 3:(2 + 2 * nv)], rescale01)), "euclidean", na.rm=T)
}

# Árbol de la clasificación de los segmentos según CLASS_METHOD (ver config.py) y distancias entre ellos
# (para la correlación cofenética). Con 'ward_area' y 'ward_area_kmeans', Ward ponderado por el
# nº de píxeles de cada segmento sobre su estadístico de cada variable tipificado con la desviación típica
# de los píxeles (dentro + entre segmentos): con la distancia de Ward entre grupos y members = nº de
# píxeles, hclust da el mismo árbol que el Ward de todos los píxeles partiendo de los segmentos
class_tree <- function(hri, method) {
	if (method == 'ward') {
		d <- segm_dist(hri)
		return(list(tree=hclust(d, "ward.D2"), dist=d))
	}
	nv <- (ncol(hri) - 3) / 3
	m <- as.matrix(hri[, 3:(2 + nv)])
	v <- as.matrix(hri[, (3 + nv):(2 + 2 * nv)])
	n <- hri$sumpamask
	z <- sapply(seq_len(nv), function(i) {
		ok <- !is.na(m[, i]) & !is.na(v[, i])
		mu <- sum(n[ok] * m[ok, i]) / sum(n[ok])
		s <- sqrt((sum(n[ok] * v[ok, i]) + sum(n[ok] * (m[ok, i] - mu)^2)) / sum(n[ok]))
		x <- (m[, i] - mu) / (if (is.finite(s) && s > 0) s else 1) # variable constante: no influye
		x[is.na(x)] <- 0 # sin datos en el segmento: la media de la variable
		x
	})
	d <- dist(z)
	w <- sqrt(2 * outer(n, n) / outer(n, n, '+')) * as.matrix(d)
	list(tree=hclust(as.dist(w), "ward.D2", members=n), dist=d, z=z, n=n)
}

# Consolidación con k-medias ponderadas por el nº de píxeles ('ward_area_kmeans'): cada segmento pasa al
# HFT de media más próxima, hasta que no cambia ninguno; los HFTs conservan los números de cl
wkmeans <- function(z, n, cl, iter=100) {
	for (it in seq_len(iter)) {
		lev <- sort(unique(cl))
		cen <- rowsum(z * n, cl) / as.vector(rowsum(n, cl))
		d2 <- sapply(seq_len(nrow(cen)), function(j) colSums((t(z) - cen[j, ])^2))
		new <- lev[max.col(-d2, ties.method='first')]
		if (all(new == cl)) break
		cl <- new
	}
	cl
}

# Clasificación de los segmentos en k HFTs con el árbol de class_tree()
classify_k <- function(ct, k, method) {
	cl <- cutree(ct$tree, k=k)
	if (method == 'ward_area_kmeans') cl <- wkmeans(ct$z, ct$n, cl)
	cl
}

# CLIP_PERCENTILE: cada variable recortada a sus percentiles p y 100 - p en los píxeles de los segmentos (inseg),
# para la clasificación: unos pocos píxeles extremos no deciden la desviación típica con que se tipifica
clip_pixels <- function(pix, inseg, p) {
	for (i in seq_len(ncol(pix))) {
		q <- quantile(pix[inseg, i], c(p, 100 - p) / 100, na.rm=TRUE, names=FALSE)
		pix[, i] <- pmin(pmax(pix[, i], q[1]), q[2])
	}
	pix
}

# hri con el estadístico (media, o mediana con STAT_AGG = 'median') y la varianza de cada variable en cada
# segmento calculados con los píxeles recortados (seg: fila de hri de cada píxel), para class_tree()
clip_hri <- function(hri, seg, pixc, median_stat) {
	nv <- (ncol(hri) - 3) / 3
	ss <- segm_sums(seg, pixc, nrow(hri))
	for (i in seq_len(nv)) {
		ok <- !is.na(seg) & !is.na(pixc[, i])
		m <- ss[[i]]$s / ss[[i]]$n # media de cada segmento, centrada en la de todos sus píxeles
		hri[, 2 + i] <- m + mean(pixc[ok, i])
		if (median_stat) {
			med <- tapply(pixc[ok, i], seg[ok], median)
			hri[, 2 + i] <- NA
			hri[as.integer(names(med)), 2 + i] <- med
		}
		hri[, 2 + nv + i] <- ss[[i]]$q / ss[[i]]$n - m^2
	}
	hri
}

# Sumas de los píxeles de cada segmento para cada variable (centrada): nº de píxeles, suma y suma de
# cuadrados. Con ellas la varianza explicada (R2) de cualquier agrupación de los segmentos es exacta,
# tanto con STAT_AGG = 'mean' como 'median'. seg: segmento (1 a nseg) de cada píxel, NA fuera de ellos
segm_sums <- function(seg, pix, nseg) {
	lapply(seq_len(ncol(pix)), function(i) {
		ok <- !is.na(seg) & !is.na(pix[, i])
		x <- pix[ok, i] - mean(pix[ok, i])
		g <- seg[ok]
		tot <- function(v) { out <- numeric(nseg); r <- rowsum(v, g); out[as.integer(rownames(r))] <- r; out }
		list(n=tot(rep(1, length(x))), s=tot(x), q=tot(x^2))
	})
}

# Segmento (fila de hri) de cada píxel de la malla de las variables, del umbral tk (NA fuera de los segmentos)
segm_cells <- function(park, tk, hri, grid) {
	v <- terra::vect(file.path('shp', paste0('park_segm_', park, '_', tk, '_diss.shp')))
	v$row <- match(as.character(v$segm_id), format(hri$segm_id, scientific=FALSE, trim=TRUE))
	terra::values(terra::rasterize(v, grid, field='row'), mat=FALSE)
}

# R2 de cada variable: varianza de los píxeles entre los grupos (cl: grupo de cada segmento) / varianza total
r2_groups <- function(ss, cl) {
	sapply(ss, function(v) {
		n <- rowsum(v$n, cl)
		s <- rowsum(v$s, cl)
		corr <- sum(v$s)^2 / sum(v$n)
		(sum(s^2 / n, na.rm=TRUE) - corr) / (sum(v$q) - corr)
	})
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

simil <- read_config_num('SIMIL', 2) # ver SIMIL en config.py
upper <- 2 # 1: default, 2: límite superior incrementado

if (!simil %in% 1:3) stop("SIMIL en config.py debe ser 1, 2 o 3")

# Método con el que moranvar.py eligió el umbral (THRESHOLD_METHOD en config.py), para los títulos
thr_method <- read_config('THRESHOLD_METHOD')
if (is.na(thr_method) || thr_method == '') thr_method <- 'f'
method_label <- switch(thr_method,
	f = sprintf('F of Johnson et al. (2015), alpha = %g', read_config_num('THRESHOLD_ALPHA', 0.5)),
	gs = 'Global Score of Espindola et al. (2006)',
	gs_fixed = 'Global Score with fixed bounds',
	crossing = 'first crossing of the curves - 1',
	crossing_first = 'first crossing of the curves',
	paste('fixed threshold', thr_method))

# Método de clasificación de los segmentos en HFTs (CLASS_METHOD en config.py)
class_method <- read_config('CLASS_METHOD')
if (is.na(class_method) || class_method == '') class_method <- 'ward_area'
if (!class_method %in% c('ward_area', 'ward_area_kmeans', 'ward')) {
	stop(paste0("CLASS_METHOD en config.py debe ser 'ward_area', 'ward_area_kmeans' o 'ward', no '", class_method, "'"))
}
class_label <- switch(class_method,
	ward_area = 'Ward weighted by area on standardized means',
	ward_area_kmeans = 'Ward weighted by area on standardized means + weighted k-means',
	ward = 'Ward on means and variances rescaled 0-1')

# Recorte de los valores extremos para la clasificación (CLIP_PERCENTILE en config.py; 0 = sin recorte)
clip_pct <- read_config_num('CLIP_PERCENTILE', 0)
if (clip_pct < 0 || clip_pct >= 50) stop(paste('CLIP_PERCENTILE en config.py debe ir de 0 a menos de 50, no', clip_pct))
if (clip_pct > 0) class_label <- paste0(class_label, ', variables clipped to percentiles ', clip_pct, '-', 100 - clip_pct)
median_stat <- identical(read_config('STAT_AGG'), 'median')

# Resultados en results/<área de estudio>/<escenario>: el área de los csv (anotada por segmentation_pca_par.py)
# o STUDY_AREA de config.py, y SCENARIO de config.py; si falta o está vacío, el nombre de la subcarpeta de las
# variables (INVARS_FOLDER, con _ en vez de los caracteres que no son letras, números, - o _), como
# scenario_name() de hft_state.py, o 'base' si las variables están directamente en ../inVars
area <- if (file.exists('csv/study_area.csv')) readLines('csv/study_area.csv', n=1, warn=FALSE) else read_config('STUDY_AREA')
if (is.na(area) || area == '') stop("No se encuentra el área de estudio: falta STUDY_AREA en config.py")
scenario <- read_config('SCENARIO')
if (is.na(scenario) || scenario == '') {
	folder <- read_config('INVARS_FOLDER')
	Encoding(folder) <- 'UTF-8' # config.py está en UTF-8: una ñ es un carácter, como en Python
	folder <- if (is.na(folder)) '' else gsub('\\\\', '/', gsub('^[/\\\\]+|[/\\\\]+$', '', trimws(folder)))
	scenario <- if (folder == '') 'base' else gsub('[^A-Za-z0-9_-]', '_', folder, perl=TRUE)
}
if (!grepl('^[A-Za-z0-9_-]+$', scenario) || grepl('_previous$', scenario)) {
	stop(paste0("SCENARIO en config.py ('", scenario, "') solo puede tener letras, números, - y _, y no puede acabar en _previous"))
}
outroot <- file.path('results', area, scenario)

# Cortes del árbol: simil = 3 genera una clasificación por corte (10, 25, 50, 75 y 90 % de la altura máxima),
# cada una en su subcarpeta, de la más detallada (p10, más HFTs) a la más general (p90, menos HFTs)
if (simil == 3) {
	cuts <- c(p10=0.10, p25=0.25, p50=0.50, p75=0.75, p90=0.90)
} else {
	cuts <- c(results=NA)
}

# Metadatos de la ejecución: opciones de config.py (sin comentarios) y de este script
cfg <- if (file.exists('config.py')) readLines('config.py', warn=FALSE) else character(0)
cfg <- sub('\\s+#.*$', '', cfg)
cfg <- cfg[!grepl('^\\s*(#|$)', cfg)]
# Los resultados anteriores de este escenario se mueven a results/<área>/<escenario>_previous (sustituyendo
# los que hubiera), para no mezclar ficheros de ejecuciones con otras opciones (p. ej. otro umbral en el nombre)
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
	paste('Scenario:', scenario),
	paste('IDs processed (csv/segm_done.csv):', paste(ecox_list, collapse=', ')),
	'',
	'## config.py',
	cfg,
	'',
	'## moranvar_plots.R',
	paste('upper =', upper, '# max. HFTs (SIMIL 1 and 2): 1: log10(n. pixels), 2: log6(n. pixels)')
), file.path(outroot, 'processing_metadata.txt'))

for (pmx in 1:mx) {
	park = ecox_list[pmx]
	print(paste("Procesando park ID:", park))
	
	# Umbral elegido por moranvar.py (THRESHOLD_METHOD en config.py) y puntuaciones de cada umbral
	namet <- paste('csv/', park, '_threshold.csv', sep='')
	if (!file.exists(namet)) { print(paste("WARNING: no threshold chosen by moranvar.py, park skipped:", namet)); next }
	thr_tab <- read.table(namet, sep=' ', header=TRUE)
	thr_all <- thr_tab$threshold # thresholds actually evaluated (some may be missing)
	res <- thr_tab$threshold[thr_tab$selected == 1][1]

	name2 <- paste('csv/', park, '_movar_thresholds.csv', sep='')
	df02 <- read.table(name2, sep=' ', header=F)
	if (upper == 1) { upplim <- ceiling(log10(df02[1,3])) }
	if (upper == 2) { upplim <- ceiling(log(df02[1,3], 6)) }
	print(paste('Upper limit is ', upplim, sep=''))

	# Curvas reescaladas (media de las variables) y puntuación del método, con el umbral elegido
	x1 <- thr_tab$moran_norm
	x2 <- thr_tab$var_norm
	png(file.path(outroot, paste(park, '_moran_var_mean.png', sep='')), width=640, height=540)
	par(mar=c(4.5, 4.5, 4.5, 1))
	plot(thr_all, x1, ylim=c(0, 1.35), col=1, typ='o', pch=19, xaxt='n', xlab='Similarity threshold (x 0.1)',
	     ylab="Rescaled Moran's I, variance and score", main=paste0(park, ': selected threshold 0.', res))
	axis(1, at=thr_all)
	mtext(method_label, side=3, line=0.4)
	lines(thr_all, x2, col=3, typ='o', pch=19)
	leg <- c("M.I: Moran's I (low: neighbours differ)", 'Var: 1 - variance (high: homogeneous segments)')
	legcol <- c(1, 3)
	if (!all(is.na(thr_tab$score))) {
		sc <- rescale01(thr_tab$score)
		if (thr_method %in% c('gs', 'gs_fixed')) sc <- 1 - sc # 1 = mejor umbral
		lines(thr_all, sc, col=4, typ='o', pch=19, lwd=2)
		leg <- c(leg, 'Score of the method, rescaled (1: best)')
		legcol <- c(legcol, 4)
	}
	if (thr_method %in% c('crossing', 'crossing_first')) {
		# Cruces de las curvas, interpolados sobre los umbrales (como en moranvar.py)
		ip <- which(diff(x1 > x2) != 0)
		tp <- (x2[ip] - x1[ip]) / ((x1[ip + 1] - x1[ip]) - (x2[ip + 1] - x2[ip]))
		points(thr_all[ip] + tp * (thr_all[ip + 1] - thr_all[ip]), x1[ip] + (x1[ip + 1] - x1[ip]) * tp, col='red', pch=19)
		leg <- c(leg, 'Crossings')
		legcol <- c(legcol, 'red')
	}
	abline(v=res, lty=2, col='grey40')
	legend('top', legend=leg, col=legcol, lty=1, pch=19, bty='n', cex=0.9)
	dev.off()

	k <- paste('0.', res, sep='')
	write.table(thr_tab, file.path(outroot, paste(park, '_optim_thresholds.csv', sep='')), sep=' ', row.names=F, quote=F)

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
	# Superficie de un píxel (ha), para la de los segmentos
	cell_ha <- if (!is.null(pixvals)) median(pixarea, na.rm=TRUE) else read_config_num('RESOLUTION', 1)^2 / 1e4

	# Sumas de los píxeles de los segmentos del umbral elegido, para el R2 de los HFTs de cada corte
	ss_sel <- NULL
	seg_sel <- NULL
	cut_r2 <- NULL
	if (!is.null(pixvals) && file.exists(file.path('shp', paste0('park_segm_', park, '_', res, '_diss.shp')))) {
		seg_sel <- segm_cells(park, res, hri, pixgrid)
		ss_sel <- segm_sums(seg_sel, pixvals, nrow(hri))
	}

	# Clasificación de los segmentos (CLASS_METHOD) y, con CLIP_PERCENTILE, con los estadísticos de los píxeles
	# recortados (solo para el árbol: las gráficas y los estadísticos de los HFTs usan todos los valores)
	pixclip <- NULL
	if (clip_pct > 0) {
		if (is.null(seg_sel)) {
			print('WARNING: CLIP_PERCENTILE needs the pixel values of ../inVars: the segments are classified without clipping')
		} else {
			pixclip <- clip_pixels(pixvals, !is.na(seg_sel), clip_pct)
		}
	}
	ct <- class_tree(if (is.null(pixclip)) hri else clip_hri(hri, seg_sel, pixclip, median_stat), class_method)
	dmh <- ct$dist
	hclust_mh <- ct$tree
	cophval <- 0
	try(chcl <- cophenetic(hclust_mh))
	try(cophval <- cor(chcl, dmh))

	# --- CORRELACIÓN ENTRE LAS VARIABLES ---
	# Correlación de Pearson entre cada par de variables en los píxeles de los segmentos (con datos en las dos),
	# en el orden de ENV_VARS. Las variables muy correlacionadas repiten el mismo gradiente, que cuenta más de una
	# vez en la clasificación (cada variable tipificada pesa lo mismo); |r| >= 0.7 en negrita
	if (!is.null(seg_sel) && num_vars > 1) {
		cm <- suppressWarnings(cor(pixvals[!is.na(seg_sel), , drop=FALSE], use='pairwise.complete.obs'))
		dimnames(cm) <- list(var_names, var_names)
		write.table(data.frame(var=var_names, cm, check.names=FALSE), file.path(outroot, paste0('cor_', park, '_vars.csv')),
		            sep=' ', row.names=F, quote=F)
		cor_labels <- short_labels(var_names)
		low <- which(lower.tri(cm), arr.ind=TRUE)
		cd <- data.frame(x=factor(cor_labels[low[, 2]], levels=cor_labels[-num_vars]),
		                 y=factor(cor_labels[low[, 1]], levels=rev(cor_labels[-1])), r=cm[low])
		corp <- ggplot(cd, aes(x=x, y=y)) +
			geom_tile(aes(fill=r), colour='white', linewidth=0.8) +
			geom_text(aes(label=sprintf('%.2f', r), colour=ifelse(abs(r) > 0.6, 'white', 'black'),
			              fontface=ifelse(abs(r) >= 0.7, 'bold', 'plain')), size=4) +
			scale_colour_identity() +
			scale_fill_gradient2(low='#2166AC', mid='white', high='#B2182B', limits=c(-1, 1), na.value='grey90', name='r') +
			scale_x_discrete(position='top') +
			coord_equal() +
			labs(title=area, subtitle=paste0('ID ', park, ' | Pearson correlation between the variables in the pixels of the segments\n',
			                                 '|r| >= 0.7 in bold: such variables repeat the same gradient, which counts more than once in the classification'),
			     x=NULL, y=NULL) +
			theme_minimal(base_size=13) +
			theme(panel.grid=element_blank(), axis.text.x.top=element_text(size=10), axis.text.y=element_text(size=10))
		ggsave(file.path(outroot, paste0('cor_', park, '_vars.png')), corp, width=3 + 1.6 * num_vars, height=2.5 + 1.3 * num_vars, dpi=300)
	}

	# SIMIL = 3: colores jerárquicos entre los cortes (ver hier_colours), solo para cortes con más de un HFT
	hier_cols <- list()
	tree_cls <- tree_cols <- list() # HFTs y colores de cada corte, para el árbol de los cortes (hft_tree)
	if (simil == 3) {
		cls_list <- lapply(cuts, function(f) cutree(hclust_mh, h=f * max(hclust_mh$height)))
		nk <- sapply(cls_list, function(cl) length(unique(cl)))
		cls_list <- cls_list[nk > 1][order(nk[nk > 1])]
		if (length(cls_list) > 0) hier_cols <- hier_colours(hclust_mh, cls_list)
	}

	# Corte sugerido: el codo de la curva de R2 del umbral elegido con 2 a 30 HFTs (Kneedle, Satopää et al. 2011:
	# el punto más alejado por encima de la recta entre sus extremos, con los dos ejes reescalados a 0-1), a
	# partir del cual cada HFT más añade poca varianza explicada, y, con SIMIL = 3, el corte con el número de HFTs
	# más próximo a él (en escala logarítmica; a igual distancia, el de menos HFTs). Sus archivos llevan _optim
	knee <- NA
	suggested <- NA
	if (!is.null(ss_sel) && nrow(hri) >= 5) {
		nh <- 2:min(nrow(hri) - 1, 30)
		r2k <- sapply(nh, function(h) mean(r2_groups(ss_sel, classify_k(ct, h, class_method)), na.rm=TRUE))
		knee <- nh[which.max(rescale01(r2k) - rescale01(nh))]
		if (simil == 3) {
			nk_cut <- sapply(cuts, function(f) length(unique(cutree(hclust_mh, h=f * max(hclust_mh$height)))))
			o <- order(nk_cut)
			suggested <- names(cuts)[o[which.min(abs(log(nk_cut[o]) - log(knee)))]]
		}
	}

	for (cn in names(cuts)) {
		outdir <- if (simil == 3) file.path(outroot, cn) else outroot
		# el corte también va en el nombre de cada fichero, con _optim si es el sugerido
		tag <- if (simil == 3) paste0('_', cn, if (identical(cn, suggested)) '_optim') else ''
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

		# 'ward_area_kmeans': los segmentos se reasignan al HFT de media más próxima (los colores y el orden
		# siguen los del árbol, y los recuadros del dendrograma muestran los HFTs antes de reasignarlos)
		if (class_method == 'ward_area_kmeans' && ncl > 1 && ncl < nrow(hri)) {
			hclust_mean <- wkmeans(ct$z, ct$n, hclust_mean)
			ncl <- length(unique(hclust_mean))
		}

		# Varianza de los píxeles explicada por los HFTs de este corte (media y por variable)
		if (!is.null(ss_sel)) {
			r2v <- r2_groups(ss_sel, hclust_mean)
			cut_r2 <- rbind(cut_r2, data.frame(cut=cn, hfts=ncl, r2=mean(r2v, na.rm=TRUE),
			                                   as.list(setNames(r2v, var_names)), check.names=FALSE))
		}

		# Colores del mapa: Set1 hasta 8 HFTs, interpolados si hay más; con SIMIL = 3, jerárquicos
		kolsn <- if (max(hclust_mean) <= length(kols)) kols else colorRampPalette(kols)(max(hclust_mean))
		if (!is.null(hier_cols[[cn]])) kolsn <- hier_cols[[cn]]
		# Orden de los HFTs en leyendas y ejes: con colores jerárquicos el del dendrograma (los tonos de un
		# mismo color quedan juntos), si no el numérico
		hft_order <- if (!is.null(hier_cols[[cn]])) unique(hclust_mean[hclust_mh$order]) else sort(unique(hclust_mean))
		# HFTs y colores de este corte en el árbol de los cortes: un solo HFT, gris, si el corte no separa ninguno
		# (los únicos sin colores jerárquicos con SIMIL = 3)
		if (simil == 3) {
			tree_cls[[cn]] <- if (is.null(hier_cols[[cn]])) rep(1, nrow(hri)) else hclust_mean
			tree_cols[[cn]] <- if (is.null(hier_cols[[cn]])) 'grey60' else kolsn
		}

		if (ncl > 1) {
			# Título de las figuras: área de estudio, umbral, corte del árbol y nº de HFTs
			plot_subtitle <- paste0('ID ', park, ' | threshold ', k, if (simil == 3) paste0(' | cut ', cn, if (identical(cn, suggested)) ' (suggested)'),
			                        ' | ', ncl, ' HFTs')

			png(file.path(outdir, paste('hclust_', park, '_', res, tag, '_segms_mean.png', sep='')))
			try(plot(hclust_mh, hang=-1, main=park, sub=cophval, xlab=class_label)); try(rect.hclust(hclust_mh, k=ncl))
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

	# --- ÁRBOL DE LOS HFTs DE LOS CORTES (SIMIL = 3) ---
	# En qué HFT de los cortes más generales está cada HFT y en qué HFTs de los más detallados se divide (ver
	# hft_tree). La tabla tiene una fila por HFT del corte más detallado (o por combinación de HFTs, si los cortes
	# no están anidados) con el HFT que la contiene en cada corte, su nº de segmentos y su superficie
	if (simil == 3) {
		rows <- rev(names(cuts)) # del corte más general al más detallado
		nk_tree <- sapply(tree_cls[rows], function(cl) length(unique(cl)))
		r2_tree <- if (!is.null(cut_r2)) ifelse(nk_tree == 1, 0, cut_r2$r2[match(rows, cut_r2$cut)])
		tree_labels <- paste0(rows, '\n', nk_tree, ifelse(nk_tree == 1, ' HFT', ' HFTs'), if (!is.null(r2_tree)) sprintf('\nR2 %.2f', r2_tree),
		                      ifelse(rows %in% suggested, '\n(suggested)', ''))
		ht <- hft_tree(hclust_mh, tree_cls[rows], tree_cols[rows], hri$sumpamask * cell_ha, tree_labels, suggested)
		treep <- ht$plot + labs(title=area, subtitle=paste0('ID ', park, ' | threshold ', k,
			' | HFTs of the cuts of the tree, from broad (top) to detailed (bottom)\n',
			'circles: HFTs with their number and colour of the maps, sized by area | lines: area shared by the HFTs of consecutive cuts'))
		ggsave(file.path(outroot, paste0('hft_tree_', park, '_', res, '.png')), treep, width=ht$width, height=ht$height, dpi=300)
		write.table(ht$paths, file.path(outroot, paste0('hft_tree_', park, '_', res, '.csv')), sep=' ', row.names=F, quote=F)
	}

	# --- VARIANZA DE LOS PÍXELES EXPLICADA POR LOS HFTs (R2) SEGÚN EL NÚMERO DE HFTs, EN CADA UMBRAL ---
	# Los segmentos de cada umbral clasificados (CLASS_METHOD) en 2, 3, ... HFTs; los puntos son los HFTs de
	# los cortes del umbral elegido. Muestra si, con el mismo número de HFTs, otro umbral recogería más varianza
	if (!is.null(ss_sel)) {
		kmax <- max(30, cut_r2$hfts + 5)
		curves <- NULL
		for (tk in thr_all) {
			hf <- paste('csv/park_', park, '_hri_results', tk, '.csv', sep='')
			if (!file.exists(hf) || !file.exists(file.path('shp', paste0('park_segm_', park, '_', tk, '_diss.shp')))) next
			hri_t <- if (tk == res) hri else read.table(hf, sep=' ', header=T)
			if (nrow(hri_t) < 3) next
			seg_t <- if (tk == res) seg_sel else segm_cells(park, tk, hri_t, pixgrid)
			ct_t <- if (tk == res) ct else class_tree(if (is.null(pixclip)) hri_t else clip_hri(hri_t, seg_t, pixclip, median_stat), class_method)
			ss_t <- if (tk == res) ss_sel else segm_sums(seg_t, pixvals, nrow(hri_t))
			nh <- 2:min(nrow(hri_t) - 1, kmax)
			r2m <- do.call(rbind, lapply(nh, function(h) r2_groups(ss_t, classify_k(ct_t, h, class_method))))
			colnames(r2m) <- var_names
			curves <- rbind(curves, data.frame(threshold=tk / 10, segments=nrow(hri_t),
			                                   r2_segments=mean(r2_groups(ss_t, seq_len(nrow(hri_t))), na.rm=TRUE),
			                                   hfts=nh, r2=rowMeans(r2m, na.rm=TRUE), r2m, check.names=FALSE))
		}
		# Corte sugerido (calculado antes de los cortes) y parte del R2 de los segmentos que conserva cada corte
		cut_r2$r2_share <- cut_r2$r2 / curves$r2_segments[curves$threshold == res / 10][1]
		cut_r2$suggested <- as.integer(cut_r2$cut %in% suggested)
		cut_r2 <- cut_r2[, c('cut', 'hfts', 'r2', 'r2_share', 'suggested', var_names)]
		write.table(curves, file.path(outroot, paste0('r2_', park, '_hfts.csv')), sep=' ', row.names=F, quote=F)
		write.table(cut_r2, file.path(outroot, paste0('r2_', park, '_cuts.csv')), sep=' ', row.names=F, quote=F)
		if (!is.na(knee)) {
			s <- cut_r2[cut_r2$suggested == 1, ]
			msg <- paste0('ID ', park, ': the R2 curve of the threshold 0.', res, ' bends at ', knee, ' HFTs',
			              if (nrow(s) == 1) sprintf('; suggested cut %s: %d HFTs, R2 %.2f (%.0f %% of the R2 of the segments)',
			                                        s$cut, s$hfts, s$r2, 100 * s$r2_share))
			print(msg)
			write(paste('Suggested cut:', msg), file.path(outroot, 'processing_metadata.txt'), append=TRUE)
		}

		# Umbrales en una rampa de azules (claro: 0.1, oscuro: 0.9) y el elegido en naranja y más grueso
		curves$thr <- factor(sprintf('%.1f', curves$threshold))
		sel <- sprintf('%.1f', res / 10)
		thr_cols <- setNames(colorRampPalette(c('#86b6ef', '#0d366b'))(nlevels(curves$thr)), levels(curves$thr))
		thr_cols[sel] <- '#eb6834'
		seg_info <- unique(curves[, c('thr', 'segments', 'r2_segments')])
		thr_labels <- setNames(paste0(seg_info$thr, ': ', seg_info$segments, ' segments, R2 ', sprintf('%.2f', seg_info$r2_segments),
		                              ifelse(seg_info$thr == sel, ' (selected)', '')), seg_info$thr)
		cut_pts <- data.frame(hfts=cut_r2$hfts, r2=cut_r2$r2,
		                      label=paste0(if (simil == 3) paste0(cut_r2$cut, ': ') else '', cut_r2$hfts, ' HFTs, R2 ', sprintf('%.2f', cut_r2$r2),
		                                   ifelse(cut_r2$suggested == 1, ' (suggested)', '')))
		r2p <- ggplot(curves, aes(x=hfts, y=r2, colour=thr)) +
			(if (!is.na(knee)) geom_vline(xintercept=knee, colour='grey45', linetype='dashed')) +
			(if (!is.na(knee)) annotate('text', x=knee, y=min(curves$r2), label=paste0(' bend of the selected curve: ', knee, ' HFTs'),
			                            hjust=0, vjust=0, size=3.5, colour='grey30')) +
			geom_line(data=curves[curves$thr != sel, ], linewidth=0.6) +
			geom_line(data=curves[curves$thr == sel, ], linewidth=1.5) +
			geom_point(data=cut_pts, aes(x=hfts, y=r2), inherit.aes=FALSE, colour='#eb6834', size=3.5) +
			ggrepel::geom_label_repel(data=cut_pts, aes(x=hfts, y=r2, label=label), inherit.aes=FALSE, size=3.5, seed=1,
			                          nudge_y=-0.04, min.segment.length=0, fill=alpha('white', 0.9)) +
			scale_colour_manual(values=thr_cols, labels=thr_labels, breaks=levels(curves$thr), name='Threshold (segments, R2 of the segments)') +
			scale_x_continuous(breaks=seq(0, max(curves$hfts), 5), minor_breaks=seq_len(max(curves$hfts))) +
			labs(title=area, subtitle=paste0('ID ', park, ' | pixel variance explained by the HFTs (R2, mean of the variables) when the segments of\n',
			                                 'each threshold are classified into 2 to ', max(curves$hfts), ' HFTs (', class_label, ')\n',
			                                 'selected threshold ', sel, ': ', method_label),
			     x='Number of HFTs', y='R2 of the HFTs') +
			theme_bw(base_size=13)
		ggsave(file.path(outroot, paste0('r2_', park, '_hfts.png')), r2p, width=11, height=7, dpi=300)
	}
}

run_state('done')
