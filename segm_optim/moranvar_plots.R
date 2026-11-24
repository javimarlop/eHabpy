# Modernized for R 4.x & Refactored for Dynamic Variables and PA Outline Overlays
library(vegan)
library(ade4)
library(ggplot2)
library(reshape2)
library(sf)
library(RColorBrewer)

source('plot_utils.R') # rescale01, short_labels, write_plot_csv, csv_label (shared with hft_compare.R)
source('subclasses.R') # hft_subclasses, subclass_colours (HFT_SUBCLASSES)

# Campos de las variables en el shapefile de los HFTs, de 10 caracteres como máximo, únicos y distintos de
# 'taken': los de abbreviate, como los acortaba sf (p. ej. G_1_NDV), y si alguno no cabe (nombres que solo se
# distinguen al final) V<n>_ y el principio de la parte del nombre que no comparten todas (n: orden en ENV_VARS)
dbf_names <- function(x, taken = character(0)) {
	a <- gsub('.', '_', abbreviate(x, minlength = 7, named = FALSE), fixed = TRUE)
	if (all(nchar(a) <= 10) && !anyDuplicated(c(a, taken))) return(a)
	p <- paste0('V', seq_along(x), '_')
	d <- abbreviate(gsub('\\s+', '_', short_labels(x, width = 1000)), minlength = 7, named = FALSE)
	paste0(p, substr(gsub('.', '_', d, fixed = TRUE), 1, 10 - nchar(p)))
}

# Fields of the shapefile of the HFTs (park_segm_<id>_<k>[_<cut>][_optim]_class.shp) besides the variables, and what
# each one is (in park_segm_..._class_fields.csv, next to it). Up to eHabpy 1.6.6 they were hclst_m (hft), cat
# (segm_cat), cat_pa and wdpaid (area_id) and label (always empty)
CLASS_FIELDS <- c(
	segm_id='ID of the segment: the ID of the study area, the position of the segmentation threshold and the number of the segment',
	segm_cat='Number of the segment in the segmentation (cat of shp/park_segm_<id>_<k>.shp)',
	area_id='ID of the study area (palist.csv)',
	hft='HFT of the segment in this cut',
	subclass='Subclass of the pixels of the HFT (HFT_SUBCLASSES, k-means): hft_<n>_sub_<y> is the subclass y of the HFT n; empty for the pixels without data in any variable',
	sub_code='Code of the subclass in the GeoTIFF hft_subclasses_*.tif and the table hft_subclasses_*.csv: 100 x HFT + subclass',
	area_ha='Area of the polygon (ha)',
	silhouette='Silhouette width of the segment between the HFTs (-1 to 1; NA with one HFT per segment)')

# Fields of the variables in the shapefile of the HFTs: the part of their names that differs (short_labels), without
# the number in front (as in gee_<area>_1_NDVI) and with _ for the other characters, cut to the 10 characters of a
# field of a .dbf (e.g. NDVI_S2_p9, CLre_lates); dbf_names if two of them would be equal
var_fields <- function(x, taken = names(CLASS_FIELDS)) {
	s <- gsub('_+', '_', gsub('[^A-Za-z0-9]', '_', short_labels(x, width = 1000)))
	s <- sub('^[0-9]+_(?=[A-Za-z])', '', sub('^_+', '', s), perl = TRUE)
	s <- sub('_+$', '', substr(s, 1, 10))
	if (all(nzchar(s)) && !any(grepl('^[0-9]', s)) && !anyDuplicated(tolower(c(s, taken)))) return(s)
	dbf_names(x, taken)
}

# The polygons of the segments with the fields of the shapefile of the HFTs: segs, the segments (sf, with cat and
# segm_id) with their HFT (hclust_mean); vars, the value of each variable for each segment (columns named by the
# fields of var_fields); sil, the silhouette width of each segment
class_layer <- function(segs, park, vars, sil) {
	ha <- st_area(st_geometry(segs))
	ha <- if (inherits(ha, 'units')) as.numeric(units::set_units(ha, 'ha', mode = 'standard')) else as.numeric(ha) / 1e4
	d <- data.frame(segm_id = segs$segm_id, segm_cat = as.integer(segs$cat), area_id = as.character(park),
	                hft = as.integer(segs$hclust_mean), area_ha = ha, silhouette = sil, stringsAsFactors = FALSE)
	st_sf(cbind(d, as.data.frame(vars, row.names = NULL, check.names = FALSE)), geometry = st_geometry(segs))
}

# Writes the shapefile of the HFTs and, next to it, park_segm_..._class_fields.csv: what each field is (vars: the
# name of the variable of each field of a variable)
write_class_layer <- function(x, path, vars) {
	st_write(x, path, delete_layer = TRUE, quiet = TRUE)
	f <- setdiff(names(x), attr(x, 'sf_column'))
	what <- ifelse(f %in% names(CLASS_FIELDS), CLASS_FIELDS[f],
	               paste0('Mean of ', vars[f], ' in the segments of the HFT, rescaled from 0 (the lowest HFT) to 1 (the ',
	                      'highest), as in the radar plot of the HFTs'))
	write.table(data.frame(field = f, meaning = unname(what)), sub('\\.shp$', '_fields.csv', path), sep = ' ',
	            row.names = FALSE, quote = TRUE)
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
	# Also what the plot draws for each HFT: its colour and horizontal position in each cut (row y of the cut: 1 at
	# the bottom, the most detailed one, as many as cuts at the top)
	for (lv in seq_len(nl)) paths[[paste0('colour_', cuts[lv])]] <- cols_list[[lv]][as.integer(as.character(paths[[cuts[lv]]]))]
	for (lv in seq_len(nl)) paths[[paste0('x_', cuts[lv])]] <- signif(as.vector(xs[[lv]][as.character(paths[[cuts[lv]]])]), 6)
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
# (n: values with data; outliers: values beyond the whiskers, not drawn; both only for the CSV of the plot)
box_stats <- function(x) {
	x <- x[!is.na(x)]
	if (length(x) == 0) return(c(ymin = NA, lower = NA, middle = NA, upper = NA, ymax = NA, n = 0, outliers = 0))
	q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
	iqr <- q[3] - q[1]
	lo <- q[1] - 1.5 * iqr
	hi <- q[3] + 1.5 * iqr
	c(ymin = min(x[x >= lo]), lower = q[1], middle = q[2], upper = q[3], ymax = max(x[x <= hi]),
	  n = length(x), outliers = sum(x < lo | x > hi))
}

# Statistics of the pixels of an HFT in a variable (original units), for the tables of the HFTs (hft_stats_*.csv):
# n: pixels with data; quartiles as in box_stats
pixel_stats <- function(x) {
	x <- x[!is.na(x)]
	if (length(x) == 0) return(c(n = 0, mean = NA, sd = NA, min = NA, q1 = NA, median = NA, q3 = NA, max = NA))
	q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
	c(n = length(x), mean = mean(x), sd = if (length(x) > 1) sd(x) else NA, min = min(x), q1 = q[1], median = q[2],
	  q3 = q[3], max = max(x))
}


# Distancias entre los segmentos del método original ('ward'): estadístico (media o mediana) y varianza
# de cada variable, reescalados a 0-1 entre los segmentos
segm_dist <- function(hri) {
	nv <- (ncol(hri) - 3) / 3
	vegdist(as.data.frame(lapply(hri[, 3:(2 + 2 * nv)], rescale01)), "euclidean", na.rm=T)
}

# Statistic (mean, or median with STAT_AGG = 'median') of each variable in each segment standardized with the
# standard deviation of the pixels (within + between segments), z, and the number of pixels of each segment, n
std_stats <- function(hri) {
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
	list(z=z, n=n)
}

# Árbol de la clasificación de los segmentos según CLASS_METHOD (ver config.py) y distancias entre ellos
# (para la correlación cofenética). Con 'ward_area' y 'ward_area_kmeans', Ward ponderado por el
# nº de píxeles de cada segmento sobre su estadístico de cada variable tipificado con la desviación típica
# de los píxeles (dentro + entre segmentos): con la distancia de Ward entre grupos y members = nº de
# píxeles, hclust da el mismo árbol que el Ward de todos los píxeles partiendo de los segmentos.
# z and n (std_stats) also with 'ward', for the silhouette of the HFTs (sil_width)
# With more than lowmem_segments() segments, the same tree from ward_lowmem.py (lowmem_tree), without the matrix of the
# distances: then dist is NULL, coph is the cophenetic correlation and lowmem is TRUE
class_tree <- function(hri, method, lowmem=nrow(hri) > lowmem_segments()) {
	s <- std_stats(hri)
	if (lowmem) {
		if (method == 'ward') {
			# the variables of segm_dist, without the distances: a segment without data in one takes its mean
			x <- sapply(hri[, 3:(2 + 2 * (ncol(hri) - 3) / 3)], function(v) {
				v <- rescale01(v)
				m <- mean(v, na.rm=TRUE)
				v[is.na(v)] <- if (is.finite(m)) m else 0
				v
			})
			lt <- lowmem_tree(x, rep(1, nrow(hri)))
		} else {
			lt <- lowmem_tree(s$z, s$n)
		}
		return(list(tree=lt$tree, dist=NULL, coph=lt$coph, z=s$z, n=s$n, lowmem=TRUE))
	}
	if (method == 'ward') {
		d <- segm_dist(hri)
		return(list(tree=hclust(d, "ward.D2"), dist=d, z=s$z, n=s$n, lowmem=FALSE))
	}
	n <- s$n
	d <- dist(s$z)
	w <- sqrt(2 * outer(n, n) / outer(n, n, '+')) * as.matrix(d)
	list(tree=hclust(as.dist(w), "ward.D2", members=n), dist=d, z=s$z, n=n, lowmem=FALSE)
}

# Clasificación con poca memoria. hclust, cophenetic y sil_width usan la matriz de las distancias entre cada par de
# segmentos (unos 70 bytes x segmentos², 200 GB con 55 000 segmentos); con más de LOWMEM_SEGMENTS segmentos en un
# umbral, ward_lowmem.py (Python con numba, el de EHAB_PYTHON) da los mismos resultados (salvo redondeo) calculando
# las distancias cuando las necesita: el árbol de Ward ponderado con la cadena de vecinos más próximos, la correlación
# cofenética y las siluetas. Por debajo, los de siempre. LOWMEM_SEGMENTS es el de hft_state.py (con él las matrices
# necesitan 3.9 GB como mucho); EHAB_LOWMEM_SEGMENTS lo cambia, para las pruebas
LOWMEM_SEGMENTS <- 7500
lowmem_segments <- function() {
	v <- suppressWarnings(as.numeric(Sys.getenv('EHAB_LOWMEM_SEGMENTS')))
	if (is.na(v)) LOWMEM_SEGMENTS else v
}

# Runs ward_lowmem.py CMD with the files of DIR (see its description) and the numbers ARGS (-s: without the packages of
# the user's site-packages, which could hide the numba of the environment, as in the app: PYTHONNOUSERSITE)
lowmem_run <- function(cmd, dir, args) {
	st <- system2(Sys.getenv('EHAB_PYTHON', 'python3'), c('-s', 'ward_lowmem.py', cmd, shQuote(dir), args))
	if (!identical(as.integer(st), 0L)) {
		stop(paste0('ward_lowmem.py ', cmd, ' failed (exit status ', st, '): see the messages above'))
	}
}
write_bin <- function(x, f, size=8) writeBin(if (size == 8) as.double(x) else as.integer(x), f, size=size, endian='little')
read_bin <- function(f, what, n, size=8) {
	v <- readBin(f, what, n + 1, size=size, endian='little')
	if (length(v) != n) stop(paste('ward_lowmem.py wrote', length(v), 'values in', basename(f), 'instead of', n))
	v
}

# Weighted Ward of the rows of x (segments) with the weights w (pixels) by ward_lowmem.py: the tree of
# hclust(method='ward.D2', members=w) on the distances sqrt(2 w_i w_j / (w_i + w_j)) ||x_i - x_j|| (as class_tree), as an
# hclust object with the labels of class_tree, and the cophenetic correlation with the Euclidean distances (coph)
lowmem_tree <- function(x, w) {
	x <- as.matrix(x)
	n <- nrow(x)
	dir <- tempfile('ward_')
	dir.create(dir)
	on.exit(unlink(dir, recursive=TRUE))
	write_bin(x, file.path(dir, 'x.f8'))
	write_bin(w, file.path(dir, 'w.f8'))
	t0 <- Sys.time()
	print(paste0('Classification of ', n, ' segments with little memory (ward_lowmem.py, more than ', lowmem_segments(),
	             ' segments)'))
	lowmem_run('tree', dir, c(n, ncol(x)))
	print(paste('Tree and cophenetic correlation in', format(round(Sys.time() - t0, 1))))
	tree <- structure(list(merge=matrix(read_bin(file.path(dir, 'merge.i4'), 'integer', 2 * (n - 1), 4), ncol=2),
	                       height=read_bin(file.path(dir, 'height.f8'), 'double', n - 1),
	                       order=read_bin(file.path(dir, 'order.i4'), 'integer', n, 4),
	                       labels=as.character(seq_len(n)), method='ward.D2', call=quote(ward_lowmem.py(x, members=w)),
	                       dist.method=NULL),
	                  class='hclust')
	list(tree=tree, coph=read_bin(file.path(dir, 'coph.f8'), 'double', 1))
}

# Silhouette widths (sil_width) of the segments in each classification of the list cls: with the matrix of the
# distances D, or, if ct$lowmem, by ward_lowmem.py with the distances between the rows of ct$z (those of D), all at once
sil_widths <- function(ct, D, cls) {
	if (!ct$lowmem) return(lapply(cls, function(cl) sil_width(D, cl)))
	ns <- nrow(ct$z)
	out <- lapply(cls, function(cl) list(s=rep(NA_real_, ns), neighbour=rep(NA, ns)))
	levs <- lapply(cls, function(cl) sort(unique(cl)))
	ok <- which(sapply(levs, function(l) length(l) >= 2 && length(l) < ns))
	if (length(ok) == 0) return(out)
	dir <- tempfile('sil_')
	dir.create(dir)
	on.exit(unlink(dir, recursive=TRUE))
	write_bin(ct$z, file.path(dir, 'x.f8'))
	write_bin(sapply(ok, function(i) match(cls[[i]], levs[[i]]) - 1L), file.path(dir, 'groups.i4'), 4)
	lowmem_run('silhouette', dir, c(ns, ncol(ct$z), length(ok)))
	s <- matrix(read_bin(file.path(dir, 's.f8'), 'double', ns * length(ok)), ncol=length(ok))
	nb <- matrix(read_bin(file.path(dir, 'neighbour.i4'), 'integer', ns * length(ok), 4), ncol=length(ok))
	for (j in seq_along(ok)) out[[ok[j]]] <- list(s=s[, j], neighbour=levs[[ok[j]]][nb[, j] + 1])
	out
}

# Silhouette width of each segment (Rousseeuw 1987, https://doi.org/10.1016/0377-0427(87)90125-7) in the HFTs cl,
# with D the matrix of distances between the segments: a = mean distance to the other segments of its HFT, b = the
# smallest mean distance to the segments of another HFT (its neighbour), s = (b - a) / max(a, b), from -1 to 1
# (negative: the segment is on average closer to the segments of its neighbour). A segment alone in its HFT gets 0,
# as in cluster::silhouette. All NA with one HFT or one HFT per segment. Returns s and the neighbour of each segment
sil_width <- function(D, cl) {
	lev <- sort(unique(cl))
	k <- length(lev)
	ns <- length(cl)
	if (k < 2 || k >= ns) return(list(s=rep(NA_real_, ns), neighbour=rep(NA, ns)))
	g <- match(cl, lev)
	own <- cbind(g, seq_len(ns))
	nh <- tabulate(g, k)
	md <- rowsum(D, g, reorder=TRUE) / nh # mean distance from each segment (column) to the segments of each HFT (row)
	a <- md[own] * nh[g] / pmax(nh[g] - 1, 1) # without the distance 0 to itself
	md[own] <- Inf
	nb <- max.col(-t(md), ties.method='first')
	b <- md[cbind(nb, seq_len(ns))]
	s <- ifelse(pmax(a, b) > 0, (b - a) / pmax(a, b), 0)
	s[nh[g] == 1] <- 0
	list(s=s, neighbour=lev[nb])
}

# Mean silhouette width of the segments weighted by their area (n: pixels), of all of them or of those of each HFT
sil_mean <- function(s, n) sum(s * n) / sum(n)
sil_hft <- function(s, n, cl) tapply(s * n, cl, sum) / tapply(n, cl, sum)

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

# segm_id as text, to match those of the shapefiles (text fields) and write them: as a number, R writes those of
# more than 15 digits rounded, and some others in scientific notation (as.character(1011000000) is "1.011e+09")
sid_txt <- function(x) if (is.numeric(x)) format(x, scientific=FALSE, trim=TRUE) else as.character(x)

# Statistics of the segments of a threshold (csv/park_<id>_hri_results<k>.csv of ehab_segm_core.py), with segm_id as
# text (sid_txt): its rows are sorted by segm_id as a number
read_hri <- function(f) read.table(f, sep=' ', header=TRUE, colClasses=c(segm_id='character'))

# Segmento (fila de hri) de cada píxel de la malla de las variables, del umbral tk (NA fuera de los segmentos)
segm_cells <- function(park, tk, hri, grid) {
	v <- terra::vect(file.path('shp', paste0('park_segm_', park, '_', tk, '_diss.shp')))
	v$row <- match(sid_txt(v$segm_id), sid_txt(hri$segm_id))
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
# 3: a step before is out of date, or an option of config.py is not valid (hft_state.py says which above)
if (st == 3) stop('moranvar_plots.R cannot run (see the message above)', call. = FALSE)
if (st != 0) warning('could not check the steps run before with hft_state.py (is python3 in the PATH?)')
progress <- function(kind, value) if (nzchar(Sys.getenv('EHAB_PROGRESS'))) cat('##', kind, ' ', value, '\n', sep='')  # for the desktop app (hft_state.progress)

run_state('start')

df0 <- read.table('csv/segm_done.csv', sep=' ', header=F)
ecox_list <- unique(df0[,1])
mx <- length(ecox_list)

# Thresholds of the segmentation (SEGM_THRESHOLDS): the value of every k, the index that names the files and is the
# column threshold of csv/<id>_threshold.csv, from csv/segm_thresholds.csv of segmentation_pca_par.py, as text for the
# labels (thr_lab) and as a number (thr_val); without it (a segmentation of eHabpy 1.5 or before) the nine of 0.1-0.9
thr_file <- 'csv/segm_thresholds.csv'
thr_text <- if (file.exists(thr_file)) {
	tt <- read.table(thr_file, sep=' ', header=TRUE, colClasses='character')
	setNames(tt$value, tt$k)
} else setNames(as.character((1:9) / 10), 1:9)
thr_lab <- function(k) unname(thr_text[as.character(k)])
thr_val <- function(k) as.numeric(thr_lab(k))
thr_of_value <- function(v) unname(thr_text[match(v, as.numeric(thr_text))]) # label of a threshold given by its value

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

thr_fixed <- !is.na(suppressWarnings(as.numeric(thr_method)))

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
# Subclases de los píxeles de cada HFT (HFT_SUBCLASSES en config.py, subclasses.R)
subclasses_on <- identical(read_config('HFT_SUBCLASSES'), 'True')

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

# Cortes del árbol: simil = 3 genera una clasificación por corte (SIMIL_CUTS, por defecto 10, 25, 50, 75 y 90 % de la
# altura máxima), cada una en su subcarpeta, de la más detallada (la de menor %, más HFTs) a la más general
if (simil == 3) {
	# SIMIL_CUTS en config.py: porcentajes de la altura máxima (1 a 99), de menos a más (más HFTs a menos)
	cut_txt <- read_config('SIMIL_CUTS')
	cut_pct <- if (is.na(cut_txt)) c(10, 25, 50, 75, 90) else sort(unique(as.integer(strsplit(gsub('[^0-9]+', ' ', cut_txt), ' ')[[1]][nzchar(strsplit(gsub('[^0-9]+', ' ', cut_txt), ' ')[[1]])])))
	if (length(cut_pct) == 0 || any(cut_pct < 1 | cut_pct > 99)) {
		stop("SIMIL_CUTS en config.py debe ser una lista de porcentajes entre 1 y 99, p. ej. [10, 25, 50, 75, 90]")
	}
	cuts <- setNames(cut_pct / 100, paste0('p', cut_pct))
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
	paste('Segmentation thresholds (k: value, csv/segm_thresholds.csv):', paste0(names(thr_text), ': ', thr_text, collapse=', ')),
	'',
	'## config.py',
	cfg,
	'',
	'## moranvar_plots.R',
	paste('upper =', upper, '# max. HFTs (SIMIL 1 and 2): 1: log10(n. pixels), 2: log6(n. pixels)')
), file.path(outroot, 'processing_metadata.txt'))
# The value of every threshold k of the files of the results (the desktop app shows it)
write.table(data.frame(k=names(thr_text), value=thr_text), file.path(outroot, 'segm_thresholds.csv'), sep=' ', row.names=F, quote=F)

progress('UNITS', mx * (length(cuts) + 3)) # per ID: thresholds, each cut, tree of the cuts and R2 curves
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

	# Con un umbral fijo solo se calculó ese umbral: no hay curvas entre umbrales que dibujar
	if (!thr_fixed) {
	# Curvas reescaladas (media de las variables) y puntuación del método, con el umbral elegido
	x1 <- thr_tab$moran_norm
	x2 <- thr_tab$var_norm
	thr_png <- file.path(outroot, paste(park, '_moran_var_mean.png', sep=''))
	png(thr_png, width=640, height=540)
	par(mar=c(4.5, 4.5, 4.5, 1))
	# one point per threshold k evaluated, equally spaced (the crossings are interpolated between consecutive ones),
	# labelled with its value
	plot(thr_all, x1, ylim=c(0, 1.35), col=1, typ='o', pch=19, xaxt='n', xlab='Similarity threshold',
	     ylab="Rescaled Moran's I, variance and score", main=paste0(park, ': selected threshold ', thr_lab(res)))
	axis(1, at=thr_all, labels=thr_lab(thr_all))
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
	# CSV of the plot: the columns of csv/<id>_threshold.csv (moranvar.py: threshold, the index k of each threshold,
	# its position on the x axis; moran_norm, black; var_norm, green; selected, dashed line), the value of each threshold
	# (threshold_value, the labels of the x axis) and what the plot computes from them: the score rescaled 0-1 with
	# 1 = best (blue) and, with THRESHOLD_METHOD 'crossing' or 'crossing_first', the crossings of the curves (red points)
	# between each threshold and the next one, as positions k on the x axis
	thr_csv <- data.frame(threshold=thr_all, threshold_value=thr_val(thr_all), thr_tab[names(thr_tab) != 'threshold'], check.names=FALSE)
	thr_csv$score_rescaled <- if (!all(is.na(thr_tab$score))) sc else NA
	thr_csv$crossing_threshold <- NA
	thr_csv$crossing_value <- NA
	if (thr_method %in% c('crossing', 'crossing_first') && length(ip) > 0) {
		thr_csv$crossing_threshold[ip] <- thr_all[ip] + tp * (thr_all[ip + 1] - thr_all[ip])
		thr_csv$crossing_value[ip] <- x1[ip] + (x1[ip + 1] - x1[ip]) * tp
	}
	write_plot_csv(thr_csv, thr_png)
	}

	progress('DONE', 1)
	k <- thr_lab(res) # value of the selected threshold (res: its index k)
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
	
	hri <- read_hri(namef)

	num_vars <- (ncol(hri) - 3) / 3
	maxnclas <- num_vars - 1

	# Valores de las variables en cada píxel (../inVars), para los boxplots por HFT
	var_names <- sub("pamean$", "", names(hri)[3:(2 + num_vars)])
	ev <- read_env_vars()
	ev <- ev[match(var_names, ev$var), ]
	# Name of each variable in the tables of the HFTs: that of its file without the folder (INVARS_FOLDER) and the
	# extension, with _ instead of other characters (as setup.py names the maps, without the folder), so that it is
	# the same in every study area whose variables have the same file names (e.g. those of get_vars_gee.py)
	var_keys <- gsub('[^A-Za-z0-9_]', '_', tools::file_path_sans_ext(basename(ev$file)))
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
	if (ct$lowmem) {
		if (is.finite(ct$coph)) cophval <- ct$coph
	} else {
		try(chcl <- cophenetic(hclust_mh))
		try(cophval <- cor(chcl, dmh))
	}
	# Distances for the silhouette of the HFTs: Euclidean between the standardized statistics z of the segments that
	# the classification used (with CLIP_PERCENTILE and STAT_AGG), with every CLASS_METHOD, also with 'ward', whose
	# tree uses other distances, so that the silhouettes of the methods are comparable (with little memory,
	# sil_widths computes them when needed)
	dsil <- if (ct$lowmem) NULL else as.matrix(dist(ct$z))
	sil_stat <- paste0(if (median_stat) 'medians' else 'means',
	                   if (!is.null(pixclip)) paste0('\n(variables clipped to percentiles ', clip_pct, '-', 100 - clip_pct, ', as for the classification)'))
	cut_sil <- c() # mean silhouette width of each cut

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

	# HFTs of each cut (cl) and the number of HFTs of the cut of the tree (n_cut, for the message), before the cuts so
	# that the silhouettes of all of them are computed at once (sil_widths)
	cut_hfts <- function(cn) {
		if (simil == 1) { hcut <- quantile(hclust_mh$height)[2] }
		if (simil == 2) { hcut <- min(hclust_mh$height) - 0.1 }
		if (simil == 3) { hcut <- cuts[[cn]] * max(hclust_mh$height) } # fracción de la altura máxima del árbol

		cutree(hclust_mh, h=hcut) -> hclust_mean
		ncl <- n_cut <- length(unique(hclust_mean))

		if (ncl == 1) {
			as.numeric(match(hri$segm_id, unique(hri$segm_id))) -> hclust_mean # one HFT per segment, in the order of the rows
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
		}
		list(cl=hclust_mean, n_cut=n_cut)
	}
	cut_cls <- lapply(setNames(names(cuts), names(cuts)), cut_hfts)
	cut_sils <- sil_widths(ct, dsil, lapply(cut_cls, `[[`, 'cl'))

	hft_rows <- list() # HFTs of each cut and statistics of their pixels, for the tables of the HFTs
	cut_sub <- cut_r2sub <- c() # subclasses of each cut (HFT_SUBCLASSES) and the part of the pixel variance they explain
	for (cn in names(cuts)) {
		progress('DONE', 1)
		outdir <- if (simil == 3) file.path(outroot, cn) else outroot
		# el corte también va en el nombre de cada fichero, con _optim si es el sugerido
		tag <- if (simil == 3) paste0('_', cn, if (identical(cn, suggested)) '_optim') else ''
		dir.create(outdir, showWarnings = FALSE)

		print(paste('Number of HFTs is ', cut_cls[[cn]]$n_cut, if (simil == 3) paste0(' (', cn, ')'), sep=''))
		hclust_mean <- cut_cls[[cn]]$cl
		ncl <- length(unique(hclust_mean))

		# Silhouette width of the segments in the HFTs of this cut (sil_width), its mean weighted by area and that of
		# each HFT
		sil <- cut_sils[[cn]]
		cut_sil[cn] <- sil_mean(sil$s, ct$n)
		sil_h <- sil_hft(sil$s, ct$n, hclust_mean)
		if (!is.na(cut_sil[cn])) print(sprintf('Mean silhouette width of the segments in the %d HFTs%s (weighted by area): %.2f',
		                                       ncl, if (simil == 3) paste0(' of ', cn) else '', cut_sil[[cn]]))

		# Varianza de los píxeles explicada por los HFTs de este corte (media y por variable)
		if (!is.null(ss_sel)) {
			r2v <- r2_groups(ss_sel, hclust_mean)
			cut_r2 <- rbind(cut_r2, data.frame(cut=cn, hfts=ncl, r2=mean(r2v, na.rm=TRUE), silhouette=cut_sil[[cn]],
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

			hclust_png <- file.path(outdir, paste('hclust_', park, '_', res, tag, '_segms_mean.png', sep=''))
			png(hclust_png)
			try(plot(hclust_mh, hang=-1, main=park, sub=cophval, xlab=class_label)); try(rect.hclust(hclust_mh, k=ncl))
			dev.off()
			# CSV of the dendrogram: rows 'leaf', the segments from left to right (leaf: their label in the plot, the
			# row of the segment in csv/park_<id>_hri_results<k>.csv; HFT; box: the rectangle that encloses them, the
			# tree cut into as many groups as HFTs, which can differ from the HFTs with 'ward_area_kmeans'), and rows
			# 'merge', the merges of hclust
			# (merged_a, merged_b: a negative number -j is the leaf j, a positive one the merge of that step) with their
			# height; cophenetic_r: the cophenetic correlation of the subtitle
			nseg <- length(hclust_mh$order)
			leaf_lab <- if (is.null(hclust_mh$labels)) seq_len(nseg) else hclust_mh$labels
			ho <- hclust_mh$order
			write_plot_csv(rbind(
				data.frame(element='leaf', step=NA, merged_a=NA, merged_b=NA, height=NA, leaf=leaf_lab[ho], segm_id=hri$segm_id[ho],
				           position=seq_len(nseg), HFT=hclust_mean[ho], box=cutree(hclust_mh, k=ncl)[ho], cophenetic_r=cophval),
				data.frame(element='merge', step=seq_len(nseg - 1), merged_a=hclust_mh$merge[, 1], merged_b=hclust_mh$merge[, 2],
				           height=hclust_mh$height, leaf=NA, segm_id=NA, position=NA, HFT=NA, box=NA, cophenetic_r=cophval)),
			           hclust_png)

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
			# CSV of the radar plot: per HFT (in the order of the legend) and variable (axis 1 at the top, then
			# clockwise), value: the mean over the segments of the HFT of their statistic (STAT_AGG: mean or median
			# of the pixels); relative: that value rescaled 0-1 among the HFTs, the one drawn
			rad_o <- match(hft_order, hri4$hclust_mean)
			rad_o <- rad_o[!is.na(rad_o)]
			rad_h <- hri4$hclust_mean[rad_o]
			write_plot_csv(data.frame(HFT=rep(rad_h, each=num_vars), colour=rep(kolsn[rad_h], each=num_vars),
			                          var=rep(clean_names, length(rad_h)), var_label=rep(csv_label(short_labels(clean_names)), length(rad_h)),
			                          axis=rep(seq_len(num_vars), length(rad_h)),
			                          value=as.vector(t(as.matrix(hri4[rad_o, 2:(num_vars + 1), drop=FALSE]))),
			                          relative=as.vector(t(as.matrix(scaled[rad_o, 1:num_vars, drop=FALSE])))), rpn)

			# --- SILHOUETTE OF THE SEGMENTS BY HFT ---
			# Silhouette width of each segment (sil_width) as a bar from 0, as thick as the area of the segment, so that
			# the mean weighted by area (dashed line) is the visual mean of the bars: the HFTs from top to bottom in the
			# order of the legends, and in each one its segments from the widest to the narrowest. Each HFT has a band as
			# tall as its bars or as its label if they are thinner (lab_in inches, with 5 in for all the bars) plus a gap,
			# with the bars and the label in its middle, and the labels left of -1, so that a label is always level with
			# its bars and does not cover any bar
			if (!is.na(cut_sil[cn])) {
				sil_ha <- ct$n * cell_ha
				hs <- hft_order[hft_order %in% hclust_mean]
				so <- order(match(hclust_mean, hs), -sil$s, seq_len(nrow(hri))) # rows sorted by segm_id
				g <- match(hclust_mean[so], hs)
				hft_ha <- as.vector(rowsum(sil_ha[so], g))
				ha_in <- 5 / sum(sil_ha) # inches per ha
				lab_in <- 0.24
				band <- pmax(hft_ha, lab_in / ha_in) + 0.06 / ha_in
				mid <- cumsum(band) - band / 2
				ymin <- mid[g] - hft_ha[g] / 2 + cumsum(sil_ha[so]) - sil_ha[so] - c(0, cumsum(hft_ha))[g]
				sil_df <- data.frame(segm_id=hri$segm_id[so], HFT=hclust_mean[so], colour=kolsn[hclust_mean[so]], area_ha=sil_ha[so],
				                     silhouette=sil$s[so], neighbour=sil$neighbour[so], ymin=ymin, ymax=ymin + sil_ha[so], label_y=mid[g])
				sil_n <- table(factor(hclust_mean, levels=hs))
				# Label of each HFT next to its bars, with its colour (number in white on the dark ones)
				sil_lab <- data.frame(y=mid, label=sprintf('HFT %d: %.2f (%d %s)', hs, sil_h[as.character(hs)], sil_n, ifelse(sil_n == 1, 'segment', 'segments')),
				                      fill=kolsn[hs])
				sil_lab$txt <- ifelse(colSums(col2rgb(sil_lab$fill) / 255 * c(0.299, 0.587, 0.114)) < 0.55, 'white', 'black')
				# Labels without border (linewidth in recent versions of ggplot2, label.size in older ones)
				no_border <- if ('linewidth' %in% names(GeomLabel$default_aes)) list(linewidth=0) else list(label.size=0)
				sil_sub <- paste0(plot_subtitle, '\n',
				                  'silhouette width of the segments (Rousseeuw 1987) on the distances between their standardized ', sil_stat, '\n',
				                  sprintf('mean %.2f (dashed line, segments weighted by area) | bars as thick as the area of the segment\n', cut_sil[[cn]]),
				                  'negative: the segment is on average closer to the segments of another HFT (its neighbour)')
				silp <- ggplot(sil_df) +
					geom_rect(aes(xmin=0, xmax=silhouette, ymin=ymin, ymax=ymax, fill=colour)) +
					geom_vline(xintercept=0, colour='grey30', linewidth=0.3) +
					geom_vline(xintercept=cut_sil[[cn]], colour='grey15', linetype='dashed') +
					annotate('text', x=cut_sil[[cn]], y=-0.5 * lab_in / ha_in, label=sprintf(' mean %.2f', cut_sil[[cn]]), hjust=0, size=3.8) +
					do.call(geom_label, c(list(data=sil_lab, mapping=aes(x=-1.03, y=y, label=label, fill=fill, colour=txt), hjust=1, size=3.5), no_border)) +
					scale_fill_identity() + scale_colour_identity() +
					scale_x_continuous(limits=c(-1.6, 1), breaks=seq(-1, 1, 0.25), minor_breaks=seq(-1, 1, 0.125), expand=expansion(mult=c(0, 0.05))) +
					scale_y_reverse(expand=expansion(add=0.15 / ha_in)) +
					labs(title=area, subtitle=sil_sub, x='Silhouette width', y='Segments by HFT') +
					theme_bw(base_size=13) +
					theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank(),
					      panel.grid.minor.y=element_blank())
				spn <- file.path(outdir, paste('silhouette_', park, '_', res, tag, '_segms.png', sep=''))
				# Height: that of the panel (the bands, the mean label and the margins) plus title, subtitle and x axis
				ggsave(spn, silp, width=11, height=(sum(band) + 0.5 * lab_in / ha_in) * ha_in + 0.3 + 0.9 + 0.22 * length(strsplit(sil_sub, '\n')[[1]]), dpi=300)
				# CSV of the plot: one row per segment from top to bottom, with its HFT and its colour, area, silhouette
				# width, neighbour (the other HFT whose segments are closest on average), bar (ymin, ymax: from the top, in
				# ha, the bands of the HFTs above it plus the segments above it in its band) and label_y, the y of the label
				# of its HFT (middle of the band and of its bars); hft_silhouette: the mean of its HFT (label),
				# mean_silhouette: that of all the segments (dashed line), both weighted by area
				write_plot_csv(cbind(sil_df, hft_silhouette=as.vector(sil_h[as.character(sil_df$HFT)]), mean_silhouette=cut_sil[[cn]]), spn)
			}

			# Integración en Shapefile y Mapa final
			shp_layer <- paste('park_segm_', park, '_', res, '_diss', sep='')
			if (file.exists(file.path('shp', paste0(shp_layer, '.shp')))) {
				segm_pa <- st_read('shp', layer=shp_layer, quiet=TRUE)
				segm_pa$segm_id <- sid_txt(segm_pa$segm_id)

				merge(segm_pa, data.frame(segm_id=hri$segm_id, hclust_mean=hclust_mean), by='segm_id') -> segm_pa_class

				# The shapefile of the HFTs: one polygon per segment, with its HFT, the mean of each variable in the HFT and
				# the silhouette width of the segment (class_layer); with HFT_SUBCLASSES it is written again below, split
				# into the subclasses of the pixels of each segment
				class_shp <- file.path(outdir, paste('park_segm_', park, '_', res, tag, '_class.shp', sep=''))
				class_var_names <- setNames(names(scaled2)[-1], var_fields(names(scaled2)[-1]))
				class_vars <- setNames(scaled2[match(segm_pa_class$hclust_mean, scaled2$group), -1, drop=FALSE], names(class_var_names))
				segm_pa_class2 <- class_layer(segm_pa_class, park, class_vars, sil$s[match(segm_pa_class$segm_id, hri$segm_id)])
				write_class_layer(segm_pa_class2, class_shp, class_var_names)
				# Color de cada HFT en los mapas y gráficas, en el orden de las leyendas (mapa interactivo de la app)
				write.table(data.frame(HFT=hft_order, colour=kolsn[hft_order]),
				            file.path(outdir, paste('hft_colours_', park, '_', res, tag, '.csv', sep='')),
				            sep=' ', row.names=F, quote=F)

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
				# CSV of the map: one row per polygon of the segments (their geometry is the shapefile _class.shp of
				# this folder; the outline is the study area of 'pa'), in the order of the legend, with its HFT, area,
				# the colour of its HFT and its centroid (x, y in the CRS of the map, that of the shapefile)
				map_g <- st_geometry(segm_pa_class)
				map_ha <- st_area(map_g)
				map_ha <- if (inherits(map_ha, 'units')) as.numeric(units::set_units(map_ha, 'ha', mode='standard')) else as.numeric(map_ha) / 1e4
				map_xy <- st_coordinates(suppressWarnings(st_centroid(map_g)))
				map_csv <- data.frame(segm_id=segm_pa_class$segm_id, HFT=segm_pa_class$hclust_mean, area_ha=map_ha,
				                      colour=kolsn[segm_pa_class$hclust_mean], x=map_xy[, 'X'], y=map_xy[, 'Y'])
				write_plot_csv(map_csv[order(match(map_csv$HFT, hfts), map_csv$segm_id), ], rpn2, precise=c('x', 'y'))

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
					# CSV of the boxplots: per variable (panel) and HFT, the statistics of its pixels that the box draws
					# (n: pixels with data; whiskers: the furthest values within 1.5 x IQR of the box; outliers: pixels
					# beyond them, not drawn) and the area of the HFT of the legend
					box_h <- as.integer(as.character(bstats$HFT))
					write_plot_csv(data.frame(var=var_names[match(as.character(bstats$var), var_labels)], var_label=csv_label(as.character(bstats$var)),
					                          HFT=box_h, area_ha=as.vector(ha[as.character(box_h)]),
					                          area_pct=as.vector(100 * ha[as.character(box_h)] / sum(ha)), colour=kolsn[box_h],
					                          n=bstats$n, lower_whisker=bstats$ymin, q1=bstats$lower, median=bstats$middle, q3=bstats$upper,
					                          upper_whisker=bstats$ymax, outliers=bstats$outliers), rpn3)

					# HFTs of this cut (segments, pixels and area) and statistics of their pixels in each variable,
					# written after the last cut to the tables of the HFTs (hft_means_ and hft_stats_<id>_<k>.csv)
					hft_cut <- if (simil == 3) cn else NA
					hft_sug <- if (simil == 3) as.integer(identical(cn, suggested)) else NA
					hft_rows[[cn]] <- list(
						info=data.frame(cut=hft_cut, suggested=hft_sug, hfts=ncl, HFT=hfts, colour=kolsn[hfts],
						                segments=as.integer(table(factor(hclust_mean, levels=hfts))), pixels=as.integer(table(hftcls)),
						                silhouette=as.vector(sil_h[as.character(hfts)]), area_ha=as.vector(ha), area_pct=as.vector(100 * ha / sum(ha))),
						stats=do.call(rbind, lapply(seq_along(var_names), function(i) {
							s <- t(sapply(split(pixvals[inhft, i], hftcls), pixel_stats))
							data.frame(cut=hft_cut, suggested=hft_sug, HFT=hfts, area_ha=as.vector(ha), variable=var_keys[i],
							           map=var_names[i], s, row.names=NULL)
						})))

					# --- SUBCLASES DE LOS HFTs (HFT_SUBCLASSES, subclasses.R) ---
					# Los píxeles de cada HFT divididos por k-means en dos, una y otra vez, mientras una división los
					# separe al menos tanto como están separados los dos HFTs más parecidos: GeoTIFF con la subclase de
					# cada píxel (código 100 x HFT + subclase), su tabla y su mapa
					if (subclasses_on && ncl > 1 && ncl < nrow(hri)) { # not when every segment is its own HFT (a cut of one HFT)
						sc <- hft_subclasses(pixvals[inhft, , drop=FALSE], cls[inhft])
						scode <- rep(NA_integer_, length(cls))
						scode[which(inhft)] <- sc$code
						codes <- sort(unique(sc$code[!is.na(sc$code)]))
						codes <- codes[order(match(codes %/% 100L, hfts), codes %% 100L)] # in the order of the legend of the HFTs
						sh <- codes %/% 100L
						sub_lab <- paste0(sh, '.', codes %% 100L)
						scl <- factor(scode[inhft], levels=codes)
						s_ha <- tapply(pixarea[inhft], scl, sum)
						s_ha[is.na(s_ha)] <- 0
						s_px <- as.vector(table(scl))
						s_col <- subclass_colours(codes, s_px, setNames(kolsn[hfts], hfts))
						h_ha <- tapply(as.vector(s_ha), sh, sum)
						cut_sub[cn] <- length(codes)
						cut_r2sub[cn] <- sc$r2_subclasses
						print(sprintf('%d subclasses in the %d HFTs%s: they explain %.0f %% of the variance of the pixels (the HFTs %.0f %%); a division had to reduce the sum of squares of the pixels by %.4g or more',
						              length(codes), ncl, if (simil == 3) paste0(' of ', cn) else '', 100 * sc$r2_subclasses, 100 * sc$r2_hfts, sc$threshold))
						# GeoTIFF of the subclasses, on the grid of the variables, cropped to the pixels of the HFTs
						sr <- terra::rast(pixgrid, nlyrs=1)
						terra::values(sr) <- scode
						sr <- terra::trim(sr)
						names(sr) <- 'subclass'
						terra::writeRaster(sr, file.path(outdir, paste0('hft_subclasses_', park, '_', res, tag, '.tif')), overwrite=TRUE,
						                   datatype='INT4S', NAflag=0, gdal='COMPRESS=LZW')
						# The shapefile of the HFTs again, its segments split into the subclasses of their pixels (the parts of a
						# segment in the same subclass are one polygon, even if they do not touch), and the colours of the
						# subclasses (code colour), for the similarity of the subclasses
						write_class_layer(subclass_pieces(segm_pa_class2, pixgrid, scode), class_shp, class_var_names)
						write.table(data.frame(code=codes, colour=s_col), file.path(outdir, paste0('subclass_colours_', park, '_', res, tag, '.csv')),
						            sep=' ', row.names=F, quote=F)
						# Table of the subclasses: as hft_means_ (label, area and the mean of their pixels in each variable,
						# after area_pct), with the HFT, its code in the GeoTIFF, its part of the area of its HFT and the
						# least reduction of the sum of squares of a division (division_threshold)
						sub_tab <- data.frame(label=paste0(area, '/', scenario, '/', park, '/', if (simil == 3) paste0(cn, '/'), 'HFT', sub_lab),
						                      study_area=area, scenario=scenario, id=park, resolution=terra::res(pixgrid)[1], threshold=thr_val(res),
						                      cut=hft_cut, suggested=hft_sug, HFT=sh, subclass=sub_lab, code=codes, colour=s_col, pixels=s_px,
						                      pct_of_HFT=100 * as.vector(s_ha) / as.vector(h_ha[as.character(sh)]), division_threshold=sc$threshold,
						                      area_ha=as.vector(s_ha), area_pct=100 * as.vector(s_ha) / sum(s_ha), check.names=FALSE)
						for (i in seq_along(var_names)) {
							sub_tab[[var_keys[i]]] <- as.vector(tapply(pixvals[inhft, i], scl, mean, na.rm=TRUE))
						}
						write_plot_csv(sub_tab, file.path(outdir, paste0('hft_subclasses_', gsub('[^A-Za-z0-9_-]', '_', area), '_', scenario, '_',
						                                                   park, '_', res, tag, '.csv')))
						# Map of the subclasses (at most 2000 pixels wide: the most frequent subclass of each block of
						# pixels) with the outline of the study area
						fact <- max(1, ceiling(max(dim(sr)[1:2]) / 2000))
						spl <- if (fact > 1) terra::aggregate(sr, fact, fun='modal', na.rm=TRUE) else sr
						sdf <- terra::as.data.frame(spl, xy=TRUE, na.rm=TRUE)
						sdf$sub <- factor(sdf$subclass, levels=codes)
						# Size: 12 inches wide (and the legend), the height of the area (5 to 16 inches), taller if the legend
						# needs it to fit in 3 columns (rows of about 0.27 inches); without legend if it would need more than 4
						ex <- terra::ext(sr)
						s_asp <- as.numeric((ex$ymax - ex$ymin) / (ex$xmax - ex$xmin))
						s_h <- min(max(10 * s_asp + 1.5, 5, 1.3 + 0.27 * ceiling(length(codes) / 3)), 16)
						s_ncol <- ceiling(length(codes) / max(1, floor((s_h - 1.3) / 0.27)))
						s_leg <- s_ncol <= 4
						submap <- ggplot() +
							geom_raster(data=sdf, aes(x=x, y=y, fill=sub)) +
							scale_fill_manual(values=setNames(s_col, codes), labels=setNames(paste0(sub_lab, ': ', format(round(as.vector(s_ha)), big.mark=',', trim=TRUE), ' ha'), codes),
							                  name='Subclass (area)', drop=FALSE, guide=if (s_leg) guide_legend(ncol=s_ncol) else 'none') +
							labs(title=area, subtitle=paste0(plot_subtitle, ' | ', length(codes), ' subclasses of the pixels of the HFTs (k-means), in shades of the colour of their HFT',
							                                 if (!s_leg) ' (legend in the .csv)')) +
							theme_bw(base_size=14) +
							theme(axis.title=element_blank())
						if (!is.null(pa_outline)) {
							po <- pa_outline
							if (!is.na(st_crs(po)) && nzchar(terra::crs(pixgrid)) && st_crs(po) != st_crs(terra::crs(pixgrid))) po <- st_transform(po, terra::crs(pixgrid))
							crd <- as.data.frame(st_coordinates(st_geometry(po)))
							crd$grp <- interaction(crd[grep('^L', names(crd))], drop=TRUE)
							submap <- submap + geom_path(data=crd, aes(x=X, y=Y, group=grp), colour='black', linewidth=0.8)
						}
						submap <- submap + coord_equal(expand=FALSE)
						rpn6 <- file.path(outdir, paste('map_', park, '_', res, tag, '_subclasses.png', sep=''))
						ggsave(rpn6, submap, width=12 + if (s_leg) 2.5 * s_ncol else 0, height=s_h, dpi=300)
						# CSV of the map: one row per subclass (its pixels are those of the GeoTIFF hft_subclasses_ of this
						# folder), in the order of the legend, with its HFT, code, area, colour and the centroid of its pixels
						# (x, y in the CRS of the variables)
						sxy <- terra::xyFromCell(pixgrid, which(inhft))
						write_plot_csv(data.frame(subclass=sub_lab, HFT=sh, code=codes, area_ha=as.vector(s_ha), colour=s_col,
						                          x=as.vector(tapply(sxy[, 1], scl, mean)), y=as.vector(tapply(sxy[, 2], scl, mean))),
						               rpn6, precise=c('x', 'y'))
					}

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
					# CSV of the heatmap: per HFT (row 1 at the top, with its area and colour) and variable (column 1 on
					# the left, sorted by eta2), the median of its pixels (the number in the cell), the mean and sd of the
					# variable in the pixels of all the HFTs, z = (median - mean) / sd (the colour, on a scale from
					# -z_limit to z_limit) and eta2: part of the variance of the variable explained by the HFTs
					heat_h <- as.integer(as.character(bstats$HFT))
					heat_csv <- data.frame(HFT=heat_h, row=match(heat_h, hord), area_ha=as.vector(ha[as.character(heat_h)]),
					                       area_pct=as.vector(100 * ha[as.character(heat_h)] / sum(ha)), colour=kolsn[heat_h],
					                       var=var_names[vi], var_label=csv_label(var_labels[vi]), column=heat$x, eta2=eta2[vi],
					                       median=heat$median, mean=as.vector(pmean[vi]), sd=as.vector(psd[vi]), z=heat$z, z_limit=zlim)
					write_plot_csv(heat_csv[order(heat_csv$row, heat_csv$column), ], rpn4)

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
							# (pc_y, n, ell_*: for the CSV of the plot) semi-axes and angle of the ellipse
							cents <- rbind(cents, data.frame(panel=panels[a], HFT=h, x=mu[1], y=mu[2], ha=ha[as.character(h)], colour=kolsn[h],
							                                 pc_y=axes2[a], n=nrow(s), ell_a=sqrt(max(e$values[1], 0) * qchisq(0.5, 2)),
							                                 ell_b=sqrt(max(e$values[2], 0) * qchisq(0.5, 2)),
							                                 ell_angle=(atan2(e$vectors[2, 1], e$vectors[1, 1]) * 180 / pi) %% 180))
							ells <- rbind(ells, data.frame(panel=panels[a], HFT=h, x=xy[, 1], y=xy[, 2], colour=kolsn[h]))
						}
						arrs <- rbind(arrs, data.frame(panel=panels[a], var=var_labels, x=loads[, 1], y=loads[, axes2[a]],
						                               pc_y=axes2[a], var_name=var_names, loading_x=loads[, 1], loading_y=loads[, axes2[a]]))
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
					# CSV of the biplot, one row per point or vector of each panel (axes pc_x and pc_y, with the % of
					# the variance that they explain): 'centroid' rows, the HFTs (x, y: mean of the scores of their
					# n_pixels pixels; area_ha: size of the point; ellipse_a, ellipse_b: semi-axes, ellipse_angle:
					# degrees from the x axis, of the ellipse with half of the pixels, NA when not drawn, with more than
					# 10 HFTs); 'variable' rows, the vectors (x, y: tip of the arrow as drawn; loading_x, loading_y: the
					# correlations of the variable with the axes, before scaling them to the cloud of the HFTs;
					# label_x, label_y: where the label points). The pixels themselves are not drawn
					pca_na <- rep(NA, nrow(arrs))
					pca_csv <- rbind(
						data.frame(element='centroid', pc_x=1, pc_y=cents$pc_y, variance_x=ve[1], variance_y=ve[cents$pc_y],
						           HFT=cents$HFT, colour=cents$colour, n_pixels=cents$n, area_ha=cents$ha, x=cents$x, y=cents$y,
						           ellipse_a=if (show_ell) cents$ell_a else NA, ellipse_b=if (show_ell) cents$ell_b else NA,
						           ellipse_angle=if (show_ell) cents$ell_angle else NA,
						           var=NA, var_label=NA, loading_x=NA, loading_y=NA, label_x=NA, label_y=NA),
						data.frame(element='variable', pc_x=1, pc_y=arrs$pc_y, variance_x=ve[1], variance_y=ve[arrs$pc_y],
						           HFT=pca_na, colour=pca_na, n_pixels=pca_na, area_ha=pca_na, x=arrs$x, y=arrs$y,
						           ellipse_a=pca_na, ellipse_b=pca_na, ellipse_angle=pca_na,
						           var=arrs$var_name, var_label=csv_label(arrs$var), loading_x=arrs$loading_x, loading_y=arrs$loading_y,
						           label_x=arrs$lx, label_y=arrs$ly))
					write_plot_csv(pca_csv[order(pca_csv$pc_y), ], rpn5)
				}
			}
		}
	}

	# --- TABLES OF THE HFTs, TO COMPARE THEM WITH THOSE OF OTHER STUDY AREAS ---
	# The HFTs of every cut, one table per ID for the whole scenario. Those of several study areas (or scenarios)
	# made with the same variables can be joined by rows for multivariate comparisons (an NMDS or a cluster analysis
	# of the HFTs of all of them): the variables are named after their files (var_keys), the values are in the
	# original units of the variables (standardize them after joining the tables, not in each area) and 'label'
	# identifies each HFT among all the areas (study area/scenario/ID/cut/HFT<n>). Compare only HFTs made with the
	# same pixel size (column resolution).
	# - hft_means_<area>_<scenario>_<id>_<k>.csv (and, with SIMIL = 3, one per cut in its folder, with _<cut>[_optim]
	#   added to the name): one row per HFT and cut with its segments, pixels, silhouette (mean silhouette width of
	#   its segments weighted by area, see sil_width) and area, and the mean of its pixels in each variable (a column
	#   per variable after area_pct): the matrix of HFTs x variables
	# - hft_stats_<area>_<scenario>_<id>_<k>.csv (and one per cut): one row per HFT, cut and variable with the statistics of its pixels (n: pixels with
	#   data, mean, sd, min, quartiles, max), e.g. for distances between HFTs that take their spread into account
	if (length(hft_rows)) {
		# resolution: pixel size of the variables whose pixels give the statistics (units of their CRS, projected: m);
		# HFTs made with other pixel sizes are not directly comparable (the variables and the segments change with it)
		hft_id <- function(d) data.frame(label=paste0(area, '/', scenario, '/', park, '/', if (simil == 3) paste0(d$cut, '/'), 'HFT', d$HFT),
		                                 study_area=area, scenario=scenario, id=park, resolution=terra::res(pixgrid)[1],
		                                 threshold=thr_val(res), d, check.names=FALSE)
		info <- do.call(rbind, lapply(hft_rows, `[[`, 'info'))
		info <- info[order(match(info$cut, names(cuts)), info$HFT), ]
		hst <- do.call(rbind, lapply(hft_rows, `[[`, 'stats'))
		hst <- hst[order(match(hst$cut, names(cuts)), hst$HFT, match(hst$map, var_names)), ]
		hst$n <- as.integer(hst$n)
		means <- hft_id(info)
		for (j in seq_along(var_names)) {
			s <- hst[hst$map == var_names[j], ]
			means[[var_keys[j]]] <- s$mean[match(paste(info$cut, info$HFT), paste(s$cut, s$HFT))]
		}
		# The names carry the study area, the scenario, the ID, the threshold and the cut, so that the tables of different
		# runs never have the same name and can be put together in a folder (the comparison folder of the app) as they are
		tname <- function(kind, tag='') paste0(kind, '_', gsub('[^A-Za-z0-9_-]', '_', area), '_', scenario, '_', park, '_', res, tag, '.csv')
		hst_all <- hft_id(hst)
		write_plot_csv(means, file.path(outroot, tname('hft_means')))
		write_plot_csv(hst_all, file.path(outroot, tname('hft_stats')))
		# With SIMIL = 3 also one table per cut, in the folder of the cut, with the cut in its name
		if (simil == 3) for (cn in unique(means$cut)) {
			ctag <- paste0('_', cn, if (identical(cn, suggested)) '_optim')
			write_plot_csv(means[means$cut %in% cn, ], file.path(outroot, cn, tname('hft_means', ctag)))
			write_plot_csv(hst_all[hst_all$cut %in% cn, ], file.path(outroot, cn, tname('hft_stats', ctag)))
		}
	}

	progress('DONE', 1)
	# --- ÁRBOL DE LOS HFTs DE LOS CORTES (SIMIL = 3) ---
	# En qué HFT de los cortes más generales está cada HFT y en qué HFTs de los más detallados se divide (ver
	# hft_tree). La tabla tiene una fila por HFT del corte más detallado (o por combinación de HFTs, si los cortes
	# no están anidados) con el HFT que la contiene en cada corte, su nº de segmentos y su superficie
	if (simil == 3) {
		rows <- rev(names(cuts)) # del corte más general al más detallado
		nk_tree <- sapply(tree_cls[rows], function(cl) length(unique(cl)))
		r2_tree <- if (!is.null(cut_r2)) ifelse(nk_tree == 1, 0, cut_r2$r2[match(rows, cut_r2$cut)])
		sil_tree <- ifelse(nk_tree == 1, NA, cut_sil[rows])
		tree_labels <- paste0(rows, '\n', nk_tree, ifelse(nk_tree == 1, ' HFT', ' HFTs'), if (!is.null(r2_tree)) sprintf('\nR2 %.2f', r2_tree),
		                      ifelse(is.na(sil_tree), '', sprintf('\nsilhouette %.2f', sil_tree)), ifelse(rows %in% suggested, '\n(suggested)', ''))
		ht <- hft_tree(hclust_mh, tree_cls[rows], tree_cols[rows], hri$sumpamask * cell_ha, tree_labels, suggested)
		treep <- ht$plot + labs(title=area, subtitle=paste0('ID ', park, ' | threshold ', k,
			' | HFTs of the cuts of the tree, from broad (top) to detailed (bottom)\n',
			'circles: HFTs with their number and colour of the maps, sized by area | lines: area shared by the HFTs of consecutive cuts'))
		tree_png <- file.path(outroot, paste0('hft_tree_', park, '_', res, '.png'))
		ggsave(tree_png, treep, width=ht$width, height=ht$height, dpi=300)
		# With the R2 and the mean silhouette width of each cut and the suggested cut of the labels (the area of the
		# circles and lines are sums of 'ha'); the colours in quotes
		for (j in seq_along(rows)) {
			if (!is.null(r2_tree)) ht$paths[[paste0('r2_', rows[j])]] <- signif(r2_tree[j], 6)
			ht$paths[[paste0('silhouette_', rows[j])]] <- signif(sil_tree[[j]], 6)
		}
		ht$paths$suggested_cut <- if (is.na(suggested)) NA else suggested
		write_plot_csv(ht$paths, tree_png, digits=NA)
	}

	progress('DONE', 1)
	# --- VARIANZA DE LOS PÍXELES EXPLICADA POR LOS HFTs (R2) SEGÚN EL NÚMERO DE HFTs, EN CADA UMBRAL ---
	# Los segmentos de cada umbral clasificados (CLASS_METHOD) en 2, 3, ... HFTs; los puntos son los HFTs de
	# los cortes del umbral elegido. Muestra si, con el mismo número de HFTs, otro umbral recogería más varianza
	if (!is.null(ss_sel)) {
		kmax <- max(30, cut_r2$hfts + 5)
		curves <- NULL
		for (tk in thr_all) {
			hf <- paste('csv/park_', park, '_hri_results', tk, '.csv', sep='')
			if (!file.exists(hf) || !file.exists(file.path('shp', paste0('park_segm_', park, '_', tk, '_diss.shp')))) next
			hri_t <- if (tk == res) hri else read_hri(hf)
			if (nrow(hri_t) < 3) next
			seg_t <- if (tk == res) seg_sel else segm_cells(park, tk, hri_t, pixgrid)
			ct_t <- if (tk == res) ct else class_tree(if (is.null(pixclip)) hri_t else clip_hri(hri_t, seg_t, pixclip, median_stat), class_method)
			ss_t <- if (tk == res) ss_sel else segm_sums(seg_t, pixvals, nrow(hri_t))
			nh <- 2:min(nrow(hri_t) - 1, kmax)
			cls_t <- lapply(nh, function(h) classify_k(ct_t, h, class_method))
			r2m <- do.call(rbind, lapply(cls_t, function(cl) r2_groups(ss_t, cl)))
			colnames(r2m) <- var_names
			# Mean silhouette width of the segments of the same HFTs (distances computed once per threshold)
			d_t <- if (tk == res) dsil else if (!ct_t$lowmem) as.matrix(dist(ct_t$z))
			sil_t <- sapply(sil_widths(ct_t, d_t, cls_t), function(sl) sil_mean(sl$s, ct_t$n))
			rm(d_t)
			curves <- rbind(curves, data.frame(threshold=thr_val(tk), segments=nrow(hri_t),
			                                   r2_segments=mean(r2_groups(ss_t, seq_len(nrow(hri_t))), na.rm=TRUE),
			                                   hfts=nh, r2=rowMeans(r2m, na.rm=TRUE), silhouette=sil_t, r2m, check.names=FALSE))
		}
		# Corte sugerido (calculado antes de los cortes) y parte del R2 de los segmentos que conserva cada corte
		cut_r2$r2_share <- cut_r2$r2 / curves$r2_segments[curves$threshold == thr_val(res)][1]
		cut_r2$suggested <- as.integer(cut_r2$cut %in% suggested)
		# silhouette: mean silhouette width of the segments weighted by area (sil_width; NA with one HFT per segment)
		cut_r2 <- cut_r2[, c('cut', 'hfts', 'r2', 'r2_share', 'silhouette', 'suggested', var_names)]
		# Also what the plot marks on the curves (r2_<id>_hfts.png): selected, the selected threshold (thick orange
		# line); bend, the bend of its curve (dashed line); cut, the cuts drawn as an orange point at that number of HFTs
		# (several separated by commas if they have the same HFTs), with their R2 (cut_r2) and mean silhouette width
		# (cut_silhouette) and suggested = 1 if one of them is the suggested cut. A cut beyond the curve (one HFT per
		# segment) gets a row of its own. silhouette: that of the HFTs of each threshold and number of HFTs (not drawn)
		r2_csv <- curves
		r2_csv$selected <- as.integer(r2_csv$threshold == thr_val(res))
		r2_csv$bend <- as.integer(r2_csv$selected == 1 & r2_csv$hfts %in% knee)
		r2_csv$cut <- NA
		r2_csv$cut_r2 <- NA
		r2_csv$cut_silhouette <- NA
		r2_csv$suggested <- 0L
		for (h in unique(cut_r2$hfts)) {
			w <- cut_r2$hfts == h
			i <- which(r2_csv$selected == 1 & r2_csv$hfts == h)
			if (length(i) == 0) {
				new <- r2_csv[NA_integer_, ]
				new[, c('threshold', 'hfts', 'selected', 'bend')] <- list(thr_val(res), h, 1L, 0L)
				new[, c('segments', 'r2_segments')] <- r2_csv[r2_csv$selected == 1, c('segments', 'r2_segments')][1, ]
				r2_csv <- rbind(r2_csv, new)
				i <- nrow(r2_csv)
			}
			r2_csv$cut[i] <- paste(cut_r2$cut[w], collapse=',')
			r2_csv$cut_r2[i] <- cut_r2$r2[w][1]
			r2_csv$cut_silhouette[i] <- cut_r2$silhouette[w][1]
			r2_csv$suggested[i] <- as.integer(any(cut_r2$suggested[w] == 1))
		}
		write.table(r2_csv, file.path(outroot, paste0('r2_', park, '_hfts.csv')), sep=' ', row.names=F, quote=F)
		# with HFT_SUBCLASSES, the subclasses of each cut and the part of the variance of the pixels that they explain
		if (length(cut_sub)) {
			cut_r2$subclasses <- as.integer(cut_sub[cut_r2$cut])
			cut_r2$r2_subclasses <- as.numeric(cut_r2sub[cut_r2$cut])
		}
		write.table(cut_r2, file.path(outroot, paste0('r2_', park, '_cuts.csv')), sep=' ', row.names=F, quote=F)
		if (!is.na(knee)) {
			s <- cut_r2[cut_r2$suggested == 1, ]
			msg <- paste0('ID ', park, ': the R2 curve of the threshold ', k, ' bends at ', knee, ' HFTs',
			              if (nrow(s) == 1) paste0(sprintf('; suggested cut %s: %d HFTs, R2 %.2f (%.0f %% of the R2 of the segments)',
			                                               s$cut, s$hfts, s$r2, 100 * s$r2_share),
			                                       if (!is.na(s$silhouette)) sprintf(', silhouette %.2f', s$silhouette)))
			print(msg)
			write(paste('Suggested cut:', msg), file.path(outroot, 'processing_metadata.txt'), append=TRUE)
		}

		# Umbrales en una rampa de azules (claro: el más bajo, oscuro: el más alto) y el elegido en naranja y más grueso
		curves$thr <- factor(thr_of_value(curves$threshold), levels=thr_of_value(sort(unique(curves$threshold))))
		sel <- k
		thr_cols <- setNames(colorRampPalette(c('#86b6ef', '#0d366b'))(nlevels(curves$thr)), levels(curves$thr))
		thr_cols[sel] <- '#eb6834'
		seg_info <- unique(curves[, c('thr', 'segments', 'r2_segments')])
		thr_labels <- setNames(paste0(seg_info$thr, ': ', seg_info$segments, ' segments, R2 ', sprintf('%.2f', seg_info$r2_segments),
		                              ifelse(seg_info$thr == sel, ' (selected)', '')), seg_info$thr)
		# The mean silhouette width in a second line, so that the labels are not wider (they would cover the bend)
		cut_pts <- data.frame(hfts=cut_r2$hfts, r2=cut_r2$r2,
		                      label=paste0(if (simil == 3) paste0(cut_r2$cut, ': ') else '', cut_r2$hfts, ' HFTs, R2 ', sprintf('%.2f', cut_r2$r2),
		                                   ifelse(is.na(cut_r2$silhouette), '', sprintf(',\nsilhouette %.2f', cut_r2$silhouette)),
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
