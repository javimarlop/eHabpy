# Exporta las sumas de los píxeles de cada segmento que necesita compare_algorithms.py.
#
#   Rscript segment_sums.R [SEGM_OPTIM] [UMBRALES] [SALIDA] [ESCENARIO]
#
# SEGM_OPTIM: carpeta segm_optim con csv/ y shp/ de una ejecución (por defecto ../segm_optim) y
# ../inVars junto a ella; UMBRALES: índices de los umbrales separados por comas (por defecto 3,4,5,6,7);
# SALIDA: carpeta de los resultados (por defecto data); ESCENARIO: carpeta de resultados de un escenario
# (results/<STUDY_AREA>/<SCENARIO>), de la que se exporta además el HFT de cada segmento en cada corte
# (hfts_<corte>.csv: segm_id, hft), la tabla de los umbrales (thresholds.csv) y las curvas de R2 (r2_hfts.csv).
# Los píxeles con datos en todas las variables se tipifican con su media y su desviación típica (como
# compare_classification.R). En data/, para cada umbral k:
#   sums_0k.csv  una fila por segmento (en el orden de csv/park_<id>_hri_results<k>.csv): segm_id, n
#                (nº de píxeles), s1..sV (sumas) y q_i_j (sumas de productos cruzados, i <= j)
#   seg_0k.bin   fila (1 a nº de segmentos) del segmento de cada píxel de la malla, 0 sin datos (int32)
#   pca_0k.csv   PCA de los píxeles tipificados (como i.pca -n de la segmentación): por eje, valor propio,
#                desviación típica y rango (máx. - mín.) de las puntuaciones de los píxeles, y cargas
# y grid.csv (filas y columnas de la malla) y vars.csv (nombres cortos de las variables).

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub('--file=', '', grep('--file=', commandArgs(), value = TRUE))))
segm <- normalizePath(if (length(args) >= 1) args[1] else file.path(here, '..', 'segm_optim'))
thr <- if (length(args) >= 2) as.integer(strsplit(args[2], ',')[[1]]) else 3:7
out <- if (length(args) >= 3) args[3] else file.path(here, 'data')
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out)
scen <- if (length(args) >= 4) normalizePath(args[4]) else NA
setwd(segm)
park <- read.table('csv/segm_done.csv')[1, 1]

# Variables (archivo y NoData) de ENV_VARS en config.py, como moranvar_plots.R
lines <- readLines('config.py', warn = FALSE)
m <- regmatches(lines, regexec("'([^']+)'\\s*:\\s*\\{\\s*'file'\\s*:\\s*'([^']+)'\\s*,\\s*'nodata'\\s*:\\s*([^,}]+)", lines))
m <- m[lengths(m) == 4]
ev <- data.frame(var = sapply(m, `[`, 2), file = sapply(m, `[`, 3), nodata = as.numeric(trimws(sapply(m, `[`, 4))))
nv <- nrow(ev)
write.csv(data.frame(var = ev$var, short = sapply(strsplit(ev$var, '_'), `[`, 3)), file.path(out, 'vars.csv'), row.names = FALSE)

rs <- lapply(file.path('..', 'inVars', ev$file), terra::rast, lyrs = 1)
grid <- rs[[1]]
write.csv(data.frame(nrow = terra::nrow(grid), ncol = terra::ncol(grid)), file.path(out, 'grid.csv'), row.names = FALSE)
P <- sapply(seq_len(nv), function(i) { x <- terra::values(rs[[i]], mat = FALSE); x[x == ev$nodata[i]] <- NA; x })
pairs <- which(upper.tri(diag(nv), diag = TRUE), arr.ind = TRUE)

for (res in thr) {
	hri <- read.table(sprintf('csv/park_%s_hri_results%d.csv', park, res), header = TRUE)
	nseg <- nrow(hri)
	v <- terra::vect(sprintf('shp/park_segm_%s_%d_diss.shp', park, res))
	v$row <- match(as.character(v$segm_id), format(hri$segm_id, scientific = FALSE, trim = TRUE))
	seg <- terra::values(terra::rasterize(v, grid, field = 'row'), mat = FALSE)
	ok <- !is.na(seg) & complete.cases(P)
	Z <- scale(P[ok, ])
	g <- seg[ok]
	tot <- function(x) { o <- numeric(nseg); r <- rowsum(x, g); o[as.integer(rownames(r))] <- r; o }
	S <- sapply(seq_len(nv), function(i) tot(Z[, i]))
	Q <- sapply(seq_len(nrow(pairs)), function(p) tot(Z[, pairs[p, 1]] * Z[, pairs[p, 2]]))
	tab <- data.frame(segm_id = format(hri$segm_id, scientific = FALSE, trim = TRUE), n = tot(rep(1, length(g))))
	tab <- cbind(tab, setNames(as.data.frame(S), paste0('s', seq_len(nv))),
	             setNames(as.data.frame(Q), sprintf('q_%d_%d', pairs[, 1], pairs[, 2])))
	write.csv(tab, file.path(out, sprintf('sums_0%d.csv', res)), row.names = FALSE)
	cells <- seg
	cells[!ok] <- 0
	writeBin(as.integer(cells), file.path(out, sprintf('seg_0%d.bin', res)), size = 4)
	eg <- eigen(cor(Z), symmetric = TRUE)
	sc <- Z %*% eg$vectors
	write.csv(data.frame(axis = seq_len(nv), eigenvalue = eg$values, sd = apply(sc, 2, sd),
	                     range = apply(sc, 2, function(x) diff(range(x))),
	                     setNames(as.data.frame(t(eg$vectors)), paste0('load', seq_len(nv)))),
	          file.path(out, sprintf('pca_0%d.csv', res)), row.names = FALSE)
	cat(sprintf('Threshold 0.%d: %d segments, %d pixels\n', res, nseg, sum(ok)))
}

# HFT de cada segmento en cada corte del escenario (park_segm_<id>_<k>_<corte>[_optim]_class.shp)
if (!is.na(scen)) {
	file.copy(sprintf('csv/%s_threshold.csv', park), file.path(out, 'thresholds.csv'), overwrite = TRUE)
	file.copy(file.path(scen, sprintf('r2_%s_hfts.csv', park)), file.path(out, 'r2_hfts.csv'), overwrite = TRUE)
	for (f in list.files(scen, pattern = '_class\\.shp$', recursive = TRUE, full.names = TRUE)) {
		cut <- sub('.*_(p[0-9]+|results)(_optim)?_class\\.shp$', '\\1', basename(f))
		if (cut == basename(f)) cut <- 'results'
		h <- sf::st_drop_geometry(sf::st_read(f, quiet = TRUE))
		h <- unique(data.frame(segm_id = format(h$segm_id, scientific = FALSE, trim = TRUE), hft = h$hclst_m))
		write.csv(h, file.path(out, sprintf('hfts_%s.csv', cut)), row.names = FALSE)
		cat(sprintf('%s: %d segments, %d HFTs\n', basename(f), nrow(h), length(unique(h$hft))))
	}
}
