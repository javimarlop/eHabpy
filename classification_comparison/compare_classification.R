# Compara métodos de clasificación de los segmentos en HFTs por la separabilidad de las clases.
#
#   Rscript compare_classification.R [SEGM_OPTIM] [UMBRAL]
#
# SEGM_OPTIM: carpeta segm_optim con csv/ y shp/ de una ejecución (por defecto ../segm_optim) y
# ../inVars junto a ella; UMBRAL: índice del umbral (por defecto 5, es decir 0.5).
# Métodos, para 2 a KMAX HFTs:
#   A  actual: Ward (ward.D2) sin ponderar sobre medias y varianzas de los segmentos reescaladas 0-1
#   B  Ward sin ponderar sobre las medias reescaladas 0-1 (sin las varianzas)
#   C  Ward ponderado por el nº de píxeles sobre las medias tipificadas con los píxeles
#   D  C + consolidación con k-medias ponderadas (los segmentos se reasignan al centro más próximo)
# Medidas (con los píxeles de ../inVars, tipificados):
#   R2   varianza de los píxeles entre HFTs / total, media de las variables
#   JM   distancia de Jeffries-Matusita entre cada par de HFTs (0-2, normales multivariantes): media y mínima
#   Sil  anchura de silueta media de los segmentos (distancias entre medias tipificadas), ponderada por área

suppressMessages({library(vegan); library(cluster); library(sf); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub('--file=', '', grep('--file=', commandArgs(), value = TRUE))))
segm <- normalizePath(if (length(args) >= 1) args[1] else file.path(here, '..', 'segm_optim'))
res <- if (length(args) >= 2) as.integer(args[2]) else 5
KMAX <- 20
setwd(segm)
park <- read.table('csv/segm_done.csv')[1, 1]

rescale01 <- function(x) { r <- range(x, na.rm = TRUE); if (diff(r) == 0) 0 * x else (x - r[1]) / diff(r) }

# Variables (archivo y NoData) de ENV_VARS en config.py, como moranvar_plots.R
lines <- readLines('config.py', warn = FALSE)
m <- regmatches(lines, regexec("'([^']+)'\\s*:\\s*\\{\\s*'file'\\s*:\\s*'([^']+)'\\s*,\\s*'nodata'\\s*:\\s*([^,}]+)", lines))
m <- m[lengths(m) == 4]
ev <- data.frame(var = sapply(m, `[`, 2), file = sapply(m, `[`, 3), nodata = as.numeric(trimws(sapply(m, `[`, 4))))
short <- sapply(strsplit(ev$var, '_'), `[`, 3)
nv <- nrow(ev)

hri <- read.table(sprintf('csv/park_%s_hri_results%d.csv', park, res), header = TRUE)
nseg <- nrow(hri)

# Píxeles: segmento de cada píxel y valores tipificados con la media y la desviación de los píxeles
rs <- lapply(file.path('..', 'inVars', ev$file), terra::rast, lyrs = 1)
grid <- rs[[1]]
v <- terra::vect(sprintf('shp/park_segm_%s_%d_diss.shp', park, res))
v$row <- match(as.character(v$segm_id), format(hri$segm_id, scientific = FALSE, trim = TRUE))
seg <- terra::values(terra::rasterize(v, grid, field = 'row'), mat = FALSE)
P <- sapply(seq_len(nv), function(i) { x <- terra::values(rs[[i]], mat = FALSE); x[x == ev$nodata[i]] <- NA; x })
ok <- !is.na(seg) & complete.cases(P)
P <- scale(P[ok, ])
g <- seg[ok]
cat(sprintf('%d segments, %d pixels with data in all the variables\n', nseg, nrow(P)))

# Sumas de cada segmento: n, suma y productos cruzados de los valores tipificados
tot <- function(x) { out <- numeric(nseg); r <- rowsum(x, g); out[as.integer(rownames(r))] <- r; out }
n <- tot(rep(1, length(g)))
S <- sapply(seq_len(nv), function(i) tot(P[, i]))
pairs <- which(upper.tri(diag(nv), diag = TRUE), arr.ind = TRUE)
Q <- sapply(seq_len(nrow(pairs)), function(p) tot(P[, pairs[p, 1]] * P[, pairs[p, 2]]))
M <- S / pmax(n, 1)  # media tipificada de cada segmento
keep <- n > 0

# --- Medidas de una partición (cl: HFT de cada segmento) ---
measures <- function(cl) {
	cl <- cl[keep]; nn <- n[keep]
	Nh <- rowsum(nn, cl); Sh <- rowsum(S[keep, ], cl); Qh <- rowsum(Q[keep, ], cl)
	N <- sum(nn)
	r2 <- mean((colSums(Sh^2 / as.vector(Nh)) - colSums(S[keep, ])^2 / N) / (colSums(Q[keep, pairs[, 1] == pairs[, 2]]) - colSums(S[keep, ])^2 / N))
	mu <- Sh / as.vector(Nh)
	cov_h <- lapply(seq_len(nrow(Sh)), function(h) {
		C <- matrix(0, nv, nv)
		C[pairs] <- Qh[h, ]
		C[pairs[, 2:1]] <- Qh[h, ]
		(C - Nh[h] * tcrossprod(mu[h, ])) / max(Nh[h] - 1, 1) + diag(1e-6, nv)
	})
	jm <- c()
	for (a in 1:(nrow(Sh) - 1)) for (b in (a + 1):nrow(Sh)) {
		Cm <- (cov_h[[a]] + cov_h[[b]]) / 2
		d <- mu[a, ] - mu[b, ]
		B <- as.numeric(t(d) %*% solve(Cm, d)) / 8 +
			0.5 * (determinant(Cm)$modulus - 0.5 * (determinant(cov_h[[a]])$modulus + determinant(cov_h[[b]])$modulus))
		jm <- c(jm, 2 * (1 - exp(-B)))
	}
	si <- silhouette(as.integer(factor(cl)), dist(M[keep, ]))[, 'sil_width']
	c(hfts = nrow(Sh), r2 = r2, jm_mean = mean(jm), jm_min = min(jm), sil = sum(si * nn) / N)
}

# Ward ponderado: distancia de Ward entre segmentos (tamaños n) y hclust con members
wward <- function(X, w) {
	D <- as.matrix(dist(X))
	hclust(as.dist(sqrt(2 * outer(w, w) / outer(w, w, '+')) * D), 'ward.D2', members = w)
}

# k-medias ponderadas desde una partición: cada segmento al centro (media ponderada) más próximo
wkmeans <- function(X, w, cl, iter = 100) {
	for (it in seq_len(iter)) {
		lev <- sort(unique(cl))
		C <- rowsum(X * w, cl) / as.vector(rowsum(w, cl))
		d2 <- sapply(seq_len(nrow(C)), function(j) colSums((t(X) - C[j, ])^2))
		new <- lev[max.col(-d2, ties.method = 'first')]
		if (all(new == cl)) break
		cl <- new
	}
	cl
}

feat_a <- as.data.frame(lapply(hri[, 3:(2 + 2 * nv)], rescale01))
feat_b <- as.data.frame(lapply(hri[, 3:(2 + nv)], rescale01))
trees <- list(A = hclust(vegdist(feat_a, 'euclidean', na.rm = TRUE), 'ward.D2'),
              B = hclust(vegdist(feat_b, 'euclidean', na.rm = TRUE), 'ward.D2'))
w <- pmax(n, 1)
trees$C <- wward(M, w)
labels <- c(A = 'A: actual (Ward, medias + varianzas 0-1, sin ponderar)',
            B = 'B: Ward, solo medias 0-1, sin ponderar',
            C = 'C: Ward ponderado por área, medias tipificadas',
            D = 'D: C + k-medias ponderadas')

parts <- list()
rows <- list()
for (k in 2:KMAX) {
	for (mth in names(labels)) {
		cl <- if (mth == 'D') wkmeans(M, w, parts[[paste('C', k)]]) else cutree(trees[[mth]], k = k)
		parts[[paste(mth, k)]] <- cl
		rows[[length(rows) + 1]] <- data.frame(method = mth, k = k, t(measures(cl)))
	}
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(here, sprintf('medidas_por_metodo_0%d.csv', res)), row.names = FALSE)

# Corte p50 del método actual con este umbral (como en moranvar_plots.R)
k50 <- length(unique(cutree(trees$A, h = 0.5 * max(trees$A$height))))
ari <- function(a, b) {
	tb <- table(a, b) ; cb <- function(x) sum(x * (x - 1) / 2)
	e <- cb(rowSums(tb)) * cb(colSums(tb)) / cb(sum(tb))
	(cb(tb) - e) / ((cb(rowSums(tb)) + cb(colSums(tb))) / 2 - e)
}
cat(sprintf('\np50 cut of the current method at threshold 0.%d: %d HFTs\n', res, k50))
for (mth in c('A', 'C')) cat(sprintf('HFTs of the cuts p25, p50, p75 of tree %s: %s\n', mth,
    paste(sapply(c(0.25, 0.5, 0.75), function(f) length(unique(cutree(trees[[mth]], h = f * max(trees[[mth]]$height))))), collapse = ', ')))
print(format(tab[tab$k %in% c(3, k50, 12), ], digits = 3), row.names = FALSE)
cat('\nAgreement (adjusted Rand index, by area) with A at k =', k50, ':\n')
for (mth in c('B', 'C', 'D')) cat(' ', mth, round(ari(rep(parts[[paste('A', k50)]], n), rep(parts[[paste(mth, k50)]], n)), 3), '\n')

# Las varianzas de los segmentos dependen de su tamaño: correlación de Spearman con el nº de píxeles
cat('\nSpearman correlation of the segment variances with their number of pixels:\n')
print(round(sapply(seq_len(nv), function(i) cor(hri[, 2 + nv + i], hri$sumpamask, method = 'spearman')), 2))

# --- Figura 1: medidas frente al número de HFTs ---
cols <- c(A = '#2a78d6', B = '#eb6834', C = '#1baf7a', D = '#6250d6')
long <- rbind(data.frame(tab[, c('method', 'k')], medida = 'R2 de los píxeles (media de las variables)', valor = tab$r2),
              data.frame(tab[, c('method', 'k')], medida = 'JM medio entre pares de HFTs (máx. 2)', valor = tab$jm_mean),
              data.frame(tab[, c('method', 'k')], medida = 'JM mínimo: el par de HFTs menos separable', valor = tab$jm_min),
              data.frame(tab[, c('method', 'k')], medida = 'Silueta media de los segmentos (por área)', valor = tab$sil))
long$medida <- factor(long$medida, levels = unique(long$medida))
p <- ggplot(long, aes(k, valor, colour = method)) +
	geom_vline(xintercept = k50, colour = 'grey60', linetype = 'dashed') +
	geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
	facet_wrap(~ medida, scales = 'free_y', ncol = 2) +
	scale_colour_manual(values = cols, labels = labels, name = NULL) +
	scale_x_continuous(breaks = seq(2, KMAX, 2)) +
	labs(title = sprintf('Separabilidad de los HFTs según el método de clasificación (umbral 0.%d, %d segmentos)', res, nseg),
	     subtitle = sprintf('Línea discontinua: %d HFTs, el corte p50 del método actual', k50), x = 'Número de HFTs', y = NULL) +
	theme_bw(base_size = 12) + theme(legend.position = 'bottom', legend.direction = 'vertical')
ggsave(file.path(here, sprintf('fig1_separabilidad_0%d.png', res)), p, width = 10, height = 8.5, dpi = 150)

# --- Figura 2: mapas de HFTs con k50 HFTs, métodos A y D ---
segs <- st_read(sprintf('shp/park_segm_%s_%d_diss.shp', park, res), quiet = TRUE)
segs$row <- match(as.character(segs$segm_id), format(hri$segm_id, scientific = FALSE, trim = TRUE))
maps <- do.call(rbind, lapply(c('A', 'D'), function(mth) {
	s <- segs
	s$HFT <- factor(parts[[paste(mth, k50)]][s$row])
	s$method <- labels[mth]
	s[, c('HFT', 'method')]
}))
ms <- tab[tab$k == k50 & tab$method %in% c('A', 'D'), ]
maps$method <- factor(maps$method, levels = labels[c('A', 'D')],
                      labels = sprintf('%s\nR2 %.2f | JM medio %.2f, mínimo %.2f | silueta %.2f', labels[c('A', 'D')],
                                       ms$r2, ms$jm_mean, ms$jm_min, ms$sil))
pm <- ggplot(maps) + geom_sf(aes(fill = HFT), colour = NA) + facet_wrap(~ method, ncol = 1) +
	scale_fill_brewer(palette = 'Set1') + theme_bw(base_size = 11) +
	labs(title = sprintf('%d HFTs con el método actual (A) y con Ward ponderado + k-medias (D), umbral 0.%d', k50, res),
	     subtitle = 'Los números de los HFTs no se corresponden entre los dos mapas')
ggsave(file.path(here, sprintf('fig2_mapas_0%d.png', res)), pm, width = 10, height = 10, dpi = 130)
cat('\nFiguras y tablas en', here, '\n')
