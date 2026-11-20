# HFT comparator: puts together the HFTs of several runs of the model (other study areas, other scenarios, other
# thresholds or cuts, other options) and shows how similar they are.
#
#   Rscript hft_compare.R [--dir=comparison] [--out=comparison/results]
#
# Input: every CSV of the folder (not of its subfolders) in the format of hft_means_<...>.csv, which moranvar_plots.R
# writes in the folder of each scenario (and, with SIMIL = 3, of each cut): one row per HFT with its `label`, and the mean of
# its pixels in each variable after the column `area_pct`. Other CSVs (hft_stats_..., ...) are listed and skipped. The
# cut is in the name of the files and in the label of every HFT, so it stays in all the tables.
# Only the variables that all the tables have are compared, and they are standardized over all the HFTs together.
# Similarity of two HFTs: 1 - d / max(d), d being the Euclidean distance between their standardized means.
#
# Output (--out): hfts_all.csv (the common table, with a unique id: the label study area/scenario/ID/cut/HFT<n>, and a
# short code of the source: A1 = HFT 1 of the source A), nmds_hfts, tree_hfts, heatmap_hfts, pca_biplot_hfts,
# similarity_hfts, similarity_sources, radarplot_hfts_<source>, best_matches.csv and comparison_metadata.txt. Every
# plot X.png has its data in X.csv.

suppressMessages({
	library(vegan)
	library(ggplot2)
})
source('plot_utils.R') # rescale01, short_labels, write_plot_csv

progress <- function(kind, value) if (nzchar(Sys.getenv('EHAB_PROGRESS'))) cat('##', kind, ' ', value, '\n', sep='') # for the desktop app

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default) {
	hit <- sub(paste0('^--', name, '='), '', grep(paste0('^--', name, '='), args, value = TRUE))
	if (length(hit)) hit[1] else default
}
dir_in <- arg('dir', 'comparison')
dir_out <- arg('out', file.path(dir_in, 'results'))
log_lines <- character(0)
say <- function(...) { m <- paste0(...); log_lines <<- c(log_lines, m); message(m) }

# --- Read the tables ---------------------------------------------------------------------------------------------
progress('STEPS', 'compare')
progress('STEP', 'compare')
files <- sort(list.files(dir_in, pattern = '\\.csv$', ignore.case = TRUE, full.names = TRUE))
if (!length(files)) stop('No CSV files in ', dir_in, ': add the hft_means tables of the runs to compare')

read_table <- function(f) {
	first <- readLines(f, n = 1, warn = FALSE)
	sep <- if (grepl(';', first)) ';' else if (grepl(',', first)) ',' else ' '
	read.table(f, header = TRUE, sep = sep, check.names = FALSE, comment.char = '', quote = '"', stringsAsFactors = FALSE,
	           na.strings = c('NA', ''), strip.white = TRUE)
}
tabs <- list()
for (f in files) {
	d <- tryCatch(read_table(f), error = function(e) NULL)
	if (is.null(d) || !all(c('label', 'area_pct') %in% names(d))) {
		say('Skipped (not a table of HFTs with the columns label and area_pct, e.g. hft_stats): ', basename(f))
		next
	}
	d$file <- basename(f)
	tabs[[f]] <- d
}
if (!length(tabs)) stop('None of the CSV files of ', dir_in, ' is a table of the HFTs (hft_means_...)')
say('Tables used: ', length(tabs), ' of ', length(files), ' CSV files')

# Variables: the columns after area_pct that all the tables have
var_of <- function(d) setdiff(names(d)[-seq_len(match('area_pct', names(d)))], 'file')
vars <- Reduce(intersect, lapply(tabs, var_of))
dropped <- setdiff(unique(unlist(lapply(tabs, var_of))), vars)
if (length(dropped)) say('Variables left out (not in every table): ', paste(dropped, collapse = ', '))

# Identification columns: when a table lacks some (made by hand), they come from the label area/scenario/ID/[cut/]HFT<n>
ids <- c('study_area', 'scenario', 'id', 'resolution', 'threshold', 'cut', 'suggested', 'HFT', 'colour', 'segments', 'pixels',
         'silhouette', 'area_ha', 'area_pct')
fill_ids <- function(d) {
	p <- strsplit(d$label, '/', fixed = TRUE)
	part <- function(i, from_end = FALSE) vapply(p, function(x) { n <- length(x); j <- if (from_end) n - i else i
		if (j >= 1 && j <= n) x[j] else NA_character_ }, character(1))
	if (!'study_area' %in% names(d)) d$study_area <- part(1)
	if (!'scenario' %in% names(d)) d$scenario <- part(2)
	if (!'id' %in% names(d)) d$id <- part(3)
	if (!'HFT' %in% names(d)) d$HFT <- as.integer(sub('^HFT', '', part(0, TRUE)))
	if (!'cut' %in% names(d)) d$cut <- ifelse(vapply(p, length, 1L) >= 5, part(4), NA_character_)
	for (n in setdiff(ids, names(d))) d[[n]] <- NA
	d
}
hf <- do.call(rbind, lapply(tabs, function(d) { d <- fill_ids(d); d[, c('label', 'file', ids, vars)] }))
rownames(hf) <- NULL
n0 <- nrow(hf)
hf <- hf[!duplicated(hf$label), ] # the same run in a table of all the cuts and in one of a cut
if (nrow(hf) < n0) say('Repeated HFTs (same label in several tables) used once: ', n0 - nrow(hf))
hf$cut[is.na(hf$cut) | hf$cut == ''] <- NA

# Variables used: with values in every HFT, and not constant
ok <- vapply(vars, function(v) { x <- hf[[v]]; !anyNA(x) && length(unique(x)) > 1 }, logical(1))
if (any(!ok)) say('Variables left out (NA in some HFT, or the same value in all): ', paste(vars[!ok], collapse = ', '))
vars <- vars[ok]
if (length(vars) < 2) stop('At least 2 variables in common (with values in all the HFTs) are needed to compare the HFTs')
res_sizes <- unique(na.omit(hf$resolution))
if (length(res_sizes) > 1) {
	say("WARNING: HFTs made with different pixel sizes (", paste(res_sizes, collapse = ', '),
	    " m): the variables and the HFTs change with it, so the comparison may tell the pixel from the landscape")
}

# Sources: one per study area, scenario, ID and cut; code A, B, C... for the plots
hf$source <- paste(hf$study_area, hf$scenario, hf$id, ifelse(is.na(hf$cut), '', hf$cut), sep = '/')
hf$source <- sub('/$', '', hf$source)
src <- unique(hf$source)
codes <- if (length(src) <= 26) LETTERS[seq_along(src)] else paste0('S', seq_along(src))
hf$code <- paste0(codes[match(hf$source, src)], hf$HFT)
hf$source_code <- codes[match(hf$source, src)]
hf$uid <- hf$label
n <- nrow(hf)
say('HFTs: ', n, ' of ', length(src), ' sources; ', length(vars), ' variables')
legend_txt <- paste0(codes, ': ', src)
for (l in legend_txt) say('  ', l)

# Clean output folder
if (dir.exists(dir_out)) unlink(list.files(dir_out, full.names = TRUE), recursive = TRUE) else dir.create(dir_out, recursive = TRUE)
progress('UNITS', 9)

# The common table
common <- hf[, c('uid', 'code', 'source', 'file', 'study_area', 'scenario', 'id', 'resolution', 'threshold', 'cut', 'suggested',
                 'HFT', 'colour', 'segments', 'pixels', 'silhouette', 'area_ha', 'area_pct', vars)]
write_plot_csv(common, file.path(dir_out, 'hfts_all.png'), digits = NA)
progress('DONE', 1)

if (n < 4) {
	say('Only ', n, ' HFTs: at least 4 are needed for the NMDS, the tree and the similarity; the common table is in hfts_all.csv')
	writeLines(log_lines, file.path(dir_out, 'comparison_metadata.txt'))
	quit(save = 'no')
}

# --- Standardize over all the HFTs together, and distances ---------------------------------------------------
x <- scale(as.matrix(hf[, vars]))
rownames(x) <- hf$uid
d <- dist(x)
sim <- 1 - as.matrix(d) / max(d)
vlab <- gsub('\n', ' ', short_labels(vars, width = 100))
pal <- if (length(codes) <= 8) setNames(RColorBrewer::brewer.pal(max(3, length(codes)), 'Set1')[seq_along(codes)], codes) else
	setNames(scales::hue_pal()(length(codes)), codes)
col_scale <- scale_colour_manual(values = pal, breaks = codes, labels = legend_txt, name = 'Source')
base_theme <- theme_bw(base_size = 13) + theme(legend.position = 'bottom', legend.direction = 'vertical')
src_h <- 0.22 * length(codes) # inches that the legend of the sources needs
progress('DONE', 1)

# --- NMDS --------------------------------------------------------------------------------------------------------
set.seed(1)
nmds <- tryCatch(metaMDS(d, k = 2, trymax = 50, trace = 0, autotransform = FALSE), error = function(e) { say('NMDS failed: ', conditionMessage(e)); NULL })
if (!is.null(nmds)) {
	pts <- data.frame(uid = hf$uid, code = hf$code, source = hf$source, source_code = hf$source_code, HFT = hf$HFT, NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2])
	fit <- envfit(nmds, as.data.frame(x), permutations = 0)
	arr <- fit$vectors$arrows
	mul <- 0.85 * max(abs(nmds$points)) / max(abs(arr))
	vec <- data.frame(variable = vars, label = vlab, NMDS1 = arr[, 1] * mul, NMDS2 = arr[, 2] * mul, r2 = fit$vectors$r)
	p <- ggplot(pts, aes(NMDS1, NMDS2)) +
		geom_hline(yintercept = 0, colour = 'grey85') + geom_vline(xintercept = 0, colour = 'grey85') +
		geom_segment(data = vec, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), inherit.aes = FALSE, colour = 'grey45',
		             arrow = arrow(length = unit(0.12, 'cm'))) +
		ggrepel::geom_text_repel(data = vec, aes(NMDS1, NMDS2, label = label), inherit.aes = FALSE, colour = 'grey30', size = 3, seed = 2, max.overlaps = Inf) +
		scale_x_continuous(expand = expansion(mult = 0.15)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
		geom_point(aes(colour = source_code), size = 3) +
		ggrepel::geom_text_repel(aes(label = code, colour = source_code), size = 3.3,
		                         show.legend = FALSE, seed = 1, max.overlaps = Inf) +
		col_scale + coord_equal() + base_theme +
		labs(title = sprintf('NMDS of the HFTs (stress %.2f)', nmds$stress),
		     subtitle = 'Euclidean distances between the standardized means of the variables; arrows: variables (envfit)')
	ggsave(file.path(dir_out, 'nmds_hfts.png'), p, width = 10, height = 8 + src_h, dpi = 200)
	write_plot_csv(pts, file.path(dir_out, 'nmds_hfts.png'))
	write_plot_csv(vec, file.path(dir_out, 'nmds_hfts_variables.png'))
}
progress('DONE', 1)

# --- Tree (Ward) -------------------------------------------------------------------------------------------------
hc <- hclust(d, method = 'ward.D2')
cols_leaf <- pal[hf$source_code]
png(file.path(dir_out, 'tree_hfts.png'), width = 1400, height = 700 + 14 * length(codes) + 8 * min(n, 60), res = 130)
par(mar = c(9, 4, 3, 1))
dend <- as.dendrogram(hc)
col_by_label <- setNames(cols_leaf, hf$uid)
code_by_label <- setNames(hf$code, hf$uid)
dend <- dendrapply(dend, function(nd) {
	if (is.leaf(nd)) {
		lab <- attr(nd, 'label')
		attr(nd, 'nodePar') <- list(lab.col = col_by_label[[lab]], pch = NA, lab.cex = 0.8)
		attr(nd, 'label') <- code_by_label[[lab]]
	}
	nd
})
plot(dend, main = 'Tree of the HFTs of all the sources (Ward, standardized means)', ylab = 'Height')
legend('topright', legend = legend_txt, col = pal[codes], pch = 15, bty = 'n', cex = 0.7)
invisible(dev.off())
lab_of <- function(i) ifelse(i < 0, hf$uid[abs(i)], paste('cluster', i))
tree_df <- data.frame(step = seq_len(nrow(hc$merge)), left = lab_of(hc$merge[, 1]), right = lab_of(hc$merge[, 2]), height = hc$height)
write_plot_csv(tree_df, file.path(dir_out, 'tree_hfts.png'))
progress('DONE', 1)

# --- Heatmap of the standardized means, in the order of the tree ---------------------------------------------
ord <- hc$order
zl <- data.frame(uid = rep(hf$uid[ord], times = length(vars)), code = rep(hf$code[ord], times = length(vars)),
                 variable = rep(vars, each = n), z = as.vector(x[ord, ]), mean = as.vector(as.matrix(hf[ord, vars])))
zl$y <- factor(zl$code, levels = rev(hf$code[ord]))
zl$v <- factor(vlab[match(zl$variable, vars)], levels = vlab)
ph <- ggplot(zl, aes(v, y, fill = z)) + geom_tile() +
	scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b', limits = c(-3, 3), oob = scales::squish, name = 'z') +
	labs(x = NULL, y = NULL, title = 'Profiles of the HFTs of all the sources (standardized means)',
	     subtitle = 'rows in the order of the tree; z: standard deviations from the mean of all the HFTs') +
	theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 0))
ggsave(file.path(dir_out, 'heatmap_hfts.png'), ph, width = 4 + 0.7 * length(vars), height = 2.5 + 0.22 * n, dpi = 200, limitsize = FALSE)
write_plot_csv(zl[, c('uid', 'code', 'variable', 'z', 'mean')], file.path(dir_out, 'heatmap_hfts.png'))
progress('DONE', 1)

# --- PCA biplot (as the model: eigen of the correlations of the standardized variables) ---------------------------
ev <- eigen(cor(x))
npc <- min(3, length(vars))
sc <- x %*% ev$vectors[, seq_len(npc), drop = FALSE]
colnames(sc) <- paste0('PC', seq_len(npc))
pct <- 100 * ev$values / sum(ev$values)
load <- ev$vectors[, seq_len(npc), drop = FALSE]
colnames(load) <- colnames(sc)
panels <- if (npc >= 3) list(c(1, 2), c(1, 3)) else list(c(1, 2))
pc_pts <- NULL; pc_arr <- NULL
for (pn in panels) {
	nm <- sprintf('PC%d (%.0f %%) vs PC%d (%.0f %%)', pn[1], pct[pn[1]], pn[2], pct[pn[2]])
	mulp <- 0.85 * max(abs(sc[, pn])) / max(abs(load[, pn]))
	pc_pts <- rbind(pc_pts, data.frame(panel = nm, uid = hf$uid, code = hf$code, source_code = hf$source_code, X = sc[, pn[1]], Y = sc[, pn[2]]))
	pc_arr <- rbind(pc_arr, data.frame(panel = nm, variable = vars, label = vlab, X = load[, pn[1]] * mulp, Y = load[, pn[2]] * mulp))
}
pc_pts$panel <- factor(pc_pts$panel, levels = unique(pc_pts$panel))
pc_arr$panel <- factor(pc_arr$panel, levels = levels(pc_pts$panel))
pp <- ggplot(pc_pts, aes(X, Y)) +
	geom_hline(yintercept = 0, colour = 'grey85') + geom_vline(xintercept = 0, colour = 'grey85') +
	geom_segment(data = pc_arr, aes(x = 0, y = 0, xend = X, yend = Y), inherit.aes = FALSE, colour = 'grey45', arrow = arrow(length = unit(0.12, 'cm'))) +
	ggrepel::geom_text_repel(data = pc_arr, aes(X, Y, label = label), inherit.aes = FALSE, colour = 'grey30', size = 3, seed = 2, max.overlaps = Inf) +
	scale_x_continuous(expand = expansion(mult = 0.15)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
	geom_point(aes(colour = source_code), size = 3) +
	ggrepel::geom_text_repel(aes(label = code, colour = source_code), size = 3.2, show.legend = FALSE, seed = 1, max.overlaps = Inf) +
	facet_wrap(~panel, scales = 'free') + col_scale + base_theme +
	labs(title = 'PCA of the HFTs of all the sources', x = NULL, y = NULL,
	     subtitle = 'standardized means of the variables; arrows: loadings of the variables')
ggsave(file.path(dir_out, 'pca_biplot_hfts.png'), pp, width = 6.5 * length(panels) + 1, height = 6.5 + src_h, dpi = 200)
write_plot_csv(cbind(uid = hf$uid, code = hf$code, as.data.frame(sc), stringsAsFactors = FALSE), file.path(dir_out, 'pca_biplot_hfts.png'))
write_plot_csv(data.frame(variable = vars, as.data.frame(load), stringsAsFactors = FALSE), file.path(dir_out, 'pca_biplot_hfts_loadings.png'))
write_plot_csv(data.frame(PC = paste0('PC', seq_along(pct)), variance_pct = pct), file.path(dir_out, 'pca_biplot_hfts_variance.png'))
progress('DONE', 1)

# --- Similarity between HFTs, in the order of the tree ---------------------------------------------------------------
sl <- expand.grid(i = seq_len(n), j = seq_len(n))
sl <- data.frame(uid_row = hf$uid[sl$i], uid_col = hf$uid[sl$j], code_row = hf$code[sl$i], code_col = hf$code[sl$j], similarity = sim[cbind(sl$i, sl$j)])
sl$r <- factor(sl$code_row, levels = rev(hf$code[ord]))
sl$cc <- factor(sl$code_col, levels = hf$code[ord])
ps <- ggplot(sl, aes(cc, r, fill = similarity)) + geom_tile() +
	scale_fill_gradientn(colours = c('#f7fbff', '#6baed6', '#08306b'), limits = c(0, 1), name = 'Similarity') +
	labs(x = NULL, y = NULL, title = 'Similarity between the HFTs of all the sources',
	     subtitle = '1 - distance / maximum distance between the standardized means; in the order of the tree') +
	theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + coord_equal()
ggsave(file.path(dir_out, 'similarity_hfts.png'), ps, width = 3 + 0.2 * n, height = 2.5 + 0.2 * n, dpi = 200, limitsize = FALSE)
write_plot_csv(sl[, c('uid_row', 'uid_col', 'code_row', 'code_col', 'similarity')], file.path(dir_out, 'similarity_hfts.png'))
progress('DONE', 1)

# --- The most similar HFT of every other source -------------------------------------------------------------------
bm <- NULL
for (i in seq_len(n)) for (s in setdiff(codes, hf$source_code[i])) {
	js <- which(hf$source_code == s)
	j <- js[which.max(sim[i, js])]
	bm <- rbind(bm, data.frame(uid = hf$uid[i], code = hf$code[i], source = hf$source[i], other_source = src[match(s, codes)], other_code = s,
	                           best_match = hf$uid[j], best_match_code = hf$code[j], similarity = sim[i, j]))
}
if (!is.null(bm)) {
	write_plot_csv(bm, file.path(dir_out, 'best_matches.png'))
	# per pair of sources: the mean over the HFTs of the first of the similarity with their best match in the second
	pairs <- aggregate(similarity ~ code_from + code_to, data = data.frame(code_from = hf$source_code[match(bm$uid, hf$uid)], code_to = bm$other_code, similarity = bm$similarity), FUN = mean)
	pairs$code_from <- factor(pairs$code_from, levels = codes); pairs$code_to <- factor(pairs$code_to, levels = rev(codes))
	pq <- ggplot(pairs, aes(code_from, code_to, fill = similarity)) + geom_tile() +
		geom_text(aes(label = sprintf('%.2f', similarity)), size = 3.5) +
		scale_fill_gradientn(colours = c('#f7fbff', '#6baed6', '#08306b'), limits = c(0, 1), name = 'Mean similarity') +
		scale_colour_identity() +
		labs(x = 'HFTs of the source', y = 'compared with the HFTs of the source',
		     title = 'Mean similarity of every HFT with its most similar HFT of another source',
		     subtitle = paste(legend_txt, collapse = '\n')) +
		theme_minimal(base_size = 12) + theme(plot.subtitle = element_text(size = 8)) + coord_equal()
	ggsave(file.path(dir_out, 'similarity_sources.png'), pq, width = 8, height = 6 + 0.17 * length(codes), dpi = 200)
	pairs$source_from <- src[match(as.character(pairs$code_from), codes)]; pairs$source_to <- src[match(as.character(pairs$code_to), codes)]
	write_plot_csv(pairs[, c('code_from', 'code_to', 'source_from', 'source_to', 'similarity')], file.path(dir_out, 'similarity_sources.png'))
}
progress('DONE', 1)

# --- Radar plots of the HFTs of each source (variables rescaled 0-1 among all the HFTs) ---------------------------
source('CreateRadialPlot.R')
rs <- as.data.frame(lapply(as.data.frame(hf[, vars]), rescale01))
max_radar <- 12
for (s in head(codes, max_radar)) {
	i <- which(hf$source_code == s)
	dat <- cbind(group = hf$HFT[i], rs[i, , drop = FALSE])
	kol <- hf$colour[i]
	if (anyNA(kol) || anyDuplicated(kol)) kol <- scales::hue_pal()(length(i))
	radar <- tryCatch(CreateRadialPlot(dat, axis.labels = short_labels(vars), axis.label.size = 3.5, axis.label.offset = 1.25,
	                                   plot.extent.x.sf = 1.9, plot.extent.y.sf = 1.35,
	                                   group.colours = setNames(kol, dat$group), group.order = dat$group),
	                  error = function(e) { say('Radar plot of ', s, ' failed: ', conditionMessage(e)); NULL })
	if (is.null(radar)) next
	radar <- radar + labs(title = paste0(s, ': ', src[match(s, codes)]), subtitle = 'variables rescaled 0-1 among the HFTs of all the sources') +
		theme(plot.title = element_text(size = 13), plot.subtitle = element_text(size = 11))
	fn <- file.path(dir_out, paste0('radarplot_hfts_', s, '.png'))
	ggsave(fn, radar, width = 10, height = 7.5, dpi = 200)
	long <- data.frame(uid = rep(hf$uid[i], times = length(vars)), code = rep(hf$code[i], times = length(vars)),
	                   variable = rep(vars, each = length(i)), value = as.vector(as.matrix(hf[i, vars])), relative = as.vector(as.matrix(rs[i, , drop = FALSE])))
	write_plot_csv(long, fn)
}
if (length(codes) > max_radar) say('Radar plots of the first ', max_radar, ' sources only')
progress('DONE', 1)

# --- Metadata ----------------------------------------------------------------------------------------------------------
writeLines(c(paste('HFT comparator,', format(Sys.time(), '%Y-%m-%d %H:%M')), paste('Folder:', dir_in), '',
             'Tables read:', paste0('  ', basename(files)), '', log_lines, '',
             'Variables compared:', paste0('  ', vars)), file.path(dir_out, 'comparison_metadata.txt'))
say('Done: the results are in ', dir_out)
