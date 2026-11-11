#### eHabitat+ HFT module: maps and charts of the similarity of the HFTs in the ecoregion (run by similarity.py)
#### License: CC BY-SA 3.0
#### Rscript similarity_plots.R --dir=<similarity run folder> [--colours=<id>=<hft_colours csv> ...] [--ecoregions=<file>]
####                            [--study=<id>=<class.shp> ...]
#### Reads hri_results.csv, the similarity maps sim_<id>_hft<h>.tif and thdi_res.csv (optional) of the folder and
#### writes there, for each study area (wdpaid) of hri_results.csv:
####	sim_maps_<id>.png: one map of the similarity per HFT, with the ecoregion outline and the study area (red)
####	sim_metrics_<id>.png: bar charts of the similarity metrics per HFT, with the THDI of the study area
#### Without --colours or --study, the hft_colours_<id>_*.csv and park_segm_<id>_*_class.shp of the folder above
#### (the cut) are used when there is only one of them.

usage <- paste(
	'Usage: Rscript similarity_plots.R --dir=<similarity run folder> [--colours=<id>=<hft_colours csv> ...]',
	'                                  [--ecoregions=<file>] [--study=<id>=<class.shp> ...] [--area=<name>]',
	'Writes sim_maps_<id>.png and sim_metrics_<id>.png in the folder of a run of similarity.py.', sep='\n')

# Arguments: one argv element each, --name=value (the value may contain spaces and '='); --colours and --study
# are <id>=<file> and can be repeated
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	cat(usage, '\n', sep='')
	quit(save='no', status=1)
}
if (any(args %in% c('-h', '--help'))) {
	cat(usage, '\n', sep='')
	quit(save='no', status=0)
}
opt <- list(dir='', ecoregions='', area='the ecoregion')  # area: what was compared, for the titles
per_id <- list(colours=character(0), study=character(0))
for (a in args) {
	m <- regmatches(a, regexec('^--([A-Za-z-]+)=(.*)$', a))[[1]]
	if (length(m) != 3 || !(m[2] %in% c(names(opt), names(per_id)))) {
		stop(paste0("unknown argument '", a, "'\n", usage), call.=FALSE)
	}
	if (m[2] %in% names(opt)) {
		opt[[m[2]]] <- m[3]
		next
	}
	p <- regexpr('=', m[3], fixed=TRUE)
	if (p < 2 || p == nchar(m[3])) stop(paste0('--', m[2], " must be <id>=<file>, not '", m[3], "'"), call.=FALSE)
	per_id[[m[2]]][substr(m[3], 1, p - 1)] <- substr(m[3], p + 1, nchar(m[3]))
}
if (opt$dir == '') stop(paste0('--dir is required\n', usage), call.=FALSE)
if (!dir.exists(opt$dir)) stop(paste('folder not found:', opt$dir), call.=FALSE)

suppressPackageStartupMessages({
	library(terra)
	library(sf)
	library(ggplot2)
	library(reshape2)
	library(RColorBrewer)
})
suppressMessages(sf_use_s2(FALSE))
options(warn=1) # warnings when they happen, not at the end

# Symbols of the texts: UTF-8 when the locale has it (or can be switched to it), ASCII otherwise (in a C locale the
# PNG device draws each byte of a non-ASCII character as a dot)
utf8 <- isTRUE(l10n_info()[['UTF-8']])
for (loc in c('C.UTF-8', 'en_US.UTF-8', 'en_GB.UTF-8')) {
	if (utf8) break
	utf8 <- nzchar(suppressWarnings(Sys.setlocale('LC_CTYPE', loc))) && isTRUE(l10n_info()[['UTF-8']])
}
SEP <- if (utf8) paste0(' ', intToUtf8(0xb7), ' ') else ' | '
KM2 <- if (utf8) paste0('km', intToUtf8(0xb2)) else 'km2'
LEQ <- if (utf8) intToUtf8(0x2264) else '<='
GEQ <- if (utf8) intToUtf8(0x2265) else '>='

rd <- opt$dir
read_tab <- function(f) {
	hdr <- scan(f, what='', sep=' ', nlines=1, quiet=TRUE)
	idcols <- intersect(c('wdpaid', 'wdpaid_hft', 'Status'), hdr)
	read.table(f, sep=' ', header=TRUE, na.strings='NA', stringsAsFactors=FALSE, check.names=FALSE,
	           colClasses=setNames(rep('character', length(idcols)), idcols))
}
hri_file <- file.path(rd, 'hri_results.csv')
if (!file.exists(hri_file)) stop(paste('table of the similarity not found:', hri_file), call.=FALSE)
hri <- read_tab(hri_file)
if (!'wdpaid' %in% names(hri)) stop(paste('no column wdpaid in', hri_file), call.=FALSE)
if (!'hft' %in% names(hri)) hri$hft <- suppressWarnings(as.integer(sub('.*_', '', hri$wdpaid_hft)))
if (!'Status' %in% names(hri)) hri$Status <- NA_character_
thdi <- if (file.exists(file.path(rd, 'thdi_res.csv'))) read_tab(file.path(rd, 'thdi_res.csv')) else NULL
stat_col <- intersect(c('AveHFTSim', 'MedianHFTSim'), names(hri))[1]
stat_lab <- if (!is.na(stat_col) && stat_col == 'MedianHFTSim') 'median' else 'mean'
num <- function(df, col) {
	if (!is.na(col) && col %in% names(df)) suppressWarnings(as.numeric(df[[col]])) else rep(NA_real_, nrow(df))
}

# Title of the run: the folders below 'results' (study area / scenario / cut)
run_title <- function(d) {
	p <- strsplit(normalizePath(d, winslash='/'), '/', fixed=TRUE)[[1]]
	i <- tail(which(p == 'results'), 1)
	if (length(i) == 1 && i < length(p) - 1) return(paste(p[(i + 1):(length(p) - 1)], collapse=' / '))
	basename(dirname(normalizePath(d)))
}
title_run <- run_title(rd)

# The only file of the cut (folder above the run) that matches 'pattern', or ''
find_one <- function(pattern) {
	f <- list.files(dirname(normalizePath(rd)), pattern=pattern, full.names=TRUE)
	if (length(f) == 1) f else ''
}

# Numbers of the texts: 3 significant digits, no exponent, thousands with a comma
fmt <- Vectorize(function(x, digits=3) {
	if (is.na(x)) return('NA')
	if (abs(x) >= 10^digits) return(format(round(x), big.mark=',', trim=TRUE))
	format(signif(x, digits), trim=TRUE, scientific=FALSE, drop0trailing=TRUE)
})
status_text <- c(ok='', few_pixels='few reference pixels', too_few_pixels='not computed: too few pixels',
                 no_reference_similarity='no HFT pixel with similarity > 0', no_similar=paste('no similar pixel in', opt$area))
status_short <- c(ok='NA', few_pixels='NA', too_few_pixels='too few pixels', no_reference_similarity='no ref. similarity',
                  no_similar='no similar pixels')
status_label <- function(s, tab=status_text) {
	ifelse(is.na(s) | s == '', '', ifelse(s %in% names(tab), tab[s], gsub('_', ' ', s)))
}

# HFT colours: the hft_colours csv (space-separated 'HFT colour', in the order of the legends), or else the
# colours of moranvar_plots.R without one (RColorBrewer Set1, interpolated beyond 8 HFTs)
hft_colours <- function(id, hfts) {
	f <- if (id %in% names(per_id$colours)) per_id$colours[[id]] else find_one(paste0('^hft_colours_', id, '_.*\\.csv$'))
	tab <- NULL
	if (f != '' && !file.exists(f)) warning(paste('HFT colours not found:', f), call.=FALSE)
	if (f != '' && file.exists(f)) {
		tab <- tryCatch({
			x <- read.table(f, sep=' ', header=TRUE, comment.char='', stringsAsFactors=FALSE, colClasses='character')
			x[!is.na(suppressWarnings(as.integer(x[[1]]))), 1:2]
		}, error=function(err) {
			warning(paste('could not read the HFT colours', f, ':', conditionMessage(err)), call.=FALSE)
			NULL
		})
	}
	kols <- brewer.pal(8, 'Set1')
	mx <- max(c(hfts, 1), na.rm=TRUE)
	kolsn <- if (mx <= length(kols)) kols else colorRampPalette(kols)(mx)
	cols <- setNames(kolsn[pmax(hfts, 1)], hfts)
	ord <- sort(hfts)
	if (!is.null(tab) && nrow(tab) > 0) {
		h <- as.integer(tab[[1]])
		ok <- h %in% hfts & vapply(tab[[2]], function(cl) !inherits(try(col2rgb(cl), silent=TRUE), 'try-error'), TRUE)
		cols[as.character(h[ok])] <- tab[[2]][ok]
		ord <- c(h[h %in% hfts], sort(setdiff(hfts, h)))
	}
	list(order=unique(ord), colours=cols)
}

# Polygons of 'file' near the polygon 'x' (sfc), in the CRS of 'x' and valid; read with a filter on the bounding
# box of 'x' (in the CRS of the file, with a margin) and cropped to it, so that global files are quick
read_near <- function(file, x) {
	lay <- st_layers(file)
	crs <- lay$crs[[1]]
	if (is.na(crs) || is.na(st_crs(x))) {
		e <- st_read(file, layer=lay$name[1], quiet=TRUE)
		return(st_set_crs(st_make_valid(e[!st_is_empty(e), ]), st_crs(x)))
	}
	bb <- st_bbox(x)
	box <- st_segmentize(st_as_sfc(bb), max(bb['xmax'] - bb['xmin'], bb['ymax'] - bb['ymin'], 1) / 50)
	be <- st_bbox(st_transform(box, crs))
	pad <- 0.05 * max(be['xmax'] - be['xmin'], be['ymax'] - be['ymin']) + 1e-6
	be[c('xmin', 'ymin')] <- be[c('xmin', 'ymin')] - pad
	be[c('xmax', 'ymax')] <- be[c('xmax', 'ymax')] + pad
	e <- st_read(file, layer=lay$name[1], wkt_filter=st_as_text(st_as_sfc(be)), quiet=TRUE)
	e <- e[!st_is_empty(e), ]
	if (nrow(e) == 0) return(e)
	e <- suppressWarnings(suppressMessages(st_crop(st_make_valid(e), be)))
	e <- e[!st_is_empty(e), ]
	if (nrow(e) == 0) return(e)
	st_make_valid(st_transform(e, st_crs(x)))
}

# Ecoregion outline in the CRS of r (read once: all the maps share the grid of the ecoregion variables)
eco_cache <- list()
eco_outline <- function(r) {
	key <- crs(r)
	if (key %in% names(eco_cache)) return(eco_cache[[key]])
	out <- NULL
	if (opt$ecoregions != '' && !file.exists(opt$ecoregions)) {
		warning(paste('ecoregions not found:', opt$ecoregions), call.=FALSE)
	}
	if (opt$ecoregions != '' && file.exists(opt$ecoregions)) {
		e <- as.vector(ext(r))
		box <- st_as_sfc(st_bbox(c(xmin=e[[1]], ymin=e[[3]], xmax=e[[2]], ymax=e[[4]]), crs=st_crs(crs(r))))
		out <- tryCatch(st_geometry(read_near(opt$ecoregions, box)), error=function(err) {
			warning(paste('could not read the ecoregions', opt$ecoregions, ':', conditionMessage(err)), call.=FALSE)
			NULL
		})
		if (!is.null(out) && length(out) == 0) out <- NULL
	}
	eco_cache[key] <<- list(out)
	out
}

# Outline of the study area (union of the HFT polygons) in the CRS of r
study_outline <- function(id, r) {
	f <- if (id %in% names(per_id$study)) per_id$study[[id]] else find_one(paste0('^park_segm_', id, '_.*_class\\.shp$'))
	if (f == '') return(NULL)
	if (!file.exists(f)) {
		warning(paste('HFT polygons of', id, 'not found:', f), call.=FALSE)
		return(NULL)
	}
	tryCatch({
		s <- st_make_valid(st_read(f, quiet=TRUE))
		s <- st_union(st_geometry(s[!st_is_empty(s), ]))
		if (is.na(st_crs(s))) st_set_crs(s, st_crs(crs(r))) else st_transform(s, st_crs(crs(r)))
	}, error=function(err) {
		warning(paste('could not read the HFT polygons', f, ':', conditionMessage(err)), call.=FALSE)
		NULL
	})
}

pal <- hcl.colors(100, 'viridis')
zero_col <- '#e2e2e2'
eco_col <- 'grey25'
study_col <- '#e8000b'

# Raster of colours of a similarity matrix: NoData transparent, 0 light grey, (0, 1] viridis
sim_colours <- function(m) {
	cl <- rep('transparent', length(m))
	z <- !is.na(m) & m <= 0
	v <- !is.na(m) & m > 0
	cl[z] <- zero_col
	cl[v] <- pal[pmin(100, pmax(1, ceiling(m[v] * 100)))]
	as.raster(matrix(cl, nrow=nrow(m)))
}

# Title of a panel, above the plot region: a chip with the HFT colour and line1 (bold), then line2 (grey)
panel_title <- function(colour, line1, line2) {
	usr <- par('usr')
	hh <- strheight('H', cex=1)
	chip <- 1.3 * hh
	gap <- 0.6 * hh
	x0 <- mean(usr[1:2]) - (chip + gap + strwidth(line1, cex=1, font=2)) / 2
	y1 <- usr[4] + 3.2 * hh
	rect(x0, y1 - chip / 2, x0 + chip, y1 + chip / 2, col=colour, border='grey20', lwd=0.6, xpd=NA)
	text(x0 + chip + gap, y1, line1, adj=c(0, 0.5), font=2, xpd=NA)
	if (line2 != '') text(mean(usr[1:2]), usr[4] + 1.3 * hh, line2, adj=c(0.5, 0.5), cex=0.8, col='grey30', xpd=NA)
}

# Legend of the maps, in a plot region in inches: similarity ramp (title above, ticks below) and, to its right,
# the 0 class and the outlines drawn
map_legend <- function(eco, study, marker) {
	par(mar=c(0, 0, 0, 0))
	plot.new()
	pin <- par('pin')
	plot.window(c(0, pin[1]), c(0, pin[2]), xaxs='i', yaxs='i')
	items <- c(paste0('0: p ', LEQ, ' 0.001'), if (!is.null(eco)) (if (opt$area == 'the ecoregion') 'ecoregion' else 'area compared'), if (!is.null(study)) 'study area')
	iw <- 0.3 + 0.07 + max(strwidth(items, units='user', cex=0.75))
	rw <- min(3.4, 0.5 * pin[1])
	x0 <- (pin[1] - (rw + 0.45 + iw)) / 2
	x1 <- x0 + rw
	yb <- 0.3 * pin[2]
	yt <- 0.52 * pin[2]
	xs <- seq(x0, x1, length.out=101)
	rect(xs[-101], yb, xs[-1], yt, col=pal, border=NA)
	rect(x0, yb, x1, yt, border='grey30', lwd=0.6)
	ticks <- c(0, 0.25, 0.5, 0.75, 1)
	tx <- x0 + ticks * rw
	segments(tx, yb, tx, yb - 0.05, col='grey30', lwd=0.6)
	text(tx, yb - 0.13, format(ticks), cex=0.72)
	text((x0 + x1) / 2, yt + 0.13, 'Similarity: chi-square p of the Mahalanobis distance to the HFT', cex=0.78)
	lx <- x1 + 0.45
	ys <- pin[2] * c(0.74, 0.5, 0.26)[seq_along(items)]
	rect(lx + 0.06, ys[1] - 0.07, lx + 0.24, ys[1] + 0.07, col=zero_col, border='grey30', lwd=0.6)
	if (!is.null(eco)) segments(lx, ys[2], lx + 0.3, ys[2], col=eco_col, lwd=1.2)
	if (!is.null(study)) {
		y <- ys[length(items)]
		segments(lx, y, lx + 0.3, y, col=study_col, lwd=1.6)
		if (marker) points(lx + 0.15, y, pch=1, cex=1.6, lwd=1.4, col=study_col)
	}
	text(lx + 0.37, ys, items, adj=c(0, 0.5), cex=0.75)
}

# Similarity maps of a study area: one panel per HFT
plot_maps <- function(id, df, hc, out) {
	hfts <- hc$order
	files <- setNames(file.path(rd, paste0('sim_', id, '_hft', hfts, '.tif')), hfts)
	have <- file.exists(files)
	if (!any(have)) {
		warning(paste0('no similarity maps (sim_', id, '_hft<h>.tif) in ', rd, ': ', basename(out), ' not written'),
		        call.=FALSE)
		return(FALSE)
	}
	n <- length(hfts)
	ncol <- if (n <= 3) n else if (n <= 12) 3 else 4
	nrow <- ceiling(n / ncol)
	width_px <- max(1300, ncol * 760 + 120)

	# Regular sample of each map (about the pixels of the panel, at most 1e6 cells): GDAL reads it from the
	# overviews, so it is quick with 60 M-pixel rasters
	r0 <- rast(files[have][1])
	maxcell <- min(1e6, ceiling((2 * width_px / ncol)^2))
	mats <- list()
	for (h in hfts[have]) {
		r <- rast(files[[as.character(h)]])[[1]]
		s <- if (ncell(r) > maxcell) spatSample(r, maxcell, method='regular', as.raster=TRUE) else r
		mats[[as.character(h)]] <- list(m=as.matrix(s, wide=TRUE), e=as.vector(ext(s)))
	}

	# Extent of the data (rows and columns with values in any map)
	m1 <- mats[[1]]
	rows <- Reduce(`|`, lapply(mats, function(x) rowSums(!is.na(x$m)) > 0))
	cols <- Reduce(`|`, lapply(mats, function(x) colSums(!is.na(x$m)) > 0))
	e <- m1$e
	dx <- (e[2] - e[1]) / ncol(m1$m)
	dy <- (e[4] - e[3]) / nrow(m1$m)
	xl <- e[1:2]
	yl <- e[3:4]
	if (any(rows) && any(cols)) {
		xl <- e[1] + dx * (range(which(cols)) + c(-1, 0))
		yl <- e[4] - dy * (rev(range(which(rows))) + c(0, -1))
	}
	pad <- 0.02 * max(diff(xl), diff(yl))
	xl <- xl + c(-pad, pad)
	yl <- yl + c(-pad, pad)

	eco <- eco_outline(r0)
	study <- study_outline(id, r0)
	sbb <- if (!is.null(study)) st_bbox(study) else NULL
	# A comparison area may not contain the study area: then it is not drawn (nor in the legend), and the subtitle
	# says so
	outside <- !is.null(sbb) && !(sbb['xmax'] > xl[1] && sbb['xmin'] < xl[2] && sbb['ymax'] > yl[1] && sbb['ymin'] < yl[2])
	if (outside) {
		study <- NULL
		sbb <- NULL
	}
	# A circle marks the study area when its outline would be too small to be seen
	marker <- !is.null(sbb) && max(sbb['xmax'] - sbb['xmin'], sbb['ymax'] - sbb['ymin']) < 0.04 * diff(xl)

	# Sizes in inches at 200 dpi: each panel keeps the aspect of the extent, plus two lines of title
	head_in <- 0.85
	leg_in <- 0.85
	pw <- (width_px / 200 - 0.3) / ncol
	ph <- min(pw * diff(yl) / diff(xl), 6) + 0.62
	height_px <- round(200 * (head_in + nrow * ph + leg_in))
	png(out, width=width_px, height=height_px, res=200,
	    type=if (capabilities('cairo')) 'cairo' else getOption('bitmapType'))
	on.exit(if (dev.cur() > 1) dev.off())
	lay <- matrix(c(seq_len(n), rep(0, nrow * ncol - n)), nrow=nrow, byrow=TRUE)
	layout(rbind(lay, rep(n + 1, ncol)), heights=c(rep(ph, nrow), leg_in))
	par(omi=c(0, 0.15, head_in, 0.15))

	for (h in hfts) {
		k <- as.character(h)
		row <- df[df$hft == h, , drop=FALSE][1, ]
		st <- row$Status
		par(mar=c(0.25, 0.25, 2.9, 0.25))
		plot.new()
		plot.window(xl, yl, asp=1, xaxs='i', yaxs='i')
		box(col='grey70', lwd=0.6)
		line1 <- paste0('HFT ', h)
		line2 <- ''
		note <- ''
		if (!is.null(mats[[k]])) {
			x <- mats[[k]]
			rasterImage(sim_colours(x$m), x$e[1], x$e[3], x$e[2], x$e[4], interpolate=FALSE)
			line1 <- paste0(line1, SEP, 'T = ', fmt(num(row, stat_col)), SEP, 'HSR = ', fmt(num(row, 'HSR')))
			sa <- num(row, 'SimAreaEco_km2')
			np <- num(row, 'NrSimLcpPatEco')
			if (!is.na(sa)) {
				line2 <- paste0('similar: ', fmt(sa), ' ', KM2, ' in ', fmt(np), ' patch', if (is.na(np) || np != 1) 'es')
			}
			if (!is.na(st) && st != 'ok') {
				line2 <- if (st == 'few_pixels' && line2 != '') paste0(status_label(st), SEP, line2) else status_label(st)
			}
		} else {
			note <- if (!is.na(st) && st != '') status_label(st) else 'no similarity map'
		}
		if (!is.null(eco)) plot(eco, add=TRUE, border=eco_col, col=NA, lwd=0.6)
		if (!is.null(study)) {
			plot(study, add=TRUE, border=study_col, col=NA, lwd=1.1)
			if (marker) {
				points(mean(sbb[c('xmin', 'xmax')]), mean(sbb[c('ymin', 'ymax')]), pch=1, cex=2.6, lwd=1.6, col=study_col)
			}
		}
		if (note != '') {
			w <- strwidth(note, cex=0.95) / 2 + strwidth('m', cex=0.95)
			hgt <- strheight(note, cex=0.95)
			rect(mean(xl) - w, mean(yl) - hgt, mean(xl) + w, mean(yl) + hgt, col='white', border='grey70', lwd=0.6)
			text(mean(xl), mean(yl), note, col='grey35', cex=0.95)
		}
		panel_title(hc$colours[[k]], line1, line2)
	}
	map_legend(eco, study, marker)

	unit <- if (isTRUE(linearUnits(r0) == 1)) ' m' else ''
	mtext(paste0('Similarity of the HFTs of study area ', id, ' in ', opt$area), side=3, outer=TRUE, line=2.3,
	      cex=1.15, font=2)
	mtext(paste0(title_run, SEP, 'pixel ', fmt(res(r0)[1], 6), unit, if (outside) paste0(SEP, 'study area outside the map')),
	      side=3, outer=TRUE, line=1.15, cex=0.8, col='grey25')
	mtext(paste0('T: ', stat_lab, ' similarity of the HFT pixels', SEP, 'HSR: pixels of the similar patches ', GEQ,
	             ' HFT (or of the largest one) / HFT pixels'), side=3, outer=TRUE, line=0.2, cex=0.72, col='grey25')

	# The data of the plot, for the users: one row per panel, with the values drawn (a regular sample of the map,
	# summarised) and the GeoTIFF with all its pixels
	rows <- lapply(seq_along(hfts), function(i) {
		h <- hfts[i]
		k <- as.character(h)
		row <- df[df$hft == h, , drop=FALSE][1, ]
		v <- if (!is.null(mats[[k]])) as.vector(mats[[k]]$m) else numeric(0)
		v <- v[!is.na(v)]
		q <- if (length(v)) quantile(v, c(0.5, 0.9), names=FALSE) else c(NA, NA)
		data.frame(id=id, panel=i, hft=h, colour=hc$colours[[k]], map=if (have[i]) basename(files[[k]]) else NA,
		           threshold=num(row, stat_col), HSR=num(row, 'HSR'), SimAreaEco_km2=num(row, 'SimAreaEco_km2'),
		           NrSimLcpPatEco=num(row, 'NrSimLcpPatEco'), Status=if (is.na(row$Status)) NA else row$Status,
		           cells_drawn=length(v), cells_zero=sum(v == 0), similarity_mean=if (length(v)) mean(v) else NA,
		           similarity_p50=q[1], similarity_p90=q[2], similarity_max=if (length(v)) max(v) else NA,
		           check.names=FALSE)
	})
	write_data(do.call(rbind, rows), out)
	TRUE
}

# CSV next to a plot (X.png -> X.csv), as the other tables of eHabitat+: separated by spaces, NA, 6 digits
# (texts with spaces, quotes or # quoted, as the colours, so that read.table reads them: # starts a comment)
write_data <- function(d, png_file) {
	for (cl in names(d)) if (is.numeric(d[[cl]])) d[[cl]] <- signif(d[[cl]], 6)
	f <- sub('\\.png$', '.csv', png_file)
	q <- which(vapply(d, function(x) is.character(x) && any(grepl('[[:space:]#"\']', x)), logical(1)))
	writeLines(paste(names(d), collapse=' '), f)
	write.table(d, f, sep=' ', row.names=FALSE, col.names=FALSE, quote=if (length(q)) q else FALSE, na='NA',
	            append=TRUE)
}

metric_labels <- c(AveHFTSim='T: mean similarity of the HFT pixels',
                   MedianHFTSim='T: median similarity of the HFT pixels',
                   HSR='HSR (dashed: patch as large as the HFT)', NrSimLcpPatEco='Patches of similar pixels',
                   SimAreaEco_km2=paste0('Similar area in ', opt$area, ' (', KM2, ')'),
                   NrPixLargLcpPat_NrPixHFT_ratio='Largest similar patch / HFT',
                   Aggregation='Aggregation (similar pixels per patch)')

# Bar charts of the metrics of a study area: one panel per metric, one bar per HFT in its colour
plot_metrics <- function(id, df, hc, out) {
	mets <- intersect(c(stat_col, 'HSR', 'NrSimLcpPatEco', 'SimAreaEco_km2', 'NrPixLargLcpPat_NrPixHFT_ratio',
	                    'Aggregation'), names(df))
	if (length(mets) == 0) {
		warning(paste('no similarity metrics in', hri_file, ':', basename(out), 'not written'), call.=FALSE)
		return(FALSE)
	}
	hfts <- hc$order
	d <- data.frame(hft=df$hft, Status=df$Status, setNames(lapply(mets, function(mt) num(df, mt)), mets),
	                check.names=FALSE)
	d <- d[!duplicated(d$hft) & d$hft %in% hfts, , drop=FALSE]
	long <- melt(d, id.vars=c('hft', 'Status'), measure.vars=mets, variable.name='metric', value.name='value')
	long$metric <- factor(long$metric, levels=mets)
	long$HFT <- factor(long$hft, levels=hfts)
	long$lab <- if (length(hfts) > 12) '' else ifelse(is.na(long$value), '', fmt(long$value)) # values on the bars
	na <- long[is.na(long$value), ]
	na$txt <- status_label(na$Status, status_short)
	na$txt[na$txt == ''] <- 'NA'
	# Panels without a bar above 0 get a y range of 0-1
	top <- tapply(long$value, long$metric, function(v) if (all(is.na(v)) || max(v, na.rm=TRUE) <= 0) 1 else NA)
	blank <- data.frame(metric=factor(names(top), levels=mets), value=as.numeric(top))
	blank <- blank[!is.na(blank$value), ]
	ref1 <- data.frame(metric=factor(intersect(c('HSR', 'NrPixLargLcpPat_NrPixHFT_ratio'), mets), levels=mets), y=1)

	# 12 inches wide (16 with more than 12 HFTs), narrower with few HFTs; panel titles wrapped to the panel width
	nc <- min(3, length(mets))
	width <- max(8, min(if (length(hfts) > 12) 16 else 12, 1.2 + nc * (1.3 + 0.38 * length(hfts))))
	wrapped <- sapply(metric_labels, function(x) paste(strwrap(x, width=floor(13 * (width - 0.8) / nc)), collapse='\n'))

	p <- ggplot(long, aes(x=HFT, y=value)) +
		geom_hline(yintercept=0, colour='grey60', linewidth=0.3) +
		(if (nrow(ref1) > 0) geom_hline(data=ref1, aes(yintercept=y), colour='grey40', linetype='dashed',
		                                linewidth=0.35)) +
		geom_col(aes(fill=HFT), width=0.72, colour='grey25', linewidth=0.2, na.rm=TRUE) +
		geom_text(aes(label=lab), vjust=-0.35, size=2.5, colour='grey20', na.rm=TRUE) +
		(if (nrow(na) > 0) geom_text(data=na, aes(x=HFT, y=0, label=txt), angle=90, hjust=-0.08, size=2.4,
		                             colour='grey45', fontface='italic')) +
		(if (nrow(blank) > 0) geom_blank(data=blank, aes(x=NULL, y=value))) +
		facet_wrap(~ metric, ncol=nc, scales='free_y', labeller=as_labeller(wrapped)) +
		scale_fill_manual(values=hc$colours[as.character(hfts)], guide='none', drop=FALSE) +
		scale_x_discrete(drop=FALSE) +
		scale_y_continuous(expand=expansion(mult=c(0.02, 0.14))) +
		labs(title=paste0('Similarity metrics of the HFTs of study area ', id),
		     subtitle=paste0(title_run, if (any(!is.na(num(df, 'PixelSize_m'))))
		                     paste0(SEP, 'pixel ', fmt(max(num(df, 'PixelSize_m'), na.rm=TRUE), 6), ' m')),
		     x='HFT', y=NULL) +
		theme_bw(base_size=12) +
		theme(strip.text=element_text(size=10), strip.background=element_rect(fill='grey95'),
		      panel.grid.major.x=element_blank(), plot.caption=element_text(hjust=0, size=9, colour='grey25'),
		      plot.title=element_text(face='bold'))

	# Caption: status of the HFTs that are not 'ok' and the THDI of the study area
	cap <- character(0)
	bad <- d[!is.na(d$Status) & d$Status != 'ok', ]
	if (nrow(bad) > 0) cap <- c(cap, paste0('Status: ', paste0('HFT ', bad$hft, ' ', bad$Status, collapse='; ')))
	if (!is.null(thdi) && id %in% thdi$wdpaid) {
		tr <- thdi[thdi$wdpaid == id, ][1, ]
		v <- intersect(c('TerrEco', 'MarEco', 'NrHFT', 'SIH', 'THD', 'THR', 'THDI'), names(tr))
		cap <- c(cap, paste0('THDI of the study area: ',
		                     paste0(v, ' ', sapply(v, function(cl) fmt(num(tr, cl))), collapse=SEP)))
	}
	if (length(cap) > 0) p <- p + labs(caption=paste(cap, collapse='\n'))
	ggsave(out, p, width=width, height=1.6 + 3 * ceiling(length(mets) / 3) + 0.25 * length(cap), dpi=200, bg='white')
	# The data of the plot: one row per HFT (in the order of the bars) with the metrics drawn and its colour
	dd <- d[order(match(d$hft, hfts)), , drop=FALSE]
	write_data(data.frame(id=id, hft=dd$hft, colour=unname(hc$colours[as.character(dd$hft)]), dd[, mets, drop=FALSE],
	                      Status=dd$Status, check.names=FALSE), out)
	TRUE
}

ids <- unique(hri$wdpaid[!is.na(hri$wdpaid)])
if (length(ids) == 0) stop(paste('no rows in', hri_file), call.=FALSE)
written <- character(0)
for (id in ids) {
	df <- hri[!is.na(hri$wdpaid) & hri$wdpaid == id & !is.na(hri$hft), , drop=FALSE]
	if (nrow(df) == 0) next
	hc <- hft_colours(id, sort(unique(df$hft)))
	for (f in file.path(rd, paste0(c('sim_maps_', 'sim_metrics_'), id, '.png'))) {
		ok <- tryCatch(if (grepl('^sim_maps_', basename(f))) plot_maps(id, df, hc, f) else plot_metrics(id, df, hc, f),
		               error=function(err) {
			warning(paste('could not plot', basename(f), ':', conditionMessage(err)), call.=FALSE)
			FALSE
		})
		if (isTRUE(ok)) written <- c(written, basename(f))
	}
}
cat(paste0('Plots: ', paste(written, collapse=' '), '\n'))
