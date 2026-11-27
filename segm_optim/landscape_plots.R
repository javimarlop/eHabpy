#### eHabitat+: plots of the landscape metrics of several maps of patches (run by landscape_metrics.py)
#### License: CC BY-SA 3.0
#### Rscript landscape_plots.R --dir=<folder of a run> [--mode=habitat|classes]
#### Reads metrics_class.csv (and metrics_landscape.csv) of the folder and writes there:
####	metrics_trend.png: one panel per class metric, its value in every map (in their order: the dates), one line per
####	                   class (the habitat in mode habitat)
####	landscape_trend.png: the same with the landscape metrics (mode classes)
#### and, next to each plot, its data (X.png -> X.csv, separated by spaces), as the other plots of eHabitat+.

suppressPackageStartupMessages({library(ggplot2); library(reshape2)})
options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
opt <- list(dir='', mode='habitat')
for (a in args) {
	m <- regmatches(a, regexec('^--([a-z]+)=(.*)$', a))[[1]]
	if (length(m) != 3 || !(m[2] %in% names(opt))) stop(paste0("unknown argument '", a, "'"), call.=FALSE)
	opt[[m[2]]] <- m[3]
}
if (opt$dir == '' || !dir.exists(opt$dir)) stop('--dir must be the folder of a run', call.=FALSE)

# Names of the metrics (FRAGSTATS) in the panels, with their units
labels <- c(ca='CA: class area (ha)', pland='PLAND: % of the landscape', np='NP: number of patches',
            pd='PD: patches per 100 ha', lpi='LPI: largest patch (% of landscape)', area_mn='AREA_MN: mean patch area (ha)',
            area_cv='AREA_CV: CV of patch area (%)', te='TE: total edge (m)', ed='ED: edge density (m/ha)',
            lsi='LSI: landscape shape index', shape_mn='SHAPE_MN: mean shape index', frac_mn='FRAC_MN: mean fractal dimension',
            gyrate_mn='GYRATE_MN: mean radius of gyration (m)', enn_mn='ENN_MN: mean nearest-neighbour distance (m)',
            enn_cv='ENN_CV: CV of nearest-neighbour distance (%)', tca='TCA: total core area (ha)',
            cpland='CPLAND: core area (% of landscape)', cai_mn='CAI_MN: mean core area index (%)', ai='AI: aggregation index (%)',
            clumpy='CLUMPY: clumpiness', pladj='PLADJ: like adjacencies (%)', cohesion='COHESION: patch cohesion',
            division='DIVISION: division index', mesh='MESH: effective mesh size (ha)', split='SPLIT: splitting index',
            contag='CONTAG: contagion (%)', shdi="SHDI: Shannon's diversity", shei="SHEI: Shannon's evenness",
            sidi="SIDI: Simpson's diversity", pr='PR: number of classes', iji='IJI: interspersion (%)',
            ent='ENT: entropy (bits)', condent='CONDENT: conditional entropy (bits)', mutinf='MUTINF: mutual information (bits)',
            relmutinf='RELMUTINF: relative mutual information')

read_tab <- function(f) {
	hdr <- scan(f, what='', sep=' ', nlines=1, quiet=TRUE, quote='"')
	txt <- intersect(c('map', 'class_name'), hdr)
	read.table(f, sep=' ', header=TRUE, na.strings='NA', quote='"', stringsAsFactors=FALSE, check.names=FALSE,
	           colClasses=setNames(rep('character', length(txt)), txt), comment.char='')
}

# The maps on the x axis: numbers (years) when every label is one, else in the order of the run
x_axis <- function(maps) {
	n <- suppressWarnings(as.numeric(maps))
	if (all(!is.na(n))) n else factor(maps, levels=unique(maps))
}

write_data <- function(d, png_file) {
	write.table(d, sub('\\.png$', '.csv', png_file), sep=' ', row.names=FALSE, quote=which(vapply(d, is.character, TRUE)),
	            na='NA')
}

trend <- function(d, id_cols, out, title) {
	mets <- intersect(names(labels), setdiff(names(d), id_cols))
	mets <- mets[vapply(mets, function(m) any(is.finite(d[[m]])), TRUE)]
	if (length(mets) == 0) {
		warning(paste('no metrics to plot in', out), call.=FALSE)
		return(invisible(NULL))
	}
	long <- melt(d[, c(id_cols, mets)], id.vars=id_cols, variable.name='metric', value.name='value')
	long$metric <- factor(long$metric, levels=mets, labels=labels[mets])
	long$x <- x_axis(long$map)
	numeric_x <- is.numeric(long$x)
	grp <- if ('class_name' %in% id_cols) 'class_name' else NULL
	nc <- min(4, length(mets))
	p <- ggplot(long, aes(x=x, y=value, group=if (is.null(grp)) 1 else .data[[grp]],
	                      colour=if (is.null(grp)) NULL else .data[[grp]])) +
		geom_line(linewidth=0.7, na.rm=TRUE) + geom_point(size=1.8, na.rm=TRUE) +
		facet_wrap(~ metric, ncol=nc, scales='free_y', labeller=label_wrap_gen(width=28)) +
		labs(title=title, x=if (numeric_x) 'Map (date)' else 'Map', y=NULL, colour='Class') +
		theme_bw(base_size=9) + theme(legend.position='bottom', strip.text=element_text(size=7.5))
	if (is.null(grp) || length(unique(long[[grp]])) == 1) p <- p + guides(colour='none')
	ggsave(out, p, width=2.6 * nc + 0.6, height=2.1 * ceiling(length(mets) / nc) + 1.1, dpi=200, bg='white')
	write_data(d[, c(id_cols, mets)], out)
}

cls <- file.path(opt$dir, 'metrics_class.csv')
if (file.exists(cls)) {
	d <- read_tab(cls)
	if (opt$mode == 'habitat') d <- d[d$class == 1, ]
	trend(d, c('map', 'class', 'class_name'), file.path(opt$dir, 'metrics_trend.png'),
	      if (opt$mode == 'habitat') 'Landscape metrics of the habitat in every map' else 'Class metrics in every map')
}
lan <- file.path(opt$dir, 'metrics_landscape.csv')
if (file.exists(lan)) {
	d <- read_tab(lan)
	trend(d, 'map', file.path(opt$dir, 'landscape_trend.png'), 'Landscape-level metrics in every map')
}
