#### eHabitat+: Total Habitat Diversity Index (THDI) of each study area, from the similarity of its HFTs
#### License: CC BY-SA 3.0
#### Run by segm_optim/similarity.py after the similarity of the HFTs in the ecoregion (hri_results.csv):
####	Rscript thdi.R --hri=<hri_results.csv> --out=<thdi_res.csv> --hfts=<id>=<class.shp> [--hfts=<id2>=<shp2> ...]
####	               [--ecoregions=<file>] [--eco-field=ECO_ID] [--marine=<file>] [--marine-field=ECO_CODE]
#### Indices of each study area (wdpaid, the ID of a --hfts), the columns of thdi_res.csv:
####	TerrEco: number of terrestrial ecoregions (distinct --eco-field values of --ecoregions) that intersect the HFT
####	         polygons of the study area; 1 without --ecoregions (the ecoregion variables are taken as one ecoregion)
####	MarEco: the same with the marine ecoregions of --marine (e.g. MEOW, field ECO_CODE); 0 without them
####	TotEco: TerrEco + MarEco
####	NrHFT: number of HFTs of the study area (its rows of hri_results.csv)
####	SIH: Shannon index of the HFTs, -sum(p ln p), p = share of each HFT in the area of all of them (column
####	     HFTArea_ha, or NrPixHFT in older tables); 1 when it is 0 (a single HFT)
####	THD: Total Habitat Diversity = TotEco x SIH: many ecoregions and many HFTs of even size give a high THD
####	THR: Total Habitat Replaceability = median of the HSR of the HFTs / TerrEco. The HSR of an HFT is the number of
####	     pixels of the patches of similar pixels at least as large as the HFT (or of the largest patch, if none is)
####	     divided by the pixels of the HFT, in its ecoregion: HFTs with a lot of similar land around have a high THR
####	THDI: Total Habitat Diversity Index = THD / THR: high in study areas with a diverse habitat that is hard to
####	      find elsewhere in their ecoregions
#### History: written in 2021 for the global run of eHabitat+ (Martinez-Lopez et al. 2023, Ecological Informatics,
#### https://www.sciencedirect.com/science/article/pii/S157495412300119X), with 1 km Mollweide pixels, the WWF
#### terrestrial ecoregions, the marine ecoregions of the world (MEOW) and one class shapefile per protected area.
#### SIH, THD and THR were first called SIHFT, EcoHab and MedianSim (THR was THI for a while); THR is divided by
#### TerrEco since 2021-11 and an SIH of 0 is set to 1 since then too. Rewritten in 2026 for the similarity step of
#### the HFT module (segm_optim/similarity.py): same formulas, inputs as arguments, ecoregions in any CRS
#### (transformed to the CRS of the HFTs), invalid polygons repaired, NA values allowed and SIH from the HFT areas.

usage <- paste(
	'Usage: Rscript thdi.R --hri=<hri_results.csv> --out=<thdi_res.csv> --hfts=<id>=<class.shp> [--hfts=<id2>=<shp2> ...]',
	'                      [--ecoregions=<file>] [--eco-field=ECO_ID] [--marine=<file>] [--marine-field=ECO_CODE]',
	'Writes the Total Habitat Diversity Index (THDI) of each study area (the IDs of --hfts) from the hri_results.csv',
	'of segm_optim/similarity.py, the HFT polygons (class shapefile) and the polygons of the ecoregions.', sep='\n')

# Arguments: one argv element each, --name=value (the value may contain spaces and '='); --hfts can be repeated
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	cat(usage, '\n', sep='')
	quit(save='no', status=1)
}
if (any(args %in% c('-h', '--help'))) {
	cat(usage, '\n', sep='')
	quit(save='no', status=0)
}
opt <- list(hri='', out='', ecoregions='', 'eco-field'='ECO_ID', marine='', 'marine-field'='ECO_CODE')
hfts <- character(0)
for (a in args) {
	m <- regmatches(a, regexec('^--([A-Za-z-]+)=(.*)$', a))[[1]]
	if (length(m) != 3 || !(m[2] %in% c(names(opt), 'hfts'))) {
		stop(paste0("unknown argument '", a, "'\n", usage), call.=FALSE)
	}
	if (m[2] != 'hfts') {
		opt[[m[2]]] <- m[3]
		next
	}
	p <- regexpr('=', m[3], fixed=TRUE)
	if (p < 2 || p == nchar(m[3])) stop(paste0("--hfts must be <id>=<class.shp>, not '", m[3], "'"), call.=FALSE)
	id <- substr(m[3], 1, p - 1)
	if (id %in% names(hfts)) stop(paste0("study area '", id, "' given twice in --hfts"), call.=FALSE)
	hfts[id] <- substr(m[3], p + 1, nchar(m[3]))
}

suppressPackageStartupMessages(library(sf))
suppressMessages(sf_use_s2(FALSE)) # planar intersects, once the ecoregions are in the (projected) CRS of the HFTs
options(warn=1) # warnings when they happen, not at the end

if (opt$hri == '' || opt$out == '' || length(hfts) == 0) {
	stop(paste0('--hri, --out and --hfts are required\n', usage), call.=FALSE)
}
if (!file.exists(opt$hri)) stop(paste('table of the similarity not found:', opt$hri), call.=FALSE)

# hri_results.csv: space-separated, NA for missing; the IDs as text (so that '01' stays '01')
hdr <- scan(opt$hri, what='', sep=' ', nlines=1, quiet=TRUE)
idcols <- intersect(c('wdpaid', 'wdpaid_hft'), hdr)
hri <- read.table(opt$hri, sep=' ', header=TRUE, na.strings='NA', stringsAsFactors=FALSE, check.names=FALSE,
                  colClasses=setNames(rep('character', length(idcols)), idcols))
if (!'wdpaid' %in% names(hri)) stop(paste('no column wdpaid in', opt$hri), call.=FALSE)
num <- function(df, col) if (col %in% names(df)) suppressWarnings(as.numeric(df[[col]])) else rep(NA_real_, nrow(df))

# Polygons of 'file' near the polygons 'x', in the CRS of 'x' and valid. Only the features that intersect the
# bounding box of 'x' (with a margin, in the CRS of the file) are read, and they are cropped to it before the
# transformation, so that a global file (e.g. RESOLVE or MEOW) is quick and its far polygons are not projected
read_near <- function(file, x) {
	lay <- st_layers(file)
	if (length(lay$name) == 0) stop('no layers')
	crs <- lay$crs[[1]]
	if (is.na(crs) || is.na(st_crs(x))) {
		warning(paste0('no CRS in ', if (is.na(crs)) file else 'the HFT polygons', ': taken as the same CRS'),
		        call.=FALSE)
		e <- st_read(file, layer=lay$name[1], quiet=TRUE)
		e <- st_make_valid(e[!st_is_empty(e), ])
		return(st_set_crs(e, st_crs(x)))
	}
	bb <- st_bbox(x)
	side <- max(bb['xmax'] - bb['xmin'], bb['ymax'] - bb['ymin'], 1)
	box <- st_segmentize(st_as_sfc(bb), side / 50)
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

# ID column of the ecoregions: 'field', or else the first integer column (NULL: each polygon is an ecoregion)
id_field <- function(e, field, file) {
	if (field %in% names(e)) return(field)
	# FID: each polygon is an ecoregion, as in similarity.py when the file has no integer ID field
	if (identical(field, 'FID')) return(NULL)
	cols <- setdiff(names(e), attr(e, 'sf_column'))
	whole <- function(v) is.integer(v) || (is.numeric(v) && all(v %% 1 == 0, na.rm=TRUE))
	int <- cols[vapply(cols, function(cl) whole(e[[cl]]), TRUE)]
	alt <- if (length(int) > 0) int[1] else NULL
	warning(paste0("no field '", field, "' in ", file, ': ', if (is.null(alt)) 'each polygon counted as an ecoregion'
	               else paste0("field '", alt, "' used instead")), call.=FALSE)
	alt
}

# Number of ecoregions of 'file' (distinct values of 'field') that intersect the polygons 'x'; NA if it fails
count_ecoregions <- function(file, field, x) {
	tryCatch({
		e <- read_near(file, x)
		if (nrow(e) == 0) return(0L)
		hit <- lengths(st_intersects(e, x)) > 0
		f <- id_field(e, field, file)
		if (is.null(f)) sum(hit) else length(unique(e[[f]][hit]))
	}, error=function(err) {
		warning(paste0('could not intersect ', file, ' with the HFTs: ', conditionMessage(err)), call.=FALSE)
		NA_integer_
	})
}

has_eco <- opt$ecoregions != '' && file.exists(opt$ecoregions)
if (!has_eco) {
	warning(paste0(if (opt$ecoregions == '') 'no --ecoregions' else paste('ecoregions not found:', opt$ecoregions),
	               ': TerrEco = 1 (the ecoregion variables are taken as one ecoregion)'), call.=FALSE)
}
has_mar <- opt$marine != '' && file.exists(opt$marine)
if (opt$marine != '' && !has_mar) warning(paste('marine ecoregions not found:', opt$marine, '(MarEco = 0)'), call.=FALSE)

ids <- names(hfts)
other <- setdiff(unique(hri$wdpaid[!is.na(hri$wdpaid)]), ids)
if (length(other) > 0) warning(paste('study areas of', opt$hri, 'without --hfts, left out:', paste(other, collapse=', ')),
                               call.=FALSE)

res <- data.frame(wdpaid=ids, TerrEco=NA_integer_, MarEco=NA_integer_, TotEco=NA_integer_, NrHFT=NA_integer_,
                  SIH=NA_real_, THD=NA_real_, THR=NA_real_, THDI=NA_real_, stringsAsFactors=FALSE)
for (i in seq_along(ids)) {
	shp <- hfts[[i]]
	pa <- NULL
	if (!file.exists(shp)) {
		warning(paste('HFT polygons of', ids[i], 'not found:', shp), call.=FALSE)
	} else {
		pa <- tryCatch(st_read(shp, quiet=TRUE), error=function(err) {
			warning(paste('could not read', shp, ':', conditionMessage(err)), call.=FALSE)
			NULL
		})
		if (!is.null(pa)) {
			pa <- st_make_valid(pa[!st_is_empty(pa), ])
			if (nrow(pa) == 0) {
				warning(paste('no HFT polygons in', shp), call.=FALSE)
				pa <- NULL
			}
		}
	}

	terr <- 1L
	if (has_eco) terr <- if (is.null(pa)) NA_integer_ else count_ecoregions(opt$ecoregions, opt[['eco-field']], pa)
	mar <- 0L
	if (has_mar) mar <- if (is.null(pa)) NA_integer_ else count_ecoregions(opt$marine, opt[['marine-field']], pa)
	tot <- terr + mar

	df <- hri[!is.na(hri$wdpaid) & hri$wdpaid == ids[i], , drop=FALSE]
	nhft <- nrow(df)

	# Shares of the HFTs in the area of all of them (HFTs without area left out): their area in ha (10 m map) or,
	# without it, their pixels
	area <- num(df, 'HFTArea_ha')
	if (!any(!is.na(area) & area > 0)) area <- num(df, 'NrPixHFT')
	area <- area[!is.na(area) & area > 0]
	sih <- NA_real_
	if (length(area) > 0) {
		p <- area / sum(area)
		sih <- -sum(p * log(p))
		if (sih == 0) sih <- 1
	}
	thd <- tot * sih

	hsr <- num(df, 'HSR')
	thr <- if (is.na(terr) || terr == 0 || all(is.na(hsr))) NA_real_ else median(hsr, na.rm=TRUE) / terr
	thdi <- if (is.na(thr) || thr == 0) NA_real_ else thd / thr

	res[i, -1] <- list(terr, mar, tot, nhft, sih, thd, thr, thdi)
}

numc <- sapply(res, is.numeric)
res[numc] <- lapply(res[numc], signif, 6)
dir.create(dirname(opt$out), showWarnings=FALSE, recursive=TRUE)
write.table(res, opt$out, sep=' ', row.names=FALSE, quote=FALSE, na='NA')

f <- function(x) ifelse(is.na(x), 'NA', format(signif(x, 4), trim=TRUE))
cat(paste0('THDI of ', nrow(res), ' study area', if (nrow(res) != 1) 's', ': ', opt$out, '\n'))
for (i in seq_len(nrow(res))) {
	r <- res[i, ]
	cat(paste0('  ', r$wdpaid, ': TerrEco ', f(r$TerrEco), ', MarEco ', f(r$MarEco), ', TotEco ', f(r$TotEco), ', NrHFT ',
	           f(r$NrHFT), ', SIH ', f(r$SIH), ', THD ', f(r$THD), ', THR ', f(r$THR), ', THDI ', f(r$THDI), '\n'))
}
