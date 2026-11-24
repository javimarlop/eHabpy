#### eHabitat+ HFT model: subclasses of the HFTs (HFT_SUBCLASSES in config.py), used by moranvar_plots.R
#### License: CC BY-SA 3.0
#
# The HFTs group whole segments, by the mean of their pixels, so they cannot separate what differs inside a segment
# (e.g. a wet meadow that is only part of a segment). The subclasses divide the PIXELS of each HFT: the subclass of
# an HFT whose division most reduces the sum of squares of its pixels (SSE; variables standardized with the mean and
# the standard deviation of the pixels of all the HFTs) is divided in two by k-means, again and again (bisection),
# while a division reduces the SSE at least as much as merging the two most similar HFTs of the cut would increase
# it (the cost of Ward of that merge, n_a n_b / (n_a + n_b) |m_a - m_b|^2): a subclass separates its pixels at
# least as much as the HFTs are separated from each other, so the number of subclasses follows from the cut, without
# another parameter. The subclasses stay inside their HFT: code = 100 x HFT + subclass (1 = the largest).

# Cost of Ward of merging the two most similar HFTs: the increase of the SSE of the pixels (z: standardized values of
# the pixels with data in every variable, cls: their HFT). Only the HFTs that could be divided (2 x MIN_PX pixels or
# more) count, if there are two: an HFT of a few pixels would give a cost near 0 (and subclasses in every HFT). Inf
# with a single HFT (no subclasses then).
subclass_threshold <- function(z, cls, min_px=10) {
	g <- factor(cls)
	if (nlevels(g) < 2) return(Inf)
	n <- tabulate(g)
	mu <- rowsum(z, g) / n
	cost <- outer(n, n) / outer(n, n, '+') * as.matrix(dist(mu))^2
	diag(cost) <- Inf
	big <- n >= 2 * min_px
	if (sum(big) >= 2) cost <- cost[big, big, drop=FALSE]
	min(cost)
}

# Best division in two of the pixels X by k-means: k-means with 10 random starts on at most NSAMPLE pixels, then
# Lloyd with all of them from those centres. list(gain = reduction of the SSE, second = the pixels of the second
# half), or NULL when X cannot be divided (too few pixels, a half with fewer than MIN_PX, all values equal).
split_in_two <- function(x, seed, nsample=50000, min_px=10) {
	if (nrow(x) < 2 * min_px) return(NULL)
	set.seed(seed)
	s <- if (nrow(x) > nsample) x[sample.int(nrow(x), nsample), , drop=FALSE] else x
	k0 <- tryCatch(suppressWarnings(kmeans(s, 2, nstart=10, iter.max=50)), error=function(e) NULL)
	if (is.null(k0)) return(NULL)
	k1 <- tryCatch(suppressWarnings(kmeans(x, k0$centers, iter.max=100, algorithm='Lloyd')), error=function(e) NULL)
	if (is.null(k1) || length(k1$size) < 2 || min(k1$size) < min_px) return(NULL)
	list(gain=k1$totss - k1$tot.withinss, second=k1$cluster == 2)
}

# Subclasses of the HFTs by bisection. z: standardized values of the pixels with data in every variable (rows), cls:
# their HFT, threshold: least reduction of the SSE of a division (subclass_threshold). Returns the subclass of every
# row as code = 100 x HFT + subclass, the subclasses of each HFT numbered from the largest (1); at most 99 per HFT.
bisect_hfts <- function(z, cls, threshold, max_sub=99) {
	hft <- sort(unique(cls))
	cl <- match(cls, hft)                   # class of each pixel: 1..number of classes
	parent <- seq_along(hft)                # HFT (index in hft) of each class
	seed <- 0
	cand <- lapply(seq_along(hft), function(c) { seed <<- seed + 1; split_in_two(z[cl == c, , drop=FALSE], seed) })
	repeat {
		gains <- vapply(cand, function(r) if (is.null(r)) -Inf else r$gain, numeric(1))
		best <- which.max(gains)
		if (!length(best) || !is.finite(gains[best]) || gains[best] < threshold) break
		if (sum(parent == parent[best]) >= max_sub) { cand[best] <- list(NULL); next }
		rows <- which(cl == best)
		new <- length(parent) + 1
		cl[rows[cand[[best]]$second]] <- new
		parent[new] <- parent[best]
		for (c in c(best, new)) {
			seed <- seed + 1
			cand[c] <- list(split_in_two(z[cl == c, , drop=FALSE], seed))
		}
	}
	# subclasses of each HFT numbered from the largest
	size <- tabulate(cl, nbins=length(parent))
	num <- integer(length(parent))
	for (p in unique(parent)) {
		cs <- which(parent == p)
		num[cs[order(-size[cs], cs)]] <- seq_along(cs)
	}
	100L * as.integer(hft[parent[cl]]) + num[cl]
}

# Subclass of the pixels with data in only some of the variables: that of the nearest centre (mean distance in the
# variables with data) among the subclasses of their HFT; NA without data in any variable.
subclass_partial <- function(zp, hft_p, cent, cent_code) {
	out <- rep(NA_integer_, nrow(zp))
	cent_hft <- cent_code %/% 100L
	for (h in unique(hft_p)) {
		r <- which(hft_p == h)
		cc <- which(cent_hft == h)
		nv <- rowSums(!is.na(zp[r, , drop=FALSE]))
		if (!length(cc)) {
			out[r] <- ifelse(nv > 0, 100L * as.integer(h) + 1L, NA_integer_)
			next
		}
		d <- vapply(cc, function(j) rowSums(sweep(zp[r, , drop=FALSE], 2, cent[j, ])^2, na.rm=TRUE) / nv, numeric(length(r)))
		d <- matrix(d, nrow=length(r))
		pick <- cent_code[cc[max.col(-d, ties.method='first')]]
		pick[nv == 0] <- NA_integer_
		out[r] <- pick
	}
	out
}

# Subclasses of the pixels of the HFTs of a cut. pix: values of the pixels in the HFTs (rows; NA without data), cls:
# their HFT. list(code = subclass of each row (NA without data), threshold, r2_hfts and r2_subclasses: part of the
# variance of the pixels with data in every variable explained by the HFTs and by the subclasses, mean of the
# variables).
hft_subclasses <- function(pix, cls) {
	mu <- colMeans(pix, na.rm=TRUE)
	sdv <- apply(pix, 2, sd, na.rm=TRUE)
	sdv[!is.finite(sdv) | sdv == 0] <- 1
	z <- sweep(sweep(pix, 2, mu), 2, sdv, '/')
	full <- stats::complete.cases(z)
	thr <- subclass_threshold(z[full, , drop=FALSE], cls[full])
	code <- rep(NA_integer_, nrow(z))
	code[full] <- if (is.finite(thr)) bisect_hfts(z[full, , drop=FALSE], cls[full], thr) else 100L * as.integer(cls[full]) + 1L
	part <- which(!full)
	if (length(part)) {
		g <- factor(code[full])
		cent <- rowsum(z[full, , drop=FALSE], g) / as.vector(table(g))
		code[part] <- subclass_partial(z[part, , drop=FALSE], cls[part], cent, as.integer(levels(g)))
		# the subclasses of each HFT numbered again from the largest, with the pixels with partial data
		ok <- !is.na(code)
		tab <- table(code[ok])
		cds <- as.integer(names(tab))
		h <- cds %/% 100L
		new <- integer(length(cds))
		for (hh in unique(h)) {
			i <- which(h == hh)
			i <- i[order(-as.vector(tab)[i], cds[i])]
			new[i] <- 100L * hh + seq_along(i)
		}
		code[ok] <- new[match(code[ok], cds)]
	}
	# part of the variance explained, mean of the variables that vary (a constant one has none to explain)
	r2 <- function(g) {
		zf <- z[full, , drop=FALSE]
		tot <- colSums(sweep(zf, 2, colMeans(zf))^2)
		within <- colSums((zf - (rowsum(zf, g) / as.vector(table(g)))[as.character(g), , drop=FALSE])^2)
		mean((1 - within / tot)[tot > 0])
	}
	list(code=code, threshold=thr, r2_hfts=r2(factor(cls[full])), r2_subclasses=r2(factor(code[full])))
}

# Colour of each subclass: a shade of the colour of its HFT by area, from the largest (the darkest) to the smallest
# (the lightest); an HFT with a single subclass keeps its colour.
subclass_colours <- function(codes, size, hft_cols) {
	out <- character(length(codes))
	for (h in unique(codes %/% 100L)) {
		i <- which(codes %/% 100L == h)
		i <- i[order(-size[i])]
		t <- if (length(i) > 1) seq(-0.42, 0.42, length.out=length(i)) else 0
		base <- grDevices::col2rgb(hft_cols[[as.character(h)]])[, 1] / 255
		for (j in seq_along(i)) {
			v <- if (t[j] < 0) base * (1 + t[j]) else base + (1 - base) * t[j]
			out[i[j]] <- grDevices::rgb(v[1], v[2], v[3])
		}
	}
	out
}

# The shapefile of the HFTs with its segments split into the subclasses of their pixels: one polygon per segment and
# subclass (with several parts when its pixels do not touch), with the fields of its segment (layer: class_layer of
# moranvar_plots.R), its subclass (hft_<n>_sub_<y>) and code (100 x HFT + subclass), and its area (that of the
# polygon, as that of the segments). The pixels of a segment without a subclass (no data in any variable) are one
# more polygon, without subclass. pixgrid: the grid of the pixels; scode: the subclass code of each of its pixels (NA
# none)
subclass_pieces <- function(layer, pixgrid, scode) {
	layer$.row <- seq_len(nrow(layer))
	row <- terra::values(terra::rasterize(terra::vect(layer), pixgrid, field='.row'), mat=FALSE)
	# the HFT of the subclass of every pixel is that of its segment (both from the same polygons, drawn alike)
	both <- !is.na(row) & !is.na(scode)
	if (any(scode[both] %/% 100L != layer$hft[row[both]])) stop('subclass_pieces: the subclasses of some pixels are not those of the HFT of their segment')
	key <- row * 100 + ifelse(is.na(scode), 0L, scode %% 100L)
	kr <- terra::rast(pixgrid, nlyrs=1)
	terra::values(kr) <- key
	names(kr) <- 'key'
	pol <- sf::st_as_sf(terra::as.polygons(kr, dissolve=TRUE))
	pol <- pol[order(pol$key), ]
	k <- pol$key
	sb <- as.integer(k %% 100)
	a <- sf::st_drop_geometry(layer)[k %/% 100, setdiff(names(layer), c('.row', attr(layer, 'sf_column'))), drop=FALSE]
	ha <- sf::st_area(sf::st_geometry(pol))
	a$area_ha <- if (inherits(ha, 'units')) as.numeric(units::set_units(ha, 'ha', mode='standard')) else as.numeric(ha) / 1e4
	a$subclass <- ifelse(sb > 0, sprintf('hft_%d_sub_%d', a$hft, sb), NA_character_)
	a$sub_code <- ifelse(sb > 0, 100L * a$hft + sb, NA_integer_)
	first <- c('segm_id', 'segm_cat', 'area_id', 'hft', 'subclass', 'sub_code', 'area_ha', 'silhouette')
	a <- a[, c(first, setdiff(names(a), first))]
	rownames(a) <- NULL
	sf::st_sf(a, geometry=sf::st_geometry(pol))
}
