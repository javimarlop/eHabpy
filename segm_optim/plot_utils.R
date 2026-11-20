# Functions shared by moranvar_plots.R (the HFTs of a run) and hft_compare.R (the HFTs of several runs)

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

# Data of a plot, next to it with the same name (X.png -> X.csv), so that it can be reused or the plot redrawn
# elsewhere (or a table given as X.csv). Same format as the other tables of the results: separated by spaces, a header, NA for missing values,
# no quotes except around the texts that need them (the colours #RRGGBB, which read.table would take as a comment).
# Numbers with decimals rounded to 'digits' significant digits (NA: as they are), those of the columns 'precise'
# (coordinates) to 10; whole numbers (IDs, counts) are written as they are
write_plot_csv <- function(df, png, digits = 6, precise = character(0)) {
	df <- as.data.frame(df, stringsAsFactors = FALSE)
	rownames(df) <- NULL
	for (i in seq_along(df)) {
		x <- df[[i]]
		if (is.factor(x)) x <- as.character(x)
		if (!is.na(digits) && is.numeric(x) && any(x != round(x), na.rm = TRUE)) {
			x <- signif(x, if (names(df)[i] %in% precise) 10 else digits)
		}
		df[[i]] <- x
	}
	q <- which(vapply(df, function(x) is.character(x) && any(grepl('[[:space:]#"\']', x)), logical(1)))
	f <- sub('\\.png$', '.csv', png)
	writeLines(paste(names(df), collapse = ' '), f)
	write.table(df, f, sep = ' ', row.names = FALSE, col.names = FALSE, quote = if (length(q)) q else FALSE, append = TRUE)
	invisible(f)
}

# Short label of the variables in the plots (short_labels) without spaces, for the CSVs of the plots
csv_label <- function(x) gsub('[[:space:]]+', '_', x)
