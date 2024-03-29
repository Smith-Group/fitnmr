#' Make a multiplet matrix with weights and scalar coupling coefficients
#'
#' @param sc_names character vector with names of scalar couplings
make_coupling_mat <- function(sc_names) {

	grid_list <- rep(list(c(0.5,-0.5)), length(sc_names))
	names(grid_list) <- sc_names
	
	offset_grid <- as.matrix(expand.grid(grid_list))
	
	offsets <- sapply(unique(sc_names), function(sc_name) rowSums(offset_grid[,colnames(offset_grid) == sc_name,drop=FALSE]))
	
	if (length(offsets) == 0) {
		offsets <- matrix(nrow=1, ncol=0)
	}
	
	offset_char <- apply(offsets, 1, paste, collapse=" ")
	
	offset_weights <- tapply(offset_char, offset_char, length)/length(offset_char)
	
	cbind(unname(offset_weights), 0, offsets[match(names(offset_weights), offset_char),,drop=FALSE])
}

#' Split string of scalar coupling names
#'
#' The only currently implemented way of splitting string is using whitespace
#'
#' @param name_char single string (character vector of length 1) with scalar coupling names
split_coupling_names <- function(name_char) {

	strsplit(name_char, "\\s+")[[1]]
}

#' Convert data frame of resonances into a parameter list
#'
#' @param spec single spectrum
#' @param resonance data frame with resonances
#' @param nuclei data frame with nuclei chemical shifts and R2 rates
#' @param data frame with scalar couplings
resonance_to_param_list <- function(spec, resonance, nuclei, couplings) {

	num_spec <- 1
	num_dim <- ncol(spec$fheader)
	
	dim_names <- c("x", "y", "z", "a")[seq_len(num_dim)]
	
	comb_list_coupling <- lapply(seq_along(dim_names), function(dim_idx) {
	
		# get scalar couplings for this dimension
		res_colname <- dim_names[dim_idx]
		sc_colname <- paste(res_colname, "_sc", sep="")
		sc <- ""
		if (sc_colname %in% names(resonance)) {
			sc <- resonance[[sc_colname]]
		}
		sc_names <- split_coupling_names(sc)
		
		make_coupling_mat(sc_names)
	})
	
	num_peaks <- 1
	
	sc_names <- unique(unlist(lapply(comb_list_coupling, function(x) colnames(x)[-(1:2)]), use.names=FALSE))
	nuclei_names <- unlist(resonance[dim_names], use.names=FALSE)
	
	# prepare start_list components
	omega0_start <- nuclei[nuclei_names,"omega0_ppm"]
	r2_start <- nuclei[nuclei_names,"r2_hz"]
	m0_start <- resonance[["m0"]]/length(nuclei_names)
	omega0_comb_start <- couplings[sc_names,"hz"]
	names(omega0_comb_start) <- sc_names
	
	# prepare group_list components
	omega0_group <- NA_integer_
	if ("omega0_group" %in% colnames(nuclei)) {
		omega0_group <- nuclei[nuclei_names,"omega0_group"]
	}
	r2_group <- NA_integer_
	if ("r2_group" %in% colnames(nuclei)) {
		r2_group <- nuclei[nuclei_names,"r2_group"]
	}
	m0_group <- NA_integer_
	if ("m0_group" %in% names(resonance)) {
		m0_group <- nuclei[["m0_group"]]
	}
	omega0_comb_group <- rep(NA_integer_, length(omega0_comb_start))
	names(omega0_comb_group) <- names(omega0_comb_start)
	if ("group" %in% colnames(couplings)) {
		omega0_comb_group[sc_names] <- couplings[sc_names,"group"]	
	}
	
	param_list <- list(
		start_list=list(
			omega0=array(omega0_start, dim=c(num_dim, num_peaks, num_spec)),
			r2=array(r2_start, dim=c(num_dim, num_peaks, num_spec)),
			m0=array(m0_start, dim=c(num_peaks, num_spec)),
			omega0_comb=omega0_comb_start
		),
		group_list=list(
			omega0=array(omega0_group, dim=c(num_dim, num_peaks, num_spec)),
			r2=array(r2_group, dim=c(num_dim, num_peaks, num_spec)),
			m0=array(m0_group, dim=c(num_peaks, num_spec)),
			omega0_comb=omega0_comb_group
		),
		comb_list=list(
			omega0=array(list(NULL), dim=c(num_dim, num_peaks, num_spec)),
			coupling=array(comb_list_coupling, dim=c(num_dim, num_peaks, num_spec))
		),
		resonance_names=rownames(resonance),
		nucleus_names=array(nuclei_names, dim=c(num_dim, num_peaks))
	)
	
	param_list
}

#' Combine multi-dimensional arrays with lists
#'
#' @param ... any number of lists with dimensions or a single list of list arrrays
abind_list <- function(..., rev.along=NULL) {

	list_list <- list(...)
	if (length(list_list) == 1) list_list <- list_list[[1]]
	
	#print(lapply(list_list, function(x) dim(x)))
	
	idx_array_list <- lapply(list_list, function(x) array(seq_along(x), dim(x)))
	max_idx <- 0
	for (i in seq_along(idx_array_list)) {
		idx_array_list[[i]] <- idx_array_list[[i]] + max_idx
		max_idx <- max(idx_array_list[[i]])
	}
	
	idx_array_bind <- abind::abind(idx_array_list, rev.along=rev.along)
	
	array(do.call(c, list_list)[as.vector(idx_array_bind)], dim=dim(idx_array_bind))
}

#' Combine parameter lists referring to different peaks
#'
#' @param ... any number of parameter lists or a single list of parameter lists
peak_bind <- function(...) {

	param_lists <- list(...)
	if (length(param_lists) == 1 && is.list(param_lists[[1]]) && is.list(param_lists[[1]][[1]]) && "start_list" %in% names(param_lists[[1]][[1]])) {
		param_lists <- param_lists[[1]]
	}

	start_list_omega0 <- lapply(param_lists, function(x) x[["start_list"]][["omega0"]])
	start_list_r2 <- lapply(param_lists, function(x) x[["start_list"]][["r2"]])
	start_list_m0 <- lapply(param_lists, function(x) x[["start_list"]][["m0"]])
	start_list_omega0_comb <- do.call(c, lapply(param_lists, function(x) x[["start_list"]][["omega0_comb"]]))

	group_list_omega0 <- lapply(param_lists, function(x) x[["group_list"]][["omega0"]])
	group_list_r2 <- lapply(param_lists, function(x) x[["group_list"]][["r2"]])
	group_list_m0 <- lapply(param_lists, function(x) x[["group_list"]][["m0"]])
	group_list_omega0_comb <- do.call(c, lapply(param_lists, function(x) x[["group_list"]][["omega0_comb"]]))

	comb_list_omega0 <- lapply(param_lists, function(x) x[["comb_list"]][["omega0"]])
	
	resonance_names <- do.call(c, lapply(param_lists, function(x) x[["resonance_names"]]))

	param_list <- list(
		start_list=list(
			omega0=abind::abind(start_list_omega0, rev.along=2),
			r2=abind::abind(start_list_r2, rev.along=2),
			m0=abind::abind(start_list_m0, rev.along=2),
			omega0_comb=start_list_omega0_comb[!duplicated(names(start_list_omega0_comb))]
		),
		group_list=list(
			omega0=abind::abind(group_list_omega0, rev.along=2),
			r2=abind::abind(group_list_r2, rev.along=2),
			m0=abind::abind(group_list_m0, rev.along=2),
			omega0_comb=group_list_omega0_comb[!duplicated(names(group_list_omega0_comb))]
		),
		comb_list=list(
			omega0=abind_list(comb_list_omega0, rev.along=2)
		),
		resonance_names=resonance_names
	)
	
	if ("coupling" %in% names(param_lists[[1]][["comb_list"]])) {
		comb_list_coupling <- lapply(param_lists, function(x) x[["comb_list"]][["coupling"]])
		param_list[["comb_list"]][["coupling"]] <- abind_list(comb_list_coupling, rev.along=2)
	}
	
	if ("nucleus_names" %in% names(param_lists[[1]])) {
		nucleus_names <- lapply(param_lists, function(x) x[["nucleus_names"]])
		param_list[["nucleus_names"]] <- abind::abind(nucleus_names, rev.along=1)
	}
	
	param_list
}

#' Combine parameter lists referring to different spectra
#'
#' @param ... any number of parameter lists or a single list of parameter lists
spec_bind <- function(...) {

	param_lists <- list(...)
	if (length(param_lists) == 1 && is.list(param_lists[[1]]) && is.list(param_lists[[1]][[1]]) && "start_list" %in% names(param_lists[[1]][[1]])) {
		param_lists <- param_lists[[1]]
	}

	start_list_omega0 <- lapply(param_lists, function(x) x[["start_list"]][["omega0"]])
	start_list_r2 <- lapply(param_lists, function(x) x[["start_list"]][["r2"]])
	start_list_m0 <- lapply(param_lists, function(x) x[["start_list"]][["m0"]])
	start_list_omega0_comb <- do.call(c, lapply(param_lists, function(x) x[["start_list"]][["omega0_comb"]]))

	group_list_omega0 <- lapply(param_lists, function(x) x[["group_list"]][["omega0"]])
	group_list_r2 <- lapply(param_lists, function(x) x[["group_list"]][["r2"]])
	group_list_m0 <- lapply(param_lists, function(x) x[["group_list"]][["m0"]])
	group_list_omega0_comb <- do.call(c, lapply(param_lists, function(x) x[["group_list"]][["omega0_comb"]]))

	comb_list_omega0 <- lapply(param_lists, function(x) x[["comb_list"]][["omega0"]])
	
	resonance_names <- param_lists[[1]][["resonance_names"]]

	param_list <- list(
		start_list=list(
			omega0=abind::abind(start_list_omega0, rev.along=1),
			r2=abind::abind(start_list_r2, rev.along=1),
			m0=abind::abind(start_list_m0, rev.along=1),
			omega0_comb=start_list_omega0_comb[!duplicated(names(start_list_omega0_comb))]
		),
		group_list=list(
			omega0=abind::abind(group_list_omega0, rev.along=1),
			r2=abind::abind(group_list_r2, rev.along=1),
			m0=abind::abind(group_list_m0, rev.along=1),
			omega0_comb=group_list_omega0_comb[!duplicated(names(group_list_omega0_comb))]
		),
		comb_list=list(
			omega0=abind_list(comb_list_omega0, rev.along=1)
		),
		resonance_names=resonance_names
	)
	
	if ("coupling" %in% names(param_lists[[1]][["comb_list"]])) {
		comb_list_coupling <- lapply(param_lists, function(x) x[["comb_list"]][["coupling"]])
		param_list[["comb_list"]][["coupling"]] <- abind_list(comb_list_coupling, rev.along=1)
	}
	
	if ("nucleus_names" %in% names(param_lists[[1]])) {
		nucleus_names <- param_lists[[1]][["nucleus_names"]]
		param_list[["nucleus_names"]] <- nucleus_names
	}
	
	param_list
}

#' Convert tables with resonance/nuclei/couplings to a parameter list
#'
#' @param spec_list list of spectra
#' @param tables list with resonance/nuclei/couplings tables
#'
#' @export
tables_to_param_list <- function(spec_list, tables) {

	resonances <- tables[["resonances"]]
	nuclei <- tables[["nuclei"]]
	couplings <- tables[["couplings"]]

	resonance_param_list <- lapply(seq_along(spec_list), function(spec_i) {
	
		resonance_param_list <- lapply(seq_len(nrow(resonances)), function(resonance_i) {
		
			resonance_to_param_list(spec_list[[spec_i]], resonances[resonance_i,,drop=FALSE], nuclei, couplings)
		})
	
		param_list <- peak_bind(resonance_param_list)
		
		param_list
	})
	
	param_list <- spec_bind(resonance_param_list)
	
	# find unique nucleus names
	nucleus_names <- param_list[["nucleus_names"]]
	nucleus_unique_idx <- which(!duplicated(as.vector(nucleus_names)))
	
	# create groups for nuclei with NA omega0 values
	nucleus_omega0_group <- param_list[["group_list"]][["omega0"]][nucleus_unique_idx]
	names(nucleus_omega0_group) <- nucleus_names[nucleus_unique_idx]
	na_idx <- is.na(nucleus_omega0_group)
	nucleus_omega0_group[na_idx] <- utils::head(setdiff(seq_along(nucleus_omega0_group), nucleus_omega0_group[!na_idx]), sum(na_idx))
	param_list[["group_list"]][["omega0"]][] <- nucleus_omega0_group[nucleus_names]
	
	# create groups for nuclei with NA r2 values
	nucleus_r2_group <- param_list[["group_list"]][["r2"]][nucleus_unique_idx]
	names(nucleus_r2_group) <- nucleus_names[nucleus_unique_idx]
	na_idx <- is.na(nucleus_r2_group)
	nucleus_r2_group[na_idx] <- utils::head(setdiff(seq_along(nucleus_r2_group), nucleus_r2_group[!na_idx]), sum(na_idx))
	param_list[["group_list"]][["r2"]][] <- nucleus_r2_group[nucleus_names]
	
	spec_names <- paste(seq_along(spec_list), "_m0", sep="")
	resonance_idx <- match(param_list[["resonance_names"]], rownames(resonances))
	for (i in seq_along(spec_names)) {
		if (spec_names[i] %in% colnames(resonances)) {
			param_list[["start_list"]][["m0"]][,i] <- resonances[resonance_idx,spec_names[i]]
		}
	}
	
	param_list
}

#' Convert a parameter list into a set of tables with resonance/nuclei/couplings
#'
#' @param param_list param_list, fit_input, or fit_output structure
#' @param tables list with resonance/nuclei/couplings tables to update
#'
#' @export
param_list_to_tables <- function(param_list, tables) {

	if ("fit_list" %in% names(param_list)) {
		params <- param_list[["fit_list"]]
	} else {
		params <- param_list[["start_list"]]
	}
	
	resonances <- tables[["resonances"]]
	nuclei <- tables[["nuclei"]]
	couplings <- tables[["couplings"]]
	
	resonances_standard_columns <- c("x", "x_sc", "y", "y_sc", "z", "z_sc", "a", "a_sc")
	
	resonances <- resonances[,intersect(colnames(resonances), resonances_standard_columns),drop=FALSE]

	m0_mat <- sapply(seq_len(ncol(params$m0)), function(i) {
		tapply(params$m0[,i,drop=FALSE], param_list$resonance_names, sum)
	})
	if (!is.matrix(m0_mat)) {
		m0_mat <- t(t(m0_mat))
	}
	m0_mat <- m0_mat[rownames(resonances),,drop=FALSE]
	colnames(m0_mat) <- paste(seq_len(ncol(m0_mat)), "_m0", sep="")
	
	resonances <- cbind(resonances, m0_mat)
	
	nuclei_names <- intersect(rownames(nuclei), as.vector(param_list$nucleus_names))
	nuclei_idx <- match(nuclei_names, as.vector(param_list$nucleus_names))
	nuclei[nuclei_names,"omega0_ppm"] <- params$omega0[nuclei_idx]
	nuclei[nuclei_names,"r2_hz"] <- params$r2[nuclei_idx]
	
	coupling_names <- intersect(rownames(couplings), names(params$omega0_comb))
	couplings[coupling_names,"hz"] <- params$omega0_comb[coupling_names]
	
	list(
		resonances=resonances,
		nuclei=nuclei,
		couplings=couplings
	)
}

#' Collapse strings of repeated NAs in a vector with numeric names
#'
#' @param x vector with names having numeric positional values
collapse_na <- function(x) {

	if (is.null(x)) {
		return(NULL)
	}

	x <- structure(as.vector(x), .Names=dimnames(x)[[1]])

	rle_list <- rle(as.integer(is.na(x)))
	
	#print(rle_list)
	
	starts <- c(1L, cumsum(rle_list$lengths)[-length(rle_list$lengths)]+1)
	ends <- cumsum(rle_list$lengths)
	
	#print(cbind(starts, ends, rle_list$value))
	
	x_new <- vector(class(x), 0)
	
	for (i in seq_along(starts)) {
	
		if (rle_list$value[i] == 0) {
			x_new <- c(x_new, x[starts[i]:ends[i]])
		} else if (i != 1 || i != length(starts)) {
			x_new <- c(x_new, structure(NA_real_, .Names=mean(as.numeric(names(x)[starts[i]:ends[i]]))))
		}
	}
	
	rle_list <- rle(as.integer(is.na(x_new)))
	
	starts <- c(1L, cumsum(rle_list$lengths)[-length(rle_list$lengths)]+1)
	ends <- cumsum(rle_list$lengths)
	
	structure(x_new, starts=starts[rle_list$value == 0], ends=ends[rle_list$value == 0])
}

#' Create a sparse axis
#'
#' @export
make_map <- function(x, max_spacing=0.125, new_spacing=0.025) {

	x <- collapse_na(drop(x))

	x_starts <- y_starts <- as.numeric(names(x)[attr(x, "starts")])
	x_ends <- y_ends <- as.numeric(names(x)[attr(x, "ends")])
	#print(x_starts)
	#print(x_ends)
	
	dec <- 0
	
	delete_start <- logical(length(x_starts))
	delete_end <- logical(length(x_starts))
	
	for (i in seq_along(x_starts)[-1]) {
		if (x_ends[i-1] - x_starts[i] > max_spacing) {
			dec <- dec + new_spacing - (x_ends[i-1] - x_starts[i])
		} else {
			delete_start[i] <- TRUE
			delete_end[i-1] <- TRUE
		}
		y_starts[i] <- y_starts[i]-dec
		y_ends[i] <- y_ends[i]-dec
	}
	
	y_ends <- y_ends-y_starts[1]
	y_starts <- y_starts-y_starts[1]
	
	map <- cbind(as.vector(rbind(x_starts, x_ends)), as.vector(rbind(y_starts, y_ends)))
	
	structure(map[order(map[,1]),], starts=x_starts[!delete_start], ends=x_ends[!delete_end])
}

overlap_groups <- function(interval_mat) {

	overlap_fn <- function(i, j, x1, x2) {

		max1 <- pmax(x1[i], x1[j])
		min2 <- pmin(x2[i], x2[j])
	
		pmax(0, min2-max1)
	}

	idx <- seq_len(nrow(interval_mat))
	overlap_mat <- outer(idx, idx, overlap_fn, x1=interval_mat[,1], x2=interval_mat[,2])
	#diag(overlap_mat) <- 0
	#print(overlap_mat)

	overlap_clust <- stats::hclust(stats::as.dist(-overlap_mat), method="single")

	stats::cutree(overlap_clust, h=-1e-8)
}

remove_overlaps <- function(interval_mat) {

	stopifnot(interval_mat[,2] > interval_mat[,1])

	if (nrow(interval_mat) == 1) {
		return(interval_mat)
	}

	group_idx <- overlap_groups(interval_mat)
	
	while (any(duplicated(group_idx))) {
	
		interval_centers <- rowMeans(interval_mat)
		
		interval_centers <- rowMeans(interval_mat)
		interval_lengths <- interval_mat[,2]-interval_mat[,1]

		cluster_list <- tapply(seq_along(group_idx), group_idx, c)
		cluster_size <- sapply(cluster_list, length)

		for (cluster_idx in cluster_list[cluster_size > 1]) {

			# reorder cluster
# 			print(cluster_idx)
			cluster_idx <- cluster_idx[order(interval_centers[cluster_idx])]
	
			cluster_cumsum <- cumsum(interval_lengths[cluster_idx])
	
			new_interval_mat <- cbind(
				c(0, cluster_cumsum[-length(cluster_cumsum)]),
				cluster_cumsum
			)
	
			cluster_com <- mean(interval_centers[cluster_idx])
	
# 			print(mean(new_interval_mat))
			new_interval_mat <- new_interval_mat - (mean(new_interval_mat) - cluster_com)
	
# 			print(mean(new_interval_mat))
# 			print(cluster_com)
# 			print(interval_mat[cluster_idx,])
# 			print(new_interval_mat)
			interval_mat[cluster_idx,] <- new_interval_mat
			#stop()
		}
	
		group_idx <- overlap_groups(interval_mat)
	}

	interval_mat
}

#' Plot spectrum from 1D fit
#'
#' @param fit_data fit_input or fit_output structure
#' @param tables list with resonances, nuclei, and couplings tables
#' @param spec_idx index of spectrum to plot
#' @param col_model color for modeled peak shapes
#' @param col_resonances colors for resonances
#'
#' @export
plot_sparse_1d <- function(fit_data, tables=NULL, spec_idx=1, col_model=2, col_resonance=NULL, lwd=1, tick_spacing=0.02, coupling_spacing=0.01, coupling_marks=0.009, xaxs="i", yaxt="n", bty="n", always_show_start=FALSE, add=FALSE, ppm_map=make_map(get_spec_int(fit_data, "input", spec_idx)[[1]])) {

	stopifnot(length(spec_idx) == 1)

	input_spec_int <- collapse_na(get_spec_int(fit_data, "input", spec_idx)[[1]])
	start_spec_int <- collapse_na(get_spec_int(fit_data, "start", spec_idx)[[1]])
	if ("fit_list" %in% names(fit_data)) {
		fit_spec_int <- collapse_na(get_spec_int(fit_data, "fit", spec_idx)[[1]])
	} else {
		fit_spec_int <- NULL
	}
	
	# determine which model spectra to be plotting
	dashed_int <- NULL
	if (is.null(fit_spec_int)) {
		solid_int <- start_spec_int
		solid_list <- fit_data[["start_list"]]
	} else {
		solid_int <- fit_spec_int
		solid_list <- fit_data[["fit_list"]]
		if (always_show_start) {
			dashed_int <- start_spec_int
		}
	}
	
	if (!is.null(tables[["resonances"]])) {
	
		solid_omega0_weights_list <- lapply(seq_len(nrow(tables[["resonances"]])), function(idx) {
			coupling_omega0_weights(
				solid_list[["omega0"]][,idx,spec_idx],
				fit_data[["comb_list"]][["coupling"]][[1,idx,spec_idx]], 
				solid_list[["omega0_comb"]], 
				fit_data[["spec_data"]][[spec_idx]][["ref_freq"]]
			)
		})
		names(solid_omega0_weights_list) <- fit_data$resonance_names
		#print(solid_omega0_weights_list)
	}
	
	# create sparse ppm map
	ppm <- as.numeric(names(input_spec_int))
	ppm_map_fn <- stats::approxfun(ppm_map[,1], ppm_map[,2])
	
	# determine initial limits
	xlim <- rev(range(ppm_map[,2]))
	ylim <- range(input_spec_int, solid_int, na.rm=TRUE)
	
	if (!is.null(tables[["couplings"]])) {
	
		resonance_couplings <- strsplit(trimws(tables$resonances[,"x_sc"]), " +")
		names(resonance_couplings) <- rownames(tables$resonances)
		unique_couplings <- unique(unlist(resonance_couplings))
	
	} else {
	
		unique_couplings <- character()
	}
	
	if (length(unique_couplings) && !is.null(tables[["couplings"]])) {
	
		# determine which resonances correspond to each coupling
		coupling_resonances <- lapply(unique_couplings, function(x) names(resonance_couplings)[sapply(resonance_couplings, function(y) any(x %in% y))])
		names(coupling_resonances) <- unique_couplings
		
		# sort the couplings by the size of the ppm range they cover
		coupling_range <- sapply(coupling_resonances, function(x) range(tables$nuclei[tables$resonances[x,"x"],"omega0_ppm"]))
		coupling_order <- order(apply(coupling_range, 2, diff))
		
		# determine y coordinates for each couplign and adjust y-limits
		ywidth <- diff(ylim)/(1-length(unique_couplings)*coupling_spacing-coupling_marks/2)
		coupling_y <- ylim[2]+seq_along(unique_couplings)*coupling_spacing*ywidth
		ylim[2] <- ylim[1]+ywidth
	}
	
	plot(ppm_map_fn(ppm), input_spec_int, type="l", lwd=lwd, xlim=xlim, ylim=ylim, xaxs=xaxs, xlab=NA, ylab=NA, xaxt="n", yaxt=yaxt, bty=bty)
	
	starts <- attr(ppm_map, "starts")
	ends <- attr(ppm_map, "ends")
	
	all_ticks <- numeric()
	
	for (j in seq_along(starts)) {
		
		ticks <- seq(floor(starts[j]/tick_spacing), ceiling(ends[j]/tick_spacing))*tick_spacing
		graphics::axis(1, ppm_map_fn(c(starts[j], ends[j])), labels=FALSE, lwd.ticks=0)
		graphics::axis(1, ppm_map_fn(ticks), ticks, labels=FALSE)
		all_ticks <- c(all_ticks, ticks)
	}
	graphics::axis(1, ppm_map_fn(all_ticks), all_ticks, tick=FALSE)
	
	graphics::points(ppm_map_fn(ppm), solid_int, type="l", lwd=lwd, col=col_model)
	
	if (!is.null(dashed_int)) {
		graphics::points(ppm_map_fn(ppm), dashed_int, type="l", lwd=lwd, col=col_model, lty="dashed")
	}
	
	if (is.null(tables)) {
	
		graphics::title(xlab=paste(names(fit_data$spec_data[[spec_idx]]$omega_contigous), "(ppm)"))
	
	} else {
	
		fit_res_int <- NULL
		start_res_int <- lapply(seq_len(nrow(tables$resonances)), function(i) get_spec_int(fit_data, "start", spec_idx=spec_idx, peak_idx=i)[[1]])
		if ("fit_list" %in% names(fit_data)) {
			fit_res_int <- lapply(seq_len(nrow(tables$resonances)), function(i) get_spec_int(fit_data, "fit", spec_idx=spec_idx, peak_idx=i)[[1]])
		}
		
		dashed_res_int <- NULL
		if (is.null(fit_spec_int)) {
			solid_res_int <- start_res_int
		} else {
			solid_res_int <- fit_res_int
			if (always_show_start) {
				dashed_res_int <- start_res_int
			}
		}
		
		res_label <- rownames(tables$resonances)
		res_label_x <- ppm_map_fn(tables$nuclei[tables$resonances[,"x"],"omega0_ppm"])
		
		if (is.null(col_resonance)) {
			col_resonance <- rep(grDevices::palette()[-c(1,2)], length.out=length(res_label))[rank(-res_label_x, ties.method="first")]
		} else if (!is.null(names(col_resonance))) {
			col_resonance <- col_resonance[res_label]
		}
		names(col_resonance) <- rownames(tables$resonances)
		
		res_label_widths <- abs(graphics::strwidth(res_label))+0.001
		res_label_x <- remove_overlaps(cbind(res_label_x-res_label_widths/2, res_label_x+res_label_widths/2))[,1] + res_label_widths/2
		#res_label_col <- cols[res_label]
		
		for (i in seq_len(nrow(tables$resonances))) {
		
			res_int <- solid_res_int[[i]]
			res_int[abs(res_int) < max(abs(res_int), na.rm=TRUE)*5e-3] <- NA
			
			graphics::points(ppm_map_fn(as.numeric(names(res_int))), res_int, type="l", lwd=lwd, col=col_resonance[i])
			
			if (!is.null(dashed_res_int)) {
				res_int_dashed <- dashed_res_int[[i]]
				res_int_dashed[abs(res_int_dashed) < max(abs(res_int_dashed), na.rm=TRUE)*5e-3] <- NA
				
				graphics::points(ppm_map_fn(as.numeric(names(res_int_dashed))), res_int_dashed, type="l", lwd=lwd, col=col_resonance[i], lty="dashed")
			}
			
			graphics::axis(1, res_label_x[i], res_label[i], tick=FALSE, mgp=graphics::par("mgp")[c(1L,1L,3L)], col.axis=col_resonance[i])
		}
	}
	
	if (length(unique_couplings) && !is.null(tables[["couplings"]])) {
	
		coupling_resonance_count <- sapply(coupling_resonances, length)
	
		for (i in seq_along(coupling_order)) {
		
			coupling_idx <- coupling_order[i]
			coupling_name <- unique_couplings[coupling_idx]
			#print(coupling_name)
			y <- coupling_y[i]
			#print(coupling_range[,coupling_idx])
			
			if (coupling_resonance_count[coupling_idx] == 2) {
				
				# get list of resonances ordered by the chemical shift
				res_names <- coupling_resonances[[coupling_idx]]
				res_names <- res_names[order(tables$nuclei[tables$resonances[res_names,"x"],"omega0_ppm"])]
				
				# use the start and end of the multiplet patterns
				x1 <- min(solid_omega0_weights_list[[res_names[1]]][,1])
				x2 <- max(solid_omega0_weights_list[[res_names[2]]][,1])
				x <- ppm_map_fn(c(x1, x2))
				
				# draw line with two segments colored by the opposite resonance
				xm <- mean(x)
				graphics::segments(x[1], y, xm, y, col=col_resonance[res_names[2]], lwd=lwd)
				graphics::segments(xm, y, x[2], y, col=col_resonance[res_names[1]], lwd=lwd)
				
				# draw vertical segments indicating the multiplet pattern created by the coupling
				count1 <- sum(resonance_couplings[[res_names[1]]] == coupling_name)
				count2 <- sum(resonance_couplings[[res_names[2]]] == coupling_name)
				x1multi <- x[1] + seq(0, by=1, length.out=count1+1)*abs(tables[["couplings"]][coupling_name,"hz"])/fit_data$spec_data[[spec_idx]]$ref_freq
				x2multi <- x[2] - seq(0, by=1, length.out=count2+1)*abs(tables[["couplings"]][coupling_name,"hz"])/fit_data$spec_data[[spec_idx]]$ref_freq
				ym <- y-coupling_marks*ywidth/2
				yp <- y+coupling_marks*ywidth/2
				graphics::segments(x1multi, rep(ym, length(x1multi)), x1multi, rep(yp, length(x1multi)), col=col_resonance[res_names[1]], lwd=2*lwd)
				graphics::segments(x2multi, rep(ym, length(x2multi)), x2multi, rep(yp, length(x2multi)), col=col_resonance[res_names[2]], lwd=2*lwd)
			
			} else {
			
				x <- ppm_map_fn(coupling_range[,coupling_idx])
				graphics::segments(x[1], y, x[2], y)
			}
		}
	}
}

#' Plot resonances from 1D fit
#'
#' @param fit_data fit_input or fit_output structure
#' @param omega0_plus length 3 vector giving ppm range for each dimension
#' @param always_show_start show start parameters even if fit already done
#'
#' @export
plot_resonances_1d <- function(fit_data, always_show_start=FALSE, omega0_plus=0.05) {

	original_int <- unlist(lapply(fit_data$spec_data, function(spec_data) spec_data$spec_int))
	
	start_int <- if (!"fit_list" %in% names(fit_data) || always_show_start) {
		start_par <- pack_fit_params(fit_data$start_list, fit_data$group_list)
		#print(start_par)
		#print(group_param_idx(names(start_par), fit_data$group_list, fit_data$start_list))
		fit_fn(start_par, fit_data, return_resid=FALSE)
	} else {
		NULL
	}
	
	fit_int <- if ("fit_list" %in% names(fit_data)) {
		fit_par <- pack_fit_params(fit_data$fit_list, fit_data$group_list)
		#print(fit_par)
		#print(group_param_idx(names(fit_par), fit_data$group_list, fit_data$start_list))
		fit_fn(fit_par, fit_data, return_resid=FALSE)
	} else {
		NULL
	}
	
	for (spec_i in seq_along(fit_data$spec_data)) {
	
		spec_data <- fit_data$spec_data[[spec_i]]
	
		omega_ppm <- spec_data$omega_eval[[1]]
		
		omega_ppm_seg_ends <- which(abs(diff(omega_ppm)) > abs(stats::median(diff(omega_ppm)))*2)
		omega_ppm_seg_starts <- c(1, omega_ppm_seg_ends+1)
		omega_ppm_seg_ends <- c(omega_ppm_seg_ends, length(omega_ppm))
		plot_idx <- unlist(lapply(seq_along(omega_ppm_seg_starts), function(i) c(NA, seq(omega_ppm_seg_starts[i], omega_ppm_seg_ends[i]))))[-1]
		
		spec_eval_idx <- seq_along(spec_data$spec_int)+spec_data$spec_offset
		
		spec_original_int <- original_int[spec_eval_idx]
		spec_start_int <- start_int[spec_eval_idx]
		spec_fit_int <- fit_int[spec_eval_idx]
		
		for (resonance in unique(fit_data$resonance_names)) {
		
			resonance_idx <- resonance == fit_data$resonance_names
		
			start_omega0_weights_list <- lapply(which(resonance_idx), function(idx) {
				coupling_omega0_weights(
					fit_data[["start_list"]][["omega0"]][,idx,spec_i],
					fit_data[["comb_list"]][["coupling"]][[1,idx,spec_i]], 
					fit_data[["start_list"]][["omega0_comb"]], 
					fit_data[["spec_data"]][[spec_i]][["ref_freq"]]
				)
			})
			#print(start_omega0_weights_list)
			start_omega0_weights <- do.call(rbind, start_omega0_weights_list)
			start_omega0_weights[,2] <- start_omega0_weights[,2]/sum(start_omega0_weights[,2])
			xlim <- rev(range(start_omega0_weights[,1]))+c(1,-1)*omega0_plus
			spec_idx <- omega_ppm < xlim[1] & omega_ppm > xlim[2]
			ylim <- range(0, spec_original_int[spec_idx], spec_start_int[spec_idx], spec_fit_int[spec_idx])
			
			graphics::plot(omega_ppm[plot_idx], spec_original_int[plot_idx], type="l", xlim=xlim, ylim=ylim, xlab=expression(delta (ppm)), ylab="Intensity", main=resonance)
			
			if (!is.null(spec_start_int)) {
				graphics::points(omega_ppm[plot_idx], spec_start_int[plot_idx], type="l", col="blue")
			}
		
			if (!is.null(spec_fit_int)) {
				graphics::points(omega_ppm[plot_idx], spec_fit_int[plot_idx], type="l", col="red")
			}
			
			#graphics::abline(v=fit_data$fit_list$omega0[,resonance_idx,spec_i], col=grDevices::rgb(0, 0, 1, 1/sqrt(sum(resonance_idx))))
			if (!is.null(spec_fit_int)) {
				fit_omega0_weights_list <- lapply(which(resonance_idx), function(idx) {
					coupling_omega0_weights(
						fit_data[["fit_list"]][["omega0"]][,idx,spec_i],
						fit_data[["comb_list"]][["coupling"]][[1,idx,spec_i]], 
						fit_data[["fit_list"]][["omega0_comb"]], 
						fit_data[["spec_data"]][[spec_i]][["ref_freq"]]
					)
				})
				fit_omega0_weights <- do.call(rbind, fit_omega0_weights_list)
				fit_omega0_weights[,2] <- fit_omega0_weights[,2]/sum(fit_omega0_weights[,2])
				graphics::abline(v=fit_omega0_weights[,1], col=grDevices::rgb(0, 0, 1, sqrt(start_omega0_weights[,2])))
			}
			
			if (!is.null(spec_start_int)) {
				graphics::abline(v=start_omega0_weights[,1], col=grDevices::rgb(0, 0, 1, sqrt(start_omega0_weights[,2])), lty="dashed")
			}
		}
		
# 		ylim <- range(0, spec_original_int, spec_start_int, spec_fit_int)
# 		
# 		graphics::plot(omega_ppm[plot_idx], spec_original_int[plot_idx], type="l", xlim=rev(range(omega_ppm)), ylim=ylim, xlab=expression(delta (ppm)), ylab="Intensity")
# 
# 		if (!is.null(spec_start_int)) {
# 			graphics::points(omega_ppm[plot_idx], spec_start_int[plot_idx], type="l", col="blue")
# 		}
# 		
# 		if (!is.null(spec_fit_int)) {
# 			graphics::points(omega_ppm[plot_idx], spec_fit_int[plot_idx], type="l", col="red")
# 		}
	}
}

#' Plot resonances from 2D fit
#'
#' @param fit_data fit_input or fit_output structure
#' @param omega0_plus length 2 vector giving ppm range for each dimension
#'
#' @export
plot_resonances_2d <- function(fit_data, omega0_plus, resonances=unique(fit_data$resonance_names)) {

	low_frac <- 0.05

	spec_i <- 1

	input_spec_int <- get_spec_int(fit_data, "input", spec_i)
	start_spec_int <- get_spec_int(fit_data, "start", spec_i)
	fit_spec_int <- get_spec_int(fit_data, "fit", spec_i)
	
	for (resonance in resonances) {
	
		peak_idx <- which(resonance==fit_data$resonance_names)
		resonance_spec_int <- get_spec_int(fit_data, "fit", peak_idx=peak_idx)
		
		limits <- t(apply(fit_data$start_list$omega0[,peak_idx,spec_i,drop=FALSE], 1, range))+c(-omega0_plus, omega0_plus)
	
		limits <- t(sapply(1:2, function(i) {
			omega_contigous <- fit_data$spec_data[[spec_i]]$omega_contigous[[i]]
			range(omega_contigous[omega_contigous >= limits[i,1] & omega_contigous <= limits[i,2]])
		}))
		
		idx <- 1:2
		
		zlim <- range(input_spec_int[[1]], start_spec_int[[1]], na.rm=TRUE)
		
		plot(1, 1, type="n", xlim=rev(limits[idx[1],]), ylim=rev(limits[idx[2],]), xaxs="i", yaxs="i", xlab=paste(names(dimnames(input_spec_int[[1]]))[1], "(ppm)"), ylab=paste(names(dimnames(input_spec_int[[1]]))[2], "(ppm)"))
		contour_pipe(input_spec_int[[1]], zlim=zlim, low_frac=low_frac, col_pos="black", col_neg="gray", add=TRUE)	
		contour_pipe(resonance_spec_int[[1]], zlim=zlim, low_frac=low_frac, col_pos="purple", col_neg="mediumpurple1", add=TRUE)
		contour_pipe(fit_spec_int[[1]], zlim=zlim, low_frac=low_frac, col_pos="red", col_neg="lightpink", add=TRUE)
		#contour_pipe(start_spec_int[[1]], zlim=zlim, low_frac=low_frac, col_pos="blue", col_neg="lightblue", add=TRUE)
		#contour_pipe(start_spec_int[[1]], low_frac=low_frac, col_pos="blue", col_neg="lightblue", add=TRUE)
		
		omega0_weights_1 <- coupling_omega0_weights(
			fit_data$fit_list$omega0[idx[1],peak_idx[1],spec_i],
			fit_data$comb_list$coupling[[idx[1],peak_idx[1],spec_i]],
			fit_data$fit_list$omega0_comb,
			fit_data$spec_data[[spec_i]]$ref_freq[idx[1]]
		)
		
		omega0_weights_2 <- coupling_omega0_weights(
			fit_data$fit_list$omega0[idx[2],peak_idx[1],spec_i],
			fit_data$comb_list$coupling[[idx[2],peak_idx[1],spec_i]],
			fit_data$fit_list$omega0_comb,
			fit_data$spec_data[[spec_i]]$ref_freq[idx[2]]
		)
	
		omega0_weights_idx <- expand.grid(seq_len(nrow(omega0_weights_1)), seq_len(nrow(omega0_weights_2)))
		omega0_weights <- cbind(
			omega0_weights_1[omega0_weights_idx[,1],1],
			omega0_weights_2[omega0_weights_idx[,2],1],
			omega0_weights_1[omega0_weights_idx[,1],2]*omega0_weights_2[omega0_weights_idx[,2],2]
		)
		graphics::points(omega0_weights[,1], omega0_weights[,2], pch=16, col=grDevices::rgb(0,0,1,sqrt(omega0_weights[,3])))
	}
}


#' Plot resonances from 3D fit
#'
#' @param fit_data fit_input or fit_output structure
#' @param omega0_plus length 3 vector giving ppm range for each dimension
#'
#' @export
plot_resonances_3d <- function(fit_data, omega0_plus, resonances=unique(fit_data$resonance_names)) {

	low_frac <- 0.05

	spec_i <- 1

	input_spec_int <- get_spec_int(fit_data, "input", spec_i)
	start_spec_int <- get_spec_int(fit_data, "start", spec_i)
	fit_spec_int <- get_spec_int(fit_data, "fit", spec_i)
	
	for (resonance in resonances) {
	
		peak_idx <- which(resonance==fit_data$resonance_names)
		resonance_spec_int <- get_spec_int(fit_data, "fit", peak_idx=peak_idx)
	
		limits <- t(apply(fit_data$start_list$omega0[,peak_idx,spec_i,drop=FALSE], 1, range))+c(-omega0_plus, omega0_plus)
		
		limits <- t(sapply(1:3, function(i) {
			omega_contigous <- fit_data$spec_data[[spec_i]]$omega_contigous[[i]]
			range(omega_contigous[omega_contigous >= limits[i,1] & omega_contigous <= limits[i,2]])
		}))
	
		graphics::par(mfcol=c(2, 3))
	
		wts <- lapply(seq_along(dim(fit_spec_int[[1]])), function(i) {
			wts <- apply(resonance_spec_int[[1]], i, mean, na.rm=TRUE)
			wts[!is.finite(wts)] <- 0
			wts
		})
	
		for (idx in list(c(1,2), c(2,3), c(3,1))) {
	
			# create 2D projection using weights from the modeled line shape in the projected dimension
			w_vec <- numeric(dim(input_spec_int[[1]])[-idx])
			names(w_vec) <- dimnames(input_spec_int[[1]])[-idx][[1]]
			w_vec[names(wts[-idx][[1]])] <- wts[-idx][[1]]
			w_vec <- w_vec/sum(w_vec, na.rm=TRUE)
		
			input_spec_2d <- apply(input_spec_int[[1]], idx, function(x) sum(x*w_vec, na.rm=TRUE))
			start_spec_2d <- apply(start_spec_int[[1]], idx, function(x) sum(x*w_vec, na.rm=TRUE))
			fit_spec_2d <- apply(fit_spec_int[[1]], idx, function(x) sum(x*w_vec, na.rm=TRUE))
			resonance_spec_2d <- apply(resonance_spec_int[[1]], idx, function(x) sum(x*w_vec, na.rm=TRUE))
		
			# create 1D projection using weights from the modeled line shape in the projected dimension
			w_vec <- numeric(dim(input_spec_int[[1]])[idx[2]])
			names(w_vec) <- dimnames(input_spec_int[[1]])[idx[2]][[1]]
			w_vec[names(wts[idx[2]][[1]])] <- wts[idx[2]][[1]]
			w_vec <- w_vec/sum(w_vec, na.rm=TRUE)
		
			input_spec_1d <- colSums(t(input_spec_2d)*w_vec, na.rm=TRUE)
			fit_spec_1d <- colSums(t(fit_spec_2d)*w_vec, na.rm=TRUE)
			resonance_spec_1d <- colSums(t(resonance_spec_2d)*w_vec, na.rm=TRUE)
			
			graphics::par(mar=c(2.9, 2.9, 1.5, 1), mgp=c(1.7, 0.6, 0))
			
			ppm <- as.numeric(names(input_spec_1d))
			ylim <- range(input_spec_1d, resonance_spec_1d, fit_spec_1d, na.rm=TRUE)
			ylim <- range(resonance_spec_1d, 0, na.rm=TRUE)
			nucleus_name <- fit_data$nucleus_names[idx[1],peak_idx[1]]
			omega0 <- fit_data$fit_list$omega0[idx[1],peak_idx[1],spec_i]
			plot(ppm, input_spec_1d, type="l", xlim=rev(limits[idx[1],]), ylim=ylim, xaxs="i", main=sprintf("%s = %s ppm", nucleus_name, signif(omega0, 5)), xlab=paste(names(dimnames(input_spec_2d))[1], "(ppm)"), ylab="Weighted Intensity", font.main=1)
			graphics::abline(h=0, col=grDevices::gray(0, 0.1))
			graphics::points(ppm, resonance_spec_1d, type="l", col="purple")
			graphics::points(ppm, fit_spec_1d, type="l", col="red")
			omega0_weights_1 <- coupling_omega0_weights(
				fit_data$fit_list$omega0[idx[1],peak_idx[1],spec_i],
				fit_data$comb_list$coupling[[idx[1],peak_idx[1],spec_i]],
				fit_data$fit_list$omega0_comb,
				fit_data$spec_data[[spec_i]]$ref_freq[idx[1]]
			)
			graphics::abline(v=omega0_weights_1[,1], col=grDevices::rgb(0,0,1,sqrt(omega0_weights_1[,2])))
			
			zlim <- range(input_spec_2d, fit_spec_2d, na.rm=TRUE)
			zlim <- range(resonance_spec_1d, na.rm=TRUE)
		
			if (abs(zlim[1]) < abs(zlim[2])) {
				r2_pos <- "topleft"
				sc_pos <- "topright"
			} else {
				r2_pos <- "bottomleft"
				sc_pos <- "bottomright"
			
			}
			
			r2 <- fit_data$fit_list$r2[idx[1],peak_idx[1],spec_i]
			ref_freq <- fit_data$spec_data[[spec_i]]$ref_freq[idx[1]]
			graphics::legend(r2_pos, legend=parse(text=sprintf("R[2]*\" = %0.2f Hz\"", r2)), bty="n", x.intersp=0, text.col="red")
			graphics::segments(omega0-r2/ref_freq, 0, omega0+r2/ref_freq, 0, col="red")
			sc_names <- colnames(fit_data$comb_list$coupling[[idx[1],peak_idx[1],spec_i]])[-(1:2)]
			if (length(sc_names)) {
				sc_values <- fit_data$fit_list$omega0_comb[sc_names]
				graphics::legend(sc_pos, legend=sprintf("%s = %0.2f Hz", sc_names, sc_values), bty="n", x.intersp=0, text.col="blue")
			}
			
			graphics::par(mar=c(2.9, 2.9, 0.5, 1), mgp=c(1.7, 0.6, 0))
			
			#graphics::abline(v=limits[idx[1],], col="green")
			
			plot(1, 1, type="n", xlim=rev(limits[idx[1],]), ylim=rev(limits[idx[2],]), xaxs="i", yaxs="i", xlab=paste(names(dimnames(input_spec_2d))[1], "(ppm)"), ylab=paste(names(dimnames(input_spec_2d))[2], "(ppm)"))
			contour_pipe(input_spec_2d, zlim=zlim, low_frac=low_frac, col_pos="black", col_neg="gray", add=TRUE)
			graphics::abline(0, 1, col=grDevices::gray(0, 0.1))
			
			#contour_pipe(start_spec_2d, zlim=zlim, low_frac=low_frac, col_pos="blue", col_neg="lightblue", add=TRUE)
			contour_pipe(resonance_spec_2d, zlim=zlim, low_frac=low_frac, col_pos="purple", col_neg="mediumpurple1", add=TRUE)
			contour_pipe(fit_spec_2d, zlim=zlim, low_frac=low_frac, col_pos="red", col_neg="lightpink", add=TRUE)
			#graphics::points(t(fit_data$fit_list$omega0[idx,peak_idx,spec_i]), pch=16, col=grDevices::rgb(0,0,1,1/sqrt(length(peak_idx))))
			
			omega0_weights_2 <- coupling_omega0_weights(
				fit_data$fit_list$omega0[idx[2],peak_idx[1],spec_i],
				fit_data$comb_list$coupling[[idx[2],peak_idx[1],spec_i]],
				fit_data$fit_list$omega0_comb,
				fit_data$spec_data[[spec_i]]$ref_freq[idx[2]]
			)
			
			omega0_weights_idx <- expand.grid(seq_len(nrow(omega0_weights_1)), seq_len(nrow(omega0_weights_2)))
			omega0_weights <- cbind(
				omega0_weights_1[omega0_weights_idx[,1],1],
				omega0_weights_2[omega0_weights_idx[,2],1],
				omega0_weights_1[omega0_weights_idx[,1],2]*omega0_weights_2[omega0_weights_idx[,2],2]
			)
			graphics::points(omega0_weights[,1], omega0_weights[,2], pch=16, col=grDevices::rgb(0,0,1,sqrt(omega0_weights[,3])))
			
			#graphics::abline(v=limits[idx[1],], col="green")
			#graphics::abline(h=limits[idx[2],], col="green")
		}
	}
}
