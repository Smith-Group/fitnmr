#' fitnmr Package Overview
#'
#' The functionality provided by the fitnmr package can be divided into several categories:
#' 
#' @section Spectrum reading, plotting, and analysis:
#' Spectra in the NMRPipe file format can be read with \code{\link{read_nmrpipe}}. Those spectra can be plotted with \code{\link{contour_pipe}}. The noise level within a spectrum (or any numeric vector) can be calculated with \code{\link{noise_estimate}}.
#' 
#' @section Low-level fitting and plotting:
#' The core fitting procedure consists of a two-step process: First, the fit input is created with \code{\link{make_fit_input}}. Second, the fit is executed with \code{\link{perform_fit}}. Prior to running the fit, constraints can be added to the fit with either \code{\link{update_fit_bounds}} or \code{\link{limit_omega0_by_r2}}.
#'
#' Before or after a fit has been performed, you can extract the raw, starting, or fit spectral intensities with \code{\link{get_spec_int}}. Furthermore, there are convenience plotting functions for plotting 1D (\code{\link{plot_fit_1d}}) or 2D (\code{\link{plot_fit_2d}}) fits.
#'
#' The \code{\link{make_fit_input}} function takes many different parameters. To make a new fit from an existing fit, possibly with a different set of spectra or otherwise modified fitting parameters, \code{\link{param_list_to_arg_list}} can be helpful in generating a list of parameters for \code{\link{make_fit_input}}.
#'
#' @section Modifying parameter lists for fits:
#' In fitnmr, a "parameter list" is a named list data structure that has all the information necessary to describe a set of peaks and interdependencies between the parameters for those peaks. It containins three elements: \code{start_list} (or also \code{fit_list} after a fit has been performed), \code{group_list}, and \code{comb_list}. Fitting constraints can also be stored in \code{lower_list} and \code{upper_list}. Each of those is itself a named list of arrays corresponding to different parameters, namely \code{omega0}, \code{r2}, \code{m0}, \code{p0}, \code{p1}, and \code{omega0_comb}. To help manage the values stored in those lists (particularly \code{omega0} and \code{omega0_comb}), there are several convenience functions to select particular subsets, including \code{\link{omega0_param_idx}} and \code{\link{coupling_param_idx}}. Once a subset is made, the parameters can be read or changed using \code{\link{param_values}}.
#'
#' @section High-level fitting and plotting:
#' \code{\link{fit_peak_iter}}
#'
#' \code{\link{param_list_to_peak_df}}, \code{\link{plot_peak_df}}, \code{\link{peak_df_to_param_list}}, \code{\link{peak_df_to_fit_input}}
#'
#' \code{\link{fit_peak_cluster}}, \code{\link{fit_footprint}}
#'
#' \code{\link{fit_peaks}}, \code{\link{make_param_list}}
#'
#' \code{\link{spec_overlap_mat}}
#'
#' @section NMRPipe simulation and processing:
#' \code{\link{ppm_to_pts}}, \code{\link{whz_to_pts}}, \code{\link{write_nmrdraw_peak_tab}}, \code{\link{sim_time_nd}}, \code{\link{nmr_pipe}}
#'
#' @section Postprocessing and assignment:
#' \code{\link{height_assign}}, \code{\link{read_nmrdraw_peak_tab}}
#'
#' @section Deprecated:
#' \code{\link{extract_params}}, \code{\link{get_spec_peak_int}}
#'
#' @docType package
#' @name fitnmr
NULL

#' Read NMRPipe spectrum
#'
#' This function reads 1D-4D spectra stored in the NMRPipe format.
#'
#' For three and four dimensional datasets, the spectral data is often spread across multiple files. To read those, inFormat should be a \code{\link[base]{sprintf}}-style string that describes how the files are named. For instance, if the files are named 001.ft3, 002.ft3, etc., then \code{inFormat} should be \code{"\%03i.ft3"}. If there is no zero-padding, as in this case, 0 should be omitted from the format. If there are fewer digits, then the first 3 should be changed accordingly.
#'
#' The default 2D NMRPipe scripts only have a single transpose ("TP") command, leaving the the indirect dimension as the first dimension in the resulting array. The 2D plotting functions in \code{fitnmr} usually plot this first dimension along the x-axis, which will make for generally non-standard contour plots. Furthermore, when peak fitting is employed, this will also be the first dimension. To fix this, you can change the spectral order with the \code{dim_order} parameter. The order of the dimensions should be specified in the same way that would be done for \code{\link[base]{aperm}}. Alternatively, you can specify a character argument to have \code{fitnmr} attempt to automatically detect and correct the array order. The only currently supported type is \code{"hx"}, which will put the dimension with the greatest observe frequency first.
#'
#' This function is partly based on the pipe2rnmr function from rNMR.
#'
#' @param inFormat character with file name or format for multiple files
#' @param dim_order integer vector used to reorder dimensions or character specifying
#' @param complex_data logical value indicating whether complex data should be read
#' @return a named list with four elements:
#'  \describe{
#'   \item{int}{multidimensional array with spectrum intensities (ppm values are given in the dimnames)}
#'   \item{ppm}{list of numeric vectors giving the ppm values associated with each int array dimension}
#'   \item{fheader}{matrix with dimension-specific header information}
#'   \item{header}{numeric vector with the complete header contents}
#' }
#' 
#' @export
read_nmrpipe <- function(inFormat, dim_order=NULL, complex_data=FALSE) {
	
	inFile <- sprintf(inFormat, 1)
	
	if (missing(inFile))
		stop('The inFile file path is required')	
	
	## Test file connection
	readCon <- file(inFile, 'rb')
	header_raw <- try(readBin(readCon, what='raw', n=512*4, size=1), silent=TRUE)
	if (class(header_raw) == "try-error" || length(header_raw) != 512*4){
		close(readCon)
		stop(paste('Could not read NMRPipe file:\n"', inFile, '"',	sep=''))
	}
	close(readCon)
	header_raw <- matrix(header_raw, nrow=4)
	
	header <- readBin(header_raw, what='numeric', n=512, size=4)
	
	## Check for correct endianness
	if (round(header[3], 3) != 2.345) {
		seek(readCon, where=0)
		header <- readBin(header_raw, what='numeric', n=512, size=4, endian='swap')
		endianness <- 'swap'
	} else {
		endianness <- .Platform$endian
	}
	
	names(header) <- c("FDMAGIC", "FDFLTFORMAT", "FDFLTORDER", "", "", "", "", "", "", "FDDIMCOUNT", "FDF3OBS", "FDF3SW", "FDF3ORIG", "FDF3FTFLAG", "FDPLANELOC", "FDF3SIZE", "FDF2LABEL", "", "FDF1LABEL", "", "FDF3LABEL", "", "FDF4LABEL", "", "FDDIMORDER1", "FDDIMORDER2", "FDDIMORDER3", "FDDIMORDER4", "FDF4OBS", "FDF4SW", "FDF4ORIG", "FDF4FTFLAG", "FDF4SIZE", "", "", "", "", "", "", "", "FDDMXVAL", "FDDMXFLAG", "FDDELTATR", "", "", "FDNUSDIM", "", "", "", "", "FDF3APOD", "FDF3QUADFLAG", "", "FDF4APOD", "FDF4QUADFLAG", "FDF1QUADFLAG", "FDF2QUADFLAG", "FDPIPEFLAG", "FDF3UNITS", "FDF4UNITS", "FDF3P0", "FDF3P1", "FDF4P0", "FDF4P1", "FDF2AQSIGN", "FDPARTITION", "FDF2CAR", "FDF1CAR", "FDF3CAR", "FDF4CAR", "FDUSER1", "FDUSER2", "FDUSER3", "FDUSER4", "FDUSER5", "FDPIPECOUNT", "FDUSER6", "FDFIRSTPLANE", "FDLASTPLANE", "FDF2CENTER", "FDF1CENTER", "FDF3CENTER", "FDF4CENTER", "", "", "", "", "", "", "", "", "", "", "", "", "FDF2APOD", "FDF2FTSIZE", "FDREALSIZE", "FDF1FTSIZE", "FDSIZE", "FDF2SW", "FDF2ORIG", "", "", "", "", "FDQUADFLAG", "", "FDF2ZF", "FDF2P0", "FDF2P1", "FDF2LB", "", "", "", "", "", "", "", "FDF2OBS", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDMCFLAG", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDF2UNITS", "FDNOISE", "", "", "", "FDTEMPERATURE", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDRANK", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDTAU", "FDF3FTSIZE", "FDF4FTSIZE", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDF1OBS", "FDSPECNUM", "FDF2FTFLAG", "FDTRANSPOSED", "FDF1FTFLAG", "", "", "", "", "", "", "FDF1SW", "", "", "", "", "FDF1UNITS", "", "", "", "", "", "", "", "", "FDF1LB", "", "FDF1P0", "FDF1P1", "FDMAX", "FDMIN", "FDF1ORIG", "FDSCALEFLAG", "FDDISPMAX", "FDDISPMIN", "FDPTHRESH", "FDNTHRESH", "", "FD2DPHASE", "FDF2X1", "FDF2XN", "FDF1X1", "FDF1XN", "FDF3X1", "FDF3XN", "FDF4X1", "FDF4XN", "", "FDDOMINFO", "FDMETHINFO", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDHOURS", "FDMINS", "FDSECS", "FDSRCNAME", "", "", "", "FDUSERNAME", "", "", "", "FDMONTH", "FDDAY", "FDYEAR", "FDTITLE", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDCOMMENT", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDLASTBLOCK", "FDCONTBLOCK", "FDBASEBLOCK", "FDPEAKBLOCK", "FDBMAPBLOCK", "FDHISTBLOCK", "FD1DBLOCK", "", "", "", "", "FDSCORE", "FDSCANS", "FDF3LB", "FDF4LB", "FDF2GB", "FDF1GB", "FDF3GB", "FDF4GB", "FDF2OBSMID", "FDF1OBSMID", "FDF3OBSMID", "FDF4OBSMID", "FDF2GOFF", "FDF1GOFF", "FDF3GOFF", "FDF4GOFF", "FDF2TDSIZE", "FDF1TDSIZE", "FDF3TDSIZE", "FDF4TDSIZE", "", "", "", "", "", "", "", "", "", "FD2DVIRGIN", "FDF3APODCODE", "FDF3APODQ1", "FDF3APODQ2", "FDF3APODQ3", "FDF3C1", "FDF4APODCODE", "FDF4APODQ1", "FDF4APODQ2", "FDF4APODQ3", "FDF4C1", "", "", "", "FDF2APODCODE", "FDF1APODCODE", "FDF2APODQ1", "FDF2APODQ2", "FDF2APODQ3", "FDF2C1", "", "FDF1APODQ1", "FDF1APODQ2", "FDF1APODQ3", "FDF1C1", "", "", "", "", "FDF1APOD", "", "", "", "", "", "", "", "", "FDF1ZF", "FDF3ZF", "FDF4ZF", "", "", "FDFILECOUNT", "FDSLICECOUNT0", "FDTHREADCOUNT", "FDTHREADID", "FDSLICECOUNT1", "FDCUBEFLAG", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "FDOPERNAME", "", "", "", "", "", "", "", "", "", "", "FDF1AQSIGN", "FDF3AQSIGN", "FDF4AQSIGN", "", "", "FDF2OFFPPM", "FDF1OFFPPM", "FDF3OFFPPM", "FDF4OFFPPM", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
	
	label_idx <- match(paste("FDF", 1:4, "LABEL", sep=""), names(header))
	
	flabels <- apply(matrix(rawToChar(header_raw[,label_idx], TRUE), nrow=4), 2, paste, collapse="")
	
	## Check that all data is contained within a single file
	if (header[1] != 0){
		stop(paste('Can not convert "', inFile, '",\n  File is not a valid NMRPipe', 
						' format spectrum.', sep=''))
	}
	
	fnames <- substr(names(header)[grep("FDF1", names(header))], 5, 100)
	
	forder <- header[c("FDDIMORDER1", "FDDIMORDER2", "FDDIMORDER3", "FDDIMORDER4")]
	
	fheader_names <- matrix(paste(rep(c("FDF1", "FDF2", "FDF3", "FDF4"), each=length(fnames)), rep(fnames, 4), sep=""), ncol=4)
	
	#fheader_names <- rbind(fheader_names, c("FDSPECNUM", "FDSIZE", "FDF3SIZE", "FDF4SIZE"))
	#fnames <- c(fnames, "SIZE")
	
	fheader <- matrix(header[fheader_names], ncol=4, dimnames=list(fnames, flabels))
	fheader <- fheader[setdiff(rownames(fheader), "LABEL"),]
	
	fsizes <- header[c("FDSIZE", "FDSPECNUM", "FDF3SIZE", "FDF4SIZE")][order(forder)]
	
	#print(header)
	
	fheader <- rbind(fheader, SIZE=fsizes, DMXVAL=0)
	
	fheader["DMXVAL",which.max(utils::head(forder, header["FDDIMCOUNT"]))] <- header["FDDMXVAL"]
	
	FDFILECOUNT <- header["FDFILECOUNT"]
	FDDIMCOUNT <- header["FDDIMCOUNT"]
	
	f_ppm <- apply(fheader, 2, function(x) (x["ORIG"]+x["SW"]*(1-seq_len(x["SIZE"])/x["SIZE"]))/x["OBS"])
	f_ppm <- f_ppm[utils::head(forder, header["FDDIMCOUNT"])]
	f_size <- sapply(f_ppm, length)
	n <- if (header["FDFILECOUNT"] == 1) prod(f_size) else prod(utils::head(f_size, 2))
	
	if (complex_data) {
		data_array <- array(NA_complex_, f_size, f_ppm)
	} else {
		data_array <- array(NA_real_, f_size, f_ppm)
	}
	
	for (i in seq_len(header["FDFILECOUNT"])) {
	
		inFile <- sprintf(inFormat, i)
		readCon <- file(inFile, "rb")
		seek(readCon, where=4*512)
		data_read <- readBin(readCon, size=4, what="numeric", n=n, endian=endianness)
		if (complex_data) {
			data_read_imaginary <- readBin(readCon, size=4, what="numeric", n=n, endian=endianness)
			data_read <- complex(real=data_read, imaginary=data_read_imaginary)
		}
		close(readCon)
		if (length(data_read) < n){
			stop(paste('Can not convert "', inFile, '",\n', 
							'  file size does not match data size.'), sep='')
		}
		data_array[seq(1+(i-1)*n, length.out=n)] <- data_read
	}
	
	
	fheader <- fheader[,utils::head(forder, header["FDDIMCOUNT"]),drop=FALSE]
	
	#print(apply(fheader, 2, function(x) (x["ORIG"]+x["SW"]*(1-x["CENTER"]/x["SIZE"]))/x["OBS"]))
	
	#print(f_ppm)

	if (!is.null(dim_order)) {
		if (is.character(dim_order) && dim_order == "hx") {
			h_idx <- which.max(fheader["OBS",])
			x_idx <- which.min(fheader["OBS",])
			dim_order <- c(h_idx, x_idx)
		}
		data_array <- aperm(data_array, dim_order)
		fheader <- fheader[,dim_order]
		f_ppm <- f_ppm[dim_order]
	}

	list(int=data_array, ppm=f_ppm, fheader=fheader, header=header)
}

lineshapes <- list(
	none=list(	
		func="(1i*(-1 + exp(-(1i*aq*(omega - omega0 - 1i*r2)))))/(omega - omega0 - 1i*r2)",
		domega0="(1i - 1i*exp(aq*(1i*(omega - omega0) + r2)) + aq*(-omega + omega0 + 1i*r2))/(exp(aq*(1i*(omega - omega0) + r2))*(-omega + omega0 + 1i*r2)^2)",
		dr2="(-1 + exp(aq*(1i*(omega - omega0) + r2)) - 1i*aq*(omega - omega0 - 1i*r2))/(exp(aq*(1i*(omega - omega0) + r2))*(-omega + omega0 + 1i*r2)^2)"
	),
	sp1=list(	
		func="(aq*((end - off)*pi*cos(end*pi) - exp(aq*(1i*omega - 1i*omega0 + r2))*(end - off)*pi*cos(off*pi) + aq*(1i*omega - 1i*omega0 + r2)*(sin(end*pi) - exp(aq*(1i*omega - 1i*omega0 + r2))*sin(off*pi))))/(exp(aq*(1i*omega - 1i*omega0 + r2))*((end - off)*pi + aq*(omega - omega0 - 1i*r2))*((-end + off)*pi + aq*(omega - omega0 - 1i*r2)))",
		domega0="(aq^2*(-1i*(end - off)*pi*((end - off)^2*pi^2 - aq^2*(-omega + omega0 + 1i*r2)^2 + 2*aq*(1i*(omega - omega0) + r2))*cos(end*pi) - 2*aq*exp(aq*(1i*(omega - omega0) + r2))*(end - off)*pi*(omega - omega0 - 1i*r2)*cos(off*pi) + (1i*(end - off)^2*pi^2 - aq*(-((end - off)^2*pi^2) + aq*(-1i + aq*(omega - omega0 - 1i*r2))*(omega - omega0 - 1i*r2))*(omega - omega0 - 1i*r2))*sin(end*pi) + exp(aq*(1i*(omega - omega0) + r2))*(1i*(end - off)*pi + aq*(omega - omega0 - 1i*r2))*(-1i*aq*omega + 1i*aq*omega0 - end*pi + off*pi - aq*r2)*sin(off*pi)))/(exp(aq*(1i*(omega - omega0) + r2))*((end - off)*pi + aq*(omega - omega0 - 1i*r2))^2*((end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2)",
		dr2="(aq^2*((end - off)*pi*((end - off)^2*pi^2 - aq^2*(-omega + omega0 + 1i*r2)^2 + 2*aq*(1i*(omega - omega0) + r2))*cos(end*pi) - 2i*aq*exp(aq*(1i*(omega - omega0) + r2))*(end - off)*pi*(omega - omega0 - 1i*r2)*cos(off*pi) + (-((end - off)^2*pi^2) + 1i*aq*(omega - omega0 - 1i*r2)*((end - off)^2*pi^2 - aq^2*(-omega + omega0 + 1i*r2)^2 + aq*(1i*(omega - omega0) + r2)))*sin(end*pi) + exp(aq*(1i*(omega - omega0) + r2))*((end - off)^2*pi^2 + aq^2*(-omega + omega0 + 1i*r2)^2)*sin(off*pi)))/(exp(aq*(1i*(omega - omega0) + r2))*((end - off)*pi + aq*(omega - omega0 - 1i*r2))^2*((end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2)"
	),
	sp2=list(
		func="((-1i*aq*(-1 + exp(-(1i*(2*(end - off)*pi + aq*(omega - omega0 - 1i*r2))))))/(exp(2*1i*off*pi)*(2*(end - off)*pi + aq*(omega - omega0 - 1i*r2))) - (1i*aq*exp(2*1i*off*pi)*(-1 + exp(1i*(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2)))))/(2*(-end + off)*pi + aq*(omega - omega0 - 1i*r2)) + (2i*(-1 + exp(-(1i*aq*(omega - omega0 - 1i*r2)))))/(omega - omega0 - 1i*r2))/4",
		domega0="(exp(1i*aq*(-omega + omega0 + 1i*r2))*(aq*((aq*exp(2*1i*end*pi)*(-1i - 2*end*pi + 2*off*pi + aq*(omega - omega0 - 1i*r2)))/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2 + (aq*(-1i + 2*end*pi - 2*off*pi + aq*(omega - omega0 - 1i*r2)))/(exp(2*1i*end*pi)*(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2) - 2/(omega - omega0 - 1i*r2)) + 2i/(-omega + omega0 + 1i*r2)^2) + (1i*aq^2*exp(2*1i*off*pi))/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2 + (1i*aq^2)/(exp(2*1i*off*pi)*(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2) - 2i/(-omega + omega0 + 1i*r2)^2)/4",
		dr2="((2*(1 - exp(-(1i*aq*(omega - omega0 - 1i*r2)))))/(-omega + omega0 + 1i*r2)^2 + (2*aq)/(exp(aq*(1i*(omega - omega0) + r2))*(1i*(omega - omega0) + r2)) + aq^2*(-(exp(2*1i*off*pi)/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2) - 1/(exp(2*1i*off*pi)*(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2) + (exp(1i*(2*end*pi + aq*(-omega + omega0 + 1i*r2)))*(1 - 2i*(end - off)*pi + aq*(1i*(omega - omega0) + r2)))/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2 + (exp((-2i)*end*pi - aq*(1i*(omega - omega0) + r2))*(1 + 2i*(end - off)*pi + aq*(1i*(omega - omega0) + r2)))/(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2))/4"
	)
)

#lineshapes_parsed <- lapply(lineshapes, function(flist) lapply(flist, function(fvec) parse(text=fvec, keep.source=FALSE)))

#lineshapes_simplified_func <- lapply(lineshapes, function(lineshape_list) compile(findSubexprs(eval(parse(text=paste("quote(", lineshape_list[[1]], ")"))), simplify=TRUE), options=list(optimize=3, suppressUndefined=TRUE)))

#lineshapes_simplified <- lapply(lineshapes, function(lineshape_list) compile(findSubexprs(eval(parse(text=paste("quote(list(", paste(unlist(lineshape_list), collapse=","), "))"))), simplify=TRUE), options=list(optimize=3, suppressUndefined=TRUE)))
#lineshapes_simplified <- lapply(lineshapes, function(lineshape_list) findSubexprs(eval(parse(text=paste("quote(list(", paste(unlist(lineshape_list), collapse=","), "))"))), simplify=TRUE))

fill_array <- function(x, array_dim) {

	if (is.null(x)) x <- NA
	
	array(x, array_dim)
}

fill_array_int <- function(x, array_dim) {

	if (is.null(x)) {
		x <- NA_integer_
	} else {
		x <- as.integer(x)
	}
	
	array(x, array_dim)
}

fill_array_list <- function(x, array_dim, comb=FALSE) {

	if (is.null(x)) {
		x <- vector("list", 1)
	} else {
		x <- as.list(x)
	}
	
	not_null_idx <- !sapply(x, is.null)
	x[not_null_idx] <- lapply(x[not_null_idx], function(y) {
		data.frame(as.integer(y[[1]]), as.numeric(y[[2]]), fix.empty.names=FALSE)
	})
	
	array(x, array_dim)
}

fill_group <- function(x, array_dim, na_dim1_same=FALSE) {

	x <- fill_array_int(x, array_dim)
	
	na_idx <- is.na(x)
	
	if (na_dim1_same) {
		x[na_idx] <- utils::head(setdiff(seq_along(x), x[!na_idx]), array_dim[1])
	} else {
		x[na_idx] <- utils::head(setdiff(seq_along(x), x[!na_idx]), sum(na_idx))
	}
	
	x
}

fill_comb_group <- function(x, comb_array) {

	not_null_idx <- which(!sapply(comb_array, is.null))
	
	if (length(not_null_idx)) {
	
		x_dim <- dim(x)
		x_dimnames <- dimnames(x)
	
		if (is.null(x)) {
			x <- NA_integer_
		} else {
			x <- as.integer(x)
		}
	
		all_comb <- do.call(rbind, comb_array[not_null_idx])
		max_comb_idx <- max(all_comb[,1])
		stopifnot(seq_len(max_comb_idx) %in% all_comb[,1])
	
		x <- rep_len(x, max_comb_idx)
	
		na_idx <- is.na(x)
		x[na_idx] <- utils::head(setdiff(seq_along(x), x[!na_idx]), sum(na_idx))
	
		if (prod(x_dim) == length(x)) {
			dim(x) <- x_dim
			dimnames(x) <- x_dimnames
		}
	
		x
	
	} else {
	
		integer()
	}
}

aq_times <- function(fheader, empirically_correct=TRUE) {

	aq <- unname(fheader["SIZE",]/fheader["SW",]*fheader["TDSIZE",]/fheader["FTSIZE",])
	
	# Correct for FID offset as a result of digital oversampling
	aq <- aq*unname(1-fheader["DMXVAL",]/fheader["TDSIZE",])
	
	if (empirically_correct) {
	
		correction_factor <- ifelse(fheader["APODCODE",] == 1, 1, 0.5)
		aq <- aq*unname(1-correction_factor/fheader["TDSIZE",])
	}
	
	aq
}

pack_fit_params <- function(start_list, group_list) {

	stopifnot(all(tapply(start_list[["omega0"]][group_list[["omega0"]] != 0], group_list[["omega0"]][group_list[["omega0"]] != 0], function(x) all(x==x[1]))))
	stopifnot(all(tapply(start_list[["r2"]][group_list[["r2"]] != 0], group_list[["r2"]][group_list[["r2"]] != 0], function(x) all(x==x[1]))))
	stopifnot(all(tapply(start_list[["m0"]][group_list[["m0"]] != 0], group_list[["m0"]][group_list[["m0"]] != 0], function(x) all(x==x[1]))))
	stopifnot(all(tapply(start_list[["p0"]][group_list[["p0"]] != 0], group_list[["p0"]][group_list[["p0"]] != 0], function(x) all(x==x[1]))))
	stopifnot(all(tapply(start_list[["p1"]][group_list[["p1"]] != 0], group_list[["p1"]][group_list[["p1"]] != 0], function(x) all(x==x[1]))))
	stopifnot(all(tapply(start_list[["omega0_comb"]][group_list[["omega0_comb"]] != 0], group_list[["omega0_comb"]][group_list[["omega0_comb"]] != 0], function(x) all(x==x[1]))))
	stopifnot(all(tapply(start_list[["field"]][group_list[["field"]] != 0], group_list[["field"]][group_list[["field"]] != 0], function(x) all(x==x[1]))))

	omega0_group_unique <- setdiff(unique(as.vector(group_list[["omega0"]])), 0)
	omega0_vec <- start_list[["omega0"]][match(omega0_group_unique, group_list[["omega0"]])]
	if (length(omega0_group_unique)) {
		names(omega0_vec) <- paste("omega0", omega0_group_unique, sep="_")
	}
	
	r2_group_unique <- setdiff(unique(as.vector(group_list[["r2"]])), 0)
	r2_vec <- start_list[["r2"]][match(r2_group_unique, group_list[["r2"]])]
	if (length(r2_group_unique)) {
		names(r2_vec) <- paste("r2", r2_group_unique, sep="_")
	}
	
	m0_group_unique <- setdiff(unique(as.vector(group_list[["m0"]])), 0)
	m0_vec <- start_list[["m0"]][match(m0_group_unique, group_list[["m0"]])]
	if (length(m0_group_unique)) {
		names(m0_vec) <- paste("m0", m0_group_unique, sep="_")
	}
	
	p0_group_unique <- setdiff(unique(as.vector(group_list[["p0"]])), 0)
	p0_vec <- start_list[["p0"]][match(p0_group_unique, group_list[["p0"]])]
	if (length(p0_group_unique)) {
		names(p0_vec) <- paste("p0", p0_group_unique, sep="_")
	}
	
	p1_group_unique <- setdiff(unique(as.vector(group_list[["p1"]])), 0)
	p1_vec <- start_list[["p1"]][match(p1_group_unique, group_list[["p1"]])]
	if (length(p1_group_unique)) {
		names(p1_vec) <- paste("p1", p1_group_unique, sep="_")
	}
	
	omega0_comb_group_unique <- setdiff(unique(as.vector(group_list[["omega0_comb"]])), 0)
	omega0_comb_vec <- start_list[["omega0_comb"]][match(omega0_comb_group_unique, group_list[["omega0_comb"]])]
	if (length(omega0_comb_group_unique)) {
		names(omega0_comb_vec) <- paste("omega0_comb", omega0_comb_group_unique, sep="_")
	}
	
	field_group_unique <- setdiff(unique(as.vector(group_list[["field"]])), 0)
	field_vec <- start_list[["field"]][match(field_group_unique, group_list[["field"]])]
	if (length(field_group_unique)) {
		names(field_vec) <- paste("field", field_group_unique, sep="_")
	}
	
	c(omega0_vec, r2_vec, m0_vec, p0_vec, p1_vec, omega0_comb_vec, field_vec)
}

unpack_fit_params <- function(start_vec, group_list, comb_list, default_list=group_list) {

	unpacked_list <- default_list
	
	unpacked_list[["omega0"]][group_list[["omega0"]] != 0] <- start_vec[paste("omega0", group_list[["omega0"]][group_list[["omega0"]] != 0], sep="_")]
	unpacked_list[["r2"]][group_list[["r2"]] != 0] <- start_vec[paste("r2", group_list[["r2"]][group_list[["r2"]] != 0], sep="_")]
	unpacked_list[["m0"]][group_list[["m0"]] != 0] <- start_vec[paste("m0", group_list[["m0"]][group_list[["m0"]] != 0], sep="_")]
	unpacked_list[["p0"]][group_list[["p0"]] != 0] <- start_vec[paste("p0", group_list[["p0"]][group_list[["p0"]] != 0], sep="_")]
	unpacked_list[["p1"]][group_list[["p1"]] != 0] <- start_vec[paste("p1", group_list[["p1"]][group_list[["p1"]] != 0], sep="_")]
	unpacked_list[["omega0_comb"]][group_list[["omega0_comb"]] != 0] <- start_vec[paste("omega0_comb", group_list[["omega0_comb"]][group_list[["omega0_comb"]] != 0], sep="_")]
	unpacked_list[["field"]][group_list[["field"]] != 0] <- start_vec[paste("field", group_list[["field"]][group_list[["field"]] != 0], sep="_")]
	
	unpacked_list[["omega0"]] <- comb_vec_to_param_array(unpacked_list[["omega0_comb"]], comb_list[["omega0"]], unpacked_list[["omega0"]])
	
	unpacked_list
}

group_param_idx <- function(param_names, group_list, start_list) {

	idx_list <- group_list
	
	for (i in seq_along(start_list)) {
		idx_list[[i]][] <- match(paste(names(group_list)[i], group_list[[i]], sep="_"), param_names)
	}
	
	idx_list
}

#' Determine array of destination parameters from a source vector
#'
#' @param vector of source parameters
#' @param comb_array array of 2xN data frames with mapping between vector and array
#' @param param_array starting array of destination parameters
#' @param na_only only overwrite values in param_array if they are NA
comb_vec_to_param_array <- function(comb_vec, comb_array, param_array=NULL, na_only=FALSE) {

	not_null_idx <- which(!sapply(comb_array, is.null))
	
	if (length(not_null_idx)) {
	
		if (is.null(param_array)) {
			param_array <- fill_array(NA_real_, dim(comb_array))
		}
	
		all_comb <- do.call(rbind, comb_array[not_null_idx])
	
		stopifnot(length(comb_vec) == max(all_comb[,1]))
	
		for (i in seq_along(not_null_idx)) {
			comb_df <- comb_array[[not_null_idx[i]]]
			if (!na_only || is.na(param_array[not_null_idx[i]])) {
				param_array[not_null_idx[i]] <- sum(comb_vec[comb_df[,1]]*comb_df[,2])
			}
		}
	}
	
	param_array
}

#' Determine vector of source parameters from destination array via least squares
#'
#' @param param_array array of destination parameters
#' @param comb_array array of 2xN data frames with mapping between vector and array
#' @param group_vec optional vector of group numbers for source parameters
#' @param resid_thresh all residuals from least squares fit must be less than this value
param_array_to_comb_vec <- function(param_array, comb_array, group_vec=NULL, resid_thresh=1e-8) {

	not_null_idx <- which(!sapply(comb_array, is.null))
	
	if (length(not_null_idx)) {
	
		# get all combination data
		all_comb <- do.call(rbind, comb_array[not_null_idx])
		
		if (is.null(group_vec)) {
			# default group vec
			group_vec <- seq_len(max(all_comb[,1]))
		} else {
			# replace any 0 values with unused numbers
			zero_idx <- group_vec == 0
			group_vec[zero_idx] <- utils::head(setdiff(seq_along(group_vec), group_vec[!zero_idx]), sum(zero_idx))
			# ensure group_vec ordered sequentially from 1
			group_vec <- match(group_vec, unique(group_vec))
		}
		
		# make sure all combinations from 1 to the maximum are represented
		stopifnot(seq_len(max(all_comb[,1])) %in% all_comb[,1])
		
		y <- param_array[not_null_idx]
		x <- matrix(0, nrow=length(not_null_idx), ncol=max(group_vec))
		
		empty_full_values <- numeric(length(group_vec))
		
		for (i in seq_along(not_null_idx)) {
			full_values <- empty_full_values
			idx <- comb_array[[not_null_idx[i]]][,1]
			full_values[idx] <- comb_array[[not_null_idx[i]]][,2]
			x[i,] <- tapply(full_values, group_vec, sum)
		}
		
		fit <- stats::lsfit(x, y, intercept=FALSE)
		
		stopifnot(abs(fit$residuals) < resid_thresh)
		
		unname(fit$coefficients[group_vec])
		
	} else {
	
		numeric()
	}
}

#' Prepare input data structure for peak fitting
#'
#' @export
make_fit_input <- function(spectra, omega0_start, omega0_plus, omega0_minus=omega0_plus, r2_start=NULL, m0_start=NULL, m0_region=(omega0_plus+omega0_minus)/2, p0_start=0, p1_start=0, omega0_group=NULL, r2_group=NULL, m0_group=NULL, p0_group=0, p1_group=0, omega0_comb=NULL, omega0_comb_start=NULL, omega0_comb_group=NULL, field_offsets=numeric(), field_start=numeric(), field_group=0, fheader=NULL) {

	if (is.data.frame(spectra)) {
		fheader <- fheader
		spec_int <- spectra[,2]
		dim(spec_int) <- length(spec_int)
		spec_ppm <- list(spectra[,1])
		names(spec_ppm) <- "1H"
		dimnames(spec_int) <- spec_ppm
		spectra <- list(list(int=spec_int, ppm=spec_ppm, fheader=fheader, dim_idx=1))
		#print(str(spectra))
	}
	
	n_dimensions <- ncol(spectra[[1]][["fheader"]])
	n_spectra <- length(spectra)
	
	if (is.null(dim(omega0_start)) || length(dim(omega0_start)) != 3) {
		dim(omega0_start) <- c(n_dimensions, length(omega0_start)/n_dimensions/n_spectra, n_spectra)
	}
	
	omega0_plus <- fill_array(omega0_plus, dim(omega0_start))
	omega0_minus <- fill_array(omega0_minus, dim(omega0_start))
	
	r2_start <- fill_array(r2_start, dim(omega0_start))
	m0_start <- fill_array(m0_start, dim(omega0_start)[-1])
	p0_start <- fill_array(p0_start, dim(omega0_start))
	p1_start <- fill_array(p1_start, dim(omega0_start))
	
	omega0_group <- fill_group(omega0_group, dim(omega0_start))
	r2_group <- fill_group(r2_group, dim(omega0_start))
	m0_group <- fill_group(m0_group, dim(omega0_start)[-1])
	p0_group <- fill_group(p0_group, dim(omega0_start), TRUE)
	p1_group <- fill_group(p1_group, dim(omega0_start), TRUE)
	
	omega0_comb <- fill_array_list(omega0_comb, dim(omega0_start), TRUE)
	
	if (is.null(omega0_comb_start)) {
		omega0_comb_start <- param_array_to_comb_vec(omega0_start, omega0_comb, omega0_comb_group)
	} else {
		omega0_start <- comb_vec_to_param_array(omega0_comb_start, omega0_comb, omega0_start)
	}
	
	omega0_comb_group <- fill_comb_group(omega0_comb_group, omega0_comb)
	
	if (!is.matrix(field_offsets)) {
		field_offsets <- matrix(field_offsets, nrow=length(field_offsets), ncol=n_spectra)
	}
	
	field_start <- fill_array(field_start, dim(field_offsets))
	field_group <- fill_group(field_group, dim(field_offsets))
	
	# make sure group setting for combinations is 0
	stopifnot(omega0_group[!sapply(omega0_comb, is.null)] == 0)
	
	if (FALSE) {
	stopifnot(all(dim(omega0_start) == dim(omega0_plus)))
	if (length(r2_start) > 1) stopifnot(all(dim(r2_start) == dim(omega0_start)))
	if (length(m0_start) > 1) stopifnot(all(dim(m0_start) == dim(omega0_start)[-1]))
	if (length(omega0_group) > 1) {
		stopifnot(all(dim(omega0_group) == dim(omega0_start)))
		stopifnot(all(tapply(omega0_start, omega0_group, function(x) all(x==x[1]))))
	}
	if (length(r2_group) > 1) {
		stopifnot(all(dim(r2_group) == dim(omega0_start)))
		if (length(r2_start) > 1) stopifnot(all(tapply(r2_start, r2_group, function(x) all(x==x[1]))))
	}
	if (length(m0_group) > 1) {
		stopifnot(all(dim(m0_group) == dim(omega0_start)))
		if (length(m0_start) > 1) stopifnot(all(tapply(m0_start, m0_group, function(x) all(x==x[1]))))
	}
	}
	
	peak_idx <- peak_dimnames <- array(list(), dim=dim(omega0_start)[2:3])
	
	spec_data <- vector("list", length(spectra))
	names(spec_data) <- names(spectra)
	
	spec_offset <- 0
	
	for (i in seq_along(spectra)) {
	
		peak_idx_list <- vector("list", dim(omega0_start)[2])
		
		spec_int <- spectra[[i]][["int"]]
		spec_ppm <- spectra[[i]][["ppm"]]
		
		for (j in seq_along(peak_idx_list)) {
		
			roi_idx <- roi_dimnames <- vector("list", length(spec_ppm))
			
			for (k in seq_along(spec_ppm)) {
				# take aliasing into account
				sw_ppm <- spectra[[i]][["fheader"]]["SW",k]/spectra[[i]][["fheader"]]["OBS",k]
				sw_ppm <- 0
				roi_idx[[k]] <- (
					(spec_ppm[[k]] >= omega0_start[k,j,i]-omega0_minus[k,j,i] & spec_ppm[[k]] <= omega0_start[k,j,i]+omega0_plus[k,j,i]) |
					(spec_ppm[[k]]+sw_ppm >= omega0_start[k,j,i]-omega0_minus[k,j,i] & spec_ppm[[k]]+sw_ppm <= omega0_start[k,j,i]+omega0_plus[k,j,i]) |
					(spec_ppm[[k]]-sw_ppm >= omega0_start[k,j,i]-omega0_minus[k,j,i] & spec_ppm[[k]]-sw_ppm <= omega0_start[k,j,i]+omega0_plus[k,j,i])
				)
				roi_dimnames[[k]] <- dimnames(spec_int)[[k]][roi_idx[[k]]]
			}
			
			peak_mask <- array(FALSE, dim=dim(spec_int))
			
			if (length(dim(spec_int)) == 1) {
				peak_mask[roi_idx[[1]]] <- TRUE
			} else if (length(dim(spec_int)) == 2) {
				peak_mask[roi_idx[[1]],roi_idx[[2]]] <- TRUE
			} else if (length(dim(spec_int)) == 3) {
				peak_mask[roi_idx[[1]],roi_idx[[2]],roi_idx[[3]]] <- TRUE
			} else if (length(dim(spec_int)) == 4) {
				peak_mask[roi_idx[[1]],roi_idx[[2]],roi_idx[[3]],roi_idx[[4]]] <- TRUE
			} else {
				stop()
			}
			
			peak_idx_list[[j]] <- which(peak_mask)
			
			peak_dimnames[[j,i]] <- roi_dimnames
		}
		
		peak_idx_unique <- unique(sort(unlist(peak_idx_list)))
		
		for (j in seq_along(peak_idx_list)) {
		
			peak_idx_sub <- match(peak_idx_list[[j]], peak_idx_unique)
			peak_idx[[j,i]] <- peak_idx_sub
		}
		
		spec_mask <- array(FALSE, dim=dim(spec_int))
		spec_mask[peak_idx_unique] <- TRUE
		
		spec_eval_idx <- which(spec_mask, arr.ind=TRUE)
		
		#print(str(spec_eval_idx))
		
		#print(spectra[[i]][["fheader"]])
		
		omega_eval_idx <- lapply(seq_len(ncol(spec_eval_idx)), function(j) unique(sort(spec_eval_idx[,j])))
		
		#spec_eval_idx <- sapply(seq_len(ncol(spec_eval_idx)), function(j) match(omega_eval_idx[[j]], spec_eval_idx[,j]))
		spec_eval_idx <- sapply(seq_len(ncol(spec_eval_idx)), function(j) match(spec_eval_idx[,j], omega_eval_idx[[j]]))
		
		omega_eval <- lapply(seq_len(ncol(spec_eval_idx)), function(j) spec_ppm[[j]][omega_eval_idx[[j]]])
		
		fheader <- spectra[[i]][["fheader"]]
		
		p1_frac <- lapply(seq_along(omega_eval), function(j) {
			frac <- seq(0, 1, length.out=fheader["FTSIZE",j]+1)
			frac <- frac[-length(frac)]
			if (all(fheader[c("X1","XN"),j] != 0)) {
				frac <- frac[seq(fheader["X1",j], fheader["XN",j])]
			}
			frac[omega_eval_idx[[j]]]
		})
		
		aq_time <- aq_times(fheader)
		
		fit_func <- vector("list", ncol(fheader))
		
		for (j in seq_len(ncol(fheader))) {
			
			fit_func_name <- "none"
			fit_func_data <- c(aq=aq_time[j])
	
			if (fheader["APODCODE",j] == 1) {
		
				stopifnot(fheader["APODQ3",j] %in% c(1,2))
		
				if (fheader["APODQ3",j] == 1) {
					fit_func_name <- "sp1"
					fit_func_data <- c(fit_func_data, off=unname(fheader["APODQ1",j]), end=unname(fheader["APODQ2",j]))
				} else if (fheader["APODQ3",j] == 2) {
					fit_func_name <- "sp2"
					fit_func_data <- c(fit_func_data, off=unname(fheader["APODQ1",j]), end=unname(fheader["APODQ2",j]))
				} else {
					stop()
				}
			}
			
			fit_func[[j]] <- list(
				formulas=lapply(lineshapes[[fit_func_name]], function(func_text) parse(text=func_text)),
				data=as.list(fit_func_data)
			)
		}
		
		omega_idx_ranges <- lapply(omega_eval_idx, range)
		
		omega_contigous <- lapply(seq_along(omega_idx_ranges), function(j) {
			spectra[[i]][["ppm"]][[j]][seq(omega_idx_ranges[[j]][1], omega_idx_ranges[[j]][2])]
		})
		names(omega_contigous) <- names(spec_ppm)
		
		spec_data[[i]] <- list(
			ref_freq=unname(fheader["OBS",]),
			omega_eval=omega_eval,
			omega_contigous=omega_contigous,
			p1_frac=p1_frac,
			spec_eval_idx=spec_eval_idx,
			spec_int=unname(spec_int[peak_idx_unique]),
			spec_offset=spec_offset,
			fit_func=fit_func
		)
		
		spec_offset <- spec_offset+length(peak_idx_unique)
	}
	
	#print(spec_data)
	
	start_list <- list(
		omega0=omega0_start,
		r2=r2_start,
		m0=m0_start,
		p0=p0_start,
		p1=p1_start,
		omega0_comb=omega0_comb_start,
		field=field_start
	)
	
	group_list <- list(
		omega0=omega0_group,
		r2=r2_group,
		m0=m0_group,
		p0=p0_group,
		p1=p1_group,
		omega0_comb=omega0_comb_group,
		field=field_group
	)
	
	comb_list <- list(
		omega0=omega0_comb
	)
	
	lower_list <- upper_list <- start_list
	
	for (i in seq_along(lower_list)) {
		lower_list[[i]][] <- -Inf
		upper_list[[i]][] <- Inf
	}
	
	lower_list[["omega0"]] <- omega0_start - omega0_minus
	upper_list[["omega0"]] <- omega0_start + omega0_plus
	lower_list[["r2"]][] <- 0
	lower_list[["field"]][] <- 0
	upper_list[["field"]][] <- 0.5
	
	#print(start_list)
	
	list(
		spec_data=spec_data,
		start_list=start_list,
		group_list=group_list,
		comb_list=comb_list,
		field_offsets=field_offsets,
		lower_list=lower_list,
		upper_list=upper_list,
		num_points=spec_offset
	)
}

eval_peak_1d <- function(func_list, func_data, ref_mhz, omega, omega0, r2, p0=0, p1=0, p1_frac=0) {

	func_data <- c(list(
		omega=omega*ref_mhz*2*pi, # convert ppm to rad/s
		omega0=omega0*ref_mhz*2*pi, # convert ppm to rad/s
		r2=r2*2*pi # convert Hz to rad/s
	), func_data)
	
	p0 <- p0*pi/180 # convert degrees to rad
	p1 <- p1*pi/180 # convert degrees to rad
	
	func <- eval(func_list$func, func_data)
	#func <- eval(lineshapes_simplified_func$none, func_data)
	
	pvec <- exp(1i*(p0+p1*p1_frac))
	
	Re(func*pvec)
}

eval_peak_1d_deriv <- function(func_list, func_data, ref_mhz, omega, omega0, r2, p0=0, p1=0, p1_frac=0) {

	func_data <- c(list(
		omega=omega*ref_mhz*2*pi, # convert ppm to rad/s
		omega0=omega0*ref_mhz*2*pi, # convert ppm to rad/s
		r2=r2*2*pi # convert Hz to rad/s
	), func_data)
	
	p0 <- p0*pi/180 # convert degrees to rad
	p1 <- p1*pi/180 # convert degrees to rad
	
	#foo <- eval(lineshapes_simplified$none, func_data)
	#print(foo)
	
	if (TRUE) {
		func <- eval(func_list$func, func_data)
		dfunc_domega0 <- eval(func_list$domega0, func_data)
		dfunc_dr2 <- eval(func_list$dr2, func_data)
	} else {
		eval_vec <- eval(lineshapes_simplified$none, func_data)
		func <- eval_vec[[1]]
		dfunc_domega0 <- eval_vec[[2]]
		dfunc_dr2 <- eval_vec[[3]]
	}
	
	#pvec <-          exp(1i*(p0 + p1*p1_frac))
	#dpvec_dp0 <-  1i*exp(1i*(p0 + p1*p1_frac))
	#dpvec_dp1 <- -1i*exp(1i*(p0 + p1*p1_frac))*p1_frac
	
	pvec <- exp(1i*(p0 + p1*p1_frac))
	dpvec_dp0 <- 1i*pvec
	dpvec_dp1 <- dpvec_dp0*p1_frac
	
	func_re <- Re(func*pvec)
	dfunc_domega0_re <- Re(dfunc_domega0*pvec)*ref_mhz*2*pi # convert rad/s back to ppm
	dfunc_dr2_re <- Re(dfunc_dr2*pvec)*2*pi # convert rad/s back to Hz 
	dfunc_dp0_re <- Re(func*dpvec_dp0)*pi/180 # convert rad back to degrees
	dfunc_dp1_re <- Re(func*dpvec_dp1)*pi/180 # convert rad back to degrees
	
	cbind(
		f=func_re,
		omega0=dfunc_domega0_re,
		r2=dfunc_dr2_re,
		p0=dfunc_dp0_re,
		p1=dfunc_dp1_re
	)
}

fit_fn <- function(par, fit_data, return_resid=TRUE) {

	func_eval <- numeric(fit_data$num_points)

	param_list_orig <- unpack_fit_params(par, fit_data$group_list, fit_data$comb_list, fit_data$start_list)

	for (spec_idx in seq_along(fit_data$spec_data)) {
	
		spec_data <- fit_data[["spec_data"]][[spec_idx]]
		
		spec_eval_idx <- seq_along(spec_data$spec_int)+spec_data$spec_offset
		
		for (field_idx in seq(0, nrow(fit_data$field_offsets))) {
		
			param_list <- param_list_orig
		
			if (field_idx == 0) {
				field_weight <- 1-sum(param_list[["field"]][,spec_idx])
			} else {
				field_weight <- param_list[["field"]][field_idx,spec_idx]
				param_list[["omega0"]] <- param_list[["omega0"]]+fit_data$field_offsets[field_idx,spec_idx]
			}
		
			for (peak_idx in seq_len(dim(fit_data[["start_list"]][["omega0"]])[2])) {
		
				func_1d_evals <- lapply(seq_along(spec_data$omega_eval), function(dim_idx) {
			
					eval_peak_1d(
						spec_data$fit_func[[dim_idx]]$formulas,
						spec_data$fit_func[[dim_idx]]$data,
						spec_data$ref_freq[dim_idx],
						spec_data$omega_eval[[dim_idx]],
						param_list[["omega0"]][dim_idx,peak_idx,spec_idx],
						param_list[["r2"]][dim_idx,peak_idx,spec_idx],
						param_list[["p0"]][dim_idx,peak_idx,spec_idx],
						param_list[["p1"]][dim_idx,peak_idx,spec_idx],
						spec_data$p1_frac[[dim_idx]]
					)
				})
			
				func_nd_prod <- func_1d_evals[[1]][spec_data$spec_eval_idx[,1]]*param_list[["m0"]][peak_idx,spec_idx]
				for (i in seq_len(ncol(spec_data$spec_eval_idx))[-1]) {
					func_nd_prod <- func_nd_prod*func_1d_evals[[i]][spec_data$spec_eval_idx[,i]]
				}
			
				func_eval[spec_eval_idx] <- func_eval[spec_eval_idx] + func_nd_prod*field_weight
			}
		}
		
		if (return_resid) {
			func_eval[spec_eval_idx] <- spec_data$spec_int - func_eval[spec_eval_idx]
		}
	}
	
	func_eval
}

fit_jac <- function(par, fit_data) {

	jac_eval <- matrix(0, nrow=fit_data$num_points, ncol=length(par), dimnames=list(NULL, names(par)))

	param_list_orig <- unpack_fit_params(par, fit_data$group_list, fit_data$comb_list, fit_data$start_list)

	idx_list <- group_param_idx(names(par), fit_data$group_list, fit_data$start_list)

	for (spec_idx in seq_along(fit_data$spec_data)) {
	
		spec_data <- fit_data[["spec_data"]][[spec_idx]]
		
		spec_eval_idx <- seq_along(spec_data$spec_int)+spec_data$spec_offset
		
		field_eval_idx <- !is.na(idx_list[["field"]][,spec_idx])
		
		for (field_idx in seq(0, nrow(fit_data$field_offsets))) {
		
			field_eval <- (field_idx == 0 && any(field_eval_idx)) || (field_idx != 0 && !is.na(idx_list[["field"]][field_idx,spec_idx]))
		
			param_list <- param_list_orig
		
			if (field_idx == 0) {
				field_weight <- 1-sum(param_list[["field"]][,spec_idx])
			} else {
				field_weight <- param_list[["field"]][field_idx,spec_idx]
				param_list[["omega0"]] <- param_list[["omega0"]]+fit_data$field_offsets[field_idx,spec_idx]
			}
		
			for (peak_idx in seq_len(dim(fit_data[["start_list"]][["omega0"]])[2])) {
		
				deriv_1d_evals <- lapply(seq_along(spec_data$omega_eval), function(dim_idx) {
			
					eval_peak_1d_deriv(
						spec_data$fit_func[[dim_idx]]$formulas,
						spec_data$fit_func[[dim_idx]]$data,
						spec_data$ref_freq[dim_idx],
						spec_data$omega_eval[[dim_idx]],
						param_list[["omega0"]][dim_idx,peak_idx,spec_idx],
						param_list[["r2"]][dim_idx,peak_idx,spec_idx],
						param_list[["p0"]][dim_idx,peak_idx,spec_idx],
						param_list[["p1"]][dim_idx,peak_idx,spec_idx],
						spec_data$p1_frac[[dim_idx]]
					)
				})
			
				#print(str(deriv_1d_evals))
			
				if (FALSE) {
				omega0_idx <- which(!is.na(idx_list[["omega0"]][,peak_idx,spec_idx]))
				for (idx in omega0_idx) {
					deriv_nd_prod <- deriv_1d_evals[[idx]][spec_data$spec_eval_idx[,idx],"omega0"]*param_list[["m0"]][peak_idx,spec_idx]
					for (i in seq_len(ncol(spec_data$spec_eval_idx))[-idx]) {
						deriv_nd_prod <- deriv_nd_prod*deriv_1d_evals[[i]][spec_data$spec_eval_idx[,i],"f"]
					}
					jac_eval[spec_eval_idx,idx_list[["omega0"]][idx,peak_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[["omega0"]][idx,peak_idx,spec_idx]] + deriv_nd_prod*field_weight
				}
				}
			
				for (var_name in c("omega0", "r2", "p0", "p1")) {
					var_idx <- which(!is.na(idx_list[[var_name]][,peak_idx,spec_idx]))
					for (idx in var_idx) {
						deriv_nd_prod <- deriv_1d_evals[[idx]][spec_data$spec_eval_idx[,idx],var_name]*param_list[["m0"]][peak_idx,spec_idx]
						for (i in seq_len(ncol(spec_data$spec_eval_idx))[-idx]) {
							deriv_nd_prod <- deriv_nd_prod*deriv_1d_evals[[i]][spec_data$spec_eval_idx[,i],"f"]
						}
						jac_eval[spec_eval_idx,idx_list[[var_name]][idx,peak_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[[var_name]][idx,peak_idx,spec_idx]] + deriv_nd_prod*field_weight
					}
				}
			
				if (!is.na(idx_list[["m0"]][peak_idx,spec_idx]) || field_eval) {
					func_nd_prod <- deriv_1d_evals[[1]][spec_data$spec_eval_idx[,1],"f"]
					for (i in seq_len(ncol(spec_data$spec_eval_idx))[-1]) {
						func_nd_prod <- func_nd_prod*deriv_1d_evals[[i]][spec_data$spec_eval_idx[,i],"f"]
					}
					if (!is.na(idx_list[["m0"]][peak_idx,spec_idx])) {
						jac_eval[spec_eval_idx,idx_list[["m0"]][peak_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[["m0"]][peak_idx,spec_idx]] + func_nd_prod*field_weight
					}
					if (field_eval) {
						if (field_idx == 0) {
							jac_eval[spec_eval_idx,idx_list[["field"]][field_eval_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[["field"]][field_eval_idx,spec_idx]] - func_nd_prod*param_list[["m0"]][peak_idx,spec_idx]
						} else {
							jac_eval[spec_eval_idx,idx_list[["field"]][field_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[["field"]][field_idx,spec_idx]] + func_nd_prod*param_list[["m0"]][peak_idx,spec_idx]
						}
					}
				}
			
				for (var_name in c("omega0")) {
					var_comb_name <- paste(var_name, "comb", sep="_")
					var_idx <- which(!sapply(fit_data$comb_list[[var_name]][,peak_idx,spec_idx], is.null))
					for (idx in var_idx) {
						deriv_nd_prod <- deriv_1d_evals[[idx]][spec_data$spec_eval_idx[,idx],var_name]*param_list[["m0"]][peak_idx,spec_idx]
						for (i in seq_len(ncol(spec_data$spec_eval_idx))[-idx]) {
							deriv_nd_prod <- deriv_nd_prod*deriv_1d_evals[[i]][spec_data$spec_eval_idx[,i],"f"]
						}
						comb_data <- fit_data$comb_list[[var_name]][[idx,peak_idx,spec_idx]]
						for (comb_idx in which(!is.na(idx_list[[var_comb_name]][comb_data[,1]]))) {
							jac_eval[spec_eval_idx,idx_list[[var_comb_name]][comb_data[comb_idx,1]]] <- jac_eval[spec_eval_idx,idx_list[[var_comb_name]][comb_data[comb_idx,1]]] + deriv_nd_prod*comb_data[comb_idx,2]*field_weight
						}
					}
				}
			}
		}
	}
	
	#print(jac_eval)
	
	-jac_eval
}

#' Perform a fit with an input data structure
#'
#' @export
perform_fit <- function(fit_input) {

	fit_par <- pack_fit_params(fit_input$start_list, fit_input$group_list)
	fit_lower <- pack_fit_params(fit_input$lower_list, fit_input$group_list)
	fit_upper <- pack_fit_params(fit_input$upper_list, fit_input$group_list)
		
	fit <- suppressWarnings(minpack.lm::nls.lm(fit_par, fit_lower, fit_upper, fn=fit_fn, jac=fit_jac, fit_data=fit_input, control=list(minpack.lm::nls.lm.control(nprint=0), maxiter = 1e5)))
	
	fit_input[["fit_list"]] <- unpack_fit_params(fit$par, fit_input$group_list, fit_input$comb_list, default_list=fit_input$start_list)
	fit_input[["fit_rsstrace"]] <- fit$rsstrace
		
	fit_input
}

fit_jac_nlfb <- function(par, fit_data) {

	jj <- fit_jac(par, fit_data)
	attr(jj, "gradient") <- jj
	jj
}

#perform_fit_nlfb <- function(fit_input) {
#
#	fit_par <- pack_fit_params(fit_input$start_list, fit_input$group_list)
#
#	fit <- nlfb(fit_par, resfn=fit_fn, jacfn=fit_jac_nlfb, fit_data=fit_input)
#	
#	fit_input[["fit_list"]] <- unpack_fit_params(fit$coefficients, fit_input$group_list, fit_input$comb_list, default_list=fit_input$start_list)
#	
#	fit_input
#}

#' Get arrays of spectral intensities for input, starting parameters, and fit peaks
#'
#' @export
get_spec_int <- function(fit_data, spec_type=c("input", "start", "fit"), spec_idx=seq_along(fit_data$spec_data), peak_idx=seq_len(dim(fit_data$start_list$omega0)[2])) {

	spec_type <- match.arg(spec_type)
	
	all_peak_idx <- seq_len(dim(fit_data$start_list$omega0)[2])
	
	if (spec_type == "input") {
	
		int <- do.call(c, lapply(fit_data$spec_data, function(x) x$spec_int))
		
		stopifnot(all(all_peak_idx %in% peak_idx))
	
	} else {
	
		if (spec_type == "start") {
		
			param_list <- fit_data[["start_list"]]
		
		} else if (spec_type == "fit") {
		
			param_list <- fit_data[["fit_list"]]
		}
		
		param_list$m0[!all_peak_idx %in% peak_idx,] <- 0
		
		#print(str(param_list))
		
		# make sure all m0 values are independent so they can be turned off
		fit_data$group_list$m0[] <- seq_along(fit_data$group_list$m0)
		params <- pack_fit_params(param_list, fit_data$group_list)
		
		int <- fit_fn(params, fit_data, FALSE)
	}
	
	lapply(spec_idx, function(spec_idx) {
	
		spec_data <- fit_data$spec_data[[spec_idx]]
	
		spec_int <- array(NA_real_, sapply(spec_data$omega_contigous, length), dimnames=spec_data$omega_contigous)
		
		omega_contig_idx <- lapply(seq_along(spec_data$omega_contigous), function(i) {
			match(spec_data$omega_eval[[i]], spec_data$omega_contigous[[i]])
		})
		
		spec_contig_idx <- spec_data$spec_eval_idx
		for (i in seq_len(ncol(spec_contig_idx))) {
			spec_contig_idx[,i] <- omega_contig_idx[[i]][spec_contig_idx[,i]]
		}
		
		spec_int[spec_contig_idx] <- int[seq_len(nrow(spec_contig_idx))+spec_data$spec_offset]
		
		spec_int
	})
}

#' Plot a one dimensional peak fit
#'
#' @export
plot_fit_1d <- function(fit_data, always_show_start=FALSE) {

	original_int <- unlist(lapply(fit_data$spec_data, function(spec_data) spec_data$spec_int))
	
	start_int <- if (!"fit_list" %in% names(fit_data) || always_show_start) {
		start_par <- pack_fit_params(fit_data$start_list, fit_data$group_list)
		fit_fn(start_par, fit_data, return_resid=FALSE)
	}
	
	fit_int <- if ("fit_list" %in% names(fit_data)) {
		fit_par <- pack_fit_params(fit_data$fit_list, fit_data$group_list)
		fit_fn(fit_par, fit_data, return_resid=FALSE)
	}
	
	for (spec_data in fit_data$spec_data) {
	
		omega_ppm <- spec_data$omega_eval[[1]]
		
		omega_ppm_seg_ends <- which(abs(diff(omega_ppm)) > abs(stats::median(diff(omega_ppm)))*2)
		omega_ppm_seg_starts <- c(1, omega_ppm_seg_ends+1)
		omega_ppm_seg_ends <- c(omega_ppm_seg_ends, length(omega_ppm))
		plot_idx <- unlist(lapply(seq_along(omega_ppm_seg_starts), function(i) c(NA, seq(omega_ppm_seg_starts[i], omega_ppm_seg_ends[i]))))[-1]
		
		spec_eval_idx <- seq_along(spec_data$spec_int)+spec_data$spec_offset
		
		spec_original_int <- original_int[spec_eval_idx]
		spec_start_int <- start_int[spec_eval_idx]
		spec_fit_int <- fit_int[spec_eval_idx]
		
		ylim <- range(0, spec_original_int, spec_start_int, spec_fit_int)
		
		graphics::plot(omega_ppm[plot_idx], spec_original_int[plot_idx], type="l", xlim=rev(range(omega_ppm)), ylim=ylim, xlab=expression(delta (ppm)), ylab="Intensity")

		if (!is.null(spec_start_int)) {
			graphics::points(omega_ppm[plot_idx], spec_start_int[plot_idx], type="l", col="blue")
		}
		
		if (!is.null(spec_fit_int)) {
			graphics::points(omega_ppm[plot_idx], spec_fit_int[plot_idx], type="l", col="red")
		}
	}
}

#' Plot a two dimensional peak fit
#'
#' @export
plot_fit_2d <- function(fit_output, spec_ord=seq_len(dim(fit_output$start_list$omega0)[1]), always_show_start=FALSE, main=NULL) {

	input_spec_int <- fitnmr::get_spec_int(fit_output, "input")
	if ("fit_list" %in% names(fit_output)) {
		fit_spec_int <- fitnmr::get_spec_int(fit_output, "fit")
		plot_start <- always_show_start
	} else {
		plot_start <- TRUE
	}
	
	ref_freq_mat <- sapply(fit_output$spec_data, "[[", "ref_freq")
	ref_freq_mat_exp <- ref_freq_mat[,rep(seq_len(ncol(ref_freq_mat)), each=dim(fit_output$start_list$r2)[2])]

	fit_r2_ppm <- fit_output$fit_list$r2/as.numeric(ref_freq_mat_exp)
	fit_r2_ppm_low <- fit_output$fit_list$omega0-fit_r2_ppm
	fit_r2_ppm_high <- fit_output$fit_list$omega0+fit_r2_ppm
	
	if (plot_start) {
		start_spec_int <- fitnmr::get_spec_int(fit_output, "start")
		
		start_r2_ppm <- fit_output$start_list$r2/as.numeric(ref_freq_mat_exp)
		start_r2_ppm_low <- fit_output$start_list$omega0-start_r2_ppm
		start_r2_ppm_high <- fit_output$start_list$omega0+start_r2_ppm
	}
	
	for (i in seq_along(input_spec_int)) {

		zlim <- range(input_spec_int[[i]], na.rm=TRUE)
		fitnmr::contour_pipe(aperm(input_spec_int[[i]], spec_ord), zlim=zlim, col_pos="black", col_neg="gray")
		graphics::title(main)
		if (plot_start) {
			fitnmr::contour_pipe(aperm(start_spec_int[[i]], spec_ord), zlim=zlim, col_pos="blue", col_neg="lightblue", add=TRUE)
			graphics::points(t(fit_output$start_list$omega0[spec_ord,,i]), pch=16, col="blue")
			graphics::segments(fit_output$start_list$omega0[spec_ord[1],,i], start_r2_ppm_low[spec_ord[2],,i], fit_output$start_list$omega0[spec_ord[1],,i], start_r2_ppm_high[spec_ord[2],,i], col="blue")
			graphics::segments(start_r2_ppm_low[spec_ord[1],,i], fit_output$start_list$omega0[spec_ord[2],,i], start_r2_ppm_high[spec_ord[1],,i], fit_output$start_list$omega0[spec_ord[2],,i], col="blue")
		}
		if ("fit_list" %in% names(fit_output)) {
			fitnmr::contour_pipe(aperm(fit_spec_int[[i]], spec_ord), zlim=zlim, col_pos="red", col_neg="pink", add=TRUE)
			graphics::points(t(fit_output$fit_list$omega0[spec_ord,,i]), pch=16, col="red")
			graphics::segments(fit_output$fit_list$omega0[spec_ord[1],,i], fit_r2_ppm_low[spec_ord[2],,i], fit_output$fit_list$omega0[spec_ord[1],,i], fit_r2_ppm_high[spec_ord[2],,i], col="red")
			graphics::segments(fit_r2_ppm_low[spec_ord[1],,i], fit_output$fit_list$omega0[spec_ord[2],,i], fit_r2_ppm_high[spec_ord[1],,i], fit_output$fit_list$omega0[spec_ord[2],,i], col="red")
		}
		omega0_1_idx <- fitnmr::omega0_param_idx(fit_output, spec_ord[1], specs=i)
		omega0_2_idx <- fitnmr::omega0_param_idx(fit_output, spec_ord[2], specs=i)
		graphics::rect(param_values(fit_output$upper_list, omega0_1_idx), param_values(fit_output$upper_list, omega0_2_idx), param_values(fit_output$lower_list, omega0_1_idx), param_values(fit_output$lower_list, omega0_2_idx), border="gray")
	}
}

#' Plot Spectra Contours
#'
#' Plot a two dimensional contour plot
#'
#' The first dimension of \code{data_matrix} is drawn along the x-axis and the second dimension is drawn along the y-axis.
#'
#' If \code{low_frac} is specified, then there can actually be up to \code{nlevels} total positive contour levels and/or \code{nlevels} total negative contour levels, whichever has the larger magnitude in \code{zlim}. The other dimension will have the mirror image of those up to the relevant limit in \code{zlim}. Note that the levels do not actually go up to \code{zlim}. There are \code{nlevels+1} log spaced levels calculated from the contour determined by \code{low_frac} up to the maximum absolute \code{zlim}, but the final level is not drawn because it is at the maximum absolute \code{zlim}.
#'
#' If \code{low_frac} is set to \code{NA}, then the contour levels are calculated in the same way as \code{\link[graphics]{contour}}, with there being approximately \code{nlevels} linearly spaced levels in a range close to \code{zlim}.
#'
#' Note that this function does not directly take the data returned by \code{\link{read_nmrpipe}}. You must pass the \code{int} matrix from the value returned by that function.
#'
#' @param data_mat matrix with spectral intensities.
#' @param nlevels number of contour levels.
#' @param zlim minimum and maximum values at which to show contours.
#' @param low_frac minimum absolute value (as a fraction of \code{zlim}) at which to show contours.
#' @param lwd width of contour lines.
#' @param main title of the plot.
#' @param col_pos color of positive contours.
#' @param col_neg color of negative contours, defaults to lighter version of \code{col_pos}.
#' @param add logical indicating whether to add to an existing plot (i.e. not start a new one).
#' @param xlab label for x-axis, defaults to \code{names(dimnames(datamat))[1]}.
#' @param ylab label for y-axis, defaults to \code{names(dimnames(datamat))[2]}.
#' @param frame.plot a logical indicating whether a box should be drawn around the plot.
#'
#' @export
contour_pipe <- function(data_mat, nlevels=10, zlim=range(data_mat, na.rm=TRUE), low_frac=0.05, lwd=0.25, main=NA, col_pos="black", col_neg=grDevices::rgb(t(grDevices::col2rgb(col_pos))/255*0.25+0.75), add=FALSE, xlab=NULL, ylab=NULL, frame.plot=TRUE) {

	zlim <- zlim

	usr <- graphics::par("usr")

	x <- as.numeric(colnames(data_mat))
	y <- as.numeric(rownames(data_mat))
	
	if (add) {
		x_lim <- range(which(x >= usr[4] & x <= usr[3]))
		x_idx <- seq(max(x_lim[1]-1, 1), min(x_lim[2]+1, length(x)))
	
		y_lim <- range(which(y >= usr[2] & y <= usr[1]))
		y_idx <- seq(max(y_lim[1]-1, 1), min(y_lim[2]+1, length(y)))
	
		data_mat <- data_mat[y_idx,x_idx]
		x <- x[x_idx]
		y <- y[y_idx]
	}

	data_mat_transform <- (data_mat[rev(seq_len(nrow(data_mat))),rev(seq_len(ncol(data_mat)))])

	levels <- pretty(zlim, nlevels)
	if (!is.na(low_frac)) {
		max_int <- max(abs(zlim))
		levels <- exp(seq(log(max_int*low_frac), log(max_int), length.out=nlevels+1))
		levels <- c(-rev(levels), 0, levels)
		levels <- c(-rev(levels), levels)
		levels <- levels[levels > zlim[1] & levels < zlim[2] & levels]
	}
	col <- rep("gray", length(levels))
	col[levels < 0] <- col_neg
	col[levels > 0] <- col_pos

	if (is.null(xlab)) xlab <- paste(names(dimnames(data_mat))[1], "(ppm)")
	if (is.null(ylab)) ylab <- paste(names(dimnames(data_mat))[2], "(ppm)")

	graphics::contour(rev(y), rev(x), data_mat_transform, levels=levels, xlim=rev(range(y)), ylim=rev(range(x)), drawlabels=FALSE, add=add, col=col, lwd=lwd, xaxs="i", yaxs="i", main=main, xlab=xlab, ylab=ylab, frame.plot=frame.plot)
}

read_nmrdraw_peak_tab_old <- function(filepath) {
	vars <- readLines(filepath, 16)[16]
	vars <- strsplit(substring(vars, 8), " ")[[1]]
	
	peaktab <- utils::read.table(filepath, skip=18)
	colnames(peaktab) <- vars
	
	noise_level <- strsplit(readLines(filepath, 5)[5], " ")[[1]][3]
	noise_level <- strsplit(noise_level, ",")[[1]][1]
	noise_level <- as.numeric(noise_level)
	#print(noise_level)
	
	attr(peaktab, "noise") <- noise_level
	
	peaktab
}

#' Read an NMRDraw formatted peak table
#'
#' @export
read_nmrdraw_peak_tab <- function(file_path) {

	peak_tab_lines <- readLines(file_path)
	
	vars_line_idx <- grep("^VARS", peak_tab_lines)
	format_line_idx <- grep("^FORMAT", peak_tab_lines)
	
	var_names <- strsplit(peak_tab_lines[vars_line_idx], " +")[[1]][-1]
	formats <- strsplit(peak_tab_lines[format_line_idx], " +")[[1]][-1]
	
	blank_line_idx <- grep("^$", peak_tab_lines)
	
	peak_tab <- utils::read.table(text=peak_tab_lines, skip=utils::tail(blank_line_idx, 1), stringsAsFactors=FALSE)
	colnames(peak_tab) <- var_names
	
	peak_tab
}

peak_tab_formats <- c("%5d", "%9.3f", "%9.3f", "%6.3f", "%6.3f", "%8.3f", "%8.3f", "%9.3f", "%9.3f", "%7.3f", "%7.3f", "%8.3f", "%8.3f", "%4d", "%4d", "%4d", "%4d", "%+e", "%+e", "%+e", "%.5f", "%d", "%s", "%4d", "%4d")
names(peak_tab_formats) <- c("INDEX", "X_AXIS", "Y_AXIS", "DX", "DY", "X_PPM", "Y_PPM", "X_HZ", "Y_HZ", "XW", "YW", "XW_HZ", "YW_HZ", "X1", "X3", "Y1", "Y3", "HEIGHT", "DHEIGHT", "VOL", "PCHI2", "TYPE", "ASS", "CLUSTID", "MEMCNT")

#' Write an NMRDraw formatted peak table
#'
#' @export
write_nmrdraw_peak_tab <- function(peak_tab, file_path) {

	if (!"INDEX" %in% colnames(peak_tab)) {
		peak_tab <- cbind(INDEX=seq_len(nrow(peak_tab)), peak_tab)
	}
	
	peak_tab_format <- paste(peak_tab_formats[colnames(peak_tab)], collapse=" ")
	
	peak_tab_lines <- c(
		paste("VARS   ", paste(colnames(peak_tab), collapse=" ")),
		paste("FORMAT ", peak_tab_format),
		"",
		do.call(sprintf, c(list(peak_tab_format), as.list(as.data.frame(peak_tab))))
	)
	
	writeLines(peak_tab_lines, file_path)
}

#' Convert PPM values to points
#'
#' @export
ppm_to_pts <- function(ppm_mat, fheader) {

	orig <- fheader["CAR",]*fheader["OBS",]-fheader["SW",]/2+fheader["SW",]/fheader["FTSIZE",]
	
	# orig+fheader["SW",] seems wrong to me...
	points_mat <- t((orig+fheader["SW",]-t(ppm_mat)*fheader["OBS",])/fheader["SW",]*fheader["FTSIZE",])
	
	colnames(points_mat) <- sub("_PPM", "_AXIS", colnames(points_mat))
	
	points_mat
}

#' Convert widths in Hz into points
whz_to_pts <- function(whz_mat, fheader) {

	wpoints_mat <- t(t(whz_mat)/fheader["SW",]*fheader["FTSIZE",])
	
	colnames(wpoints_mat) <- sub("_HZ", "", colnames(wpoints_mat))
	
	wpoints_mat
}

collapse_args <- function(named_args) {

	named_args <- named_args[!is.na(named_args)]

	named_args <- as.vector(rbind(paste("-", names(named_args), sep=""), named_args))
	
	named_args[named_args != ""]
}

#' Extract parameters from fit object for use with make_fit_input
extract_params <- function(fit, expand=0) {

	fit_params <- fit$fit_list
	group_params <- fit$group_list
	comb_params <- fit$comb_list
	
	if (expand) {
	
		expand_idx <- c(seq_len(dim(fit_params$omega0)[2]), rep(NA, expand))
		new_idx <- is.na(expand_idx)
		
		fit_params$omega0 <- fit_params$omega0[,expand_idx,,drop=FALSE]
		fit_params$r2 <- fit_params$r2[,expand_idx,,drop=FALSE]
		fit_params$m0 <- fit_params$m0[expand_idx,,drop=FALSE]
		fit_params$p0 <- fit_params$p0[,expand_idx,,drop=FALSE]
		fit_params$p1 <- fit_params$p1[,expand_idx,,drop=FALSE]
		
		#fit_params$m0[new_idx,] <- 0
		fit_params$p0[,new_idx,] <- 0
		fit_params$p1[,new_idx,] <- 0
		
		group_params$omega0 <- group_params$omega0[,expand_idx,,drop=FALSE]
		group_params$r2 <- group_params$r2[,expand_idx,,drop=FALSE]
		group_params$m0 <- group_params$m0[expand_idx,,drop=FALSE]
		group_params$p0 <- group_params$p0[,expand_idx,,drop=FALSE]
		group_params$p1 <- group_params$p1[,expand_idx,,drop=FALSE]
		
		group_params$p0[,new_idx,] <- 0
		group_params$p1[,new_idx,] <- 0
		
		comb_params$omega0 <- comb_params$omega0[,expand_idx,,drop=FALSE]
	}
	
	names(fit_params) <- paste(names(fit_params), "_start", sep="")
	names(group_params) <- paste(names(group_params), "_group", sep="")
	names(comb_params) <- paste(names(comb_params), "_comb", sep="")
	
	c(fit_params, group_params)
}

#' Determine the region of a spectrum containing the majority of the fit peaks
fit_footprint <- function(fit, frac=0.12) {

	fit_spec_int <- fitnmr::get_spec_int(fit, "fit")
	
	fit_footprints <- lapply(seq_along(fit_spec_int), function(i) {
	
		sort_int <- sort(abs(fit_spec_int[[i]]))
		cum_int <- cumsum(sort_int)
		cum_frac_int <- cum_int/utils::tail(cum_int, 1)
	
		thresh <- sort_int[which(cum_frac_int > frac)[1]]
	
		fit_spec_int[[i]][is.na(fit_spec_int[[i]])] <- 0
	
		abs(fit_spec_int[[i]]) > thresh
		
		#abs(fit_spec_int[[1]]) > max(abs(fit_spec_int[[1]]))*.025
	})
	
	if (length(fit_footprints) == 1) {
		fit_footprints <- fit_footprints[[1]]
	}
	
	fit_footprints
}

#' Add upper/lower limits based on the r2 value
#'
#' @export
limit_omega0_by_r2 <- function(fit_input, factor=1.5) {

	ref_freq_mat <- sapply(fit_input$spec_data, "[[", "ref_freq")
	ref_freq_mat_exp <- ref_freq_mat[,rep(seq_len(ncol(ref_freq_mat)), each=dim(fit_input$start_list$r2)[2])]

	#print(fit_input$start_list$r2)
	#print(ref_freq_mat_exp)

	#r2_ppm <- fit_input$start_list$r2/spec_list[[1]]$fheader["OBS",]
	r2_ppm <- fit_input$start_list$r2/as.numeric(ref_freq_mat_exp)
	
	omega0_idx <- omega0_param_idx(fit_input)
	source_idx <- omega0_comb_source_idx(fit_input, omega0_idx)
	
	omega0_ppm <- param_values(fit_input$start_list, omega0_idx)
	param_values(fit_input$upper_list, omega0_idx) <- omega0_ppm+r2_ppm[source_idx]*factor
	param_values(fit_input$lower_list, omega0_idx) <- omega0_ppm-r2_ppm[source_idx]*factor
	
	fit_input
}

#' Update bounds on fitting parameters
#'
#' @export
update_fit_bounds <- function(fit_input, omega0_r2_factor=NULL, r2_bounds=NULL, sc_bounds=NULL) {

	if (!is.null(omega0_r2_factor)) {
		fit_input <- limit_omega0_by_r2(fit_input, omega0_r2_factor)
	}
	
	if (!is.null(r2_bounds)) {
		fit_input$lower_list$r2[] <- r2_bounds[1]
		fit_input$upper_list$r2[] <- r2_bounds[2]
	}
	
	if (!is.null(sc_bounds)) {
		coupling_idx <- coupling_param_idx(fit_input)
		param_values(fit_input$lower_list, coupling_idx) <- sc_bounds[1]
		param_values(fit_input$upper_list, coupling_idx) <- sc_bounds[2]
	}
	
	fit_input
}

omega0_bound_distance <- function(fit_output) {

	omega0_idx <- omega0_param_idx(fit_output)
	
	fit_ppm <- param_values(fit_output$fit_list, omega0_idx)
	lower_ppm <- param_values(fit_output$lower_list, omega0_idx)
	upper_ppm <- param_values(fit_output$upper_list, omega0_idx)
	
	pmin(abs(fit_ppm-lower_ppm), abs(fit_ppm-upper_ppm))
}

#' Get spectra for individual peaks
get_spec_peak_int <- function(fit_data, spec_type=c("start", "fit"), spec_idx=seq_along(fit_data$spec_data), peak_idx=seq_len(dim(fit_data$start_list$omega0)[2])) {

	spec_type <- match.arg(spec_type)
	
	spec_peak_list <- vector("list", length(spec_idx))
	
	for (i in seq_along(peak_idx)) {
	
		spec_list <- fitnmr::get_spec_int(fit_data, spec_type, spec_idx, peak_idx[i])
		
		for (j in seq_along(spec_list)) {
			spec_peak_list[[j]] <- c(spec_peak_list[[j]], spec_list[j])
		}
	}
	
	spec_peak_list
}

#' Determine a matrix of fractional peak overlap
spec_overlap_mat <- function(peak_int_list) {

	norm_peak_int_list <- lapply(peak_int_list, function(x) abs(x)/sum(abs(x), na.rm=TRUE))

	overlap_mat <- matrix(NA_real_, nrow=length(peak_int_list), ncol=length(peak_int_list))
	
	for (i in seq_len(nrow(overlap_mat)-1)) {
		for (j in seq(i+1, nrow(overlap_mat))) {
			overlap_mat[i,j] <- overlap_mat[j,i] <- sum(pmin(norm_peak_int_list[[i]], norm_peak_int_list[[j]]))
		}
	}
	
	overlap_mat
}

#' Fit peaks from a table of chemical shifts
#'
#' @export
fit_peaks <- function(spec_list, cs_mat, fit_prev=NULL, spec_ord=1:2, omega0_plus=c(0.075, 0.75), r2_start=5, r2_bounds=c(0.5, 20), sc_start=NULL, sc_bounds=c(0, Inf), positive_only=TRUE, plot_fit_stages=TRUE) {
	
	if (FALSE) {
	sc_start <- NULL
	if (!is.null(j_bounds)) {
		sc_start <- rep(NA_real_, length(j_bounds))
		bounds_not_null_idx <- !sapply(sc_bounds, is.null)
		if (sum(bounds_not_null_idx) > 0) {
			# the code currently only supports one type of scalar coupling
			stopifnot(sum(bounds_not_null_idx) == 1)
			sc_start[bounds_not_null_idx] <- sapply(sc_bounds[bounds_not_null_idx], mean)
			sc_bounds <- sc_bounds[[which(bounds_not_null_idx)]]
		}
	}
	
	if (is.null(sc_start) || all(is.na(sc_start))) {
		sc_bounds <- c(-Inf, Inf)
	}
	}
	
	if (length(r2_start) < length(spec_ord)) {
		r2_start <- rep(r2_start, length(spec_ord))
	}
	
	param_list_orig <- make_param_list(spec_list, cs_mat[,spec_ord,drop=FALSE], fit_prev=fit_prev, r2_start=r2_start[spec_ord], sc_start=sc_start[spec_ord], same_r2=TRUE, same_coupling=TRUE) 
	param_list <- param_list_orig
	omega0_idx <- omega0_param_idx(param_list)
	coupling_idx <- coupling_param_idx(param_list)
	
	# fix omega0 values
	param_values(param_list$group_list, omega0_idx) <- 0
	# fix r2 values
	param_list$group_list$r2[] <- 0
	# fix coupling values
	param_values(param_list$group_list, coupling_idx) <- 0
	
	fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list_to_arg_list(param_list)))
	if (positive_only) {
		fit_input$lower_list$m0[] <- 0
	}
	fit_output <- fitnmr::perform_fit(fit_input)
	param_list <- fit_output[c("fit_list", "group_list", "comb_list")]
	names(param_list)[1] <- "start_list"
	#print(lapply(param_list_orig, "[[", "omega0"))
	#print(lapply(param_list_orig, "[[", "omega0_comb"))
	#print(lapply(param_list_orig, "[[", "m0"))
	
	if (plot_fit_stages) plot_fit_2d(fit_output, spec_ord, always_show_start=TRUE, "Unfixed: m0")
	
	# unfix r2 values
	param_list$group_list$r2[] <- param_list_orig$group_list$r2[]
	# unfix coupling values
	param_values(param_list$group_list, coupling_idx) <- param_values(param_list_orig$group_list, coupling_idx)
	
	fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list_to_arg_list(param_list)))
	if (positive_only) {
		fit_input$lower_list$m0[] <- 0
	}
	fit_input$lower_list$r2[] <- r2_bounds[1]
	fit_input$upper_list$r2[] <- r2_bounds[2]
	param_values(fit_input$lower_list, coupling_idx) <- sc_bounds[1]
	param_values(fit_input$upper_list, coupling_idx) <- sc_bounds[2]
	fit_output <- fitnmr::perform_fit(fit_input)
	param_list <- fit_output[c("fit_list", "group_list", "comb_list")]
	names(param_list)[1] <- "start_list"
	
	if (plot_fit_stages) plot_fit_2d(fit_output, spec_ord, always_show_start=TRUE, "Unfixed: m0, r2, sc")
	
	# unfix omega0 values
	param_values(param_list$group_list, omega0_idx) <- param_values(param_list_orig$group_list, omega0_idx)
	
	fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list_to_arg_list(param_list)))
	if (positive_only) {
		fit_input$lower_list$m0[] <- 0
	}
	fit_input$lower_list$r2[] <- r2_bounds[1]
	fit_input$upper_list$r2[] <- r2_bounds[2]
	param_values(fit_input$lower_list, coupling_idx) <- sc_bounds[1]
	param_values(fit_input$upper_list, coupling_idx) <- sc_bounds[2]
	fit_input <- limit_omega0_by_r2(fit_input)
	fit_output <- fitnmr::perform_fit(fit_input)
	
	if (plot_fit_stages) plot_fit_2d(fit_output, spec_ord, always_show_start=TRUE, "Unfixed: m0, r2, sc, omega0")
	
	refit_thresh <- c(1e-3, 1e-2)[spec_ord]*2
	
	if (any(omega0_bound_distance(fit_output) < refit_thresh)) {
	
		#print("refitting")
		param_list <- fit_output[c("fit_list", "group_list", "comb_list")]
		names(param_list)[1] <- "start_list"
		fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list_to_arg_list(param_list)))
		if (positive_only) {
			fit_input$lower_list$m0[] <- 0
		}
		fit_input <- limit_omega0_by_r2(fit_input)
		fit_input$lower_list$r2[] <- r2_bounds[1]
		fit_input$upper_list$r2[] <- r2_bounds[2]
		param_values(fit_input$lower_list, coupling_idx) <- sc_bounds[1]
		param_values(fit_input$upper_list, coupling_idx) <- sc_bounds[2]
		fit_output <- fitnmr::perform_fit(fit_input)
		
		if (plot_fit_stages) plot_fit_2d(fit_output, spec_ord, always_show_start=TRUE, "Unfixed: m0, r2, sc, omega0 (Refit 1)")
	}
	
	if (any(omega0_bound_distance(fit_output) < refit_thresh)) {
	
		#print("refitting2")
		param_list <- fit_output[c("fit_list", "group_list", "comb_list")]
		names(param_list)[1] <- "start_list"
		fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list_to_arg_list(param_list)))
		if (positive_only) {
			fit_input$lower_list$m0[] <- 0
		}
		fit_input$lower_list$r2[] <- r2_bounds[1]
		fit_input$upper_list$r2[] <- r2_bounds[2]
		param_values(fit_input$lower_list, coupling_idx) <- sc_bounds[1]
		param_values(fit_input$upper_list, coupling_idx) <- sc_bounds[2]
		fit_input <- limit_omega0_by_r2(fit_input)
		fit_output <- fitnmr::perform_fit(fit_input)
		
		if (plot_fit_stages) plot_fit_2d(fit_output, spec_ord, always_show_start=TRUE, "Unfixed: m0, r2, sc, omega0 (Refit 2)")
	}
	
	if (any(omega0_bound_distance(fit_output) < refit_thresh)) {
	
		#print("refitting")
		param_list <- fit_output[c("fit_list", "group_list", "comb_list")]
		names(param_list)[1] <- "start_list"
		fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list_to_arg_list(param_list)))
		if (positive_only) {
			fit_input$lower_list$m0[] <- 0
		}
		fit_input$lower_list$r2[] <- r2_bounds[1]
		fit_input$upper_list$r2[] <- r2_bounds[2]
		param_values(fit_input$lower_list, coupling_idx) <- sc_bounds[1]
		param_values(fit_input$upper_list, coupling_idx) <- sc_bounds[2]
		fit_input <- limit_omega0_by_r2(fit_input)
		fit_output <- fitnmr::perform_fit(fit_input)
		
		if (plot_fit_stages) plot_fit_2d(fit_output, spec_ord, always_show_start=TRUE, "Unfixed: m0, r2, sc, omega0 (Refit 3)")
	}
	
	fit_output
}

#' Make a parameter list for a set of spectra and chemical shifts
#'
#' @export
make_param_list <- function(spec_list, cs_mat, fit_prev=NULL, r2_start=5, m0_start=1, sc_start=NULL, same_r2=FALSE, same_coupling=FALSE) {

	num_spec <- length(spec_list)
	
	spec_dim <- sapply(spec_list, function(x) ncol(x$fheader))
	num_dim <- spec_dim[1]
	
	stopifnot(spec_dim == num_dim, ncol(cs_mat) == num_dim)

	r2_group_offset <- 0
	m0_group_offset <- 0
	comb_idx_offset <- 0
	comb_group_offset <- 0
	
	if (!is.null(fit_prev)) {
	
		r2_group_offset <- max(fit_prev[["group_list"]][["r2"]])
		m0_group_offset <- max(fit_prev[["group_list"]][["m0"]])
		comb_idx_offset <- length(fit_prev[["group_list"]][["omega0_comb"]])
		comb_group_offset <- max(fit_prev[["group_list"]][["omega0_comb"]], 0)
	}
	
	if (is.null(sc_start)) {
	
		num_peaks <- 1
		omega0_start <- t(cs_mat)
		omega0_group <- rep(NA_integer_, num_dim)
		r2_group <- rep(NA_integer_, num_dim)
		m0_group <- NA_integer_
		omega0_comb <- vector("list", 1)
		omega0_comb_start <- numeric()
		omega0_comb_group <- integer()
	
	} else {
	
		stopifnot(length(sc_start) == num_dim)
		num_peaks <- 2^sum(!is.na(sc_start))
		omega0_start <- NA_real_
		omega0_group <- 0L
		m0_start <- rep(m0_start, each=num_peaks)
		r2_vec <- seq_len(num_dim*nrow(cs_mat))
		if (same_r2) {
			r2_vec <- seq_len(num_dim)
		}
		r2_group <- r2_group_offset+matrix(r2_vec, nrow=num_dim, ncol=nrow(cs_mat))[rep(seq_len(num_dim), num_peaks),,drop=FALSE]
		m0_group <- m0_group_offset+rep(seq_len(nrow(cs_mat)*length(spec_list)), each=num_peaks)
		
		comb_spacing <- num_dim+sum(!is.na(sc_start))
		comb_grid_params <- c(lapply(seq_along(sc_start), function(i) if (is.na(sc_start[i])) 0 else c(-0.5, 0.5)/spec_list[[1]]$fheader["OBS",i]), list(seq(comb_idx_offset, by=comb_spacing, length.out=nrow(cs_mat))))
		comb_grid <- do.call(expand.grid, comb_grid_params)
		scalar_coupling_coef_mat <- t(as.matrix(comb_grid[,-ncol(comb_grid),drop=FALSE]))
		omega0_idx_mat <- outer(seq_len(num_dim), comb_grid[,ncol(comb_grid)], "+")
		scalar_coupling_idx_vec <- rep(NA_integer_, num_dim)
		scalar_coupling_idx_vec[!is.na(sc_start)] <- num_dim+seq_len(sum(!is.na(sc_start)))
		scalar_coupling_idx_mat <- outer(scalar_coupling_idx_vec, comb_grid[,ncol(comb_grid)], "+")
		
		omega0_comb <- lapply(seq_along(scalar_coupling_coef_mat), function(i) {
			if (is.na(scalar_coupling_idx_mat[i])) {
				data.frame(omega0_idx_mat[i], 1, fix.empty.names=FALSE)
			} else {
				data.frame(c(omega0_idx_mat[i], scalar_coupling_idx_mat[i]), c(1, scalar_coupling_coef_mat[i]), fix.empty.names=FALSE)
			}
		})
		
		if (is.null(colnames(cs_mat))) {
			colnames(cs_mat) <- paste("omega0_ppm_", seq_len(ncol(cs_mat)), sep="")
		}
		sc_mat <- matrix(rep(sc_start[!is.na(sc_start)], each=nrow(cs_mat)), ncol=nrow(cs_mat), byrow=TRUE)
		if (nrow(sc_mat)) {
			rownames(sc_mat) <- paste("sc_hz_", which(!is.na(sc_start)), sep="")
		}
		omega0_comb_start <- rbind(t(cs_mat), sc_mat)
		omega0_comb_group <- comb_group_offset+array(seq_along(omega0_comb_start), dim=dim(omega0_comb_start))
		if (same_coupling) {
			omega0_comb_group[-seq_len(ncol(cs_mat)),] <- omega0_comb_group[-seq_len(ncol(cs_mat)),1]
		}
		rownames(omega0_comb_group) <- rownames(omega0_comb_start)
	}
	
	param_list <- list(
		start_list=list(
			omega0=array(omega0_start, dim=c(num_dim, nrow(cs_mat)*num_peaks, num_spec)),
			r2=array(r2_start, dim=c(num_dim, nrow(cs_mat)*num_peaks, num_spec)),
			m0=array(m0_start, dim=c(nrow(cs_mat)*num_peaks, num_spec)),
			omega0_comb=omega0_comb_start
		),
		group_list=list(
			omega0=array(omega0_group, dim=c(num_dim, nrow(cs_mat)*num_peaks, num_spec)),
			r2=array(r2_group, dim=c(num_dim, nrow(cs_mat)*num_peaks, num_spec)),
			m0=array(m0_group, dim=c(nrow(cs_mat)*num_peaks, num_spec)),
			omega0_comb=omega0_comb_group
		),
		comb_list=list(
			omega0=array(omega0_comb, dim=c(num_dim, nrow(cs_mat)*num_peaks, num_spec))
		)
	)
	
	if (!is.null(fit_prev)) {
		
		if (same_r2) {
			param_list[["start_list"]][["r2"]][] <- fit_prev[["fit_list"]][["r2"]][,1,]
			param_list[["group_list"]][["r2"]][] <- fit_prev[["group_list"]][["r2"]][,1,]
		}
		
		if (same_coupling) {
			new_coupling_idx <- coupling_param_idx(param_list, comb_idx_offset)
			prev_coupling_idx <- coupling_param_idx(fit_prev)
			if (length(prev_coupling_idx[["omega0_comb"]])) {
				param_values(param_list[["start_list"]], new_coupling_idx) <- fit_prev[["fit_list"]][["omega0_comb"]][which(prev_coupling_idx[["omega0_comb"]][,1])]
				param_values(param_list[["group_list"]], new_coupling_idx) <- fit_prev[["group_list"]][["omega0_comb"]][which(prev_coupling_idx[["omega0_comb"]][,1])]
			}
		}
		
		param_list[["start_list"]][["omega0"]] <- abind::abind(fit_prev[["fit_list"]][["omega0"]], param_list[["start_list"]][["omega0"]], rev.along=2)
		param_list[["start_list"]][["r2"]] <- abind::abind(fit_prev[["fit_list"]][["r2"]], param_list[["start_list"]][["r2"]], rev.along=2)
		param_list[["start_list"]][["m0"]] <- abind::abind(fit_prev[["fit_list"]][["m0"]], param_list[["start_list"]][["m0"]], rev.along=2)
		param_list[["start_list"]][["omega0_comb"]] <- abind::abind(fit_prev[["fit_list"]][["omega0_comb"]], param_list[["start_list"]][["omega0_comb"]])
		
		param_list[["group_list"]][["omega0"]] <- abind::abind(fit_prev[["group_list"]][["omega0"]], param_list[["group_list"]][["omega0"]], rev.along=2)
		param_list[["group_list"]][["r2"]] <- abind::abind(fit_prev[["group_list"]][["r2"]], param_list[["group_list"]][["r2"]], rev.along=2)
		param_list[["group_list"]][["m0"]] <- abind::abind(fit_prev[["group_list"]][["m0"]], param_list[["group_list"]][["m0"]], rev.along=2)
		param_list[["group_list"]][["omega0_comb"]] <- abind::abind(fit_prev[["group_list"]][["omega0_comb"]], param_list[["group_list"]][["omega0_comb"]])
		
		# abind::abind doesn't work for arrays of lists...
		#param_list[["comb"]][["omega0"]] <- abind::(list(fit_prev[["comb_list"]][["omega0"]], param_list[["comb"]][["omega0"]]), rev.along=2)
		idx_array_1 <- array(seq_along(fit_prev[["comb_list"]][["omega0"]]), dim(fit_prev[["comb_list"]][["omega0"]]))
		idx_array_2 <- array(seq_along(param_list[["comb_list"]][["omega0"]])+length(fit_prev[["comb_list"]][["omega0"]]), dim(param_list[["comb_list"]][["omega0"]]))
		idx_array_bind <- abind::abind(idx_array_1, idx_array_2, rev.along=2)
		param_list[["comb_list"]][["omega0"]] <- array(c(fit_prev[["comb_list"]][["omega0"]], param_list[["comb_list"]][["omega0"]])[as.vector(idx_array_bind)], dim=dim(fit_prev[["comb_list"]][["omega0"]])+c(0, dim(param_list[["comb_list"]][["omega0"]])[2], 0))
	}
	
	param_list
}

#' Get a list of logical arrays indicating which parameters correspond to peak positions
#'
#' @export
omega0_param_idx <- function(param_list, dims=seq_len(dim(param_list[["group_list"]][["omega0"]])[1]), peaks=seq_len(dim(param_list[["group_list"]][["omega0"]])[2]), specs=seq_len(dim(param_list[["group_list"]][["omega0"]])[3])) {
	
	idx_list <- lapply(param_list[["group_list"]], function(x) if (is.array(x)) array(FALSE, dim(x)) else logical(length(x)))
	names(idx_list) <- names(param_list[["group_list"]])
	
	subset_idx <- idx_list[["omega0"]]
	subset_idx[dims, peaks, specs] <- TRUE
	
	idx_list[["omega0"]][subset_idx] <- sapply(param_list[["comb_list"]][["omega0"]][subset_idx], is.null)
	
	for (i in which(!idx_list[["omega0"]] & subset_idx)) {
		unity_idx <- param_list[["comb_list"]][["omega0"]][[i]][,2] == 1
		idx_list[["omega0_comb"]][ param_list[["comb_list"]][["omega0"]][[i]][unity_idx,1] ] <- TRUE
	}
	
	idx_list
}

#' Get a list of logical arrays indicating which parameters correspond to scalar couplings
#'
#' @export
coupling_param_idx <- function(param_list, comb_idx_offset=0) {
	
	idx_list <- lapply(param_list[["group_list"]], function(x) if (is.array(x)) array(FALSE, dim(x)) else logical(length(x)))
	names(idx_list) <- names(param_list[["group_list"]])
	
	not_null_idx <- which(!sapply(param_list[["comb_list"]][["omega0"]], is.null))
	
	for (i in which(!idx_list[["omega0"]] )) {
		not_unity_idx <- param_list[["comb_list"]][["omega0"]][[i]][,2] != 1
		idx_list[["omega0_comb"]][ param_list[["comb_list"]][["omega0"]][[i]][not_unity_idx,1]-comb_idx_offset ] <- TRUE
	}
	
	idx_list
}

#' Get the first index in the omega0 array corresponding to each TRUE value in omega0_idx
omega0_comb_source_idx <- function(param_list, omega0_idx=omega0_param_idx(param_list)) {

	source_idx <- param_values(param_list[["group_list"]], omega0_idx)
	source_idx[] <- NA
	
	# first get the index of omega0 values not generated through linear combination
	comb_null_idx <- sapply(param_list[["comb_list"]][["omega0"]], is.null)
	
	omega0_source_idx <- which(omega0_idx[["omega0"]])
	source_idx[seq_along(omega0_source_idx)] <- omega0_source_idx
	
	for (i in which(!comb_null_idx)) {
		unity_idx <- param_list[["comb_list"]][["omega0"]][[i]][,2] == 1
		comb_idx <- param_list[["comb_list"]][["omega0"]][[i]][unity_idx,1]
		if (omega0_idx[["omega0_comb"]][comb_idx]) {
			idx <- match(comb_idx, which(omega0_idx[["omega0_comb"]]))+length(omega0_source_idx)
			source_idx[idx] <- i
		}
	}
	
	source_idx
}

#' Get/set a subset of fitting parameters specified by a list of logical vectors
#'
#' @param params a list of fitting parameters
#' @param idx_list a list with the same structure as \code{params} but with logical
#' vectors indicating which values should be return/set
#' @return a vector of parameter values or the modified parameter list
#' @export
param_values <- function(params, idx_list) {

	params_subset <- lapply(seq_along(params), function(i) params[[i]][idx_list[[i]]])
	
	do.call(c, params_subset)
}

#' @rdname param_values
#' @export
"param_values<-" <- function(params, idx_list, value) {

	stopifnot(sapply(params, length)[names(idx_list)] == sapply(idx_list, length))
	value <- rep(value, sum(sapply(idx_list, sum)))

	idx_offset <- 0
	
	for (i in names(idx_list)) {
		idx <- which(idx_list[[i]])
		params[[i]][idx_offset+idx] <- value[idx_offset+seq_along(idx)]
		idx_offset <- idx_offset+length(idx)
	}
	
	params
}

#' Convert a list of parameters for use with make_fit_input
#'
#' @export
param_list_to_arg_list <- function(param_list) {

	for (i in seq_along(param_list)) {
		names(param_list[[i]]) <- paste(names(param_list[[i]]), sub("_list", "", names(param_list)[i]), sep="_")
	}
	names(param_list) <- NULL

	do.call(c, param_list)
}

#' Fit a cluster nearby peaks starting from a seed table of chemical shifts
#'
#' @export
fit_peak_cluster <- function(spec_list, cs_start, spec_ord, f_alpha_thresh=0.001, r2_start=5, r2_bounds=c(0.5, 20), sc_start=NULL, sc_bounds=c(0, Inf), plot_main_prefix=NULL, peak_num_offset=0, plot_fit_stages=FALSE) {

	cs_new <- cs_start
	fit_output <- NULL
	fit_residuals <- NULL
	footprint <- NULL
	num_params <- 0
	f_alpha <- 0
	f_p_trace <- numeric()
	
	while (f_alpha < f_alpha_thresh) {
	
		# perform trial fit
		trial_fit_output <- fit_peaks(spec_list, cs_new, fit_output, spec_ord=spec_ord, sc_start=sc_start, sc_bounds=sc_bounds, positive_only=TRUE, plot_fit_stages=plot_fit_stages)
		
		# get new data from trial fit
		trial_input_spec_int <- fitnmr::get_spec_int(trial_fit_output, "input")
		trial_fit_spec_int <- fitnmr::get_spec_int(trial_fit_output, "fit")
		if (length(spec_list) > 1) {
			trial_fit_residuals <- lapply(seq_along(spec_list), function(i) trial_input_spec_int[[i]]-trial_fit_spec_int[[i]])
		} else {
			trial_fit_residuals <- trial_input_spec_int[[1]]-trial_fit_spec_int[[1]]
		}
		trial_footprint <- fit_footprint(trial_fit_output)
		
		# terminate search if any peak had zero volume in every spectrum
		if (any(rowSums(trial_fit_output$fit_list$m0) == 0)) {
			cat(" Terminating search because fit produced zero volume", sep="\n")
			if (is.null(fit_output)) {
				if (length(spec_list) > 1) {
					fit_output <- trial_input_spec_int
				} else {
					fit_output <- trial_input_spec_int[[1]]
				}
			}
			break
		}
		
		# determine footprint of the peaks from the union of previous and current fits
		if (is.null(footprint)) {
			common_footprint <- trial_footprint
		} else {
			if (length(spec_list) > 1) {
				for (i in seq_along(common_footprint)) {
					common_rows <- intersect(dimnames(fit_spec_int[[i]])[[1]], dimnames(trial_fit_spec_int[[i]])[[1]])
					common_cols <- intersect(dimnames(fit_spec_int[[i]])[[2]], dimnames(trial_fit_spec_int[[i]])[[2]])
					common_footprint[[i]] <- footprint[[i]][common_rows,common_cols] | trial_footprint[[i]][common_rows,common_cols]
				}
			} else {
				common_rows <- intersect(dimnames(fit_spec_int[[1]])[[1]], dimnames(trial_fit_spec_int[[1]])[[1]])
				common_cols <- intersect(dimnames(fit_spec_int[[1]])[[2]], dimnames(trial_fit_spec_int[[1]])[[2]])
				common_footprint <- footprint[common_rows,common_cols] | trial_footprint[common_rows,common_cols]
			}
		}
		
		if (length(spec_list) > 1) {
			common_footprint_idx <- lapply(common_footprint, which, arr.ind=TRUE)
			common_footprint_idx <- lapply(seq_along(common_footprint), function(i) {
				cbind(dimnames(common_footprint[[i]])[[1]][common_footprint_idx[[i]][,1]], dimnames(common_footprint[[i]])[[2]][common_footprint_idx[[i]][,2]])
			})
		} else {
			common_footprint_idx <- which(common_footprint, arr.ind=TRUE)
			common_footprint_idx <- cbind(dimnames(common_footprint)[[1]][common_footprint_idx[,1]], dimnames(common_footprint)[[2]][common_footprint_idx[,2]])
		}
		
		if (is.null(fit_residuals)) {
		
			# initial model is having no peak (residuals = input)
			if (length(spec_list) > 1) {
				fit_residuals <- trial_input_spec_int
			} else {
				fit_residuals <- trial_input_spec_int[[1]]
			}
			
		} else {
		
			# remove footprint locations where residuals are undefined
			if (length(spec_list) > 1) {
				for (i in seq_along(common_footprint_idx)) {
					common_footprint_idx[[i]] <- common_footprint_idx[[i]][!is.na(fit_residuals[[i]][common_footprint_idx[[i]]]) & !is.na(trial_fit_residuals[[i]][common_footprint_idx[[i]]]),,drop=FALSE]
				}
			} else {
				common_footprint_idx <- common_footprint_idx[!is.na(fit_residuals[common_footprint_idx]) & !is.na(trial_fit_residuals[common_footprint_idx]),,drop=FALSE]
			}
		}
		
		# calculate residual sum of squares for both models
		if (length(spec_list) > 1) {
			rss <- sum(sapply(seq_along(spec_list), function(i) sum(fit_residuals[[i]][common_footprint_idx[[i]]]^2)))
			trial_rss <- sum(sapply(seq_along(spec_list), function(i) sum(trial_fit_residuals[[i]][common_footprint_idx[[i]]]^2)))
		} else {
			rss <- sum(fit_residuals[common_footprint_idx]^2)
			trial_rss <- sum(trial_fit_residuals[common_footprint_idx]^2)
		}
		
		# calculate new degrees of freedom
		trial_num_params <- sum(sapply(trial_fit_output$group_list, function(x) length(unique(x[x!=0]))))
		
		# scale the number of points to account for zero-filling
		if (length(spec_list) > 1) {
			num_pts <- sum(sapply(seq_along(spec_list), function(i) {
				nrow(common_footprint_idx[[i]])*prod(spec_list[[i]]$fheader["TDSIZE",]/(spec_list[[i]]$fheader["FTSIZE",]/2))
			}))
		} else {
			num_pts <- nrow(common_footprint_idx)*prod(spec_list[[1]]$fheader["TDSIZE",]/(spec_list[[1]]$fheader["FTSIZE",]/2))
		}
		
		# calculate value of F statistic and corresponding P-value
		f_val <- ((rss-trial_rss)/(trial_num_params-num_params))/(trial_rss/(num_pts-trial_num_params))
		f_alpha <- 1-stats::pf(f_val, trial_num_params-num_params, num_pts-trial_num_params)
		f_alpha <- -expm1(stats::pf(f_val, trial_num_params-num_params, num_pts-trial_num_params, log.p=TRUE))
		f_p_trace <- c(f_p_trace, f_alpha)
		
		cat(sprintf(" %2i -> %2i fit parameters: F = %0.1f (p = %g)", num_params, trial_num_params, f_val, f_alpha), sep="\n")
		
		if (f_alpha < f_alpha_thresh)  {
		
			fit_output <- trial_fit_output
			input_spec_int <- trial_input_spec_int
			fit_spec_int <- trial_fit_spec_int
			fit_residuals <- trial_fit_residuals
			footprint <- trial_footprint
			num_params <- trial_num_params
			
			if (length(spec_list) > 1) {
				norm_max_resid <- sapply(seq_along(spec_list), function(i) max(abs(fit_residuals[[i]][footprint[[i]]]))/max(abs(spec_list[[i]]$int)))
				max_spec_i <- which.max(norm_max_resid)
				resid_max_idx <- which(fit_residuals[[max_spec_i]] == max(fit_residuals[[max_spec_i]][footprint[[max_spec_i]]]), arr.ind=TRUE)
				cs_new <- matrix(as.numeric(c(dimnames(footprint[[max_spec_i]])[[1]][resid_max_idx[1]], dimnames(footprint[[max_spec_i]])[[2]][resid_max_idx[2]])), nrow=1)[,spec_ord,drop=FALSE]
			} else {
				resid_max_idx <- which(fit_residuals == max(fit_residuals[footprint]), arr.ind=TRUE)
				cs_new <- matrix(as.numeric(c(dimnames(footprint)[[1]][resid_max_idx[1]], dimnames(footprint)[[2]][resid_max_idx[2]])), nrow=1)[,spec_ord,drop=FALSE]
			}
		
		} else {
		
			cat(sprintf(" Terminating search because F-test p-value < %g", f_alpha_thresh), sep="\n")
			if (is.null(fit_output)) {
				if (length(spec_list) > 1) {
					fit_output <- trial_fit_spec_int
				} else {
					fit_output <- trial_fit_spec_int[[1]]
				}
			}
		}
	}
	
	if (is.list(fit_output)) {
		fit_output$fit_fptrace <- f_p_trace
	}
	
	if (!is.null(plot_main_prefix) && "fit_list" %in% names(fit_output))  {
		
		omega_comb_ids <- apply(array(sapply(fit_output$comb_list$omega0, "[[", 1, 1), dim=dim(fit_output$comb_list$omega0)), c(2:3), paste, collapse="-")
		omega_comb_ids_unique <- unique(as.vector(omega_comb_ids))
		
		for (i in seq_along(spec_list)) {
		
			zlim <- range(input_spec_int[[i]], na.rm=TRUE)
	
			fitnmr::contour_pipe(aperm(input_spec_int[[i]], spec_ord), zlim=zlim, col_pos="black", col_neg="gray")
			main <- paste(plot_main_prefix, "Spectrum", i)
			if (length(fit_output$fit_fptrace) > length(omega_comb_ids_unique)) {
				main <- paste(main, " (next p: ", sprintf("%.1e", fit_output$fit_fptrace[length(omega_comb_ids_unique)+1]), ")", sep="")
			} else {
				main <- paste(main, " (next 0 volume)", sep="")
			}
			graphics::title(main)
			#graphics::title(paste("Cluster", j))
			fitnmr::contour_pipe(aperm(fit_spec_int[[i]], spec_ord), zlim=zlim, col_pos="red", col_neg="pink", add=TRUE)
			#fitnmr::contour_pipe(aperm(trial_fit_spec_int[[1]], spec_ord), zlim=zlim, col_pos="purple", col_neg="plum", add=TRUE)
	
			graphics::points(cs_start, col="green")
	
			#graphics::points(common_footprint_idx[,spec_ord], col="gray")
	
			#graphics::rect(fit_output$upper_list$omega0[spec_ord[1],,1], fit_output$upper_list$omega0[spec_ord[2],,1], fit_output$lower_list$omega0[spec_ord[1],,1], fit_output$lower_list$omega0[spec_ord[2],,1], border="gray")
	
			lab_coord <- matrix(nrow=length(omega_comb_ids_unique), ncol=2)
			for (j in seq_along(omega_comb_ids_unique)) {
				id <- omega_comb_ids_unique[j]
				lab_coord[j,] <- c(
					mean(fit_output$fit_list$omega0[spec_ord[1],id==omega_comb_ids[,i],i]),
					max(fit_output$fit_list$omega0[spec_ord[2],id==omega_comb_ids[,i],i])
				)
				graphics::points(t(fit_output$fit_list$omega0[spec_ord,id==omega_comb_ids[,i],i]), type="l", col="blue")
			}
			frac <- tapply(fit_output$fit_list$m0[,i], match(omega_comb_ids[,i], omega_comb_ids_unique), sum)
			graphics::points(t(fit_output$fit_list$omega0[spec_ord,,i]), pch=16, col="blue")
			graphics::text(lab_coord, labels=sprintf("%s: %.0f%%\np: %.0e", seq_along(frac)+peak_num_offset, 100*frac/sum(frac), fit_output$fit_fptrace[seq_along(frac)]), pos=1, cex=0.6)
	
			#other_clust_idx <- peak_tab_list[[i]][,"TYPE"] == 1 & peak_tab_list[[i]][,"CLUSTID"] != j
			#graphics::points(peak_tab_list[[i]][other_clust_idx,paste(hn_name_mat[i,], "_PPM", sep=""),drop=FALSE], col="purple")
			#graphics::text(peak_tab_list[[i]][other_clust_idx,paste(hn_name_mat[i,], "_PPM", sep=""),drop=FALSE], labels=peak_tab_list[[i]][other_clust_idx,"CLUSTID"], pos=1, cex=0.6, col="purple")
		}
	}
	
	fit_output
}

#' Iterative Peak Fitting
#'
#' Iteratively fit peaks for whole spectra
#'
#' This function uses an iterative algorithm to fit all the peaks in a given list of spectra, assuming identical peak positions and shapes across the spectra. Each iteration starts by identifying the maximum value across all spectra, using that position to fit a cluster of overlapping peaks with the \code{\link{fit_peak_cluster}} function. After the cluster of peaks is fit, the modeled intensity is subtracted from all spectra and another iteration is performed. Iterations are terminated if \code{max_iter} is reached or the residual intensity in all spectra falls below \code{noise_sigma*noise_cutoff}.
#'
#' This function currently only supports fitting of 2D spectra, but will be generalized to work with spectra of any dimensionality in the near future. To reduce the number of false positives/negatives, the most important parameters to adjust are \code{noise_cutoff}, \code{f_alpha}, and \code{iter_max}. If \code{iter_max} is reached before all peaks have been identified, then you can call this function again, setting the \code{fit_list} parameter to the return value of the previous invocation. In that case, \code{iter_max} new iterations will be performed and appended to \code{fit_list}.
#'
#' For visualizing the iterative algorithm as it progresses, you can enable either the \code{plot_fit} or \code{plot_fit_stages} parameters.
#'
#' @param spectra list of spectrum objects read by \code{\link{read_nmrpipe}}.
#' @param noise_sigma numeric vector of noise levels associated with each spectrum. If \code{NULL}, it is calculated with \code{\link{noise_estimate}}.
#' @param noise_cutoff numeric value multiplied by \code{noise_sigma} to determine cutoffs for each spectrum. Peak fitting will terminate if the maximum residuals for all spectra fall below these cutoffs.
#' @param f_alpha numeric value giving a F-test p-value threshold above which a peak will not be accepted.
#' @param iter_max integer maximum number of iterations to apply.
#' @param omega0_plus numeric vector giving the window size (ppm plus or minus the starting \code{omega0} values) around which to use points from the spectra for fitting.
#' @param r2_start numeric vector giving the starting \code{r2} value(s) for the fit (in Hz).
#' @param r2_bounds numeric vector of length two giving the lower and upper bounds for  \code{r2}.
#' @param sc_start numeric vector giving the starting scalar coupling values for doublets. It should be the same length as the number of dimensions in the spectrum. Set the value to \code{NA} for a given dimension to make it a singlet.
#' @param sc_bounds numeric vector of length two giving the lower and upper bounds for  scalar couplings.
#' @param fit_list list of previous fits to which the new fits should be appended.
#' @param plot_fit logical indicating whether produce a fit cluster plot for each iteration.
#' @param plot_fit_stages logical indicating whether to plot each stage of fitting within the iterations.
#' @return List of fit objects returned by \code{\link{fit_peak_cluster}}, one for each iteration. They are appended to \code{fit_list} if supplied.
#'
#' @export
fit_peak_iter <- function(spectra, noise_sigma=NULL, noise_cutoff=15, f_alpha=1e-3, iter_max=100, omega0_plus=c(0.075, 0.75), r2_start=5, r2_bounds=c(0.5, 20), sc_start=c(6, NA), sc_bounds=c(2, 12), fit_list=list(), plot_fit=FALSE, plot_fit_stages=FALSE) {

	if (is.null(noise_sigma)) {
		noise_sigma <- sapply(spectra, function(x) fitnmr::noise_estimate(x$int, plot_distributions=FALSE))["sigma",]
	}

	spec_sub_list <- spectra
	
	fit_num <- 1
	peak_num_offset <- 0
	
	# subtract out the intensities from any inputs given
	for (iter in seq_along(fit_list)) {
	
		fit_spec_int <- fitnmr::get_spec_int(fit_list[[iter]], "fit")
		
		for (i in seq_along(spectra)) {
			idx_1 <- dimnames(fit_spec_int[[i]])[[1]]
			idx_2 <- dimnames(fit_spec_int[[i]])[[2]]
			not_na_idx <- !is.na(fit_spec_int[[i]])
			
			spec_sub_list[[i]]$int[idx_1,idx_2][not_na_idx] <- (spec_sub_list[[i]]$int[idx_1,idx_2]-fit_spec_int[[i]])[not_na_idx]
		}
		
		fit_num <- fit_num + length(unique(fit_list[[iter]]$group_list$r2[1,,]))
		peak_num_offset <- peak_num_offset + length(unique(fit_list[[iter]]$group_list$m0[,1]))
	}
	
	if (is.character(plot_fit)) {
	
		grDevices::pdf(plot_fit)
		graphics::par(mar=c(2.9, 2.9, 1.5, 1), mgp=c(1.7, 0.6, 0))
		plot_fit <- TRUE
		on.exit(grDevices::dev.off())
	
	} else {
	
		plot_main_prefix <- NULL
	}
	
	for (iter in seq_len(iter_max)) {
	
		# check if 
		spec_max_val <- sapply(spec_sub_list, function(x) max(x$int))
		if (all(spec_max_val < noise_sigma*noise_cutoff)) {
			break
		}
	
		cat(paste("Fit iteration ", iter, ":", sep=""), sep="\n")
		
		spec_max_idx <- which.max(spec_max_val)
		max_idx <- which(spec_sub_list[[spec_max_idx]]$int == spec_max_val[spec_max_idx], arr.ind=TRUE)[1,,drop=FALSE]
		max_cs <- matrix(sapply(seq_along(max_idx), function(i) spec_sub_list[[spec_max_idx]]$ppm[[i]][max_idx[i]]), nrow=1)
		
		if (plot_fit) {
			plot_main_prefix <- paste("Fit", fit_num)
		}
		fit_output <- fitnmr::fit_peak_cluster(spec_sub_list, max_cs, spec_ord=1:2, f_alpha_thresh=f_alpha, r2_start=r2_start, r2_bounds=r2_bounds, sc_start=sc_start, sc_bounds=sc_bounds, plot_main_prefix=plot_main_prefix, peak_num_offset=peak_num_offset, plot_fit_stages=plot_fit_stages)
		
		if ("fit_list" %in% names(fit_output)) {

			fit_list[[length(fit_list)+1]] <- fit_output
			peak_num_offset <- peak_num_offset+dim(fit_output$fit_list$omega0)[2]/2^sum(!is.na(sc_start))
	
			fit_spec_int <- fitnmr::get_spec_int(fit_output, "fit")
	
			for (i in seq_along(spectra)) {
				idx_1 <- dimnames(fit_spec_int[[i]])[[1]]
				idx_2 <- dimnames(fit_spec_int[[i]])[[2]]
				not_na_idx <- !is.na(fit_spec_int[[i]])
				
				spec_sub_list[[i]]$int[idx_1,idx_2][not_na_idx] <- (spec_sub_list[[i]]$int[idx_1,idx_2]-fit_spec_int[[i]])[not_na_idx]
			}
			
			fit_num <- fit_num+1
	
		} else if (is.matrix(fit_output[[1]]) || is.matrix(fit_output)) {
	
			if (is.matrix(fit_output)) {
				fit_output <- list(fit_output)
			}
			
			for (i in seq_along(spectra)) {
				idx_1 <- dimnames(fit_output[[i]])[[1]]
				idx_2 <- dimnames(fit_output[[i]])[[2]]
				not_na_idx <- !is.na(fit_output[[i]])
		
				spec_sub_list[[i]]$int[idx_1,idx_2][not_na_idx] <- (spec_sub_list[[i]]$int[idx_1,idx_2]-fit_output[[i]])[not_na_idx]
			}
		}
	}
	
	fit_list
}

#' Convert Fit to Data Frame
#'
#' Convert a parameter list into a peak data frame
#'
#' This function takes the input or (if present) output parameters from a fit and converts them into a data frame. A parameter list must contain three lists:
#' \describe{
#'   \item{\code{start_list} or \code{fit_list}}{input or output values of the respective fit parameters}
#'   \item{\code{group_list}}{group numbers for the fit parameters}
#'   \item{\code{comb_list}}{coefficients for deriving fit parameters from a linear combination other auxiliary parameters}
#' }
#'
#' This function currently assumes the fit parameters were generated by \code{\link{fit_peak_iter}}, \code{\link{fit_peak_cluster}}, or \code{\link{peak_df_to_param_list}}. These functions use a particular convention for \code{group_list} and \code{comb_list} to represent either singlets or doublets in each dimension of a 2D spectrum.
#'
#' This function can take either a single parameter list or a list of parameter lists. If the latter is given, then the results from all the parameter lists will be combined into a single table.
#'
#' @param param_list list of fit parameters (or a list of such lists)
#' @param spec_names character vector of spectrum names
#' @return A data frame with the following columns:
#' \describe{
#'   \item{\code{peak}}{peak number}
#'   \item{\code{fit}}{fit cluster number, with all peaks in the same cluster having the same \code{r2}, scalar couplings, and \code{m0}}
#'   \item{\code{f_pvalue}}{optional, p-value determined from F-test during iterative fitting}
#'   \item{\code{omega0_ppm_1}}{chemical shift of singlet/doublet center in the first dimension (ppm)}
#'   \item{\code{omega0_ppm_2}}{chemical shift of singlet/doublet center in the second dimension (ppm)}
#'   \item{\code{sc_hz_1}}{optional, scalar coupling of doublet in first dimension (Hz)}
#'   \item{\code{sc_hz_2}}{optional, scalar coupling of doublet in second dimension (Hz)}
#'   \item{\code{r2_hz_1}}{R2 in first dimension (Hz)}
#'   \item{\code{r2_hz_2}}{R2 in first dimension (Hz)}
#'   \item{...}{\code{m0} values for each spectrum}
#' }
#'
#' @export
param_list_to_peak_df <- function(param_list, spec_names=NULL) {

	if ("group_list" %in% names(param_list)) {
	
		list_name <- "fit_list"
		if (!list_name %in% names(param_list)) {
			list_name <- "start_list"
		}
		stopifnot(list_name %in% names(param_list))
		
		comb_omega0_null <- sapply(param_list$comb_list$omega0, is.null)
		stopifnot(all(comb_omega0_null) || all(!comb_omega0_null))
	
		if (all(comb_omega0_null)) {
		
			stop("No scalar couplings detected")
		
		} else {
	
			# determine which column of omega0_comb determines each omega0 value
			omega0_comb_idx <- array(sapply(param_list$comb_list$omega0, function(x) {
				if (is.null(x)) {
					NA
				} else {
					col_idx <- unique((x[,1]-1) %/% nrow(param_list[[list_name]]$omega0_comb) + 1)
					stopifnot(length(col_idx) == 1)
					col_idx
				}
			}), dim=dim(param_list$comb_list$omega0))
		
			stopifnot(apply(omega0_comb_idx, 2, function(x) all(x==x[1])))
			
			r2_mat <- param_list[[list_name]]$r2[,!duplicated(omega0_comb_idx[1,,1]),1]
			r2_mat <- t(matrix(r2_mat, nrow=dim(param_list[[list_name]]$r2)[1]))
			colnames(r2_mat) <- paste("r2_hz_", seq_len(ncol(r2_mat)), sep="")
			
			m0_mat <- simplify2array(tapply(seq_len(nrow(param_list[[list_name]]$m0)), omega0_comb_idx[1,,1], function(idx) {
				colSums(param_list[[list_name]]$m0[idx,,drop=FALSE])
			}))
			m0_mat <- t(matrix(m0_mat, nrow=ncol(param_list[[list_name]]$m0)))
			if (is.null(spec_names)) {
				colnames(m0_mat) <- names(param_list$spec_data)
			} else {
				colnames(m0_mat) <- spec_names
			}
			
			peak_df <- data.frame(
				"peak"=seq_len(ncol(param_list$group_list$omega0_comb)),
				"fit"=match(param_list$group_list$r2[1,!duplicated(omega0_comb_idx[1,,1]),1], unique(param_list$group_list$r2[1,!duplicated(omega0_comb_idx[1,,1]),1]))
			)
			
			if ("fit_fptrace" %in% names(param_list)) {
				peak_df <- cbind(peak_df, f_pvalue=param_list$fit_fptrace[seq_len(nrow(peak_df))])
			}
			
			peak_df <- cbind(peak_df, t(param_list[[list_name]]$omega0_comb), r2_mat, m0_mat)
			
			peak_df
		}
	
	} else {
	
		peak_df_list <- lapply(param_list, param_list_to_peak_df, spec_names=spec_names)
		
		cum_peaks <- cumsum(sapply(peak_df_list, function(x) length(unique(x[,"peak"]))))
		cum_fits <- cumsum(sapply(peak_df_list, function(x) length(unique(x[,"fit"]))))
		
		for (i in seq_along(peak_df_list)[-1]) {
			peak_df_list[[i]][,"peak"] <- peak_df_list[[i]][,"peak"] + cum_peaks[i-1]
			peak_df_list[[i]][,"fit"] <- peak_df_list[[i]][,"fit"] + cum_fits[i-1]
		}
		
		do.call(rbind, peak_df_list)
	}
}

#' Convert a peak data frame to a parameter list
#'
#' @export
peak_df_to_param_list <- function(peak_df, spectra) {

	omega0_cols <- grep("^omega0_ppm_", colnames(peak_df), value=TRUE)
	sc_cols <- grep("^sc_hz_", colnames(peak_df), value=TRUE)
	r2_cols <- grep("^r2_hz_", colnames(peak_df), value=TRUE)
	other_cols <- intersect(c("peak", "fit", "f_pvalue"), colnames(peak_df))
	m0_cols <- setdiff(colnames(peak_df), c(other_cols, omega0_cols, sc_cols, r2_cols))
	
	sc_means <- colMeans(peak_df[sc_cols])
	sc_start <- sc_means[match(paste("sc_hz_", seq_along(omega0_cols), sep=""), sc_cols)]

	param_list <- fitnmr::make_param_list(spectra, as.matrix(peak_df[,omega0_cols]), sc_start=sc_start)

	param_list$start_list$r2[] <- t(matrix(rep(as.vector(as.matrix(peak_df[,c("r2_hz_1", "r2_hz_2")])), each=2^sum(!is.na(sc_start))), ncol=2))
	param_list$start_list$m0[] <- rep(as.vector(as.matrix(peak_df[,m0_cols])/2^sum(!is.na(sc_start))), each=2^sum(!is.na(sc_start)))
	for (sc_col in intersect(c("sc_hz_1", "sc_hz_1"), colnames(peak_df))) {
		param_list$start_list$omega0_comb[sc_col,] <- peak_df[,sc_col]
	}

	fit_vec <- match(peak_df[,"fit"], unique(peak_df[,"fit"]))

	all_unique_group_mat <- matrix(seq_len(length(unique(peak_df[,"fit"]))*2), nrow=2)

	param_list$group_list$r2[] <- all_unique_group_mat[,rep(fit_vec,each=2^sum(!is.na(sc_start)))]
	param_list$group_list$m0[] <- rep(seq_len(length(param_list$group_list$m0)/2^sum(!is.na(sc_start))),each=2^sum(!is.na(sc_start)))
	#param_list$group_list$omega0_comb[] <- matrix(seq_len(length(unique(peak_df[,"fit"]))*3), nrow=3)[,fit_vec]
	param_list$group_list$omega0_comb[1:2,] <- matrix(seq_len(length(fit_vec)*2),nrow=2)
	for (i in seq_len(nrow(param_list$group_list$omega0_comb))[-(1:2)]) {
		param_list$group_list$omega0_comb[i,] <- fit_vec+length(fit_vec)*(i-1)
	}
	
	param_list
}

#' Convert a peak data frame to fit input
#'
#' @export
peak_df_to_fit_input <- function(peak_df, spectra, ...) {

	param_list <- peak_df_to_param_list(peak_df, spectra)
	arg_list <- param_list_to_arg_list(param_list)
	do.call(fitnmr::make_fit_input, c(list(spectra, ...), arg_list))
}

#' Plot Peaks from a Peak Table
#'
#' Plot fits for series of spectra with parameters from a peak data frame
#'
#' The raw spectral data is shown in black contours and the modeled peak intensity is shown in red. The centers of peaks are shown with semi-transparent blue dots, with the area of the dot proportional to the volume of the peak (\code{m0}). Blue lines connect peaks from modeled doublets. Singlets or doublets are labeled with the syntax <peak>:<fit>. If an F-test p-value column is present (\code{f_pvalue}), that will be given below the peak label.
#'
#' @param peak_df data frame with peak data as produced by \code{\link{param_list_to_peak_df}}.
#' @param spectra list of spectra corresponding to the volumes found in \code{peak_df}.
#' @param noise_sigma numeric vector of noise levels associated with each spectrum. If \code{NULL}, it is calculated with \code{\link{noise_estimate}}.
#' @param noise_cutoff numeric value used to calculate the lowest contour level according to \code{noise_sigma*noise_cutoff}.
#' @param cex numeric value by which to scale blue points and labels.
#' @param lwd numeric value giving width of contour lines.
#' @param label logical indicating whether to draw text labels and connecting lines.
#' @param p0 zero order phase for plotting modeled peaks.
#' @param p1 first order phase for plotting modeled peaks.
#' @param add logical indicating whether to suppress generation of a new plot and add to an existing plot.
#'
#' @export
plot_peak_df <- function(peak_df, spectra, noise_sigma=NULL, noise_cutoff=4, cex=0.2, lwd=0.25, label=TRUE, p0=NULL, p1=NULL, add=FALSE) {

	if (is.null(noise_sigma)) {
		noise_sigma <- sapply(spectra, function(x) fitnmr::noise_estimate(x$int, plot_distributions=FALSE))["sigma",]
	}

	param_list <- peak_df_to_param_list(peak_df, spectra)

	fit_input <- do.call(fitnmr::make_fit_input, c(list(spectra, omega0_plus=c(0.075, 0.75)*2), param_list_to_arg_list(param_list)))
	if (!is.null(p0)) {
		fit_input$start_list$p0[] <- p0
	}
	if (!is.null(p1)) {
		fit_input$start_list$p1[] <- p1
	}

	int_input <- fitnmr::get_spec_int(fit_input, "input")
	int_start <- fitnmr::get_spec_int(fit_input, "start")

	zlim_mat <- sapply(int_input, range, na.rm=TRUE)

	omega_comb_ids <- apply(array(sapply(fit_input$comb_list$omega0, "[[", 1, 1), dim=dim(fit_input$comb_list$omega0)), c(2:3), paste, collapse="-")
	omega_comb_ids_unique <- unique(as.vector(omega_comb_ids))

	for (spec_i in seq_along(int_input)) {

		low_frac <- noise_sigma[spec_i]/max(spectra[[spec_i]]$int)*noise_cutoff
		fitnmr::contour_pipe(spectra[[spec_i]]$int, zlim=zlim_mat[,spec_i], col_pos="black", col_neg="gray", low_frac=low_frac, lwd=lwd, add=add)
		fitnmr::contour_pipe(int_start[[spec_i]], zlim=zlim_mat[,spec_i], col_pos="red", col_neg="pink", low_frac=low_frac, lwd=lwd, add=TRUE)
		m0_vec <- fit_input$start_list$m0[,spec_i]
		graphics::points(t(fit_input$start_list$omega0[,,spec_i]), col=grDevices::rgb(0, 0, 1, 0.5), pch=16, cex=sqrt(m0_vec/max(m0_vec))*5*cex)
	
		if (label) {
			lab_coord <- matrix(nrow=length(omega_comb_ids_unique), ncol=2)
			for (j in seq_along(omega_comb_ids_unique)) {
				id <- omega_comb_ids_unique[j]
				lab_coord[j,] <- c(
					mean(fit_input$start_list$omega0[1,id==omega_comb_ids[,spec_i],spec_i]),
					max(fit_input$start_list$omega0[2,id==omega_comb_ids[,spec_i],spec_i])
				)
				graphics::points(t(fit_input$start_list$omega0[,id==omega_comb_ids[,spec_i],spec_i]), type="l", col="blue", lwd=cex*1.25)
			}
			#graphics::text(lab_coord, labels=sprintf("%s-%s\np: %.0e", peak_df$fit, peak_df$peak, peak_df$f_pvalue), pos=1, offset=0.1, cex=0.2)
			graphics::text(lab_coord, labels=sprintf("%s:%s", peak_df$peak, peak_df$fit), pos=1, offset=cex*0.5, cex=cex)
			if ("f_pvalue" %in% names(peak_df)) {
				graphics::text(lab_coord, labels=sprintf("%.0e", peak_df$f_pvalue), pos=1, offset=cex*1.25, cex=cex)
			}
		}
	}
}

#' Simulate an FID using the NMRPipe SimTimeND function
#'
#' @export
sim_time_nd <- function(peak_tab, fheader, rms=0, iseed=stats::runif(1,max=.Machine$integer.max), file_path=NULL, verbose=FALSE) {

	tab_path <- tempfile(fileext=".tab")
	write_nmrdraw_peak_tab(peak_tab, tab_path)
	on.exit(unlink(tab_path))

	tn_mat <- fheader[c("TDSIZE", "FTSIZE"),,drop=FALSE]
	rownames(tn_mat) <- c("T", "N")
	
	arg_mat <- rbind(
		tn_mat, 
		MODE=rep("Complex", ncol(fheader)),
		fheader[c("SW", "OBS", "CAR"),,drop=FALSE],
		LAB=colnames(fheader)
	)
	
	cmd_args <- as.vector(arg_mat)
	names(cmd_args) <- paste(c("x","y","z","a")[col(arg_mat)], rownames(arg_mat)[row(arg_mat)], sep="")
	
	cmd_args <- c(
		ndim=ncol(fheader), 
		"in"=tab_path, 
		rms=rms, 
		iseed=as.integer(iseed),
		scale="1.0", 
		nots="", 
		cmd_args, 
		aq2D="States"
	)
	
	if (!is.null(file_path)) {
		cmd_args <- c(cmd_args, out=file_path, ov="")
	}
	
	if (verbose) {
		cmd_args <- c(cmd_args, verb="")
	}
	
	cmd_args <- as.vector(rbind(paste("-", names(cmd_args), sep=""), cmd_args))
	cmd_args <- cmd_args[cmd_args != ""]
	
	if (verbose) {
		cat(paste(c("SimTimeND", cmd_args), collapse=" "), sep="\n")
	}

	system2("SimTimeND", cmd_args)
}

#' Process an FID with NMRPipe
#'
#' @export
nmr_pipe <- function(in_path, out_path, ndim=1, apod=NULL, sp=rbind(off=0.5, end=1.0, pow=1, c=0.5), zf=rbind(auto=""), ps=rbind(p0=0, p1=0, di=""), ext=NULL) {

	# more work is needed to support 3 and 4 dimensional data
	stopifnot(ndim >= 1, ndim <= 2)
	
	commands <- paste("nmrPipe -in ", in_path, sep="")
	
	for (i in seq_len(ndim)) {
		
		if (i == 2) {
			commands <- c(commands, paste(c("nmrPipe", "-fn", "TP"), collapse=" "))
		}
		
		if (!is.null(apod)) {
			commands <- c(commands, paste(c("nmrPipe", "-fn", "APOD", collapse_args(apod[,min(ncol(apod),i)])), collapse=" "))
		}
		
		if (!is.null(sp)) {
			commands <- c(commands, paste(c("nmrPipe", "-fn", "SP", collapse_args(sp[,min(ncol(sp),i)])), collapse=" "))
		}
		
		if (!is.null(zf)) {
			commands <- c(commands, paste(c("nmrPipe", "-fn", "ZF", collapse_args(zf[,min(ncol(zf),i)])), collapse=" "))
		}
		
		commands <- c(commands, paste(c("nmrPipe", "-fn", "FT", "-auto"), collapse=" "))
		
		if (!is.null(ps)) {
			commands <- c(commands, paste(c("nmrPipe", "-fn", "PS", collapse_args(ps[,min(ncol(ps),i)])), collapse=" "))
		}
		
		if (!is.null(ext)) {
			commands <- c(commands, paste(c("nmrPipe", "-fn", "EXT", collapse_args(ext[,min(ncol(ext),i)])), collapse=" "))
		}
		
		if (i == 2) {
			commands <- c(commands, paste(c("nmrPipe", "-fn", "TP"), collapse=" "))
		}
	}
	
	commands[length(commands)] <- paste(commands[length(commands)], " -ov -out ", out_path, sep="")
	
	system(paste(commands, collapse=" | "))
}

#' Estimate Noise
#'
#' Estimate properties of noise by fitting a Gaussian to a histogram of intensities
#'
#' This function estimates noise using iterative calculation of the mean and standard deviation, followed by fitting a Gaussian function to a histogram of the values within a range determined at the end of the iterations.
#'
#' The iterative algorithm first calculates the mean and standard deviation of the values in \code{x}. In the next iteration, only points within \code{thresh} times the standard deviation of the mean value are used for calculating a new mean and standard deviation. This is repeated 20 times. A histogram with 512 bins is then computed within a range determined from the final mean and standard deviation. The returned parameters are based on a Gaussian fit to this histogram.
#'
#' @param x numeric values for which to estimate noise.
#' @param height logical indicating whether to use Gaussian function whose area is not fixed at one because of an additional height scaling factor.
#' @param thresh numeric value specifying the factor by which to multiply the standard deviation to determine the threshold away from the mean value within which to include values for the next iteration and final histogram for fitting.
#' @param plot_distributions logical indicating whether to plot the distribution and Gaussian fit used to estimate the noise.
#' @param peak_intensities numeric values of peak intensities to determine the signal to noise.
#' @return a named numeric vector with values:
#'  \describe{
#'   \item{mean}{mean value from the Gaussian fit}
#'   \item{mu}{standard deviation from the Gaussian fit}
#'   \item{h}{height of Gaussian fit (optional depending on value of \code{height} parameter)}
#'   \item{max}{maximum value of \code{x}}
#'  }
#'
#' @export
noise_estimate <- function(x, height=TRUE, thresh=10, plot_distributions=TRUE, peak_intensities=NULL) {

	x <- x[!is.na(x)]

	origname <- NULL
	max_data <- NULL
	
	max_data <- c(max=max(x), max_data)
	
	itermat <- matrix(nrow=20, ncol=2)
	colnames(itermat) <- c("mean", "sd")
	
	itermat[1,"mean"] <- mean(x)
	itermat[1,"sd"] <- stats::sd(as.vector(x))
	
	for (i in 2:nrow(itermat)) {
	
		minthresh <- itermat[i-1,"mean"]-itermat[i-1,"sd"]*thresh
		maxthresh <- itermat[i-1,"mean"]+itermat[i-1,"sd"]*thresh
		idx <- which(x > minthresh & x < maxthresh)
		
		x_sub <- x[idx]
		itermat[i,"mean"] <- mean(x_sub)
		itermat[i,"sd"] <- stats::sd(x_sub)
	}
	
	xhist <- graphics::hist(x[idx], breaks=512, plot=FALSE)
	
	normdist_height_formula <- y ~ h * exp(-(x-mu)^2/(2*sigma^2))
	normdist_formula <- y ~ 1 / (sigma*sqrt(2*pi)) * exp(-(x-mu)^2/(2*sigma^2))

	fit_data <- data.frame(y=xhist$density, x=xhist$mids)
	fit_start <- c(mu=unname(itermat[i,"mean"]), sigma=unname(itermat[i,"sd"]), h=max(xhist$density))
	fit_lower <- fit_start
	fit_lower[] <- -Inf
	fit_lower[2] <- 0
	
	if (height) {
		fit <- stats::nls(normdist_height_formula, fit_data, fit_start, algorithm="port", lower=fit_lower)
	} else {
		fit <- stats::nls(normdist_formula, fit_data, fit_start[1:2], algorithm="port", lower=fit_lower[1:2])
	}
	
	if (plot_distributions) {
		
		fit_pred <- fit$m$predict(data.frame(x=fit_data$x))
		xlim <- c(-thresh, thresh)*fit$m$getPars()["sigma"]
		ylim <- range(xhist$density, fit_pred)
	
		graphics::plot(xhist$mids, xhist$density, type="n", xlim=xlim, ylim=ylim, col="black", xlab="Signal Intensity", ylab="", yaxt="n")
		graphics::abline(h=0, v=0, col="gray")
		graphics::points(xhist$mids, xhist$density, type="l", lwd=0.75)
		graphics::points(fit_data$x, fit_pred, type="l", col="blue", lwd=0.75)
		graphics::abline(v=fit$m$getPars()["mu"], col="blue")
		graphics::abline(v=fit$m$getPars()["mu"]+c(-1,1)*fit$m$getPars()["sigma"], col="blue", lty="dashed")
		
		if (is.null(peak_intensities)) {
			snval <- signif(max_data["max"]/abs(fit$m$getPars()["sigma"]), 3)
		} else {
			snval <- paste(round(min(peak_intensities)/abs(fit$m$getPars()["sigma"])), "-", round(max(peak_intensities)/abs(fit$m$getPars()["sigma"])))
		}
		
		legtext <- as.expression(c(
			substitute(sigma: ~ sigmaval, list(sigmaval=signif(abs(fit$m$getPars()["sigma"]), 3))),
			substitute(mu: ~ muval, list(muval=signif(fit$m$getPars()["mu"], 3)))
		))
		
		legend("topleft", legend=legtext, lwd=c(1, 1), lty=c("dashed", "solid"), col="blue", bty="n", x.intersp=0.5)
		legend("topright", legend=paste("S/N:", snval), bty="n")
	}
	
	c(fit$m$getPars()[1:2], max_data)
}

#' Calculate mapping from assigned peak list onto an unknown peak list
#'
#' This uses a greedy algorithm. It iterates over the unknown peaks in order of
#' decreasing height and assigns each unknown to the closest assigned peak within a
#' distance determined by the thresh parameter, as long as that peak wasn't already
#' assigned.
#'
#' @param assigned two column matrix with assigned peak chemical shifts
#' @param unknown three column matrix with unknown peak chemical shifts and heights
#' @param thresh maximum distance (as a fraction of the two chemical shift ranges)
#' @export
height_assign <- function(assigned, unknown, thresh=0.01) {

	# find widths peak ranges in each dimension
	wd1 <- abs(diff(range(assigned[,1], unknown[,1], na.rm=TRUE)))
	wd2 <- abs(diff(range(assigned[,2], unknown[,2], na.rm=TRUE)))
	
	# normalize assigned peaks by their ranges
	assigned[,1] <- assigned[,1]/wd1
	assigned[,2] <- assigned[,2]/wd2
	
	# normalize unknown peaks by their ranges
	unknown[,1] <- unknown[,1]/wd1
	unknown[,2] <- unknown[,2]/wd2
	
	t_assigned <- t(assigned)
	
	# initialize empty output
	assign_idx <- integer(nrow(assigned))
	
	for (i in order(-abs(unknown[,3]))) {
		
		# distances of all the peaks
		i_dist <- sqrt(colSums((t_assigned - as.numeric(unknown[i,1:2]))^2))
		# index of closest peak
		i_min <- which(i_dist == min(i_dist, na.rm=TRUE))
		
		# check threshold and previous assignment criteria
		if (i_dist[i_min] < thresh && assign_idx[i_min] == 0) {
			assign_idx[i_min] <- i
		}
	}
	
	# set unassigned peaks to NA
	assign_idx[assign_idx == 0] <- NA
	
	assign_idx
}

