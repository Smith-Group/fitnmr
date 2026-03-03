#' Infer acquisition time for each dimension
#'
#' @param fheader matrix of generalized ND parameters
#' @param empirically_correct logical indicating whether to apply empirical correction
infer_acquisition_time <- function(fheader, empirically_correct=TRUE) {

	# showhdr.c does NDAPOD/NDSW which gives a slightly different answer than below
	acquisition_time <- unname(fheader["TDSIZE",]/fheader["SW",]*ifelse(fheader["FTSIZE",] == 0, 1, fheader["SIZE",]/fheader["FTSIZE",]))
	
	# correct for FID offset as a result of digital oversampling
	acquisition_time <- acquisition_time*unname(1-fheader["DMXVAL",]/fheader["TDSIZE",])
	
	if (empirically_correct) {
	
		correction_factor <- ifelse(fheader["APODCODE",] == 1, 1, 0.5)
		acquisition_time <- acquisition_time*unname(1-correction_factor/fheader["TDSIZE",])
	}
	
	acquisition_time
}

#' Infer original sweep width for each dimension
#'
#' @param fheader matrix of generalized ND parameters
infer_sweep_width <- function(fheader) {

	fheader["SW",]/fheader["OBS",]*ifelse(fheader["FTSIZE",] == 0, 1, fheader["FTSIZE",]/fheader["SIZE",])
}

#' Infer which dimension was directly acquired
#'
#' @param fheader matrix of generalized ND parameters
infer_direct <- function(fheader) {

	direct <- integer(ncol(fheader))
	
	direct[order(fheader["TDSIZE",], decreasing=TRUE)[1]] <- 1
	
	direct
}

#' Infer which dimension was directly acquired
#'
#' @param fheader matrix of generalized ND parameters
#' @param phase_tolerance tolerance in degrees for detecting half-dwell delay
infer_aliasing <- function(fheader, phase_tolerance=10) {

	# assume everything except direct dimension aliases
	alias_vec <- 1L-infer_direct(fheader)
	
	#half_dwell_delay <- abs(fheader["P0",] - -90) < phase_tolerance & abs(fheader["P1",] - 180) < phase_tolerance
	# getFold in fdatap.c only checks for P1 == 180 degrees
	half_dwell_delay <- abs(fheader["P1",] - 180) < phase_tolerance

	alias_vec[half_dwell_delay] <- -alias_vec[half_dwell_delay]

	alias_vec
}

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
#' The \code{fheader} row definitions are as follows (taken from NMRPipe fdatp.h):
#'
#' | Row Name  | Description                            |
#' | --------- | -------------------------------------- |
#' | SIZE      | Number of points in dimension          |
#' | APOD      | Current valid time-domain size         |
#' | SW        | Sweep Width, Hz                        |
#' | ORIG      | Axis Origin (Last Point), Hz           |
#' | OBS       | Obs Freq, MHz                          |
#' | FTFLAG    | 1=Freq Domain 0=Time Domain            |
#' | QUADFLAG  | Data Type Code (See Below)             |
#' | UNITS     | Axis Units Code (See Below)            |
#' | P0        | Zero Order Phase, Degrees              |
#' | P1        | First Order Phase, Degrees             |
#' | CAR       | Carrier Position, PPM                  |
#' | CENTER    | Point Location of Zero Freq            |
#' | AQSIGN    | Sign adjustment needed for FT          |
#' | APODCODE  | Window function used                   |
#' | APODQ1    | Window parameter 1                     |
#' | APODQ2    | Window parameter 2                     |
#' | APODQ3    | Window parameter 3                     |
#' | C1        | Add 1.0 to get First Point Scale       |
#' | ZF        | Negative of Zero Fill Size             |
#' | X1        | Extract region origin, if any, pts     |
#' | XN        | Extract region endpoint, if any, pts   |
#' | OFFPPM    | Additional PPM offset (for alignment)  |
#' | FTSIZE    | Size of data when FT performed         |
#' | TDSIZE    | Original valid time-domain size        |
#' | LB        | Extra Exponential Broadening, Hz       |
#' | GB        | Extra Gaussian Broadening, Hz          |
#' | GOFF      | Offset for Gaussian Broadening, 0 to 1 |
#' | OBSMID    | Original Obs Freq before 0.0ppm adjust |
#'
#' In addition several rows contain information inferred from the header data:
#'
#' | Row Name  | Description                            |
#' | --------- | -------------------------------------- |
#' | aq_s      | Acquisition time, seconds              |
#' | sw_ppm    | Original Sweep Width, PPM              |
#' | direct    | Direct (1) or indirect (0) dimension   |
#' | alias     | Aliasing (0/1) with inversion (-1)     |
#' | mag       | Magnitude mode (direct from FDMCFLAG)  |
#'
#' @examples
#' spec_file <- system.file("extdata", "t1", "1.ft2", package="fitnmr")
#' spec <- read_nmrpipe(spec_file, dim_order="hx")
#' str(spec)
#'
#' @md
#' @export
read_nmrpipe <- function(inFormat, dim_order=NULL, complex_data=FALSE) {
	
	if (length(grep("%[0-9]+[di]", inFormat))) {
		inFile <- sprintf(inFormat, 1)
	} else {
		inFile <- inFormat
	}
	
	if (missing(inFile))
		stop('The inFile file path is required')	
	
	## Test file connection
	readCon <- file(inFile, 'rb')
	header_raw <- try(readBin(readCon, what='raw', n=512*4, size=1), silent=TRUE)
	if (inherits(header_raw, "try-error") || length(header_raw) != 512*4){
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
	
	fnames <- c("APOD", "SW", "ORIG", "OBS", "FTFLAG", "QUADFLAG", "UNITS", "P0", "P1", "CAR", "CENTER", "AQSIGN", "APODCODE", "APODQ1", "APODQ2", "APODQ3", "C1", "ZF", "X1", "XN", "OFFPPM", "FTSIZE", "TDSIZE", "LB", "GB", "GOFF", "OBSMID")
	
	forder <- header[c("FDDIMORDER1", "FDDIMORDER2", "FDDIMORDER3", "FDDIMORDER4")]
	
	fheader_names <- matrix(paste(rep(c("FDF1", "FDF2", "FDF3", "FDF4"), each=length(fnames)), rep(fnames, 4), sep=""), ncol=4)
	
	#fheader_names <- rbind(fheader_names, c("FDSPECNUM", "FDSIZE", "FDF3SIZE", "FDF4SIZE"))
	#fnames <- c(fnames, "SIZE")
	
	fheader <- matrix(header[fheader_names], ncol=4, dimnames=list(fnames, flabels))
	
	fsizes <- header[c("FDSIZE", "FDSPECNUM", "FDF3SIZE", "FDF4SIZE")][order(forder)]
	
	#print(header)
	
	fheader <- rbind(SIZE=unname(fsizes), fheader, DMXVAL=0)
	
	fheader <- rbind(
		fheader,
		aq_s=infer_acquisition_time(fheader),
		sw_ppm=infer_sweep_width(fheader),
		direct=infer_direct(fheader),
		alias=infer_aliasing(fheader),
		magnitude=0
	)
	
	# placeholder code assuming only direct dimension is magnitude mode (for HMBC spectra)
	if (header["FDMCFLAG"]) {
		fheader["magnitude",which(fheader["direct",] == 1)] <- 1
	}
	
	fheader["DMXVAL",fheader["direct",] == 1] <- header["FDDMXVAL"]
	
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
	
		if (header["FDFILECOUNT"] > 1) {
			inFile <- sprintf(inFormat, 1)
		} else {
			inFile <- inFormat
		}
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
#' @param file_path path to the NMRDraw peak table file
#' @return A `data.frame` containing the peak table, with columns defined by the
#'   `VARS` line in the NMRDraw file.
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
#' @param peak_tab data frame containing peak table data
#' @param file_path path to write the NMRDraw peak table file
#' @return No return value, called for side effects (writes `file_path`).
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
#' @param ppm_mat matrix of ppm values
#' @param fheader matrix of generalized ND parameters
#' @return A numeric matrix of point indices with the same dimensions as
#'   `ppm_mat`, with column names converted from `*_PPM` to `*_AXIS`.
#' @export
ppm_to_pts <- function(ppm_mat, fheader) {

	orig <- fheader["CAR",]*fheader["OBS",]-fheader["SW",]/2+fheader["SW",]/fheader["FTSIZE",]
	
	# orig+fheader["SW",] seems wrong to me...
	points_mat <- t((orig+fheader["SW",]-t(ppm_mat)*fheader["OBS",])/fheader["SW",]*fheader["FTSIZE",])
	
	colnames(points_mat) <- sub("_PPM", "_AXIS", colnames(points_mat))
	
	points_mat
}

#' Convert widths in Hz into points
#'
#' @param whz_mat matrix of widths in Hz
#' @param fheader matrix of generalized ND parameters
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

#' Simulate an FID using the NMRPipe SimTimeND function
#'
#' @param peak_tab NMRDraw peak table
#' @param fheader matrix of generalized ND parameters
#' @param rms RMS noise level
#' @param iseed random seed
#' @param file_path optional output path for the simulated FID
#' @param verbose logical indicating whether to print the command
#' @return The integer exit status from `system2("SimTimeND", ...)` (typically
#'   `0` on success).
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
#' @param in_path input file path
#' @param out_path output file path
#' @param ndim number of dimensions
#' @param apod optional apodization argument matrix
#' @param sp optional sine-bell argument matrix
#' @param zf optional zero-fill argument matrix
#' @param ps optional phase correction argument matrix
#' @param ext optional extraction argument matrix
#' @return The integer exit status from `system(...)` (typically `0` on success).
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

# spectrum processing functions that closely mimic the functionality in NMRPipe

shifter <- function(x, n = 1) {

	if (n == 0) {
  		x
  	} else {
  		c(utils::tail(x, -n), utils::head(x, n))
  	}
}

#' Apply sine-based window function to a 1D FID
#'
#' This function aims to numerically match the NMRPipe `SP` function.
#'
#' See the NMRPipe documentation: https://www.nmrscience.com/ref/nmrpipe/sp.html
#'
#' @param fid list with `int`, `header`, and `fheader` elements containing FID data
#' @param off equivalent to `-off` flag
#' @param end equivalent to `-end` flag
#' @param pow equivalent to `-pow` flag
#' @param cval equivalent to `-c` flag
#' @param dmx equivalent to `-dmx` flag
#'
#' @examples
#' fid_path <- system.file("extdata", "noesy1d", "11", "test.fid", package = "fitnmr")
#' fid <- read_nmrpipe(fid_path, complex_data = TRUE)
#' sp <- nmrpipe_sp(fid)
#'
#' @return A modified FID `list` with updated `int`, `header`, and `fheader`
#'   values after sine-bell apodization.
#' @export
nmrpipe_sp <- function(fid, off=0, end=1, pow=1, cval=1, dmx=FALSE) {

	x <- seq(0, 1, length.out=length(fid$int))
	
	fid$int[] <- fid$int[] * sin( pi*off + pi*(end-off)*x )^pow
	
	if (dmx) {
		fid$int[fid$header["FDDMXVAL"]+1] <- fid$int[fid$header["FDDMXVAL"]+1]*cval
	} else {
		fid$int[1] <- fid$int[1]*cval
	}
	
	fid$header["FDF2APODCODE"] <- fid$fheader["APODCODE",1] <- 1
	fid$header["FDF2APODQ1"] <- fid$fheader["APODQ1",1] <- off
	fid$header["FDF2APODQ2"] <- fid$fheader["APODQ2",1] <- end
	fid$header["FDF2APODQ3"] <- fid$fheader["APODQ3",1] <- pow
	fid$header["FDF2C1"] <- fid$fheader["C1",1] <- cval
	
	fid
}

#' Apply zero filling to a 1D FID
#'
#' This function aims to numerically match the NMRPipe `ZF` function.
#'
#' See the NMRPipe documentation: https://www.nmrscience.com/ref/nmrpipe/zf.html
#'
#' @param fid list with `int`, `header`, and `fheader` elements containing FID data
#' @param zf equivalent to `-zf` flag
#' @param pad equivalent to `-pad` flag
#' @param size equivalent to `-size` flag
#' @param auto equivalent to `-auto` flag
#'
#' @examples
#' fid_path <- system.file("extdata", "noesy1d", "11", "test.fid", package = "fitnmr")
#' fid <- read_nmrpipe(fid_path, complex_data = TRUE)
#' sp <- nmrpipe_sp(fid)
#' zf <- nmrpipe_zf(sp)
#'
#' @return A modified FID `list` with zero-filled `int` data and updated
#'   `header`/`fheader` metadata.
#' @export
nmrpipe_zf <- function(fid, zf=1, pad=NULL, size=NULL, auto=FALSE) {
	
	stopifnot(!((!is.null(pad)) && (!is.null(size))))

	if (!is.null(pad)) {
	
		size <- length(fid$int)+pad
		
	} else if (is.null(size)) {
	
		size <- length(fid$int)*2^zf
	}
			
	if (auto) {
		size <- 2^ceiling(log(size,2))
	}
	
	pad <- size-length(fid$int)
	
	fid$int <- array(c(fid$int[seq_len(min(length(fid$int),size))], complex(max(0, pad))), size)
	
	dimnames(fid$int)[[1]] <- seq(-fid$header["FDDMXVAL"]/fid$fheader["SW",], by=1/fid$fheader["SW",], length.out=size)
	
	fid$header["FDF2CENTER"] <- fid$fheader["CENTER",1] <- size/2+1
	fid$header["FDSIZE"] <- fid$fheader["SIZE",1] <- size
	fid$header["FDF2ORIG"] <- fid$fheader["ORIG",1] <- fid$fheader["CAR",]*fid$fheader["OBS",]+(fid$fheader["CENTER",]/fid$fheader["SIZE",]-1)*fid$fheader["SW",]
	fid$header["FDF2ZF"] <- fid$fheader["ZF",1] <- -size
	
	fid
}

#' Fourier transform a 1D FID
#'
#' This function aims to numerically match the NMRPipe `FT` function.
#'
#' See the NMRPipe documentation: https://www.nmrscience.com/ref/nmrpipe/ft.html
#'
#' @param fid list with `int`, `header`, and `fheader` elements containing FID data
#'
#' @examples
#' fid_path <- system.file("extdata", "noesy1d", "11", "test.fid", package = "fitnmr")
#' fid <- read_nmrpipe(fid_path, complex_data = TRUE)
#' sp <- nmrpipe_sp(fid)
#' zf <- nmrpipe_zf(sp)
#' ft <- nmrpipe_ft(zf)
#'
#' @return A transformed spectrum `list` with updated `int`, `ppm`, `header`,
#'   and `fheader` fields.
#' @export
nmrpipe_ft <- function(fid) {
	
	if (fid$header["FDDMXFLAG"] <= 0 && fid$header["FDDMXVAL"] != 0) {
		fid$int[] <- shifter(fid$int, fid$header["FDDMXVAL"])
		fid$header["FDDMXFLAG"] <- 1
	}
	
	fid$int[] <- rev(shifter(stats::fft(fid$int), length(fid$int)/2+1))
	
	fid$header["FDF2FTSIZE"] <- fid$fheader["FTSIZE",1] <- length(fid$int)
	fid$header["FDF2FTFLAG"] <- fid$fheader["FTFLAG",1] <- 1
	
	fid$header["FDMAX"] <- fid$header["FDDISPMAX"] <- max(Re(fid$int[]))
	fid$header["FDMIN"] <- fid$header["FDDISPMIN"] <- min(Re(fid$int[]))
	
	ppm <- apply(fid$fheader, 2, function(x) (x["ORIG"]+x["SW"]*(1-seq_len(x["SIZE"])/x["SIZE"]))/x["OBS"])
	if (!is.list(ppm)) {
		ppm <- list(ppm)
	}
	fid$ppm <- ppm
	dimnames(fid$int) <- ppm
	
	fid
}

#' Inverse Fourier transform a 1D spectrum
#'
#' This function aims to numerically match the NMRPipe `FT -inv` function.
#'
#' See the NMRPipe documentation: https://www.nmrscience.com/ref/nmrpipe/ft.html
#'
#' @param ft list with `int`, `header`, and `fheader` elements containing spectrum data
#'
#' @examples
#' fid_path <- system.file("extdata", "noesy1d", "11", "test.fid", package = "fitnmr")
#' fid <- read_nmrpipe(fid_path, complex_data = TRUE)
#' sp <- nmrpipe_sp(fid)
#' zf <- nmrpipe_zf(sp)
#' ft <- nmrpipe_ft(zf)
#' fti <- nmrpipe_fti(ft)
#'
#' @return A modified FID-like `list` with inverse-transformed `int` values and
#'   updated `header`/`fheader` metadata.
#' @export
nmrpipe_fti <- function(ft) {
	
	ft$int[] <- stats::fft(shifter(rev(ft$int), -(length(ft$int)/2+1)), inverse=TRUE)/length(ft$int)
	
	ft$header["FDF2FTSIZE"] <- ft$fheader["FTSIZE",1] <- length(ft$int)
	ft$header["FDF2FTFLAG"] <- ft$fheader["FTFLAG",1] <- 0
	
	ft$header["FDMAX"] <- ft$header["FDDISPMAX"] <- max(Re(ft$int[]))
	ft$header["FDMIN"] <- ft$header["FDDISPMIN"] <- min(Re(ft$int[]))
	
	if (ft$header["FDDMXFLAG"] == 1 && ft$header["FDDMXVAL"] != 0) {
		ft$int[] <- shifter(ft$int, -ft$header["FDDMXVAL"])
		ft$header["FDDMXFLAG"] <- -1
	}
	
	dimnames(ft$int)[[1]] <- seq(-ft$header["FDDMXVAL"]/ft$fheader["SW",], by=1/ft$fheader["SW",], length.out=length(ft$int))
	
	ft
}

#' Inverse Fourier transform a 1D spectrum
#'
#' This function aims to numerically match the NMRPipe `PS` function.
#'
#' See the NMRPipe documentation: https://www.nmrscience.com/ref/nmrpipe/ps.html
#'
#' @param ft list with `int`, `header`, and `fheader` elements containing spectrum data
#' @param p0 equivalent to `-p0` flag
#' @param p1 equivalent to `-p1` flag
#'
#' @examples
#' fid_path <- system.file("extdata", "noesy1d", "11", "test.fid", package = "fitnmr")
#' fid <- read_nmrpipe(fid_path, complex_data = TRUE)
#' sp <- nmrpipe_sp(fid)
#' zf <- nmrpipe_zf(sp)
#' ft <- nmrpipe_ft(zf)
#' ps <- nmrpipe_ps(ft, p0 = 129, p1 = 3)
#'
#' @return A modified spectrum `list` with phase-corrected `int` values and
#'   updated phase metadata in `header` and `fheader`.
#' @export
nmrpipe_ps <- function(ft, p0=0, p1=0) {
	
	ft$header["FDF2P0"] <- ft$fheader["P0",1] <- p0
	ft$header["FDF2P1"] <- ft$fheader["P1",1] <- p1

	p1_frac <- as.vector(sapply(1, function(j) {
		frac <- seq(0, 1, length.out=ft$fheader["FTSIZE",j]+1)
		frac <- frac[-length(frac)]
		if (all(ft$fheader[c("X1","XN"),j] != 0)) {
			frac <- frac[seq(ft$fheader["X1",j], ft$fheader["XN",j])]
		}
		frac
	}))
	
	p0 <- p0*pi/180 # convert degrees to rad
	p1 <- p1*pi/180 # convert degrees to rad
	
	pvec <- exp(1i*(p0+p1*p1_frac))
	
	ft$int[] <- ft$int[]*pvec
	
	ft
}
