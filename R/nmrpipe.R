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
	
	fheader["DMXVAL",header["FDDIMCOUNT"]] <- header["FDDMXVAL"]
	
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