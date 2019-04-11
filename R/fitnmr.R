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
	
	fheader["DMXVAL",forder==1] <- header["FDDMXVAL"]
	
	FDFILECOUNT <- header["FDFILECOUNT"]
	FDDIMCOUNT <- header["FDDIMCOUNT"]
	
	f_ppm <- apply(fheader, 2, function(x) (x["ORIG"]+x["SW"]*(1-seq_len(x["SIZE"])/x["SIZE"]))/x["OBS"])
	f_ppm <- f_ppm[head(forder, header["FDDIMCOUNT"])]
	f_size <- sapply(f_ppm, length)
	n <- if (header["FDFILECOUNT"] == 1) prod(f_size) else prod(head(f_size, 2))
	
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
	
	
	fheader <- fheader[,head(forder, header["FDDIMCOUNT"]),drop=FALSE]
	
	#print(apply(fheader, 2, function(x) (x["ORIG"]+x["SW"]*(1-x["CENTER"]/x["SIZE"]))/x["OBS"]))
	
	#print(f_ppm)

	if (!is.null(dim_order)) {
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

lineshapes_parsed <- lapply(lineshapes, function(flist) lapply(flist, function(fvec) parse(text=fvec, keep.source=FALSE)))

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

fill_group <- function(x, array_dim, na_dim1_same=FALSE) {

	x <- fill_array_int(x, array_dim)
	
	na_idx <- is.na(x)
	
	if (na_dim1_same) {
		x[na_idx] <- head(setdiff(seq_along(x), x[!na_idx]), array_dim[1])
	} else {
		x[na_idx] <- head(setdiff(seq_along(x), x[!na_idx]), sum(na_idx))
	}
	
	x
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
	
	c(omega0_vec, r2_vec, m0_vec, p0_vec, p1_vec)
}

unpack_fit_params <- function(start_vec, group_list, default_list=group_list) {

	unpacked_list <- default_list
	
	unpacked_list[["omega0"]][group_list[["omega0"]] != 0] <- start_vec[paste("omega0", group_list[["omega0"]][group_list[["omega0"]] != 0], sep="_")]
	unpacked_list[["r2"]][group_list[["r2"]] != 0] <- start_vec[paste("r2", group_list[["r2"]][group_list[["r2"]] != 0], sep="_")]
	unpacked_list[["m0"]][group_list[["m0"]] != 0] <- start_vec[paste("m0", group_list[["m0"]][group_list[["m0"]] != 0], sep="_")]
	unpacked_list[["p0"]][group_list[["p0"]] != 0] <- start_vec[paste("p0", group_list[["p0"]][group_list[["p0"]] != 0], sep="_")]
	unpacked_list[["p1"]][group_list[["p1"]] != 0] <- start_vec[paste("p1", group_list[["p1"]][group_list[["p1"]] != 0], sep="_")]
	
	unpacked_list
}

group_param_idx <- function(param_names, group_list) {

	idx_list <- group_list
	
	for (i in seq_along(group_list)) {
		idx_list[[i]][] <- match(paste(names(group_list)[i], group_list[[i]], sep="_"), param_names)
	}
	
	idx_list
}

#' @export
make_fit_input <- function(spectra, omega0_start, omega0_plus, omega0_minus=omega0_plus, r2_start=NULL, m0_start=NULL, m0_region=(omega0_plus+omega0_minus)/2, p0_start=0, p1_start=0, omega0_group=NULL, r2_group=NULL, m0_group=NULL, p0_group=0, p1_group=0, fheader=NULL) {

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
	
	if (is.null(dim(omega0_start))) {
		dim(omega0_start) <- c(n_dimensions, length(omega0_start)/n_dimensions/n_spectra, n_spectra)
	}
	
	omega0_plus <- fill_array(omega0_plus, dim(omega0_start))
	omega0_minus <- fill_array(omega0_plus, dim(omega0_minus))
	
	r2_start <- fill_array(r2_start, dim(omega0_start))
	m0_start <- fill_array(m0_start, dim(omega0_start)[-1])
	p0_start <- fill_array(p0_start, dim(omega0_start))
	p1_start <- fill_array(p1_start, dim(omega0_start))
	
	omega0_group <- fill_group(omega0_group, dim(omega0_start))
	r2_group <- fill_group(r2_group, dim(omega0_start))
	m0_group <- fill_group(m0_group, dim(omega0_start)[-1])
	p0_group <- fill_group(p0_group, dim(omega0_start), TRUE)
	p1_group <- fill_group(p1_group, dim(omega0_start), TRUE)
	
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
			( omega_eval[[j]] - min(spec_ppm[[j]]) ) / (fheader["SW",j] / fheader["OBS",j] )
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
		p1=p1_start
	)
	
	group_list <- list(
		omega0=omega0_group,
		r2=r2_group,
		m0=m0_group,
		p0=p0_group,
		p1=p1_group
	)
	
	lower_list <- upper_list <- start_list
	
	for (i in seq_along(lower_list)) {
		lower_list[[i]][] <- -Inf
		upper_list[[i]][] <- Inf
	}
	
	lower_list[["omega0"]] <- omega0_start - omega0_minus
	upper_list[["omega0"]] <- omega0_start + omega0_plus
	
	#print(start_list)
	
	list(
		spec_data=spec_data,
		start_list=start_list,
		group_list=group_list,
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
	
	pvec <- exp(1i*(p0-p1*p1_frac))
	
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
	
	#pvec <-          exp(1i*(p0 - p1*p1_frac))
	#dpvec_dp0 <-  1i*exp(1i*(p0 - p1*p1_frac))
	#dpvec_dp1 <- -1i*exp(1i*(p0 - p1*p1_frac))*p1_frac
	
	pvec <- exp(1i*(p0 - p1*p1_frac))
	dpvec_dp0 <- 1i*pvec
	dpvec_dp1 <- -dpvec_dp0*p1_frac
	
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

	param_list <- unpack_fit_params(par, fit_data$group_list, fit_data$start_list)

	for (spec_idx in seq_along(fit_data$spec_data)) {
	
		spec_data <- fit_data[["spec_data"]][[spec_idx]]
		
		spec_eval_idx <- seq_along(spec_data$spec_int)+spec_data$spec_offset
		
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
			
			func_eval[spec_eval_idx] <- func_eval[spec_eval_idx] + func_nd_prod
		}
		
		if (return_resid) {
			func_eval[spec_eval_idx] <- spec_data$spec_int - func_eval[spec_eval_idx]
		}
	}
	
	func_eval
}

fit_jac <- function(par, fit_data) {

	jac_eval <- matrix(0, nrow=fit_data$num_points, ncol=length(par), dimnames=list(NULL, names(par)))

	param_list <- unpack_fit_params(par, fit_data$group_list, fit_data$start_list)

	idx_list <- group_param_idx(names(par), fit_data$group_list)

	for (spec_idx in seq_along(fit_data$spec_data)) {
	
		spec_data <- fit_data[["spec_data"]][[spec_idx]]
		
		spec_eval_idx <- seq_along(spec_data$spec_int)+spec_data$spec_offset
		
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
				jac_eval[spec_eval_idx,idx_list[["omega0"]][idx,peak_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[["omega0"]][idx,peak_idx,spec_idx]] + deriv_nd_prod
			}
			}
			
			for (var_name in c("omega0", "r2", "p0", "p1")) {
				var_idx <- which(!is.na(idx_list[[var_name]][,peak_idx,spec_idx]))
				for (idx in var_idx) {
					deriv_nd_prod <- deriv_1d_evals[[idx]][spec_data$spec_eval_idx[,idx],var_name]*param_list[["m0"]][peak_idx,spec_idx]
					for (i in seq_len(ncol(spec_data$spec_eval_idx))[-idx]) {
						deriv_nd_prod <- deriv_nd_prod*deriv_1d_evals[[i]][spec_data$spec_eval_idx[,i],"f"]
					}
					jac_eval[spec_eval_idx,idx_list[[var_name]][idx,peak_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[[var_name]][idx,peak_idx,spec_idx]] + deriv_nd_prod
				}
			}
			
			if (!is.na(idx_list[["m0"]][peak_idx,spec_idx])) {
				func_nd_prod <- deriv_1d_evals[[1]][spec_data$spec_eval_idx[,1],"f"]
				for (i in seq_len(ncol(spec_data$spec_eval_idx))[-1]) {
					func_nd_prod <- func_nd_prod*deriv_1d_evals[[i]][spec_data$spec_eval_idx[,i],"f"]
				}
				jac_eval[spec_eval_idx,idx_list[["m0"]][peak_idx,spec_idx]] <- jac_eval[spec_eval_idx,idx_list[["m0"]][peak_idx,spec_idx]] + func_nd_prod
			}
		}
	}
	
	#print(jac_eval)
	
	-jac_eval
}

#' @export
perform_fit <- function(fit_input) {

	fit_par <- pack_fit_params(fit_input$start_list, fit_input$group_list)
	fit_lower <- pack_fit_params(fit_input$lower_list, fit_input$group_list)
	fit_upper <- pack_fit_params(fit_input$upper_list, fit_input$group_list)
		
	fit <- suppressWarnings(minpack.lm::nls.lm(fit_par, fit_lower, fit_upper, fn=fit_fn, jac=fit_jac, fit_data=fit_input, control=list(minpack.lm::nls.lm.control(nprint=0), maxiter = 1e5)))
	
	fit_input[["fit_list"]] <- unpack_fit_params(fit$par, fit_input$group_list, default_list=fit_input$start_list)
	fit_input[["fit_rsstrace"]] <- fit$rsstrace
		
	fit_input
}

fit_jac_nlfb <- function(par, fit_data) {

	jj <- fit_jac(par, fit_data)
	attr(jj, "gradient") <- jj
	jj
}

perform_fit_nlfb <- function(fit_input) {

	fit_par <- pack_fit_params(fit_input$start_list, fit_input$group_list)

	fit <- nlfb(fit_par, resfn=fit_fn, jacfn=fit_jac_nlfb, fit_data=fit_input)
	
	fit_input[["fit_list"]] <- unpack_fit_params(fit$coefficients, fit_input$group_list, default_list=fit_input$start_list)
	
	fit_input
}

#' @export
get_spec_int <- function(fit_data, spec_type=c("input", "start", "fit"), spec_idx=seq_along(fit_data$spec_data), peak_idx=seq_len(dim(fit_data$start_list$omega0)[2])) {

	spec_type <- match.arg(spec_type)
	
	all_peak_idx <- seq_len(dim(fit_data$start_list$omega0)[2])
	
	if (spec_type == "input") {
	
		int <- fit_data$spec_data[[1]]$spec_int
		
		stopifnot(all(all_peak_idx %in% peak_idx))
	
	} else {
	
		if (spec_type == "start") {
		
			param_list <- fit_data[["start_list"]]
		
		} else if (spec_type == "fit") {
		
			param_list <- fit_data[["fit_list"]]
		}
		
		param_list$m0[!all_peak_idx %in% peak_idx,] <- 0
		
		#print(str(param_list))
		
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
		
		omega_ppm_seg_ends <- which(abs(diff(omega_ppm)) > abs(median(diff(omega_ppm)))*2)
		omega_ppm_seg_starts <- c(1, omega_ppm_seg_ends+1)
		omega_ppm_seg_ends <- c(omega_ppm_seg_ends, length(omega_ppm))
		plot_idx <- unlist(lapply(seq_along(omega_ppm_seg_starts), function(i) c(NA, seq(omega_ppm_seg_starts[i], omega_ppm_seg_ends[i]))))[-1]
		
		spec_eval_idx <- seq_along(spec_data$spec_int)+spec_data$spec_offset
		
		spec_original_int <- original_int[spec_eval_idx]
		spec_start_int <- start_int[spec_eval_idx]
		spec_fit_int <- fit_int[spec_eval_idx]
		
		ylim <- range(0, spec_original_int, spec_start_int, spec_fit_int)
		
		plot(omega_ppm[plot_idx], spec_original_int[plot_idx], type="l", xlim=rev(range(omega_ppm)), ylim=ylim, xlab=expression(delta (ppm)), ylab="Intensity")

		if (!is.null(spec_start_int)) {
			points(omega_ppm[plot_idx], spec_start_int[plot_idx], type="l", col="blue")
		}
		
		if (!is.null(spec_fit_int)) {
			points(omega_ppm[plot_idx], spec_fit_int[plot_idx], type="l", col="red")
		}
	}
}

#' Plot a two dimensional peak fit
#'
#' @export
plot_fit_2d <- function(fit_output, spec_ord, plot_start=FALSE, main=NULL) {

	input_spec_int <- fitnmr::get_spec_int(fit_output, "input")
	fit_spec_int <- fitnmr::get_spec_int(fit_output, "fit")
	
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
		title(main)
		if (plot_start) {
			fitnmr::contour_pipe(aperm(input_spec_int[[i]], spec_ord), zlim=zlim, col_pos="blue", col_neg="lightblue", add=TRUE)
			points(t(fit_output$start_list$omega0[spec_ord,,i]), pch=16, col="blue")
			segments(fit_output$start_list$omega0[spec_ord[1],,i], start_r2_ppm_low[spec_ord[2],,i], fit_output$start_list$omega0[spec_ord[1],,i], start_r2_ppm_high[spec_ord[2],,i], col="blue")
			segments(start_r2_ppm_low[spec_ord[1],,i], fit_output$start_list$omega0[spec_ord[2],,i], start_r2_ppm_high[spec_ord[1],,i], fit_output$start_list$omega0[spec_ord[2],,i], col="blue")
		}
		fitnmr::contour_pipe(aperm(fit_spec_int[[i]], spec_ord), zlim=zlim, col_pos="red", col_neg="pink", add=TRUE)
		rect(fit_output$upper_list$omega0[spec_ord[1],,i], fit_output$upper_list$omega0[spec_ord[2],,i], fit_output$lower_list$omega0[spec_ord[1],,i], fit_output$lower_list$omega0[spec_ord[2],,i], border="gray")
		points(t(fit_output$fit_list$omega0[spec_ord,,i]), pch=16, col="red")
		segments(fit_output$fit_list$omega0[spec_ord[1],,i], fit_r2_ppm_low[spec_ord[2],,i], fit_output$fit_list$omega0[spec_ord[1],,i], fit_r2_ppm_high[spec_ord[2],,i], col="red")
		segments(fit_r2_ppm_low[spec_ord[1],,i], fit_output$fit_list$omega0[spec_ord[2],,i], fit_r2_ppm_high[spec_ord[1],,i], fit_output$fit_list$omega0[spec_ord[2],,i], col="red")
	}
}

#' @export
contour_pipe <- function(data_mat, nlevels=10, zlim=range(data_mat, na.rm=TRUE), low_frac=0.05, lwd=0.25, main=NA, col_pos="blue", col_neg="red", add=FALSE, xlab=NULL, ylab=NULL, frame.plot=TRUE) {

	zlim <- zlim

	usr <- par("usr")

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

	contour(rev(y), rev(x), data_mat_transform, levels=levels, xlim=rev(range(y)), ylim=rev(range(x)), drawlabels=FALSE, add=add, col=col, lwd=lwd, xaxs="i", yaxs="i", main=main, xlab=xlab, ylab=ylab, frame.plot=frame.plot)
}

read_nmrdraw_peak_tab_old <- function(filepath) {
	vars <- readLines(filepath, 16)[16]
	vars <- strsplit(substring(vars, 8), " ")[[1]]
	
	peaktab <- read.table(filepath, skip=18)
	colnames(peaktab) <- vars
	
	noise_level <- strsplit(readLines(filepath, 5)[5], " ")[[1]][3]
	noise_level <- strsplit(noise_level, ",")[[1]][1]
	noise_level <- as.numeric(noise_level)
	#print(noise_level)
	
	attr(peaktab, "noise") <- noise_level
	
	peaktab
}

#' @export
read_nmrdraw_peak_tab <- function(file_path) {

	peak_tab_lines <- readLines(file_path)
	
	vars_line_idx <- grep("^VARS", peak_tab_lines)
	format_line_idx <- grep("^FORMAT", peak_tab_lines)
	
	var_names <- strsplit(peak_tab_lines[vars_line_idx], " +")[[1]][-1]
	formats <- strsplit(peak_tab_lines[format_line_idx], " +")[[1]][-1]
	
	blank_line_idx <- grep("^$", peak_tab_lines)
	
	peak_tab <- read.table(text=peak_tab_lines, skip=tail(blank_line_idx, 1), stringsAsFactors=FALSE)
	colnames(peak_tab) <- var_names
	
	peak_tab
}

peak_tab_formats <- c("%5d", "%9.3f", "%9.3f", "%6.3f", "%6.3f", "%8.3f", "%8.3f", "%9.3f", "%9.3f", "%7.3f", "%7.3f", "%8.3f", "%8.3f", "%4d", "%4d", "%4d", "%4d", "%+e", "%+e", "%+e", "%.5f", "%d", "%s", "%4d", "%4d")
names(peak_tab_formats) <- c("INDEX", "X_AXIS", "Y_AXIS", "DX", "DY", "X_PPM", "Y_PPM", "X_HZ", "Y_HZ", "XW", "YW", "XW_HZ", "YW_HZ", "X1", "X3", "Y1", "Y3", "HEIGHT", "DHEIGHT", "VOL", "PCHI2", "TYPE", "ASS", "CLUSTID", "MEMCNT")

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

ppm_to_pts <- function(ppm_mat, fheader) {

	orig <- fheader["CAR",]*fheader["OBS",]-fheader["SW",]/2+fheader["SW",]/fheader["FTSIZE",]
	
	# orig+fheader["SW",] seems wrong to me...
	points_mat <- t((orig+fheader["SW",]-t(ppm_mat)*fheader["OBS",])/fheader["SW",]*fheader["FTSIZE",])
	
	colnames(points_mat) <- sub("_PPM", "_AXIS", colnames(points_mat))
	
	points_mat
}

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
	}
	
	names(fit_params) <- paste(names(fit_params), "_start", sep="")
	names(group_params) <- paste(names(group_params), "_group", sep="")
	
	c(fit_params, group_params)
}

#' Determine the region of a spectrum containing the majority of the fit peaks
fit_footprint <- function(fit, frac=0.12) {

	fit_spec_int <- fitnmr::get_spec_int(fit, "fit")
	
	sort_int <- sort(abs(fit_spec_int[[1]]))
	cum_int <- cumsum(sort_int)
	cum_frac_int <- cum_int/tail(cum_int, 1)
	
	thresh <- sort_int[which(cum_frac_int > frac)[1]]
	
	fit_spec_int[[1]][is.na(fit_spec_int[[1]])] <- 0
	
	abs(fit_spec_int[[1]]) > thresh
	
	#abs(fit_spec_int[[1]]) > max(abs(fit_spec_int[[1]]))*.025
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
	
	fit_input$upper_list$omega0 <- fit_input$start_list$omega0+r2_ppm*factor
	fit_input$lower_list$omega0 <- fit_input$start_list$omega0-r2_ppm*factor
	
	fit_input
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
fit_peaks <- function(spec_list, cs_mat, fit_prev=NULL, spec_ord=1:2, omega0_plus=c(0.075, 0.75), r2_start=c(5,5), positive_only=TRUE) {
	
	plot=FALSE
	
	idx_prev <- rep(FALSE, nrow(cs_mat))
	
	if (is.null(fit_prev)) {
		
		omega0_start <- t(cs_mat[,spec_ord,drop=FALSE])
		dim(omega0_start) <- c(dim(omega0_start), 1)
		
		param_list <- list(
			omega0_start=omega0_start,
			r2_start=r2_start[spec_ord],
			m0_start=1,
			omega0_group=0,
			r2_group=0,
			m0_group=NA
		)
		
		r2_group <- array(seq_len(prod(dim(cs_mat))), dim(omega0_start))
		
	} else {
	
		param_list <- extract_params(fit_prev, nrow(cs_mat))
		
		idx_prev <- c(rep(TRUE, dim(fit_prev$fit_list$omega0)[2]), idx_prev)
		
		r2_group <- param_list$r2_group
		
		peak_int_list <- get_spec_peak_int(fit_prev)[[1]]
		norm_peak_int_list <- lapply(peak_int_list, function(x) abs(x)/sum(x, na.rm=TRUE))
		
		int_ppm <- lapply(dimnames(peak_int_list[[1]]), as.numeric)
		
		cs_mat_idx <- matrix(sapply(seq_along(spec_ord), function(i) apply(abs(outer(int_ppm[[i]], cs_mat[,spec_ord[i]], "-")), 2, which.min)), ncol=ncol(cs_mat))
		
		prev_overlap_idx <- sapply(nrow(cs_mat_idx), function(i) which.max(sapply(norm_peak_int_list, "[", cs_mat_idx[i,,drop=FALSE])))
		
		r2_group[,!idx_prev,] <- r2_group[,prev_overlap_idx,]
		
		#print(r2_group)
		
		#print(idx_prev)
		
		param_list$omega0_start[,!idx_prev,] <- as.numeric(cs_mat[,spec_ord])
		param_list$r2_start[,!idx_prev,] <- param_list$r2_start[,prev_overlap_idx,]
		param_list$m0_start[!idx_prev,] <- 1
		
		param_list$omega0_group[] <- 0
		param_list$r2_group[] <- 0
		param_list$m0_group[!idx_prev,] <- 1
		param_list$m0_group[idx_prev,] <- 0
		param_list$m0_group[] <- NA
	}
	
	fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list))
	if (positive_only) {
		fit_input$lower_list$m0[] <- 0
	}
	fit_output <- fitnmr::perform_fit(fit_input)
	param_list <- extract_params(fit_output)
	
	if (plot) plot_fit_2d(fit_output, spec_ord, plot_start=TRUE, "Fit Parameters: m0")
	
	# unfix r2 values
	param_list$r2_group[] <- r2_group[]
	
	fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list))
	if (positive_only) {
		fit_input$lower_list$m0[] <- 0
	}
	fit_input$lower_list$r2[] <- 0
	fit_input$upper_list$r2[] <- 20
	fit_output <- fitnmr::perform_fit(fit_input)
	param_list <- extract_params(fit_output)
	
	if (plot) plot_fit_2d(fit_output, spec_ord, plot_start=TRUE, "Fit Parameters: m0, r2")
	
	# unfix omega0 values
	param_list$omega0_group[] <- NA
	
	fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list))
	if (positive_only) {
		fit_input$lower_list$m0[] <- 0
	}
	fit_input$lower_list$r2[] <- 0
	fit_input$upper_list$r2[] <- 20
	fit_input <- limit_omega0_by_r2(fit_input)
	fit_output <- fitnmr::perform_fit(fit_input)
	
	if (plot) plot_fit_2d(fit_output, spec_ord, plot_start=TRUE, "Fit Parameters: m0, r2, omega0")
	
	refit_thresh <- c(1e-3, 1e-2)[spec_ord]*2
	
	if (any(c(fit_output$fit_list$omega0 > fit_output$upper_list$omega0 - refit_thresh, fit_output$fit_list$omega0 < fit_output$lower_list$omega0 + refit_thresh))) {
	
		#print("refitting")
		param_list <- extract_params(fit_output)
		fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list))
		if (positive_only) {
			fit_input$lower_list$m0[] <- 0
		}
		fit_input <- limit_omega0_by_r2(fit_input)
		fit_input$lower_list$r2[] <- 0
		fit_input$upper_list$r2[] <- 20
		fit_output <- fitnmr::perform_fit(fit_input)
		
		if (plot) plot_fit_2d(fit_output, spec_ord, plot_start=TRUE, "Fit Parameters: m0, r2, omega0 (Refit 1)")
	}
	
	if (any(fit_output$fit_list$omega0 > fit_output$upper_list$omega0 - refit_thresh, fit_output$fit_list$omega0 < fit_output$lower_list$omega0 + refit_thresh)) {
	
		#print("refitting2")
		param_list <- extract_params(fit_output)
		fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list))
		if (positive_only) {
			fit_input$lower_list$m0[] <- 0
		}
		fit_input$lower_list$r2[] <- 0
		fit_input$upper_list$r2[] <- 20
		fit_input <- limit_omega0_by_r2(fit_input)
		fit_output <- fitnmr::perform_fit(fit_input)
		
		if (plot) plot_fit_2d(fit_output, spec_ord, plot_start=TRUE, "Fit Parameters: m0, r2, omega0 (Refit 1)")
	}
	
	if (any(fit_output$fit_list$omega0 > fit_output$upper_list$omega0 - refit_thresh, fit_output$fit_list$omega0 < fit_output$lower_list$omega0 + refit_thresh)) {
	
		#print("refitting")
		param_list <- extract_params(fit_output)
		fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus[spec_ord]), param_list))
		if (positive_only) {
			fit_input$lower_list$m0[] <- 0
		}
		fit_input$lower_list$r2[] <- 0
		fit_input$upper_list$r2[] <- 20
		fit_input <- limit_omega0_by_r2(fit_input)
		fit_output <- fitnmr::perform_fit(fit_input)
		
		if (plot) plot_fit_2d(fit_output, spec_ord, plot_start=TRUE, "Fit Parameters: m0, r2, omega0 (Refit 1)")
	}
	
	fit_output
}

#' Fit a cluster nearby peaks starting from a seed table of chemical shifts
#'
#' @export
fit_peak_cluster <- function(spec_list, cs_start, spec_ord, f_alpha_thresh=0.001) {

	stopifnot(length(spec_list) == 1)

	cs_new <- cs_start
	fit_output <- NULL
	fit_residuals <- NULL
	footprint <- NULL
	ndf <- 0
	f_alpha <- 0
	
	while (f_alpha < f_alpha_thresh) {
	
		# perform trial fit
		trial_fit_output <- fit_peaks(spec_list, cs_new, fit_output, spec_ord=spec_ord, positive_only=TRUE)
		
		# get new data from trial fit
		trial_input_spec_int <- fitnmr::get_spec_int(trial_fit_output, "input")
		trial_fit_spec_int <- fitnmr::get_spec_int(trial_fit_output, "fit")
		trial_fit_residuals <- trial_input_spec_int[[1]]-trial_fit_spec_int[[1]]
		trial_footprint <- fit_footprint(trial_fit_output)
		
		# terminate search if any peak had zero volume
		if (any(trial_fit_output$fit_list$m0 == 0)) {
			cat("Terminating search because fit produced zero volume", sep="\n")
			if (is.null(fit_output)) {
				fit_output <- trial_input_spec_int[[1]]
			}
			break
		}
		
		# determine footprint of the peaks from the union of previous and current fits
		if (is.null(footprint)) {
			common_footprint <- trial_footprint
		} else {
			common_rows <- intersect(dimnames(fit_spec_int[[1]])[[1]], dimnames(trial_fit_spec_int[[1]])[[1]])
			common_cols <- intersect(dimnames(fit_spec_int[[1]])[[2]], dimnames(trial_fit_spec_int[[1]])[[2]])
			common_footprint <- footprint[common_rows,common_cols] | trial_footprint[common_rows,common_cols]
		}
		
		common_footprint_idx <- which(common_footprint, arr.ind=TRUE)
		common_footprint_idx <- cbind(dimnames(common_footprint)[[1]][common_footprint_idx[,1]], dimnames(common_footprint)[[2]][common_footprint_idx[,2]])
		
		if (is.null(fit_residuals)) {
		
			# initial model is having no peak (residuals = input)
			fit_residuals <- trial_input_spec_int[[1]]
			
		} else {
		
			# remove footprint locations where residuals are undefined
			common_footprint_idx <- common_footprint_idx[!is.na(fit_residuals[common_footprint_idx]) & !is.na(trial_fit_residuals[common_footprint_idx]),,drop=FALSE]
		}
		
		# calculate residual sum of squares for both models
		rss <- sum(fit_residuals[common_footprint_idx]^2)
		trial_rss <- sum(trial_fit_residuals[common_footprint_idx]^2)
		
		# calculate new degrees of freedom
		trial_ndf <- sum(sapply(trial_fit_output$group_list, function(x) length(unique(x[x!=0]))))
		
		# scale the number of points to account for zero-filling
		num_pts <- nrow(common_footprint_idx)*prod(spec_list[[1]]$fheader["TDSIZE",]/(spec_list[[1]]$fheader["FTSIZE",]/2))
		
		# calculate value of F statistic and corresponding P-value
		f_val <- ((rss-trial_rss)/(trial_ndf-ndf))/(trial_rss/(num_pts-trial_ndf))
		f_alpha <- 1-pf(f_val, trial_ndf-ndf, num_pts-trial_ndf)
		
		cat(sprintf("%2i -> %2i degrees of freedom: F = %0.1f (p = %g)", ndf, trial_ndf, f_val, f_alpha), sep="\n")
		
		if (f_alpha < f_alpha_thresh)  {
		
			fit_output <- trial_fit_output
			input_spec_int <- trial_input_spec_int
			fit_spec_int <- trial_fit_spec_int
			fit_residuals <- trial_fit_residuals
			footprint <- trial_footprint
			ndf <- trial_ndf
			
			resid_max_idx <- which(fit_residuals == max(fit_residuals[footprint]), arr.ind=TRUE)
		
			cs_new <- matrix(as.numeric(c(dimnames(footprint)[[1]][resid_max_idx[1]], dimnames(footprint)[[2]][resid_max_idx[2]])), nrow=1)[,spec_ord,drop=FALSE]
		
		} else {
		
			cat(sprintf("Terminating search because F-test p-value < %g", f_alpha_thresh), sep="\n")
			if (is.null(fit_output)) {
				fit_output <- trial_fit_spec_int[[1]]
			}
		}
	}
	
	if (is.list(fit_output))  {
		
		zlim <- range(input_spec_int[[1]], na.rm=TRUE)
	
		fitnmr::contour_pipe(aperm(input_spec_int[[1]], spec_ord), zlim=zlim, col_pos="black", col_neg="gray")
		#title(paste("Cluster", j))
		fitnmr::contour_pipe(aperm(fit_spec_int[[1]], spec_ord), zlim=zlim, col_pos="red", col_neg="pink", add=TRUE)
		#fitnmr::contour_pipe(aperm(trial_fit_spec_int[[1]], spec_ord), zlim=zlim, col_pos="purple", col_neg="plum", add=TRUE)
	
		#points(common_footprint_idx[,spec_ord], col="gray")
	
		rect(fit_output$upper_list$omega0[spec_ord[1],,1], fit_output$upper_list$omega0[spec_ord[2],,1], fit_output$lower_list$omega0[spec_ord[1],,1], fit_output$lower_list$omega0[spec_ord[2],,1], border="gray")
	
		points(t(fit_output$fit_list$omega0[spec_ord,,1]), pch=16)
		lab_text <- sprintf("%s: %.0f%%", seq_len(dim(fit_output$fit_list$omega0)[2]), 100*fit_output$fit_list$m0/sum(abs(fit_output$fit_list$m0)))
		text(t(fit_output$fit_list$omega0[spec_ord,,1]), labels=lab_text, pos=3, cex=0.6)
	
		#other_clust_idx <- peak_tab_list[[i]][,"TYPE"] == 1 & peak_tab_list[[i]][,"CLUSTID"] != j
		#points(peak_tab_list[[i]][other_clust_idx,paste(hn_name_mat[i,], "_PPM", sep=""),drop=FALSE], col="purple")
		#text(peak_tab_list[[i]][other_clust_idx,paste(hn_name_mat[i,], "_PPM", sep=""),drop=FALSE], labels=peak_tab_list[[i]][other_clust_idx,"CLUSTID"], pos=1, cex=0.6, col="purple")
	
		points(cs_start, col="green")
	}
	
	fit_output
}

#' @export
sim_time_nd <- function(peak_tab, fheader, rms=0, iseed=runif(1,max=.Machine$integer.max), file_path=NULL, verbose=FALSE) {

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

#' @export
noise_estimate <- function(x, height=TRUE, thresh=10, plot_distributions=TRUE, peak_intensities=NULL, absolute=FALSE) {

	x <- x[!is.na(x)]

	origname <- NULL
	max_data <- NULL
	
	max_data <- c(max=max(x), max_data)
	
	itermat <- matrix(nrow=20, ncol=2)
	colnames(itermat) <- c("mean", "sd")
	
	itermat[1,"mean"] <- mean(x)
	itermat[1,"sd"] <- sd(as.vector(x))
	
	for (i in 2:nrow(itermat)) {
	
		minthresh <- itermat[i-1,"mean"]-itermat[i-1,"sd"]*thresh
		maxthresh <- itermat[i-1,"mean"]+itermat[i-1,"sd"]*thresh
		idx <- which(x > minthresh & x < maxthresh)
		
		x_sub <- x[idx]
		itermat[i,"mean"] <- mean(x_sub)
		itermat[i,"sd"] <- sd(x_sub)
	}
	
	xhist <- hist(x[idx], breaks=512, plot=FALSE)
	
	normdist_height_formula <- y ~ h * exp(-(x-mu)^2/(2*sigma^2))
	normdist_formula <- y ~ 1 / (sigma*sqrt(2*pi)) * exp(-(x-mu)^2/(2*sigma^2))

	fit_data <- data.frame(y=xhist$density, x=xhist$mids)
	fit_start <- c(mu=unname(itermat[i,"mean"]), sigma=unname(itermat[i,"sd"]), h=max(xhist$density))
	
	if (height) {
		fit <- nls(normdist_height_formula, fit_data, fit_start)
	} else {
		fit <- nls(normdist_formula, fit_data, fit_start[1:2])
	}
	
	if (plot_distributions) {
		
		fit_pred <- fit$m$predict(data.frame(x=fit_data$x))
		xlim <- c(minthresh, maxthresh)
		ylim <- range(xhist$density, fit_pred)
		
		origname <- gsub("/[^/]+$", "", origname)
	
		plot(xhist$mids, xhist$density, type="n", xlim=xlim, ylim=ylim, col="black", main=origname, xlab="Signal Intensity", ylab="", yaxt="n")
		abline(h=0, col="gray")
		points(xhist$mids, xhist$density, type="l", lwd=0.75)
		points(fit_data$x, fit_pred, type="l", col="blue", lwd=0.75)
		
		legtext <- as.expression(c(
			expression(""),
			substitute(sigma: ~ sigmaval, list(sigmaval=round(fit$m$getPars()["sigma"]))),
			substitute(mu: ~ muval, list(muval=signif(fit$m$getPars()["mu"]))),
			substitute(SD: ~ sdval, list(sdval=signif(sd(x[idx])))),
			substitute(mean: ~ meanval, list(meanval=round(mean(x[idx]))))
		))
		
		if (is.null(peak_intensities)) {
			snval <- signif(max_data["max"]/abs(fit$m$getPars()["sigma"]), 2)
		} else {
			snval <- paste(round(min(peak_intensities)/abs(fit$m$getPars()["sigma"])), "-", round(max(peak_intensities)/abs(fit$m$getPars()["sigma"])))
		}
		
		legtext <- as.expression(c(
			substitute("S/N:" ~ snval, list(snval=snval)),
			substitute(sigma: ~ sigmaval, list(sigmaval=signif(abs(fit$m$getPars()["sigma"])))),
			substitute(mu: ~ muval, list(muval=signif(fit$m$getPars()["mu"])))#,
			#substitute(SD: ~ sdval, list(sdval=round(sd(x[idx])))),
			#substitute(mean: ~ meanval, list(meanval=round(mean(x[idx]))))
		))
		
		usr <- par("usr")
		pin <- par("pin")
		topright <- c(usr[2]-diff(usr[1:2])/pin[1]*0.05, usr[4]-diff(usr[3:4])/pin[2]*0.05)
		textwidth <- strwidth(legtext)
		textheight <- max(strheight(legtext))*1.3
		
		text(topright[1] - textwidth/2, topright[2]-seq(0, by=textheight, length.out=length(legtext)), legtext, pos=1, offset=0, col="blue")
		
	}
	
	c(fit$m$getPars(), max_data)
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

