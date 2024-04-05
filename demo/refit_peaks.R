# Input Files:
#  start_resonances.csv: input resonances table
#  start_nuclei.csv: input nuclei table
#  start_couplings.csv: input couplings table
#  *.ft*: spectra with peaks to be refit

# Output Files:
#  resonances.csv: output resonances table
#  nuclei.csv: output resonances table
#  couplings.csv: output resonances table


# dimension order
dim_order <- NULL

# Fitting Parameters:

# data +/- these ppm values will be used for fitting (X, Y, Z, A nuclei, respectively)
omega0_plus <- c(0.075, 0.75, 0.075)

# omega0 values are constrained to be within this factor times R2 of the starting omega0
omega0_r2_factor <- 1.5

# R2 values are constrained to be between these two numbers
r2_bounds <- c(0.05, 20)

# scalar coupling values are constrained to be between these two numbers
sc_bounds <- c(2, 12)

# enable refitting of omega0
fit_omega0 <- TRUE

# enable refitting of R2
fit_r2 <- TRUE

# enable refitting of scalar couplings
fit_sc <- TRUE

# enable refitting of phases
fit_phases <- FALSE

# enable constraint of peak volume sign
preserve_m0_sign <- FALSE

library(fitnmr)

ft_files <- list.files(".", pattern="[.]ft[1-4]?$", full.names=TRUE, recursive=TRUE)

if (!"spec_list" %in% ls()) {
spec_list <- lapply(ft_files, read_nmrpipe, dim_order=dim_order)
# remove ./ from spectrum labels
names(spec_list) <- sub("^[.]/", "", ft_files)
}
omega0_plus <- omega0_plus[seq_len(ncol(spec_list[[1]]$fheader))]

start_resonances <- read.csv(text=readLines("start_resonances.csv", warn=FALSE), row.names=1, check.names=FALSE)
start_resonances[,"x_sc"] <- as.character(start_resonances[,"x_sc"])
start_resonances[is.na(start_resonances[,"x_sc"]),"x_sc"] <- ""
start_nuclei <- read.csv(text=readLines("start_nuclei.csv", warn=FALSE), row.names=1, check.names=FALSE)
start_couplings <- read.csv(text=readLines("start_couplings.csv", warn=FALSE), row.names=1, check.names=FALSE)

start_tables <- list(resonances=start_resonances, nuclei=start_nuclei, couplings=start_couplings)

# make sure names of nuclei and couplings are all unique
stopifnot(length(unique(c(rownames(start_nuclei), rownames(start_couplings)))) == nrow(start_nuclei) + nrow(start_couplings))

# create list of fitting parameters from input tables
#param_list <- fitnmr:::tables_to_param_list_old(spec_list, start_tables)
param_list <- tables_to_param_list(spec_list, start_tables)

# first fit only m0
fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus), param_list_to_arg_list(param_list)))
# set NA starting volumes to maximum spectrum intensity
fit_input$start_list$m0[is.na(fit_input$start_list$m0)] <- max(sapply(fit_input$spec_data, function(x) max(abs(x$spec_int))))
# turn off omega0, r2, and scalar coupling fitting at first
fit_input$group_list$omega0[] <- 0
fit_input$group_list$r2[] <- 0
fit_input$group_list$omega0_comb[] <- 0
# run the fit
fit_output <- perform_fit(fit_input)

# next fit other variables
fit_input <- do.call(fitnmr::make_fit_input, c(list(spec_list, omega0_plus=omega0_plus), param_list_to_arg_list(param_list)))
fit_input$start_list$m0 <- fit_output$fit_list$m0

# turn off omega0, r2, and scalar coupling fitting if requested
if (!fit_omega0) {
	param_values(fit_input$group_list, omega0_param_idx(fit_input)) <- 0
}
if (!fit_r2) {
	fit_input$group_list$r2[] <- 0
}
if (!fit_sc) {
	param_values(fit_input$group_list, coupling_param_idx(fit_input)) <- 0
}
if (fit_phases) {
	for (j in seq_len(dim(fit_input$group_list$p0)[3])) {
		for (i in seq_len(dim(fit_input$group_list$p0)[1])) {
			fit_input$group_list$p0[i,,j] <- i+(j-1)*dim(fit_input$group_list$p0)[3]
			fit_input$group_list$p1[i,,j] <- i+(j-1)*dim(fit_input$group_list$p0)[3]
		}
	}
}
# apply bounds
fit_input <- update_fit_bounds(fit_input, omega0_r2_factor=omega0_r2_factor, r2_bounds=r2_bounds, sc_bounds=sc_bounds)
# run the fit
fit_output <- perform_fit(fit_input)

tables <- param_list_to_tables(fit_output, start_tables)

write.csv(tables[["resonances"]], "resonances.csv", quote=FALSE)
write.csv(tables[["nuclei"]], "nuclei.csv", quote=FALSE)
write.csv(tables[["couplings"]], "couplings.csv", quote=FALSE)

# plot output
if (ncol(spec_list[[1]]$fheader) == 1) {
	pdf("sparse_1d.pdf", width=10, height=4)
	par(mar=c(5.1, 2.1, 1.1, 2.1))
	plot_sparse_1d(fit_output, tables)
	dev.off()
	pdf("resonances_1d.pdf")
	plot_resonances_1d(fit_output, omega0_plus=omega0_plus, always_show_start=FALSE)
	dev.off()
} else if (ncol(spec_list[[1]]$fheader) == 2) {
	pdf("resonances_2d.pdf")
	plot_resonances_2d(fit_output, omega0_plus=omega0_plus, low_frac=0.01)
	dev.off()
} else if (ncol(spec_list[[1]]$fheader) == 3) {
	pdf("resonances_2d.pdf")
	plot_resonances_3d(fit_output, omega0_plus=omega0_plus)
	dev.off()
}
