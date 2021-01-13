# Input Files:
#  start_volume.csv: starting peak parameters for refitting
#  *.ft2: 2D spectra with peaks to be refit

# Output Files:
#  *_volume.csv: refit peak parameters for each *.ft2 file
#  *_fit.pdf: contour plots of the peaks fit to each *.ft2 file


# Fitting Parameters:

# data +/- these ppm values will be used for fitting (1H and X nuclei, respectively)
omega0_plus <- c(0.075, 0.75)

# omega0 values are constrained to be within this factor times R2 of the starting omega0
omega0_r2_factor <- 1.5

# R2 values are constrained to be between these two numbers
r2_bounds <- c(0.5, 20)

# scalar coupling values are constrained to be between these two numbers
sc_bounds <- c(2, 12)

# enable refitting of omega0
fit_omega0 <- TRUE

# enable refitting of R2
fit_r2 <- TRUE

# enable refitting of scalar couplings
fit_sc <- TRUE

# enable constraint of peak volume sign
preserve_m0_sign <- FALSE


# Plotting Parameters:

# lowest contour in *_fit.pdf will be this number times the noise level
plot_noise_cutoff <- 4

# scaling factor for *_fit.pdf labels
cex <- 0.2

# show omega0 constraints imposed by omega0_r2_factor
plot_omega0_bounds <- TRUE


# Computing Parameters:

# number of cores to use for refitting
mc_cores <- parallel::detectCores()


library(fitnmr)

ft2_files <- list.files(".", pattern=".ft2", full.names=TRUE, recursive=TRUE)

spec_list <- lapply(ft2_files, read_nmrpipe, dim_order="hx")
# remove ./ from spectrum labels
names(spec_list) <- sub("^[.]/", "", ft2_files)

peak_df <- read.csv("start_volume.csv", check.names=FALSE)

parallel::mclapply(seq_along(spec_list), function(spec_i) {

	output_basename <- sub(".ft2$", "_", names(spec_list)[spec_i])

	fit_input <- peak_df_to_fit_input(peak_df, spec_list[spec_i], omega0_plus=omega0_plus)
	# remember original group_list to facilitate later regeneration of peak_df
	unfixed_group_list <- fit_input$group_list
	# disable fitting peak shape parameters
	fit_input$group_list$omega0[] <- 0
	fit_input$group_list$r2[] <- 0
	fit_input$group_list$omega0_comb[] <- 0
	
	# assign bounds
	if (preserve_m0_sign) {
		pos_idx <- fit_input$start_list$m0 >= 0
		neg_idx <- fit_input$start_list$m0 < 0
		fit_input$lower_list$m0[pos_idx,] <- 0
		fit_input$upper_list$m0[neg_idx,] <- 0
	}
	
	cat(paste("Fitting m0 for spectrum ", names(spec_list)[spec_i], "...", sep=""), sep="\n")
	fit_output <- perform_fit(fit_input)
	fit_output$group_list <- unfixed_group_list
	
	peak_df <- param_list_to_peak_df(fit_output)
	
	if (fit_omega0 || fit_r2 || fit_sc) {
	
		fit_input <- peak_df_to_fit_input(peak_df, spec_list[spec_i], omega0_plus=omega0_plus)
		fit_input <- update_fit_bounds(fit_input, omega0_r2_factor=omega0_r2_factor, r2_bounds=r2_bounds, sc_bounds=sc_bounds)
		# remember original group_list to facilitate later regeneration of peak_df
		unfixed_group_list <- fit_input$group_list
		
		fit_params <- character()
		
		if (fit_omega0) {
			fit_params <- c(fit_params, "omega0")
		} else {
			# determine which parameters have omega0 values
			omega0_idx <- omega0_param_idx(fit_input)
			param_values(fit_input$group_list, omega0_idx) <- 0
		}
		
		if (fit_sc) {
			fit_params <- c(fit_params, "sc")
		} else {
		    # determine which parameters have scalar coupling values
			sc_idx <- coupling_param_idx(fit_input)
			param_values(fit_input$group_list, sc_idx) <- 0
		}
		
		if (fit_r2) {
			fit_params <- c(fit_params, "r2")
		} else {
			fit_input$group_list$r2[] <- 0
		}
		
		# reassign bounds
		if (preserve_m0_sign) {
			fit_input$lower_list$m0[pos_idx,] <- 0
			fit_input$upper_list$m0[neg_idx,] <- 0
		}
		
		fit_params <- c(fit_params, "m0")
	
		cat(paste("Fitting ", paste(fit_params, collapse=","), " for spectrum ", names(spec_list)[spec_i], "...", sep=""), sep="\n")
		fit_output <- perform_fit(fit_input)
		fit_output$group_list <- unfixed_group_list
	
		peak_df <- param_list_to_peak_df(fit_output)
	}
	
	
	write.csv(peak_df, paste(output_basename, "volume.csv", sep=""), row.names=FALSE)
	
	
	pdf(paste(output_basename, "fit.pdf", sep=""), height=10, width=10, pointsize=12)

	par(mar=c(3, 3, 1.5, 1), mgp=c(2, 0.8, 0))

	plot_peak_df(
		peak_df,
		spec_list[spec_i],
		noise_cutoff=plot_noise_cutoff,
		cex=cex,
	)
	
	# plot omega0 bounds if requested
	if (plot_omega0_bounds) {
	
		# determine which parameters have omega0 values in each dimension
		omega0_1_idx <- omega0_param_idx(fit_output, 1, specs=1)
		omega0_2_idx <- omega0_param_idx(fit_output, 2, specs=1)
		
		# plot values in blue and bounds with rectangles
		points(param_values(fit_output$fit_list, omega0_1_idx), param_values(fit_output$fit_list, omega0_2_idx), pch=16, col=rgb(0, 0, 1, 1), cex=0.5*cex)
		rect(param_values(fit_output$upper_list, omega0_1_idx), param_values(fit_output$upper_list, omega0_2_idx), param_values(fit_output$lower_list, omega0_1_idx), param_values(fit_output$lower_list, omega0_2_idx), border=gray(0, 0.25), lwd=cex)
	}

	dev.off()

}, mc.cores=mc_cores)

