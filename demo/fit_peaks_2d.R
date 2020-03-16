# Fitting Parameters:

# peak height must be at least this times the noise level for a new fitting iteration
noise_cutoff <- 15

# F-test p-value must be less than this value to accept the addition of a new peak
f_alpha <- 0.001

# maximum number of peak fitting iterations to run
iter_max <- 100

# data +/- these ppm values will be used for fitting (1H and X nuclei, respectively)
omega0_plus <- c(0.075, 0.75)

# starting R2 value for the fits
r2_start <- 5

# R2 values are constrained to be between these two numbers
r2_bounds <- c(0.5, 20)

# starting doublet scalar coupling (1H and X nuclei, respectively), NA for singlet
sc_start <- c(6, NA)

# scalar coupling values are constrained to be between these two numbers
sc_bounds <- c(2, 12)


# Plotting Parameters:

# lowest contour in fit_spectra.pdf will be this number times the noise level
plot_noise_cutoff <- 4

# scaling factor for fit_spectra.pdf labels
cex <- 0.4


library(fitnmr)

ft2_files <- list.files(".", pattern=".ft2", full.names=TRUE, recursive=TRUE)

spec_list <- lapply(ft2_files, read_nmrpipe, dim_order="hx")
# remove ./ from spectrum labels
names(spec_list) <- sub("^[.]/", "", ft2_files)


pdf("noise_histograms.pdf", height=3, width=4, pointsize=12)

par(mar=c(3, 1, 1.5, 1), mgp=c(2, 0.8, 0))

noise_mat <- sapply(seq_along(spec_list), function(i) {
	
	noise_data <- noise_estimate(spec_list[[i]]$int)
	title(main=names(spec_list)[i])
	noise_data
})

dev.off()


pdf("fit_iterations.pdf", height=6, width=6, pointsize=12)

par(mar=c(3, 3, 1.5, 1), mgp=c(2, 0.8, 0))

peak_fits <- fit_peak_iter(
	spec_list,
	noise_sigma=noise_mat["sigma",],
	noise_cutoff=noise_cutoff,
	f_alpha=f_alpha,
	omega0_plus=omega0_plus,
	iter_max=iter_max,
	r2_start=r2_start,
	r2_bounds=r2_bounds,
	sc_start=sc_start,
	sc_bounds=sc_bounds,
	plot_fit=TRUE
)

dev.off()


peak_df <- param_list_to_peak_df(peak_fits)

write.csv(peak_df, "fit_volume.csv", row.names=FALSE)


pdf("fit_spectra.pdf", height=10, width=10, pointsize=12)

par(mar=c(3, 3, 1.5, 1), mgp=c(2, 0.8, 0))

plot_peak_df(
	peak_df,
	spec_list,
	noise_sigma=noise_mat["sigma",],
	noise_cutoff=plot_noise_cutoff,
	cex=cex,
)

dev.off()
