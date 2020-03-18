# Assignment Parameters:

# peak assignment stage 1 adjustment (for first and second dimensions, respectively)
cs_adj_1 <- c(0, 0)

# peak assignment stage 1 distance threshold (fraction of chemical shift range)
thresh_1 <- 0.1

# peak assignment stage 1 adjustment (for first and second dimensions, respectively)
cs_adj_2 <- c(0, 0)

# peak assignment stage 2 distance threshold (fraction of chemical shift range)
thresh_2 <- 0.025


# Plotting Parameters:

# lowest contour in *_fit.pdf will be this number times the noise level
plot_noise_cutoff <- 4

# scaling factor for fit_spectra.pdf labels
cex <- 0.2

# plot ovals showing the thresholds for peak assignment transfer
plot_thresh <- TRUE

# color for unmatched peaks/assignments
unmatched_col <- "purple"


library("fitnmr")

ovals <- function(x, y=NULL, x_radius, y_radius, border=1, lwd=1) {

	coords <- xy.coords(x, y)
	
	angles <- seq(0, 2*pi, length.out=60)
	
	oval_x <- sin(angles)
	oval_y <- cos(angles)
	
	for (i in seq_along(coords$x)) {
		polygon(oval_x*x_radius+coords$x[i], oval_y*y_radius+coords$y[i], border=border, lwd=lwd)
	}
}

# load the assignments
assignments_df <- read.csv("assignments.csv", check.names=FALSE)

# find *_volume.csv files that have corresponding *.ft2 files
peak_files <- list.files(".", pattern="_volume.csv", full.names=TRUE, recursive=TRUE)
peak_files <- peak_files[file.exists(sub("_volume.csv", ".ft2", peak_files))]

# get the *.ft2 filenames
ft2_files <- list.files(".", pattern=".ft2", full.names=TRUE, recursive=TRUE)

# read peak lists into a multidimensional array
peak_df_list <- lapply(peak_files, function(x) read.csv(x, check.names=FALSE))
peak_array <- simplify2array(lapply(peak_df_list, as.matrix))

# determine the column numbers for the chemical shifts and volume
omega0_ppm_1_idx <- match("omega0_ppm_1", dimnames(peak_array)[[2]])
omega0_ppm_2_idx <- match("omega0_ppm_2", dimnames(peak_array)[[2]])
peak_vol_idx <- match("r2_hz_2", dimnames(peak_array)[[2]])+1

# read the starting volume
start_df <- read.csv("start_volume.csv", check.names=FALSE)
start_vol_idx <- match("r2_hz_2", colnames(start_df))+1

# determine the peak list with the lowest root mean square volume compared with start
vol_rms <- colMeans((peak_array[,peak_vol_idx,]-start_df[,start_vol_idx])^2)
min_rms_idx <- which.min(vol_rms)

# read in the spectrum for that peak list
spec_list <- lapply(ft2_files[min_rms_idx], read_nmrpipe, dim_order="hx")
# remove ./ from spectrum labels
names(spec_list) <- sub("^[.]/", "", ft2_files[min_rms_idx])

# create table with unknown peak ppm values and volume
unknown_tab <- peak_df_list[[min_rms_idx]][,c(omega0_ppm_1_idx, omega0_ppm_2_idx, peak_vol_idx)]

# create table with previously known assignments
assigned_tab <- t(t(assignments_df[,2:3])+cs_adj_1)

# calculate initial assignments
assign_idx <- height_assign(assigned_tab, unknown_tab, thresh=thresh_1)

# create adjusted assigned chemical shift table
assign_adj <- apply(unknown_tab[assign_idx,1:2]-assigned_tab, 2, median, na.rm=TRUE)
assigned_tab_adj <- t(t(assigned_tab)+assign_adj+cs_adj_2)

# calculate adjusted assignments
assign_idx_adj <- height_assign(assigned_tab_adj, unknown_tab, thresh=thresh_2)


median_tab <- apply(peak_array[,-peak_vol_idx,], 1:2, median)
volume_tab <- peak_array[,peak_vol_idx,]
colnames(volume_tab) <- sub("^[.]/", "", ft2_files)

assignment_labels <- assignments_df[match(seq_len(nrow(median_tab)), assign_idx_adj),1]

assign_tab <- data.frame(assignment=assignment_labels, median_tab, volume_tab, check.names=FALSE)

write.csv(assign_tab, "assign_volume.csv", row.names=FALSE, na="")


pdf("assign_stages.pdf", width=10, height=10, pointsize=12)

par(mar=c(3, 3, 1.5, 1), mgp=c(2, 0.8, 0))

plot_peak_df(
	peak_df_list[[min_rms_idx]],
	spec_list,
	noise_cutoff=plot_noise_cutoff,
	cex=cex,
	label_col=ifelse(seq_len(nrow(peak_df_list[[min_rms_idx]])) %in% assign_idx_adj, "black", unmatched_col)
)

# plot initial assignments
segments(assigned_tab[,1], assigned_tab[,2], unknown_tab[assign_idx,1], unknown_tab[assign_idx,2], lwd=1, col="lightgray")
points(assigned_tab, pch=16, cex=2.5*cex, col="lightgray")
if (plot_thresh) {
	ovals(assigned_tab, x_radius=attr(assign_idx, "thresh")[1], y_radius=attr(assign_idx, "thresh")[2], border="lightgray", lwd=cex*2)
}
text(assigned_tab[,1], assigned_tab[,2], assignments_df[,1], cex=cex, col="darkgray")

# plot adjusted assignments
segments(assigned_tab_adj[,1], assigned_tab_adj[,2], unknown_tab[assign_idx_adj,1], unknown_tab[assign_idx_adj,2], lwd=1, col="green")
points(assigned_tab_adj, pch=16, cex=2.5*cex, col="green")
if (plot_thresh) {
	ovals(assigned_tab_adj, x_radius=attr(assign_idx_adj, "thresh")[1], y_radius=attr(assign_idx_adj, "thresh")[2], border="green", lwd=cex*2)
}
text(assigned_tab_adj[,1], assigned_tab_adj[,2], assignments_df[,1], cex=cex, col=ifelse(is.na(assign_idx_adj), unmatched_col, "black"))

# add legends
legend("topleft", legend=c("Stage 1", "Stage 2", "Stage 2 Unmatched"), pch=16, col=c("lightgray", "green", "green"), text.col=c("darkgray", "black", unmatched_col), bty="n", cex=cex)

legend("bottomleft", legend=c("Peak:Fit", "Peak:Fit Unmatched"), text.col=c("black", unmatched_col), bty="n", cex=cex, x.intersp=0)

dev.off()


peak_cols <- rep("black", nrow(peak_df_list[[min_rms_idx]]))
for (fit_num in unique(peak_df_list[[min_rms_idx]][,"fit"])) {
	idx <- which(peak_df_list[[min_rms_idx]][,"fit"] == fit_num)
	peak_cols[idx] <- rep_len(c("green3", "cyan", "magenta", "yellow", "gray", "purple"), length(idx))
}

pdf("assign_omegas.pdf", width=10, height=10, pointsize=12)

par(mar=c(3, 3, 1.5, 1), mgp=c(2, 0.8, 0))

plot_peak_df(
	peak_df_list[[min_rms_idx]],
	spec_list,
	noise_cutoff=plot_noise_cutoff,
	cex=cex,
	label_col=ifelse(seq_len(nrow(peak_df_list[[min_rms_idx]])) %in% assign_idx_adj, "black", unmatched_col)
)

for (i in seq_along(peak_df_list)) {
	points(peak_df_list[[i]][,c(omega0_ppm_1_idx, omega0_ppm_2_idx)], pch=16, col=peak_cols, cex=0.25*cex)
}

text(median_tab[,c("omega0_ppm_1","omega0_ppm_2")], labels=assignment_labels, pos=3, offset=cex*0.5, cex=cex)

dev.off()
