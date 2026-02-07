if (!requireNamespace("HamiltonianMultiplet", quietly = TRUE)) {
	source(file.path("..", "..", "R", "HamiltonianMultiplet.R"))
} else {
	library(HamiltonianMultiplet)
}

test_that("analytical derivatives match finite differences", {
	n <- 2
	shift_labels <- c("Ha", "Hb")
	coupling_labels <- matrix(NA_character_, n, n)
	coupling_labels[1, 2] <- "Jab"
	coupling_labels[2, 1] <- "Jab"
	label_matrix <- coupling_labels
	diag(label_matrix) <- shift_labels

	params <- c(Ha = -10, Hb = 20, Jab = 7.0)
	dh <- HamiltonianMultiplet$new(label_matrix, params)

	base <- dh$multiplet("I1x", coherence = NULL)
	eps <- 1e-6

	for (lbl in names(params)) {
		params_eps <- params
		params_eps[[lbl]] <- params_eps[[lbl]] + eps
		dh$set_params(params_eps)
		bumped <- dh$multiplet("I1x", coherence = NULL)

		df_num <- (bumped$frequency - base$frequency) / eps
		dI_num <- (bumped$intensity - base$intensity) / eps

		df_ana <- base[[paste0("dFreq_", lbl)]]
		dI_ana <- base[[paste0("dInt_", lbl)]]

		expect_equal(df_ana, df_num, tolerance = 1e-4)
		expect_equal(dI_ana, dI_num, tolerance = 1e-4)
	}

	# reset params
	dh$set_params(params)
})

test_that("frequency shifts are invariant to constant offsets", {
	n <- 2
	shift_labels <- c("Ha", "Hb")
	coupling_labels <- matrix(NA_character_, n, n)
	coupling_labels[1, 2] <- "Jab"
	coupling_labels[2, 1] <- "Jab"
	label_matrix <- coupling_labels
	diag(label_matrix) <- shift_labels

	params <- c(Ha = -10, Hb = 20, Jab = 7.0)
	dh <- HamiltonianMultiplet$new(label_matrix, params)
	base <- dh$multiplet("I1x", coherence = 1)
	base_freqs <- sort(base$frequency)

	offsets <- c(-60, -30, -10, 0, 10, 30, 60)
	signs <- logical(length(offsets))
	for (i in seq_along(offsets)) {
		params2 <- params
		params2[c("Ha", "Hb")] <- params2[c("Ha", "Hb")] + offsets[[i]]
		dh$set_params(params2)
		m <- dh$multiplet("I1x", coherence = 1)
		cur_freqs <- sort(m$frequency)
		delta <- cur_freqs - base_freqs
		expect_lt(max(abs(delta - mean(delta))), 1e-8)
		signs[i] <- all(cur_freqs > 0) || all(cur_freqs < 0)
	}
	expect_true(any(signs))
})

test_that("covers Ha=Hb, Δν=J, Δν<J, and Δν>J cases", {
	n <- 2
	shift_labels <- c("Ha", "Hb")
	coupling_labels <- matrix(NA_character_, n, n)
	coupling_labels[1, 2] <- "Jab"
	coupling_labels[2, 1] <- "Jab"
	label_matrix <- coupling_labels
	diag(label_matrix) <- shift_labels

	cases <- list(
		c(Ha = 10, Hb = 10, Jab = 7.0),	# Ha=Hb
		c(Ha = 10, Hb = 17, Jab = 7.0),	# Δν=J
		c(Ha = 10, Hb = 15, Jab = 7.0),	# Δν<J
		c(Ha = 10, Hb = 30, Jab = 7.0)	 # Δν>J
	)

	for (params in cases) {
		dh <- HamiltonianMultiplet$new(label_matrix, params)
		m <- dh$multiplet("I1x", coherence = 1)
		expect_true(nrow(m) > 0)
		expect_true(is.numeric(m$frequency))
	}
})
