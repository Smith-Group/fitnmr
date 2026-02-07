if (!requireNamespace("HamiltonianMultiplet", quietly = TRUE)) {
	source(file.path("..", "..", "R", "HamiltonianMultiplet.R"))
} else {
	library(HamiltonianMultiplet)
}

test_that("label ordinals resolve correctly", {
	n <- 3
	shift_labels <- c("Ha", "Ha", "Hb")
	coupling_labels <- matrix(NA_character_, n, n)
	coupling_labels[1, 2] <- "Jaa"
	coupling_labels[2, 1] <- "Jaa"
	coupling_labels[1, 3] <- "Jab"
	coupling_labels[3, 1] <- "Jab"
	coupling_labels[2, 3] <- "Jab"
	coupling_labels[3, 2] <- "Jab"
	label_matrix <- coupling_labels
	diag(label_matrix) <- shift_labels

	params <- c(Ha = -10, Hb = 20, Jaa = 2, Jab = 7)
	dh <- HamiltonianMultiplet$new(label_matrix, params)

	expect_equal(dh$operator("Ha"), dh$operator("Ha1"))
	expect_equal(dh$operator("Ha2"), dh$operator("Ha2"))
	expect_error(dh$operator("Ha3"))
})

test_that("anti-phase shorthand matches explicit operators", {
	n <- 2
	shift_labels <- c("Ha", "Hb")
	coupling_labels <- matrix(NA_character_, n, n)
	coupling_labels[1, 2] <- "Jab"
	coupling_labels[2, 1] <- "Jab"
	label_matrix <- coupling_labels
	diag(label_matrix) <- shift_labels
	params <- c(Ha = -10, Hb = 20, Jab = 10)
	dh <- HamiltonianMultiplet$new(label_matrix, params)

	m1 <- dh$multiplet("Ha~Hb")
	m2 <- dh$multiplet("2I1xI2z")
	expect_equal(m1$frequency, m2$frequency, tolerance = 1e-10)
	expect_equal(m1$intensity, m2$intensity, tolerance = 1e-10)

	m3 <- dh$multiplet("Hb~Ha")
	m4 <- dh$multiplet("2I2xI1z")
	expect_equal(m3$frequency, m4$frequency, tolerance = 1e-10)
	expect_equal(m3$intensity, m4$intensity, tolerance = 1e-10)
})

test_that("multi-partner anti-phase shorthand matches explicit", {
	n <- 3
	shift_labels <- c("Ha", "Hb", "Hc")
	coupling_labels <- matrix(NA_character_, n, n)
	coupling_labels[1, 2] <- "Jab"
	coupling_labels[2, 1] <- "Jab"
	coupling_labels[1, 3] <- "Jac"
	coupling_labels[3, 1] <- "Jac"
	coupling_labels[2, 3] <- "Jbc"
	coupling_labels[3, 2] <- "Jbc"
	label_matrix <- coupling_labels
	diag(label_matrix) <- shift_labels
	params <- c(Ha = -10, Hb = 20, Hc = 5, Jab = 10, Jac = 3, Jbc = 4)
	dh <- HamiltonianMultiplet$new(label_matrix, params)

	m1 <- dh$multiplet("Ha~Hb~Hc")
	m2 <- dh$multiplet("4I1xI2zI3z")
	expect_equal(m1$frequency, m2$frequency, tolerance = 1e-10)
	expect_equal(m1$intensity, m2$intensity, tolerance = 1e-10)
})

test_that("anti-phase shorthand handles sums and errors", {
	n <- 2
	shift_labels <- c("Ha", "Hb")
	coupling_labels <- matrix(NA_character_, n, n)
	coupling_labels[1, 2] <- "Jab"
	coupling_labels[2, 1] <- "Jab"
	label_matrix <- coupling_labels
	diag(label_matrix) <- shift_labels
	params <- c(Ha = -10, Hb = 20, Jab = 10)
	dh <- HamiltonianMultiplet$new(label_matrix, params)

	m <- dh$multiplet("Ha~Hb+Hb")
	expect_true(nrow(m) > 0)

	expect_error(dh$multiplet("Ha~"))
	expect_error(dh$multiplet("~Hb"))
	expect_error(dh$multiplet("Ha~Hx"))
})
