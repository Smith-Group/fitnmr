#' HamiltonianMultiplet R6 class
#'
#' Build and differentiate NMR Hamiltonians for fixed-size spin systems.
#'
#' @export
HamiltonianMultiplet <- R6::R6Class(
	"HamiltonianMultiplet",
	public = list(
		#' @field n Number of spins in the system.
		n = NULL,
		#' @field shift_labels Character vector of shift labels (diagonal of label_matrix).
		shift_labels = NULL,
		#' @field coupling_labels Matrix of coupling labels (off-diagonal of label_matrix).
		coupling_labels = NULL,
		#' @field params Named numeric vector of parameter values.
		params = NULL,
		#' @field ops List of spin operators.
		ops = NULL,
		#' @field dH List of Hamiltonian derivative matrices by parameter label.
		dH = NULL,
		#' @field cache Internal cache of eigendecomposition results.
		cache = NULL,
		#' @description
		#' Create a HamiltonianMultiplet from labels and parameters.
		#' @param label_matrix Square matrix or data.frame of shift (diagonal) and coupling labels.
		#' @param params Named numeric vector or list of parameter values.
		#' @return A new HamiltonianMultiplet object.
		initialize = function(label_matrix, params = NULL) {
			parsed <- private$parse_label_matrix(label_matrix)
			self$n <- parsed$n
			self$shift_labels <- parsed$shift_labels
			self$coupling_labels <- parsed$coupling_labels
			if (is.null(params)) {
				params <- setNames(rep(0, length(private$all_labels())), private$all_labels())
			}
			self$params <- private$normalize_params(params)
			self$ops <- private$build_spin_ops(self$n)
			self$dH <- private$build_dH()
			self$cache <- list(dirty = TRUE, E = NULL, V = NULL, mvals = NULL, spin_index = NULL, spin_label = NULL)
		},
		#' @description
		#' Update parameter values and mark the cache dirty if changed.
		#' @param params Named numeric vector or list of parameter values.
		set_params = function(params) {
			params <- private$normalize_params(params)
			old <- self$params[names(params)]
			self$params[names(params)] <- params
			if (any(old != params)) {
				self$cache$dirty <- TRUE
			}
		},
		#' @description
		#' Return sorted unique parameter labels.
		get_param_labels = function() {
			sort(unique(c(self$shift_labels, self$coupling_labels[upper.tri(self$coupling_labels)])))
		},
		#' @description
		#' Return the Hamiltonian matrix for current parameters.
		hamiltonian = function() {
			private$build_H(self$params)
		},
		#' @description
		#' Return the operator matrix for a label expression.
		#' @param label Operator label string.
		operator = function(label) {
			private$parse_operator_sum(label)
		},
		#' @description
		#' Return a data.frame of transitions, intensities, and derivatives.
		#' @param state_label Operator label string for the prepared state.
		#' @param detect_label Operator label string for detection. Defaults to \code{"allx"}.
		#' @param intensity_tol Minimum absolute intensity to keep in the output.
		#' @param coherence Allowed coherence order(s) (integer vector) or NULL for all.
		#' @param coherence_tol Tolerance for matching coherence orders.
		#' @param normalize Logical; if TRUE, normalize intensities.
		multiplet = function(state_label,
		                     detect_label = "allx",
		                     intensity_tol = 0.1,
		                     coherence = 1,
		                     coherence_tol = 1e-6,
		                     normalize = TRUE) {
			if (!is.null(coherence)) {
				if (!is.numeric(coherence) || length(coherence) < 1) {
					stop("coherence must be a numeric vector of integers or NULL")
				}
				if (any(abs(coherence - round(coherence)) > 1e-6)) {
					stop("coherence values must be integers")
				}
			}
			cache <- private$get_eigendecomp()
			E <- cache$E
			V <- cache$V
			O <- private$parse_operator_sum(state_label)
			D <- private$parse_operator_sum(detect_label)
			Oe <- t(Conj(V)) %*% O %*% V
			De <- t(Conj(V)) %*% D %*% V

			transitions <- private$transition_index(
				length(E),
				mode = if (is.null(coherence)) "upper" else "all"
			)
			freqs <- E[transitions$l] - E[transitions$k]
			idx <- cbind(transitions$k, transitions$l)
			intens <- Re(De[idx] * Conj(Oe[idx]))

			if (!is.null(coherence)) {
				dq <- cache$mvals[transitions$l] - cache$mvals[transitions$k]
				keep <- Reduce(`|`, lapply(coherence, function(c) abs(dq - c) < coherence_tol))
				transitions <- transitions[keep, , drop = FALSE]
				freqs <- freqs[keep]
				idx <- idx[keep, , drop = FALSE]
				intens <- intens[keep]
			}
			derivs <- private$compute_derivatives(E, V, Oe, De, transitions)

			if (isTRUE(normalize)) {
				scale <- 2^(3 - self$n)
				intens <- intens * scale
				int_cols <- grep("^dInt_", colnames(derivs))
				derivs[, int_cols] <- derivs[, int_cols, drop = FALSE] * scale
			}

			out <- data.frame(
				k = transitions$k,
				l = transitions$l,
				frequency = Re(freqs),
				intensity = Re(intens),
				spin_index = cache$spin_index[idx],
				spin_label = cache$spin_label[idx],
				derivs,
				check.names = FALSE
			)
			if (intensity_tol > 0) {
				out <- out[abs(out$intensity) > intensity_tol, , drop = FALSE]
			} else {
				out <- out[out$intensity != 0, , drop = FALSE]
			}
			rownames(out) <- NULL
			out
		}
	),
	private = list(
		normalize_params = function(params) {
			if (is.null(names(params))) {
				stop("params must be a named numeric vector or list")
			}
			params <- unlist(params)
			labels <- private$all_labels()
			missing_labels <- setdiff(labels, names(params))
			if (length(missing_labels) > 0) {
				stop("Missing params for labels: ", paste(missing_labels, collapse = ", "))
			}
			params[labels]
		},
		all_labels = function() {
			labels <- c(self$shift_labels, self$coupling_labels[upper.tri(self$coupling_labels)])
			labels <- labels[!is.na(labels) & nzchar(labels)]
			sort(unique(labels))
		},
		build_spin_ops = function(n) {
			Ix <- matrix(c(0, 1, 1, 0), 2, 2) / 2
			Iy <- matrix(c(0, -1i, 1i, 0), 2, 2) / 2
			Iz <- matrix(c(1, 0, 0, -1), 2, 2) / 2
			Id <- diag(2)

			ops <- list(Ix = vector("list", n), Iy = vector("list", n), Iz = vector("list", n))
			for (i in seq_len(n)) {
				ops$Ix[[i]] <- private$op_on_spin(Ix, i, n, Id)
				ops$Iy[[i]] <- private$op_on_spin(Iy, i, n, Id)
				ops$Iz[[i]] <- private$op_on_spin(Iz, i, n, Id)
			}
			ops
		},
		parse_label_matrix = function(label_matrix) {
			if (is.data.frame(label_matrix)) {
				label_matrix <- as.matrix(label_matrix)
			}
			if (!is.matrix(label_matrix)) {
				stop("label_matrix must be a matrix or data.frame")
			}
			if (nrow(label_matrix) != ncol(label_matrix)) {
				stop("label_matrix must be square")
			}
			n <- nrow(label_matrix)
			shift_labels <- as.character(diag(label_matrix))
			if (any(is.na(shift_labels) | !nzchar(shift_labels))) {
				stop("All diagonal entries of label_matrix must be non-empty shift labels")
			}
			coupling_labels <- matrix(NA_character_, n, n)
			for (i in seq_len(n)) {
				for (j in seq_len(n)) {
					if (i == j) next
					a <- label_matrix[i, j]
					b <- label_matrix[j, i]
					a <- if (is.na(a)) "" else as.character(a)
					b <- if (is.na(b)) "" else as.character(b)
					a <- gsub("^\\s+|\\s+$", "", a)
					b <- gsub("^\\s+|\\s+$", "", b)
					if (nzchar(a) && nzchar(b) && a != b) {
						stop("Conflicting coupling labels at (", i, ",", j, ") and (", j, ",", i, ")")
					}
					lbl <- if (nzchar(a)) a else if (nzchar(b)) b else NA_character_
					coupling_labels[i, j] <- lbl
				}
			}
			list(n = n, shift_labels = shift_labels, coupling_labels = coupling_labels)
		},
		op_on_spin = function(op, idx, n, Id) {
			mats <- vector("list", n)
			for (i in seq_len(n)) {
				mats[[i]] <- if (i == idx) op else Id
			}
			Reduce(kronecker, mats)
		},
		build_dH = function() {
			labels <- private$all_labels()
			dimH <- 2^self$n
			dH <- lapply(labels, function(x) matrix(0, dimH, dimH))
			names(dH) <- labels

			for (i in seq_len(self$n)) {
				lbl <- self$shift_labels[[i]]
				dH[[lbl]] <- dH[[lbl]] + self$ops$Iz[[i]]
			}

			for (i in seq_len(self$n - 1)) {
				for (j in (i + 1):self$n) {
					lbl <- self$coupling_labels[i, j]
					if (!is.na(lbl) && nzchar(lbl)) {
						term <- self$ops$Ix[[i]] %*% self$ops$Ix[[j]] +
							self$ops$Iy[[i]] %*% self$ops$Iy[[j]] +
							self$ops$Iz[[i]] %*% self$ops$Iz[[j]]
						dH[[lbl]] <- dH[[lbl]] + term
					}
				}
			}
			dH
		},
		build_H = function(params) {
			H <- matrix(0+0i, 2^self$n, 2^self$n)
			for (lbl in names(self$dH)) {
				H <- H + params[[lbl]] * self$dH[[lbl]]
			}
			Re(H)
		},
		get_eigendecomp = function() {
			if (isTRUE(self$cache$dirty) || is.null(self$cache$E)) {
				H <- private$build_H(self$params)
				eig <- eigen(H, symmetric = TRUE)
				Iz_total <- Reduce(`+`, self$ops$Iz)
				Iz_e <- t(Conj(eig$vectors)) %*% Iz_total %*% eig$vectors
				self$cache$E <- eig$values
				self$cache$V <- eig$vectors
				self$cache$mvals <- Re(diag(Iz_e))
				self$cache$spin_index <- private$compute_transition_spin_index(eig$vectors)
				self$cache$spin_label <- matrix(
					self$shift_labels[self$cache$spin_index],
					nrow = nrow(self$cache$spin_index),
					ncol = ncol(self$cache$spin_index)
				)
				diag(self$cache$spin_label) <- NA_character_
				self$cache$dirty <- FALSE
			}
			list(
				E = self$cache$E,
				V = self$cache$V,
				mvals = self$cache$mvals,
				spin_index = self$cache$spin_index,
				spin_label = self$cache$spin_label
			)
		},
		compute_transition_spin_index = function(V) {
			nstates <- nrow(V)
			amp <- array(0, dim = c(nstates, nstates, self$n))
			for (i in seq_len(self$n)) {
				Dei <- t(Conj(V)) %*% self$ops$Ix[[i]] %*% V
				amp[, , i] <- abs(Dei)
			}
			spin_index <- apply(amp, c(1, 2), function(x) {
				if (all(x == 0)) NA_integer_ else which.max(x)
			})
			diag(spin_index) <- NA_integer_
			spin_index
		},
		parse_operator_sum = function(label) {
			if (is.null(label) || !nzchar(label)) {
				stop("Operator label cannot be empty")
			}
			label <- gsub("\\s+", "", label)
			if (label %in% self$shift_labels) {
				idx <- which(self$shift_labels == label)[[1]]
				return(self$ops$Ix[[idx]])
			}
			# Support label with ordinal suffix, e.g., "Ha1" or "Ha2"
			m <- regexec("^(.*?)([0-9]+)$", label, perl = TRUE)
			reg <- regmatches(label, m)[[1]]
			if (length(reg) > 0) {
				base <- reg[[2]]
				ord <- as.integer(reg[[3]])
				if (base %in% self$shift_labels) {
					hits <- which(self$shift_labels == base)
					if (ord < 1 || ord > length(hits)) {
						stop("Label ordinal out of range: ", label)
					}
					return(self$ops$Ix[[hits[[ord]]]])
				}
			}
			# Anti-phase shorthand: A*B*C => 2^(n-1) I(A)x I(B)z I(C)z
			if (label %in% c("allx", "allX", "sumIx")) {
				op <- matrix(0, 2^self$n, 2^self$n)
				for (i in seq_len(self$n)) {
					op <- op + self$ops$Ix[[i]]
				}
				return(op)
			}
			# Split on + and - while preserving signs
			label <- sub("^\\+", "", label)
			label <- gsub("-", "+-", label, fixed = TRUE)
			terms <- strsplit(label, "\\+", fixed = FALSE)[[1]]
			terms <- terms[nzchar(terms)]
			op <- matrix(0, 2^self$n, 2^self$n)
			for (t in terms) {
				if (grepl("~", t, fixed = TRUE)) {
					op <- op + private$parse_antiphase(t)
				} else if (t %in% self$shift_labels || private$is_label_ordinal(t)) {
					idx <- private$resolve_label_index(t)
					op <- op + self$ops$Ix[[idx]]
				} else {
					op <- op + private$parse_operator_product(t)
				}
			}
			op
		},
		parse_antiphase = function(label) {
			label <- gsub("\\s+", "", label)
			coef <- 1
			if (substr(label, 1, 1) == "-") {
				coef <- -1
				label <- substr(label, 2, nchar(label))
			}
			m <- regexec("^([0-9.]+)", label)
			reg <- regmatches(label, m)
			if (length(reg[[1]]) > 0) {
				coef <- coef * as.numeric(reg[[1]][[2]])
				label <- sub("^[0-9.]+", "", label)
			}
			parts <- strsplit(label, "~", fixed = TRUE)[[1]]
			parts <- parts[nzchar(parts)]
			if (length(parts) < 2) {
				stop("Anti-phase shorthand requires at least two labels")
			}
			first_idx <- private$resolve_label_index(parts[[1]])
			op <- self$ops$Ix[[first_idx]]
			for (p in parts[-1]) {
				idx <- private$resolve_label_index(p)
				op <- op %*% self$ops$Iz[[idx]]
			}
			coef * (2^(length(parts) - 1)) * op
		},
		is_label_ordinal = function(label) {
			m <- regexec("^(.*?)([0-9]+)$", label, perl = TRUE)
			reg <- regmatches(label, m)[[1]]
			if (length(reg) == 0) return(FALSE)
			base <- reg[[2]]
			base %in% self$shift_labels
		},
		resolve_label_index = function(label) {
			if (label %in% self$shift_labels) {
				return(which(self$shift_labels == label)[[1]])
			}
			m <- regexec("^(.*?)([0-9]+)$", label, perl = TRUE)
			reg <- regmatches(label, m)[[1]]
			if (length(reg) > 0) {
				base <- reg[[2]]
				ord <- as.integer(reg[[3]])
				if (base %in% self$shift_labels) {
					hits <- which(self$shift_labels == base)
					if (ord < 1 || ord > length(hits)) {
						stop("Label ordinal out of range: ", label)
					}
					return(hits[[ord]])
				}
			}
			stop("Unknown label: ", label)
		},
		parse_operator_product = function(label) {
			label <- gsub("\\s+", "", label)
			coef <- 1
			if (substr(label, 1, 1) == "-") {
				coef <- -1
				label <- substr(label, 2, nchar(label))
			}
			m <- regexec("^([0-9.]+)", label)
			reg <- regmatches(label, m)
			if (length(reg[[1]]) > 0) {
				coef <- coef * as.numeric(reg[[1]][[2]])
				label <- sub("^[0-9.]+", "", label)
			}

			matches <- gregexpr("I([0-9]+)([xyz])", label, perl = TRUE)
			parts <- regmatches(label, matches)[[1]]
			if (length(parts) == 0) {
				stop("Could not parse operator label")
			}

			op <- diag(2^self$n)
			for (p in parts) {
				idx <- as.integer(sub("I([0-9]+)([xyz])", "\\1", p, perl = TRUE))
				axis <- sub("I([0-9]+)([xyz])", "\\2", p, perl = TRUE)
				if (idx < 1 || idx > self$n) {
					stop("Spin index out of range in operator label")
				}
				if (axis == "x") op <- op %*% self$ops$Ix[[idx]]
				if (axis == "y") op <- op %*% self$ops$Iy[[idx]]
				if (axis == "z") op <- op %*% self$ops$Iz[[idx]]
			}
			coef * op
		},
		transition_index = function(n, mode = c("upper", "all")) {
			mode <- match.arg(mode)
			if (mode == "upper") {
				k <- rep(seq_len(n - 1), times = (n - 1):1)
				l <- unlist(lapply(seq_len(n - 1), function(i) (i + 1):n))
				return(data.frame(k = k, l = l))
			}
			k <- rep(seq_len(n), each = n)
			l <- rep(seq_len(n), times = n)
			keep <- k != l
			data.frame(k = k[keep], l = l[keep])
		},
		# placeholder to preserve private list comma structure after removal
		dummy = NULL,
		compute_derivatives = function(E, V, Oe, De, transitions) {
			labels <- names(self$dH)
			nstates <- length(E)
			ntrans <- nrow(transitions)
			out <- matrix(0, ntrans, 2 * length(labels))
			colnames(out) <- c(
				paste0("dFreq_", labels),
				paste0("dInt_", labels)
			)

			for (p in seq_along(labels)) {
				dH <- self$dH[[labels[[p]]]]
				M <- t(Conj(V)) %*% dH %*% V
				dE <- diag(M)

				denom <- outer(E, E, FUN = function(ei, ej) ej - ei)
				K <- matrix(0+0i, nstates, nstates)
				mask <- abs(denom) > 0
				K[mask] <- M[mask] / denom[mask]
				diag(K) <- 0

				dOe <- t(Conj(K)) %*% Oe + Oe %*% K
				dDe <- t(Conj(K)) %*% De + De %*% K

				dfreq <- dE[transitions$l] - dE[transitions$k]
				idx <- cbind(transitions$k, transitions$l)
				dI <- Re(dDe[idx] * Conj(Oe[idx]) + De[idx] * Conj(dOe[idx]))

				out[, p] <- Re(dfreq)
				out[, length(labels) + p] <- Re(dI)
			}

			out
		}
	)
)
