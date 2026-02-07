testthat::test_that("fit-jac-asp-strong", {
  spec_path <- file.path("inst", "extdata", "asp", "400mhz.ft1")
  if (!file.exists(spec_path)) {
    spec_path <- file.path("..", "..", "inst", "extdata", "asp", "400mhz.ft1")
  }
  if (!file.exists(spec_path)) {
    spec_path <- file.path("..", "..", "..", "inst", "extdata", "asp", "400mhz.ft1")
  }
  testthat::expect_true(file.exists(spec_path))

  spec_list <- list(fitnmr::read_nmrpipe(spec_path, dim_order = NULL))
  names(spec_list) <- basename(spec_path)

  start_resonances <- data.frame(
    name = c("HA", "HB2", "HB3"),
    x = c("HA", "HB2", "HB3"),
    x_sc = c("HA-HB2 HA-HB3", "HA-HB2", "HA-HB3"),
    x_ss = c("", "HB", "HB"),
    stringsAsFactors = FALSE,
    row.names = c("HA", "HB2", "HB3"),
    check.names = FALSE
  )
  start_nuclei <- data.frame(
    name = c("HA", "HB2", "HB3"),
    omega0_ppm = c(3.78, 2.698, 2.558),
    r2_hz = c(0.7, 0.7, 0.7),
    stringsAsFactors = FALSE,
    row.names = c("HA", "HB2", "HB3"),
    check.names = FALSE
  )
  start_couplings <- data.frame(
    name = c("HB2-HB3", "HA-HB2", "HA-HB3"),
    hz = c(-17.313, 3.471, 8.864),
    stringsAsFactors = FALSE,
    row.names = c("HB2-HB3", "HA-HB2", "HA-HB3"),
    check.names = FALSE
  )
  spinsystems <- list(
    HB = matrix(
      c("HB2", "HB2-HB3",
        "", "HB3"),
      nrow = 2,
      byrow = TRUE
    )
  )

  start_tables <- list(
    resonances = start_resonances,
    nuclei = start_nuclei,
    couplings = start_couplings,
    spinsystems = spinsystems
  )
  param_list <- fitnmr::tables_to_param_list(spec_list, start_tables)

  fit_input <- do.call(
    fitnmr::make_fit_input,
    c(list(spec_list, omega0_plus = 0.02), fitnmr::param_list_to_arg_list(param_list))
  )
  fit_input$start_list$m0[is.na(fit_input$start_list$m0)] <- max(
    vapply(fit_input$spec_data, function(x) max(abs(x$spec_int)), numeric(1))
  )

  fit_par <- fitnmr:::pack_fit_params(fit_input$start_list, fit_input$group_list)
  jac_analytic <- fitnmr:::fit_jac(fit_par, fit_input)

  step_default <- as.numeric(Sys.getenv("FIT_JAC_STEP", "1e-6"))
  step_by_type <- c(
    omega0 = as.numeric(Sys.getenv("FIT_JAC_STEP_OMEGA0", "1e-6")),
    r2 = as.numeric(Sys.getenv("FIT_JAC_STEP_R2", "1e-4")),
    m0 = as.numeric(Sys.getenv("FIT_JAC_STEP_M0", "1e+01")),
    omega0_comb = as.numeric(Sys.getenv("FIT_JAC_STEP_OMEGA0_COMB", "1e-4"))
  )
  rel_tol <- 1e-4
  fail_names <- character()
  pass_names <- character()
  max_rel_by_param <- numeric(length(fit_par))
  names(max_rel_by_param) <- names(fit_par)

  for (col_idx in seq_along(fit_par)) {
    param_type <- if (startsWith(names(fit_par)[col_idx], "omega0_comb")) {
      "omega0_comb"
    } else {
      sub("_.*$", "", names(fit_par)[col_idx])
    }
    step <- step_by_type[[param_type]] %||% step_default
    h <- step * max(1, abs(fit_par[col_idx]))
    par_plus <- fit_par
    par_minus <- fit_par
    par_plus[col_idx] <- par_plus[col_idx] + h
    par_minus[col_idx] <- par_minus[col_idx] - h

    res_plus <- fitnmr:::fit_fn(par_plus, fit_input)
    res_minus <- fitnmr:::fit_fn(par_minus, fit_input)
    fd_col <- (res_plus - res_minus) / (2 * h)
    jac_col <- jac_analytic[, col_idx]

    diff <- abs(fd_col - jac_col)
    scale <- pmax(1, abs(fd_col), abs(jac_col))
    max_rel <- max(diff / scale)

    max_rel_by_param[col_idx] <- max_rel
    if (max_rel < rel_tol) {
      pass_names <- c(pass_names, names(fit_par)[col_idx])
    } else {
      fail_names <- c(fail_names, names(fit_par)[col_idx])
    }
  }

  if (length(fail_names) > 0) {
    rel_lines <- paste0(
      "  - ",
      names(max_rel_by_param),
      ": ",
      format(max_rel_by_param, digits = 4, scientific = TRUE)
    )
    message(
      paste(
        "Max relative error per parameter:",
        paste(rel_lines, collapse = "\n"),
        sep = "\n"
      )
    )
  }

  testthat::expect_true(
    length(fail_names) == 0,
    info = paste(
      "finite-diff mismatches:",
      paste(fail_names, collapse = ", ")
    )
  )
})
