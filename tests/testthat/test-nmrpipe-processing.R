test_that("nmrpipe_sp/zf/ft/fti/ps match NMRPipe outputs", {

  find_pkg_root <- function(path = getwd()) {
    cur <- normalizePath(path, mustWork = FALSE)
    while (!is.na(cur) && cur != dirname(cur)) {
      if (file.exists(file.path(cur, "DESCRIPTION"))) {
        return(cur)
      }
      cur <- dirname(cur)
    }
    NA_character_
  }

  env_dir <- Sys.getenv("NMRPIPE_TEST_DIR", unset = "")
  use_env_dir <- nzchar(env_dir) && dir.exists(env_dir)

  base <- ""
  if (use_env_dir) {
    base <- env_dir
  } else if (nzchar(Sys.which("nmrPipe"))) {
    root <- find_pkg_root()
    skip_if(is.na(root), "Package root not found; cannot locate test script")

    fid_override <- Sys.getenv("NMRPIPE_TEST_FID", unset = "")
    if (nzchar(fid_override)) {
      skip_if_not(file.exists(fid_override), "NMRPIPE_TEST_FID does not exist")
      fid_src <- fid_override
    } else {
      base_src <- if (!is.na(root)) {
        file.path(root, "inst", "extdata", "noesy1d/11")
      } else {
        ""
      }
      if (!dir.exists(base_src)) {
        base_src <- system.file("extdata", "noesy1d/11", package = "fitnmr")
      }
      skip_if(base_src == "" || !dir.exists(base_src), "extdata/noesy1d/11 not available")
      fid_src <- file.path(base_src, "test.fid")
    }

    workdir <- file.path(tempdir(), paste0("fitnmr-nmrpipe-", Sys.getpid()))
    dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

    file.copy(fid_src, workdir, overwrite = TRUE)
    file.copy(file.path(root, "tests", "testthat", "test-nmrpipe-processing.csh"),
      workdir,
      overwrite = TRUE
    )

    old_wd <- getwd()
    setwd(workdir)
    on.exit(setwd(old_wd), add = TRUE)
    status <- system2("csh", "test-nmrpipe-processing.csh")
    skip_if_not(status == 0, "Failed to generate NMRPipe outputs via nmrPipe")
    base <- workdir
  } else {
    skip("nmrPipe not available and NMRPIPE_TEST_DIR not set")
  }

  fid_path <- file.path(base, "test.fid")
  sp_path <- file.path(base, "SP.fid")
  zf_path <- file.path(base, "SP_ZF.fid")
  ft_path <- file.path(base, "SP_ZF_FT.ft1")
  ps_path <- file.path(base, "SP_ZF_FT_PS.ft1")
  fti_path <- file.path(base, "SP_ZF_FT_PS_FTI.ft1")

  skip_if_not(file.exists(fid_path), "Missing test.fid")
  skip_if_not(file.exists(sp_path), "Missing SP.fid")
  skip_if_not(file.exists(zf_path), "Missing SP_ZF.fid")
  skip_if_not(file.exists(ft_path), "Missing SP_ZF_FT.ft1")
  skip_if_not(file.exists(ps_path), "Missing SP_ZF_FT_PS.ft1")
  skip_if_not(file.exists(fti_path), "Missing SP_ZF_FT_PS_FTI.ft1")

  fid_ref <- fitnmr::read_nmrpipe(fid_path, complex_data = TRUE)
  sp_ref <- fitnmr::read_nmrpipe(sp_path, complex_data = TRUE)
  zf_ref <- fitnmr::read_nmrpipe(zf_path, complex_data = TRUE)
  ft_ref <- fitnmr::read_nmrpipe(ft_path, complex_data = TRUE)
  ps_ref <- fitnmr::read_nmrpipe(ps_path, complex_data = TRUE)
  fti_ref <- fitnmr::read_nmrpipe(fti_path, complex_data = TRUE)

  sp <- fitnmr::nmrpipe_sp(fid_ref)
  sp_calc <- as.vector(sp$int)
  sp_exp <- as.vector(sp_ref$int)
  expect_lt(
    max(abs(sp_calc - sp_exp)) / max(abs(sp_exp)),
    1e-6
  )

  zf <- fitnmr::nmrpipe_zf(sp_ref)
  zf_calc <- as.vector(zf$int)
  zf_exp <- as.vector(zf_ref$int)
  expect_lt(
    max(abs(zf_calc - zf_exp)) / max(abs(zf_exp)),
    1e-6
  )

  ft <- fitnmr::nmrpipe_ft(zf_ref)
  ft_calc <- as.vector(ft$int)
  ft_exp <- as.vector(ft_ref$int)
  expect_lt(
    max(abs(ft_calc - ft_exp)) / max(abs(ft_exp)),
    1e-4
  )

  ps <- fitnmr::nmrpipe_ps(ft_ref, p0 = 10, p1 = 10)
  ps_calc <- as.vector(ps$int)
  ps_exp <- as.vector(ps_ref$int)
  expect_lt(
    max(abs(ps_calc - ps_exp)) / max(abs(ps_exp)),
    1e-6
  )

  fti <- fitnmr::nmrpipe_fti(ps_ref)
  fti_calc <- as.vector(fti$int)
  fti_exp <- as.vector(fti_ref$int)
  expect_lt(
    max(abs(fti_calc - fti_exp)) / max(abs(fti_exp)),
    1e-4
  )
})
