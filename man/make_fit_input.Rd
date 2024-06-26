% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitnmr.R
\name{make_fit_input}
\alias{make_fit_input}
\title{Prepare input data structure for peak fitting}
\usage{
make_fit_input(
  spectra,
  omega0_start,
  omega0_plus,
  omega0_minus = omega0_plus,
  omega0_trunc = NULL,
  r2_start = NULL,
  m0_start = NULL,
  m0_region = (omega0_plus + omega0_minus)/2,
  p0_start = 0,
  p1_start = 0,
  omega0_group = NULL,
  r2_group = NULL,
  m0_group = NULL,
  p0_group = 0,
  p1_group = 0,
  omega0_comb = NULL,
  omega0_comb_start = NULL,
  omega0_comb_group = NULL,
  coupling_comb = NULL,
  resonance_names = NULL,
  nucleus_names = NULL,
  field_offsets = numeric(),
  field_start = numeric(),
  field_group = 0,
  fheader = NULL
)
}
\description{
Prepare input data structure for peak fitting
}
