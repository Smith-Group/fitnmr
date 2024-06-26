% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resonance.R
\name{plot_resonances_3d}
\alias{plot_resonances_3d}
\title{Plot resonances from 3D fit}
\usage{
plot_resonances_3d(
  fit_data,
  omega0_plus,
  resonances = unique(fit_data$resonance_names)
)
}
\arguments{
\item{fit_data}{fit_input or fit_output structure}

\item{omega0_plus}{length 3 vector giving ppm range for each dimension}
}
\description{
Plot resonances from 3D fit
}
