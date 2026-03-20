# FitNMR

[![CRAN version](https://www.r-pkg.org/badges/version/fitnmr)](https://CRAN.R-project.org/package=fitnmr)
[![R-universe version](https://smith-group.r-universe.dev/badges/fitnmr)](https://smith-group.r-universe.dev/fitnmr)
[![R-universe checks](https://smith-group.r-universe.dev/fitnmr/badges/checks)](https://smith-group.r-universe.dev/fitnmr#checktable)

FitNMR provides tools for fitting and analyzing 1D-4D nuclear magnetic resonance spectra in R.

## Installation

The most recent release of `fitnmr` can be installed from [CRAN](https://cran.r-project.org) using this R command:

```
install.packages("fitnmr", repos = "https://cloud.r-project.org")
```

The latest developmental version of `fitnmr` can be installed from the [Smith Lab R-universe repository](https://smith-group.r-universe.dev/) using this R command:

```
install.packages("fitnmr", repos = c("https://smith-group.r-universe.dev", "https://cloud.r-project.org"))
```

Please look at the [R-universe check status](https://smith-group.r-universe.dev/fitnmr#checktable) before using the developmental version.

## Usage

An example of how to use FitNMR for 2D peak fitting of 15N HSQC spectra can be found in the [Automated 2D Peak Fitting Scripts](https://smith-group.github.io/fitnmr/articles/peak2d_scripts.html) vignette. This explains how to use FitNMR in a scripted workflow without writing any R code.

For more background on the 2D HSQC fitting algorithms, or if you are comfortable with R code, see the [Automated 2D Peak Fitting Code](https://smith-group.github.io/fitnmr/articles/peak2d.html) vignette. 

FitNMR has an interface for fitting spectra using tables of nuclei, couplings, and resonances visible in the spectra. An application of this to fitting 1D spectra of free amino acids was described in [Syed et al. 2024](https://doi.org/10.5194/mr-5-103-2024).

Code to entirely reproduce the manuscript is [available in a GitHub repository](https://github.com/smith-group/syed2024). The text and accompanying figures/tables can be compiled from an R Markdown file run from within the browser using the Binder service. See the repository for more details.

There is also functionality for doing 1D NMR data processing using functions that emulate the behavior of NMRPipe. Those were used in [Blejec et al. 2026](https://doi.org/10.1002/pro.70508) to preprocess unlabeled 1D protein NMR spectra for quantitative analysis.

See the [1D Time Series Preprocessing/Two-State Fitting](https://smith-group.github.io/fitnmr/articles/timeseries1d.html) vignette for more information.

For a roadmap describing the different categories of functions available, see the [FitNMR Package Overview](https://smith-group.github.io/fitnmr/reference/fitnmr.html).

If you are new to R, chapters 1, 2, 5, and 6 of [An Introduction to R](https://cran.r-project.org/doc/manuals/r-release/R-intro.html) are highly recommended.


## License

FitNMR is released under the GNU Public License version 3.0.
