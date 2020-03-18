# FitNMR

## Installation

The best way to install FitNMR is downloading and installing one of the [released packages](https://smith-group.github.io/fitnmr_releases/). These include prebuilt vignettes, which you can also view online as described below.

The latest developmental version of FitNMR can be installed directly from GitHub by running the following command from within R. Before doing so, you should make sure that you have a writable library directory set using the `.libPaths()` command.

```
source("https://raw.githubusercontent.com/r-lib/remotes/master/install-github.R")$value("Smith-Group/fitnmr")
```

Installing `fitnmr` in this way will not include the vignettes.

## Usage

An example of how to use FitNMR for 2D peak fitting can be found in the [Automated 2D Peak Fitting Scripts](https://smith-group.github.io/fitnmr/articles/peak2d_scripts.html) vignette. For more background on the algorithms, or if you are comfortable with R code, see the [Automated 2D Peak Fitting Code](https://smith-group.github.io/fitnmr/articles/peak2d.html) vignette.

For a roadmap describing the different categories of functions available, see the [FitNMR Package Overview](https://smith-group.github.io/fitnmr/reference/fitnmr.html).

To get documentation for individual functions from within R, use the following commands:

```
library(fitnmr)
?read_nmrpipe
```

To get an HTML listing of all the function-level documentation, you can run:

```
help.start()
```

After doing so, navigate to the "Packages" page, then click the "fitnmr" package.

If you are new to R, chapters 1, 2, 5, and 6 of [An Introduction to R](https://cran.r-project.org/doc/manuals/r-release/R-intro.html) are highly recommended.

## License

FitNMR is released under the GNU Public License version 3.0.
