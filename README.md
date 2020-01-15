# FitNMR

## Installation

The latest developmental version of FitNMR can be installed by running the following command from within R. Before doing so, you should make sure that you have a writable library directory set using the `.libPaths()` command.

```
source("https://raw.githubusercontent.com/r-lib/remotes/master/install-github.R")$value("Smith-Group/fitnmr")
```

Installing `fitnmr` in this way will not include the vignettes. To read those, you must install one of the forthcoming release packages, or view them online as described below.

## Usage

An example of how to use FitNMR for 2D peak fitting can be found in the [Automated 2D Peak Fitting](https://smith-group.github.io/fitnmr/peak2d.html) vignette. If you are not comfortable with using R at the command line, we plan to very soon add a set of prewritten scripts that handle 2D peak identification and refitting across multiple spectra. All you will need to do is adjust a few settings in each script to do automatic fitting of 2D spectra.

To get documentation for individual functions, you can use the following commands:

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