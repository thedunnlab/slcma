B## Useful information for creating packages
## https://kbroman.org/pkg_primer/

## run this script from outside the package base directory

library(roxygen2)
library(devtools)

devtools::document("slcma")

devtools::build_vignettes("slcma", clean=FALSE, install=TRUE, keep_md=TRUE)

## for some reason build_vignettes() does not keep the markdown output file
library(knitr)
setwd("slcma/doc")
knit("slcma.Rmd")
