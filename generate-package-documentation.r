B## Useful information for creating packages
## https://kbroman.org/pkg_primer/

## run this script from outside the package base directory

library(roxygen2)
library(devtools)

devtools::document("slcma")

devtools::build_vignettes("slcma", clean=FALSE, install=TRUE, keep_md=TRUE)

library(knitr)
setwd("slcma/doc")
knit("slcm.Rmd")
