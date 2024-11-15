% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slcma.R
\name{slcmaInfer}
\alias{slcmaInfer}
\title{Structured Life Course Modelling Approach: stage 1}
\usage{
slcmaInfer(
  x,
  step = 1L,
  method = "slcmaFLI",
  alpha = 0.05,
  do.maxtCI = FALSE,
  search.width = 3,
  ...
)
}
\arguments{
\item{x}{Object of class \code{slcma} from the \code{slcma()} function.}

\item{step}{Integer specifying which step of the LARS procedure produces 
the model on which inference is to be performed.}

\item{method}{Character string or vector containing the method or methods of inference
to be performed. (Default: \code{slcmaFLI} - Fixed Lasso Inference).}

\item{alpha}{Level of significance for confidence interval calculations. Confidence
intervals will have \code{(1 - alpha) * 100}\% coverage. (Default: 0.05 - 95\% coverage).}

\item{do.maxtCI}{Logical, indicating whether confidence intervals should be 
calculated for the max-|t| test. (Default: FALSE).}

\item{search.width}{Width of initial interval used in numerical estimation of the confidence intervals
related to the max-|t| test, in multiples of the width of the naive confidence interval. (Default: 3).}

\item{...}{Additional arguments to \code{fixedLassoInf()} (see the \code{selectiveinference} package).}
}
\value{
An list of class \code{slcmaInfer} with one element providing the output for each inference method.
}
\description{
Performs stage 2 of the SLCMA for user-specified methods
}
