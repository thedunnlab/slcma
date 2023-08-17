#### Draft slcma R package
#### Version 0.5
#### Andrew Smith
#### 15 December 2022


#' Accumulation hypothesis
#'
#' Function for deriving an accumulation variable
#' from longitudinally measured exposures.
#' 
#' @param ... Binary exposure variables to accumulate.
#' @return Derived accumulation variable, i.e. sum of exposures.
#' 
#' @export
Accumulation <- function(...) {
  X <- cbind(...)
  apply(X,1,sum)
}

#' Recency hypothesis
#'
#' Function for deriving a recency variable
#' from longitudinally measured exposures.
#'
#' @param weights Vector of numeric weights for combining exposures,
#' one per exposure.
#' @param ... Binary exposure variables to combine.
#' @return Derived recency variable, i.e. weighted
#' sum of exposures.
#' 
#' @export
Recency <- function(weights,...) {
  X <- cbind(...)
  if(length(weights)!=ncol(X)) {
    stop("number of weights differs from number of variables") 
  }   
  apply(t(weights*t(X)),1,sum)
}

#' Change hypothesis
#'
#' Function for deriving a change variable
#' from longitudinally measured exposures.
#' 
#' @param x1 First binary exposure measurement variable.
#' @param x2 Second binary exposure measurement variable.
#' @return Derived change variable, i.e. set to
#' 1 for all cases where x1 differs from x2, otherwise 0.
#' 
#' @export
Change <- function(x1,x2) {
  I(x2-x1)
}

#' Mobility hypothesis
#'
#' Function for deriving a mobility variable
#' from longitudinally measured exposures.
#' 
#' @param x1 First binary exposure measurement variable.
#' @param x2 Second binary exposure measurement variable.
#' @return Derived mobility variables, i.e. matrix
#' of two columns, first column indicating
#' upward mobility (`(1-x1)*x2`) and the second
#' column indicating downward mobility (`(1-x2)*x1`).
#' 
#' @export
Mobility <- function(x1,x2) {
  X <- cbind((1-x1)*x2,(1-x2)*x1)
  colnames(X) <- c("up","down")
  X
}

#' Always hypothesis
#'
#' Function for deriving an 'always' variable
#' from longitudinally measured exposures.
#' 
#' @param ... Binary exposure variables to accumulate.
#' @return Derived 'always' variable, i.e. set to
#' 1 if exposure present for all exposure variables, and 0 otherwise.
#' 
#' @export
Always <- function(...) {
  X <- cbind(...)
  apply(X,1,prod)
}

#' Ever hypothesis
#'
#' Function for deriving an 'ever' variable
#' from longitudinally measured exposures.
#' 
#' @param ... Binary exposure variables to accumulate.
#' @return Derived 'ever' variable, i.e.
#' set to 1 if exposure present for any exposure variable,
#' otherwise 0.
#' 
#' @export
Ever <- function(...) {
  X <- cbind(...)
  1-apply(1-X,1,prod)
}

#' Dummy datasets with specified means and covariances
#'
#' Creates a dummy dataset with the same samples means
#' and sample covariances as that estimated from pooled imputed datasets.
#'
#' @param formula An object of class `formula` which specifies the
#' SLCMA model including outcome variable, exposures
#' and lifecourse hypotheses.
#' @param data An optional data frame, list or environment
#' containing the variables in the model. If not found in `data`,
#' the variables are taken from `environment(formula)`, typically
#' the environment form which the function is called.
#' @param seed Integer input to `set.seed` to ensure that
#' dataset generation is reproducible (Default: 1).
#' @return Data frame of class `MI2dummy` including randomly
#' generated variables specified by the specified `formula` with
#' means and covariances obtained from the input `data`.
#'
#' @export
MI2dummy <- function(formula, data, seed=1) {

  if(!("mids" %in% class(data))) {
    stop("'data' is not of class 'mids'")
  }

  complete_data <- complete(data, action="stacked")
  n <- dim(complete_data)[1] / data$m

  trms <- terms(formula, data=complete_data, keep.order=TRUE)
  mf <- model.frame(trms, data=complete_data)
  mm <- model.matrix(trms, data=complete_data)
  assign <- attributes(mm)$assign
  if(attributes(trms)$intercept==1) {
    mm <- mm[,-1, drop=FALSE]
    assign <- assign[-1]
  }

  y <- model.response(mf)
  X <- cbind(mm,y)
  varnames <- colnames(mm)
  if(!is.null(y)) {
    assign <- c(assign,-1)
    varnames <- c(varnames,as.character(formula[[2]]))
  }

  Xt1 <- t(X) %*% rep(1,n*data$m) / data$m
  XtX <- t(X) %*% X / data$m
  qr_decomp <- qr(XtX)

  set.seed(seed)
  Z <- matrix(rnorm(n*qr_decomp$rank), nrow=n)
  Z <- scale(Z)
  chol_decomp <- chol(t(Z) %*% Z)
  Z <-  Z %*% solve(chol_decomp)

  if(qr_decomp$rank < dim(X)[2]) {
    bi <- qr_decomp$pivot[1:qr_decomp$rank]
    reconstruct <- solve(XtX[bi,bi], XtX[bi,-bi])
  } else {
    bi <- 1:dim(X)[2]
  }

  Mvec <- Xt1[bi] / n
  Vmat <- XtX[bi,bi]

  chol_decomp <- chol(Vmat - n * Mvec %*% t(Mvec))
  dummy <- Z %*% chol_decomp + rep(1,n) %*% t(Mvec)
  
  if(qr_decomp$rank < dim(X)[2]) {
    reconstructed <- dummy %*% reconstruct
    dummy <- cbind(dummy, reconstructed)
    colnames(dummy) <- c(varnames[bi],varnames[-bi])
    assign <- c(assign[bi],assign[-bi])
  }

  dummy <- data.frame(dummy, check.names=FALSE)
  class(dummy) <- c(class(dummy),"MI2dummy")
  attributes(dummy)$assign <- assign
  dummy

}

#' Structured Approach to Evaluating Life-Course Hypotheses: stage 1
#'
#' Performs stage 1 of the SLCMA
#'
#' @param formula An object of class `formula` which specifies the
#' SLCMA model including outcome variable, exposures
#' and lifecourse hypotheses.
#' @param data An optional data frame, list or environment
#' containing the variables in the model. If not found in `data`,
#' the variables are taken from `environment(formula)`, typically
#' the environment form which the function is called.
#' @param adjust MISSING
#' @param ... Additional arguments to be passed to `MI2dummy`.
#' @return An object of class `"sclma"` which inherits from
#' the `"lars"` class containing output from the least angle regression
#' analysis as well as the following items:
#' * `X` Numeric matrix with columns corresponding to SLCMA hypotheses.
#' * `y` Output variable.
#' * `ncov`
#' * `multiple`
#' * `sanity`
#'  
#' @export 
slcma <- function(formula, data=environment(formula), adjust=NULL, ...) {

  # Checks
  formula <- update(formula, . ~ . + 1)
  if("mids" %in% class(data)) {
    check_data <- complete(data, action="stacked")
  } else {
    check_data <- data
  }
  trms <- terms(formula, data=check_data, keep.order=TRUE)
  if(attributes(trms)$response==0) {
    stop("No response variable specified")
  }
  mf <- model.frame(trms, data=check_data)
  if(attributes(attributes(mf)$terms)$dataClasses[1]!="numeric") {
    stop("Response is not of class: numeric")
  }

  if("mids" %in% class(data)) {
    dummy <- MI2dummy(formula, data, ...)
    y <- dummy[,attributes(dummy)$assign==-1]
    mm <- as.matrix(cbind(`(Intercept)`=rep(1,dim(dummy)[1]), dummy[,attributes(dummy)$assign!=-1]))
    assign <- c(0, attributes(dummy)$assign[attributes(dummy)$assign!=-1])
    r <- assign %in% c(0, adjust)
  } else {
    y <- model.response(mf)
    mm <- model.matrix(trms, data=data)
    r <- attributes(mm)$assign %in% c(0,adjust)
  }

  covars  <- mm[, r, drop=FALSE]
  X_hypos <- mm[,!r, drop=FALSE]

  sanity <- data.frame(Term=colnames(mm),
                       Role=factor(r, levels=c(T,F), labels=c("Adjusted for in all models","Available for variable selection")))

  ncov <- dim(covars)[2] - 1
  y <- lm(y ~ covars)$residuals
  X_hypos <- lm(X_hypos ~ covars)$residuals

  full_model <- lm(y ~ X_hypos)
  df1 <- dim(X_hypos)[1] - full_model$df.residual - 1
  df2 <- full_model$df.residual - dim(covars)[2] + 1
  RSS <- sum(full_model$residuals^2)
  sigma <- sqrt(RSS/df2)
  r.squared <- 1 - RSS/sum(y^2)
  fstatistic <- r.squared / (1-r.squared) * df2/df1
  p <- pf(fstatistic, df1, df2, lower.tail=FALSE)
  multiple <- list(sigma = sigma, r.squared = r.squared, fstatistic = fstatistic, p = p)

  output <- lars(X_hypos, y)

  output$X <- X_hypos
  output$y <- y
  output$ncov <- ncov
  output$multiple <- multiple
  output$sanity <- sanity
  class(output) <- "slcma"
  print(sanity, row.names=FALSE)
  output
}

#' SLCMA stage 1 output
#'
#' Print method for objects of class "slcma".
#'
#' @param x Object of class "`slcma`" from the `slcma()` function.
#'
#' @export
print.slcma <- function(x) {
  print(x$sanity, row.names=FALSE)
}

#' SLCMA output summary
#' 
#' Summary method for class "slcma".
#'
#' @param x Object of class "`slcma`" from the `slcma()` function.
#'
#' @export
summary.slcma <- function(x) {
  unl <- unlist(x$actions)
  varnames <- names(unl)
  entries <- varnames
  exits <- varnames
  entries[unl < 0] <- ""
  exits[unl > 0] <- ""
  output <- data.frame(Step=seq(along=unl),
                       `Variable selected`=entries, 
                       `Variable removed`=exits,
                       Variables=x$df[-1]-1, 
                       R2=round(x$R2[-1],3), check.names=FALSE)
  cat("\nSummary of LARS procedure\n")
  print(output, row.names=FALSE)
  invisible(output)
}

#' SLCMA elbow plot
#'
#' Elbow plot (plot method for class "slcma")
#'
#' @param x Object of class "`slcma`" from the `slcma()` function.
#'
#' @export
plot.slcma <- function(slcma, relax=FALSE, show.remove=FALSE, show.labels=TRUE, labels=NULL, selection.order=FALSE,
                       xlab="Variables selected", ylab="R-squared", ...) {
  maxadd <- rle(slcma$actions>0)$lengths[1]
  vars <- slcma$df - 1
  R2 <- slcma$R2
  unl <- unlist(slcma$actions)

  if(relax) {
    for(j in seq_along(unl)) {
      selected <-   seq_len(dim(slcma$X)[2]) %in%  unl[1:j]
      selected <- !(seq_len(dim(slcma$X)[2]) %in% -unl[1:j]) & selected
      relaxed_model <- lm(slcma$y ~ slcma$X[,selected])
      R2[j+1] <- summary(relaxed_model)$r.sq
    }
  }
  maxR2 <- slcma$multiple$r.sq

  labs <- attr(unl,"names")
  if(!is.null(labels)) {
    if(selection.order) {
      if(length(labels) > maxadd) {
        stop("'labels' is longer than the number of selected variables")
      }
      labs[seq_along(labels)] <- labels
      print(data.frame(Label=labs, Term=attr(unl,"names")), row.names=FALSE)
    } else {
      if(length(labels) > dim(slcma$X)[2]) {
        stop("'labels' is longer than the maximum number of variables")
      }
      use.names <- colnames(slcma$X)
      use.names[seq_along(labels)] <- labels
      print(data.frame(Label=use.names, Term=colnames(slcma$X)), row.names=FALSE)
      labs <- use.names[unl[1:maxadd]]
    }
  }
  labs[-(1:maxadd)] <- ""

  if(!show.remove) {
    vars <- vars[1:(maxadd+1)]
    R2 <- R2[1:(maxadd+1)]
    labs <- labs[1:maxadd]
  }

  plot(c(0,max(vars)), c(0,maxR2), type="n", xaxt="n",
       xlab=xlab, ylab=ylab, ...)
  axis(1, at=vars)
  if(show.labels) {  
    text(vars[-1],0, labs, srt=90, adj=0)
  }
  lines(vars, R2, ...)
  abline(h=maxR2, lty="dotted")
}


#' Structured Approach to Evaluating Life-Course Hypotheses: stage 1
#'
#' Performs stage 2 of the SLCMA for user-specified methods
#'
#' @param slcma
#' @param step
#' @param method
#' @param alpha
#' @param do.maxtCI
#' @param ... Additional arguments to `slcmaFLI()`.
#' @return
#' 
#' @export
slcmaInfer <- function(slcma, step=1L, method="slcmaFLI", alpha=0.05, do.maxtCI=FALSE, ...) {
  if(step<1 | abs(step-round(step)) > .Machine$double.eps^0.5) {
    stop("'step' is not a positive integer")
  }
  if(step > length(slcma$action)) {
    stop("'step' is greater than the number of LARS steps in 'slcma'") 
  }

  nvars <- slcma$df[step+1] - 1
  R2.lasso <- slcma$R2[step+1]
  output <- list(step=step, nvars=nvars, R2.lasso=R2.lasso)

  if(any(c("fixedLassoInf","selectiveInference","fli","slcmaFLI") %in% method)) {
    fli <- slcmaFLI(slcma, step=step, alpha=alpha, ...)
    output$fli <- fli
  }
  if(any(c("maxt","slcmaMaxt") %in% method)) {
    if(step > 1) {
      stop("max-|t| test is only available at Step 1")
    }
    maxt <- slcmaMaxt(slcma, alpha=alpha, do.CI=do.maxtCI)
    output$maxt <- maxt
  }
  if(any(c("relaxed","relax","slcmaRelax") %in% method)) {
    relax <- slcmaRelax(slcma, step=step)
    output$relax <- relax
  }
  if(any(c("Bayes","slcmaBayes") %in% method)) {
    if(step > 1) {
      stop("Bayesian posterior probabilities are only available at Step 1")
    }
    Bayes <- slcmaBayes(slcma)
    output$Bayes <- Bayes
  } 

  class(output) <- "slcmaInfer"
  output
}

#' SLCMA stage 2 output
#'
#' Print method for objects of class `"slcmaInfer"`.
#'
#' @param x Object of class "`slcmaInfer`" from the `slcmaInfer()` function.
#'
#' @export
print.slcmaInfer <- function(x) {
  cat(sprintf("\nInference for model at Step %1.0f of LARS procedure\n",x$step))
  cat(sprintf("\nNumber of selected variables: %1.0f\n", x$nvars))
  cat(sprintf("Lasso R-squared: %.3f\n", x$R2.lasso))
  if(!is.null(x$fli)) {
    cat("\nResults from fixed lasso inference (selective inference):\n")
    print(x$fli)
  }
  if(!is.null(x$maxt)) {
    cat("\nResults from max-|t| test:\n")
    print(x$maxt)
  }
  if(!is.null(x$relax)) {
    cat("\nResults from relaxed lasso:\n")
    print(x$relax)
  }
  if(!is.null(x$Bayes)) {
    print(x$Bayes)
  }
}

#' slcmaFLI
#'
#' Short description of slcmaFLI.
#'
#' @param slcma
#' @param step
#' @param alpha
#' @param ...
#' @return
#'
#' @export
slcmaFLI <- function(slcma, step, alpha=0.05, ...) {
  if(step<1 | abs(step-round(step)) > .Machine$double.eps^0.5) {
    stop("'step' is not a positive integer")
  }
  if(step > length(slcma$action)) {
    stop("'step' is greater than the number of LARS steps in 'slcma'") 
  }
  if(step == length(slcma$action)) {
    stop("Fixed lasso inference and/or selective inference is not available at the final step of the LARS procedure")
  }
  sumsq <- slcma$normx
  X_normed <- scale(slcma$X, scale=sumsq)
  fli <- fixedLassoInf(X_normed, slcma$y, 
                       slcma$beta[step+1,], slcma$lambda[step+1], 
                       type="partial", alpha=alpha, ...)
  sumsq <- sumsq[fli$vars]
  fli$coef0 <- fli$coef0 / sumsq
  fli$sd <- fli$sd / sumsq
  fli$ci <- fli$ci / cbind(sumsq,sumsq)
  class(fli) <- "slcmaFLI"
  fli
}

#' Print slcmaFLI output
#'
#' Print method for class "slcmaFLI" (adapted from 'selectiveinference' package).
#'
#' @param x
#' @param tailarea (Default: TRUE)
#' @param ...
#' @return
#'
#' @export
print.slcmaFLI <- function (x, tailarea = TRUE, ...) 
{
  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n", 
              x$sigma))
  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n", 
              x$lambda, x$alpha))
  cat("", fill = T)
  tab = cbind(round(x$coef0, 3), round(x$pv, 3), round(x$ci, 3))
  colnames(tab) = c("Coef", "P-value", "CI.lo", 
                    "CI.up")
  rownames(tab) = attr(x$vars,"names")
  if (tailarea) {
    tab = cbind(tab, round(x$tailarea, 3))
    colnames(tab)[(ncol(tab) - 1):ncol(tab)] = c("LoTailArea", 
                                                 "UpTailArea")
  }
  print(tab)
  invisible()
}

#' Max-|t| test
#'
#' Function for max-|t| test and associated confidence intervals
#'
#' @param slcma
#' @param alpha (Default: 0.05)
#' @param do.CI (Default: TRUE)
#' @param seed (Default: 12345)
#' @param ...
#' @return
#'
#' @export
slcmaMaxt <- function(slcma, alpha=0.05, do.CI=TRUE, seed=12345, ...) {
  # assumes that X and y have mean subtracted from them, which should have happened by default
  y <- slcma$y
  X_hypos <- slcma$X
  n <- length(y)
  p <- dim(X_hypos)[2]
  d <- n-p-slcma$ncov-1
  sumsq <- slcma$normx
  X_normed <- scale(X_hypos, scale=sumsq)
  Xt <- t(X_normed)
  XtX <- Xt %*% X_normed
  Xty <- Xt %*% y
  selection <- slcma$action[[1]]
  r <- Xty[selection]
  s <- slcma$multiple$sigma
  set.seed(seed)
  p.maxt <- 1 - pmvt(lower=-rep(abs(r),p),
                     upper= rep(abs(r),p),
                     delta= rep(0,p), 
                     df= d, sigma= s^2 * XtX, ...)
  coefficients <- cbind(r / sumsq[selection], p.maxt)
  rownames(coefficients) <- names(selection)
  colnames(coefficients) <- c("Coef", "P-value")
  output <- list(coefficients = coefficients, sigma = s, df = d,
                 error = attributes(p.maxt)$error, msg = attributes(p.maxt)$msg, alpha = alpha)
  if(do.CI) {
    search_middle <- r
    search_radius <- 1.5*qnorm(1-alpha/2)*s*sqrt(XtX[selection,selection])
    lower <-  uniroot(function(beta0) {
                        P6(beta0, beta0+abs(r-beta0), p, selection, s, XtX, p, d, ...) -
                        P6(beta0, beta0-abs(r-beta0), p, selection, s, XtX, p, d, ...) -
                        (1-alpha)
                      }, lower=search_middle-search_radius, upper=search_middle)$root
    upper <-  uniroot(function(beta0) {
                        P6(beta0, beta0+abs(r-beta0), p, selection, s, XtX, p, d, ...) -
                        P6(beta0, beta0-abs(r-beta0), p, selection, s, XtX, p, d, ...) -
                        (1-alpha)
                      }, lower=search_middle, upper=search_middle+search_radius)$root
    output$coefficients <- cbind(output$coefficients, cbind(lower,upper) / sumsq[selection])
    colnames(output$coefficients)[3:4] <- c("CI.lo","CI.up")
  }
  class(output) <- "slcmaMaxt"
  output
}

#' Print slcmaMaxt
#' 
#' Print method for objects of class `"slcmaMaxt"`.
#'
#' @param x
#'
#' @export
print.slcmaMaxt <- function (x) {
  cat(sprintf("\nError standard deviation: %0.3f on %0.0f degrees of freedom", x$sigma, x$df))
  if(dim(x$coefficients)[2]>2) {
    cat(sprintf("\nConfidence interval coverage %0.1f percent\n", 100*(1-x$alpha)))
  } else{
    cat("\n")
  }
  cat("", fill=TRUE)
  print(x$coefficients)
  cat(sprintf("\nMessage regarding calculation of P-value: %s", x$msg))
  cat(sprintf("\n with estimated numeric absolute error: %e\n", x$error))
  invisible()
}

# Functions needed by maxt
# Calculates the probability in (7)
P7 <- function(r, p, selection, mu, Sigma, df, ...) {
  Cmatrix <- rbind(diag(p),-diag(p))
  Cmatrix[,selection] <- 1
  Cmatrix <- Cmatrix[-(p+selection),] 
  lower <- rep(0, 2*p-1)
  lower[selection] <- r
  pmvt(lower=lower, 
       upper=rep(Inf, 2*p-1), 
       delta=as.vector(Cmatrix %*% mu),
       df=df, sigma=Cmatrix %*% Sigma %*% t(Cmatrix), type="shifted")
} 
# Calculates the probability in (6)
P6 <- function(beta0, r, p, selection, s, XtX, df, ...) {
  upper.denom <- P7(0, p, selection, XtX[,selection]*beta0, s^2 * XtX, df, ...)
  lower.denom <- P7(0, p, selection,-XtX[,selection]*beta0, s^2 * XtX, df, ...)
  if(r >= 0) {
    numer <- P7(r, p, selection, XtX[,selection]*beta0, s^2 * XtX, df, ...)
    prob <- 1 - numer / (lower.denom + upper.denom)
  } else {
    numer <- P7(-r,p, selection,-XtX[,selection]*beta0, s^2 * XtX, df, ...)
    prob <- numer / (lower.denom + upper.denom)
  }
  prob
}


# A function for simple confidence interval calculations
slcmaRelax <- function(slcma, step, alpha=0.05) {
  if(step<1 | abs(step-round(step)) > .Machine$double.eps^0.5) {
    stop("'step' is not a positive integer")
  }
  if(step > length(slcma$action)) {
    stop("'step' is greater than the number of LARS steps in 'slcma'") 
  }
  unl <- unlist(slcma$actions)
  selected <-   seq_len(dim(slcma$X)[2]) %in%  unl[1:step]
  selected <- !(seq_len(dim(slcma$X)[2]) %in% -unl[1:step]) & selected
  relaxed_model <- lm(slcma$y ~ slcma$X[,selected])
  relaxed_summary <- summary(relaxed_model)
  coefs <- relaxed_summary$coefficients
  SEs <- coefs[-1,2] * slcma$multiple$sigma / relaxed_summary$sigma
  lower <- coefs[-1,1] - qnorm(1-alpha/2)*SEs
  upper <- coefs[-1,1] + qnorm(1-alpha/2)*SEs
  coefficients <- cbind(coefs[-1,1], SEs, lower, upper)
  rownames(coefficients) <- colnames(slcma$X)[selected]
  colnames(coefficients) <- c("Coef", "SE", "CI.lo", "CI.up")
  relax <- list(coefficients = coefficients, R2.relax = relaxed_summary$r.sq)
  class(relax) <- "slcmaRelax"
  relax
}

#' Print slcmaRelax
#'
#' Print method for objects of class `"slcmaRelax"`
#'
#' @param x
#'
#' @export
print.slcmaRelax <- function (x) {
  cat(sprintf("\nRelaxed r-squared: %0.3f\n", x$R2.relax))
  cat("", fill=TRUE)
  print(x$coefficients)
  invisible()
}

slcmaBayes <- function(slcma) {
  y <- slcma$y
  X_hypos <- slcma$X
  p <- dim(X_hypos)[2]
  sigma <- slcma$multiple$sigma
  numer <- numeric(p)
  for(j in 1:p) {
    residuals <- lm(y ~ X_hypos[,j])$residuals
    RSSj <- sum(residuals^2)
    numer[j] <- exp(-RSSj/(2*sigma^2))
  }
  probs <- numer / sum(numer)
  names(probs) <- colnames(X_hypos)
  Bayes <- list(probs=probs)
  class(Bayes) <- "slcmaBayes"
  Bayes
}

# Print method for class "slcmaBayes"
print.slcmaBayes <- function (x) {
  cat("\nPosterior probabilities:")
  cat("", fill=TRUE)
  probs <- data.frame(x$probs)
  names(probs) <- NULL
  print(round(probs,min(max(dim(probs)[1],3),8)))
  invisible()
}

