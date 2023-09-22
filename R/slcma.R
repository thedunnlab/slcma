#### Draft slcma R package
#### Version 0.7
#### Andrew Smith
#### 21 September 2023


#' Specifying hypotheses
#'
#' Functions for deriving variables corresponding
#' various hypotheses about longitudinally
#' measured exposures.
#'
#' @param weights Vector of numeric weights for combining exposures,
#' one per exposure.
#' @param ... Binary exposure variables.
#' @return Derived hypothesis variable(s) as follows:
#' \itemize{
#' \item \code{Accumulation(...)} - sum of exposures
#' \item \code{Recency(weight,...)} - weighted sum of exposures
#' \item \code{Change(x1,x2)} - variable indicating where the two exposures differ
#' \item \code{Mobility(x1,x2)} - two variables with the first indicating upward mobility (i.e. \code{x1==0 && x2==1}) and the second indicating downward mobility (i.e. \code{x1==1 && x2==0})
#' \item \code{Always(...)} - variable indicating whether all exposures present
#' \item \code{Ever(...)} - variable indicating whether at least one exposure present
#' }
#' 
#' 
#' @name SpecifyingHypotheses
NULL

#' @export
#' @rdname SpecifyingHypotheses
Accumulation <- function(...) {
  X <- cbind(...)
  apply(X,1,sum)
}

#' @rdname SpecifyingHypotheses
#' @export
Recency <- function(weights,...) {
  X <- cbind(...)
  if(length(weights)!=ncol(X)) {
    stop("number of weights differs from number of variables") 
  }   
  apply(t(weights*t(X)),1,sum)
}
 
#' @rdname SpecifyingHypotheses
#' @export
Change <- function(...) {
  X <- cbind(...)
  if (ncol(X) != 2)
    stop("requires exactly two variables")
  I(X[,2]-X[,1])
}

#' @rdname SpecifyingHypotheses
#' @export
Mobility <- function(...) {
  X <- cbind(...)
  if (ncol(X) != 2)
    stop("requires exactly two variables")
  X <- cbind((1-X[,1])*X[,2],(1-X[,2])*X[,1])
  colnames(X) <- c("up","down")
  X
}

#' @rdname SpecifyingHypotheses
#' @export
Always <- function(...) {
  X <- cbind(...)
  apply(X,1,prod)
}

#' @rdname SpecifyingHypotheses
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
#' @param formula An object of class \code{formula} that specifies all terms
#' to be included in the dataset.
#' @param data An object of class \code{mids} containing the variables
#' optional data frame, list or environment specified in \code{formula}.
#' @param seed Integer input to \code{set.seed} to ensure that
#' dataset generation is reproducible (Default: 1).
#' @return Data frame of class \code{MI2dummy} including randomly
#' generated variables specified by the specified \code{formula} with
#' means and covariances obtained from the input \code{data}.
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

#' Structured Life Course Modelling Approach: stage 1
#'
#' Performs stage 1 (variable selection) of the SLCMA
#'
#' @param formula An object of class \code{formula} which specifies the
#' SLCMA model including outcome variable, exposures
#' and lifecourse hypotheses.
#' @param data An optional data frame, list or environment
#' containing the variables in the model. If not found in \code{data},
#' the variables are taken from \code{environment(formula)}, typically
#' the environment form which the function is called.
#' @param adjust A numerical vector containing the position of terms to
#' be adjusted for in all models. (Default: NULL - no adjustment).
#' @param seed Integer input to \code{set.seed} to ensure that dataset 
#' generation, if \code{data} contains multiple imputations, is reproducible (Default: 1).
#' @param silent Logical. The function does not print any output when TRUE (Default: FALSE).
#' @param ... Additional arguments to be passed to \code{lars}.
#' @return An list containing the following items:
#' \itemize{
#' \item \code{X_hypos} Numeric matrix with columns corresponding to SLCMA hypotheses.
#' \item \code{y_resid} Output variable, adjusted for terms specified in \code{adjust}.
#' \item \code{multiple} Output from multiple regression of all terms.
#' \item \code{sanity} Data frame containing details of all terms specified by \code{formula}
#' and whether they are included in all models.
#' \item \code{fit} Output of \code{lars}.
#' \item \code{base.r.squared} R-squared for the base model containing only those terms
#' specified in \code{adjust}.
#' \item \code{adjust.df} Degrees of freedom for the terms specified in \code{adjust}.
#' }
#' @export 
slcma <- function(formula, data=environment(formula), adjust=NULL, seed=1, silent=FALSE, ...) {

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
    dummy <- MI2dummy(formula, data, seed=seed)
    y <- dummy[,attributes(dummy)$assign==-1]
    mm <- as.matrix(cbind(`(Intercept)`=rep(1,dim(dummy)[1]), dummy[,attributes(dummy)$assign!=-1]))
    assign <- c(0, attributes(dummy)$assign[attributes(dummy)$assign!=-1])
    r <- assign %in% c(0, adjust)
  } else {
    y <- model.response(mf)
    mm <- model.matrix(trms, data=data)
    r <- attributes(mm)$assign %in% c(0,adjust)
  }
  
  sanity <- data.frame(Term=colnames(mm),
                       Role=factor(r, levels=c(T,F), labels=c("Adjusted for in all models","Available for variable selection")))

  covars  <- mm[, r, drop=FALSE]
  X_hypos <- mm[,!r, drop=FALSE]
  
  adjust_model <- lm(y ~ covars)
  y_resid <- adjust_model$residuals         # even if no covariates, this makes y_resid have mean zero
  X_hypos <- lm(X_hypos ~ covars)$residuals # even if no covariates, this makes every column of X_hypos have mean zero
  base.r.squared <- 1 - sum(y_resid^2)/sum((y-mean(y))^2)
  adjust.df <- adjust_model$rank - 1 
  
  full_model <- lm(y_resid ~ X_hypos)
  df1 <- full_model$rank - 1
  df2 <- full_model$df.residual - adjust.df 
  RSS <- sum(full_model$residuals^2)
  sigma <- sqrt(RSS/df2)
  r.squared <- 1 - RSS/sum(y_resid^2)
  fstatistic <- r.squared / (1-r.squared) * df2/df1
  p <- pf(fstatistic, df1, df2, lower.tail=FALSE)
  multiple <- list(sigma = sigma, r.squared = r.squared, fstatistic = fstatistic, df1 = df1, df2 = df2, p = p)

  fit <- lars(X_hypos, y_resid, ...)

  output <- list(sanity=sanity,
                 base.r.squared=base.r.squared,
                 adjust.df=adjust.df,
                 X_hypos=X_hypos,
                 y_resid=y_resid,
                 multiple=multiple,
                 fit=fit)
  class(output) <- "slcma"
  if(!silent) print(sanity, row.names=FALSE)
  invisible(output)
}

#' SLCMA stage 1 output
#'
#' Print method for objects of class "slcma".
#'
#' @param x Object of class \code{slcma} from the \code{slcma()} function.
#'
#' @export
print.slcma <- function(x) {
  print(x$sanity, row.names=FALSE)
}

#' SLCMA output summary
#' 
#' Summary method for class "slcma".
#'
#' @param x Object of class \code{slcma} from the \code{slcma()} function.
#'
#' @export
summary.slcma <- function(x) {
  unl <- unlist(x$fit$actions)
  varnames <- names(unl)
  selected <- varnames
  removed <- varnames
  selected[unl < 0] <- ""
  removed[unl > 0] <- ""
  output <- x
  output$step <- c(0,seq(along=unl))
  output$selected <- c("",selected)
  output$removed <- c("",removed)
  output$vars <- x$fit$df - 1
  class(output) <- "summary.slcma"
  output
}

#' SLCMA summary output
#'
#' Print method for objects of class \code{summary.slcma}.
#'
#' @param x Object of class \code{summary.slcma} from the \code{summary.slcma()} function.
#'
#' @export
print.summary.slcma <- function(x) {
  cat("\nSummary of LARS procedure\n")
  summarytable <- data.frame(Step=x$step,
                             `Variable selected`=x$selected,
                             `Variable removed`=x$removed,
                             Variables=x$vars, 
                             `R-squared`=round(x$fit$R2, 3), check.names=FALSE)
  if(x$adjust.df > 0) {
    names(summarytable)[5] <- "Partial R-squared"
  }
  print(summarytable, row.names=FALSE)
}

#' SLCMA elbow plot
#'
#' Elbow plot (plot method for class "slcma")
#'
#' @param x Object of class \code{slcma} from the \code{slcma()} function.
#'
#' @export
plot.slcma <- function(x, relax=FALSE, show.remove=FALSE, show.labels=TRUE, labels=NULL, selection.order=FALSE,
                       xlab="Variables selected", ylab=ifelse(x$adjust.df>0,"Partial R-squared","R-squared"), ...) {
  maxadd <- rle(x$fit$actions>0)$lengths[1]
  vars <- x$fit$df - 1
  R2 <- x$fit$R2
  unl <- unlist(x$fit$actions)

  if(relax) {
    for(j in seq_along(unl)) {
      selected <-   seq_len(dim(x$X_hypos)[2]) %in%  unl[1:j]
      selected <- !(seq_len(dim(x$X_hypos)[2]) %in% -unl[1:j]) & selected
      relaxed_model <- lm(x$y_resid ~ x$X_hypos[,selected])
      R2[j+1] <- summary(relaxed_model)$r.sq
    }
  }
  maxR2 <- x$multiple$r.sq

  labs <- attr(unl,"names")
  if(!is.null(labels)) {
    if(selection.order) {
      if(length(labels) > maxadd) {
        stop("'labels' is longer than the number of selected variables")
      }
      labs[seq_along(labels)] <- labels
      print(data.frame(Label=labs, Term=attr(unl,"names")), row.names=FALSE)
    } else {
      if(length(labels) > dim(x$X_hypos)[2]) {
        stop("'labels' is longer than the maximum number of variables")
      }
      use.names <- colnames(x$X_hypos)
      use.names[seq_along(labels)] <- labels
      print(data.frame(Label=use.names, Term=colnames(x$X_hypos)), row.names=FALSE)
      labs <- use.names[unl[1:maxadd]]
    }
  }
  labs[-(1:maxadd)] <- ""

  if(!show.remove) {
    vars <- vars[1:(maxadd+1)]
    R2 <- R2[1:(maxadd+1)]
    labs <- labs[1:maxadd]
  }

  plot(c(0,max(vars)), c(0,maxR2), type="n", xaxt="n", # change 0 to min(vars) and min(R2) (and below in text)
       xlab=xlab, ylab=ylab, ...)
  axis(1, at=vars)
  if(show.labels) {  
    text(vars[-1],0, labs, srt=90, adj=0)
  }
  lines(vars, R2, ...)
  abline(h=maxR2, lty="dotted")
}


#' Structured Life Course Modelling Approach: stage 1
#'
#' Performs stage 2 of the SLCMA for user-specified methods
#'
#' @param x Object of class \code{slcma} from the \code{slcma()} function.
#' @param step Integer specifying which step of the LARS procedure produces 
#' the model on which inference is to be performed.
#' @param method Character string or vector containing the method or methods of inference
#' to be performed. (Default: \code{slcmaFLI} - Fixed Lasso Inference).
#' @param alpha Level of significance for confidence interval calculations. Confidence
#' intervals will have \code{(1 - alpha) * 100}\% coverage. (Default: 0.05 - 95\% coverage).
#' @param do.maxtCI Logical, indicating whether confidence intervals should be 
#' calculated for the max-|t| test. (Default: FALSE).
#' @param ... Additional arguments to \code{fixedLassoInf()} (see the \code{selectiveinference} package).
#' @return An list of class \code{slcmaInfer} with one element providing the output for each inference method.
#' 
#' @export
slcmaInfer <- function(x, step=1L, method="slcmaFLI", alpha=0.05, do.maxtCI=FALSE, ...) {
  if(step<0 | abs(step-round(step)) > .Machine$double.eps^0.5) {
    stop("'step' is not a non-negative integer")
  }
  if(step > length(x$fit$action)) {
    stop("'step' is greater than the number of LARS steps") 
  }

  output <- x
  output$inference <- list(step=step, vars=x$fit$df[step+1]-1, R2=x$fit$R2[step+1])
  if(x$adjust.df > 0) {
    output$inference$totalR2 <- x$base.r.squared + output$inference$R2 * (1 - x$base.r.squared)
  }

  if(any(c("fixedLassoInf","selectiveInference","fli","slcmaFLI") %in% method)) {
    fli <- slcmaFLI(x, step=step, alpha=alpha, ...)
    output$fli <- fli
  }
  if(any(c("maxt","slcmaMaxt") %in% method)) {
    if(step == 0 | step > 1) {
      stop("max-|t| test is only available at Step 1")
    }
    maxt <- slcmaMaxt(x, alpha=alpha, do.CI=do.maxtCI)
    output$maxt <- maxt
  }
  if(any(c("relaxed","relax","slcmaRelax") %in% method)) {
    relax <- slcmaRelax(x, step=step)
    output$relax <- relax
  }
  if(any(c("Bayes","slcmaBayes") %in% method)) {
    if(step == 0 | step > 1) {
      stop("Bayesian posterior probabilities are only available at Step 1")
    }
    Bayes <- slcmaBayes(x)
    output$Bayes <- Bayes
  } 

  class(output) <- "slcmaInfer"
  output
}

#' SLCMA stage 2 output
#'
#' Print method for objects of class \code{slcmaInfer}.
#'
#' @param x Object of class \code{slcmaInfer}
#' from the \code{slcmaInfer()} function.
#'
#' @export
print.slcmaInfer <- function(x) {
  adjust <- x$sanity[,2] == "Adjusted for in all models" & x$sanity[,1] != "(Intercept)"
  adjust.list <- x$sanity[adjust,1]

  cat(sprintf("\nInference for model at Step %1.0f of LARS procedure\n",x$inference$step))
  cat(sprintf("\nNumber of selected variables: %1.0f\n", x$inference$vars))
  if(x$adjust.df > 0) {
    cat(sprintf("Partial R-squared from lasso fit: %.3f\n", x$inference$R2))
    cat(sprintf("  Total R-squared from lasso fit: %.3f\n", x$inference$totalR2))
    cat("Adjusted for:",paste(adjust.list, collapse=", "),"\n")
  }
  else {
    cat(sprintf("R-squared from lasso fit: %.3f\n", x$inference$R2))
  }
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
#' @param x
#' @param step
#' @param alpha
#' @param ...
#' @return
#'
#' @export
slcmaFLI <- function(x, step, alpha=0.05, ...) {
  if(step<0 | abs(step-round(step)) > .Machine$double.eps^0.5) {
    stop("'step' is not a non-negative integer")
  }
  if(step > length(x$fit$action)) {
    stop("'step' is greater than the number of LARS steps") 
  }
  if(step == length(x$fit$action)) {
    stop("Fixed lasso inference and/or selective inference is not available at the final step of the LARS procedure")
  }
  sumsq <- x$fit$normx
  X_normed <- scale(x$X_hypos, scale=sumsq)
  if(step == 0) {
    fli <- list(lambda=x$fit$lambda[step+1])
  }
  else {
    fli <- fixedLassoInf(X_normed, x$y_resid, 
                         x$fit$beta[step+1,], x$fit$lambda[step+1], 
                         type="partial", alpha=alpha, ...)
    sumsq <- sumsq[fli$vars]
    fli$coef0 <- fli$coef0 / sumsq
    fli$sd <- fli$sd / sumsq
    fli$ci <- fli$ci / cbind(sumsq,sumsq)
  }
  class(fli) <- "slcmaFLI"
  fli
}

#' Print slcmaFLI output
#'
#' Print method for class \code{slcmaFLI} (adapted from \code{selectiveinference} package).
#'
#' @param x
#' @param tailarea (Default: TRUE)
#' @param ...
#' @return
#'
#' @export
print.slcmaFLI <- function (x, tailarea = TRUE, ...) 
{
  if(!is.null(x$sigma)) {
    cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n", 
                x$sigma))
  }
  if(!is.null(x$alpha)) {
    cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n", 
                x$lambda, x$alpha))
  }
  else {
    cat(sprintf("\nTesting results at lambda = %0.3f\n", x$lambda))
  }
  cat("", fill = T)
  if(!is.null(x$coef0)) {
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
  }
  else {
    cat("No variables selected\n")
  }
  invisible()
}

#' Max-|t| test
#'
#' Function for max-|t| test and associated confidence intervals
#'
#' @param x
#' @param alpha (Default: 0.05)
#' @param do.CI (Default: FALSE)
#' @param seed (Default: 12345)
#' @param ...
#' @return
#'
#' @export
slcmaMaxt <- function(x, alpha=0.05, do.CI=TRUE, seed=12345, ...) {
  # assumes that X and y have mean subtracted from them, which should have happened by default
  y_resid <- x$y_resid
  X_hypos <- x$X_hypos
  n <- length(y_resid)
  p <- dim(X_hypos)[2]
  d <- x$multiple$df2
  sumsq <- x$fit$normx
  X_normed <- scale(X_hypos, scale=sumsq)
  Xt <- t(X_normed)
  XtX <- Xt %*% X_normed
  Xty <- Xt %*% y_resid
  selection <- x$fit$action[[1]]
  r <- Xty[selection]
  s <- x$multiple$sigma
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
#' Print method for objects of class \code{slcmaMaxt}.
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
slcmaRelax <- function(x, step, alpha=0.05) {
  if(step<0 | abs(step-round(step)) > .Machine$double.eps^0.5) {
    stop("'step' is not a non-negative integer")
  }
  if(step > length(x$fit$action)) {
    stop("'step' is greater than the number of LARS steps") 
  }
  unl <- unlist(x$fit$actions)
  if(step == 0) {
    selected <- 0
    relaxed_model <- lm(x$y_resid ~ 1)
  }
  else {
    selected <-   seq_len(dim(x$X_hypos)[2]) %in%  unl[1:step]
    selected <- !(seq_len(dim(x$X_hypos)[2]) %in% -unl[1:step]) & selected
    relaxed_model <- lm(x$y_resid ~ x$X_hypos[,selected])
  }
  relaxed_summary <- summary(relaxed_model)
  coefs <- relaxed_summary$coefficients
  SEs <- coefs[-1,2] * x$multiple$sigma / relaxed_summary$sigma
  lower <- coefs[-1,1] - qnorm(1-alpha/2)*SEs
  upper <- coefs[-1,1] + qnorm(1-alpha/2)*SEs
  coefficients <- cbind(coefs[-1,1], SEs, lower, upper)
  rownames(coefficients) <- colnames(x$X_hypos)[selected]
  colnames(coefficients) <- c("Coef", "SE", "CI.lo", "CI.up")
  relax <- list(coefficients = coefficients, relaxR2 = relaxed_summary$r.sq)
  if(x$adjust.df > 0) {
    relax$totalR2 <- x$base.r.squared + relax$relaxR2 * (1 - x$base.r.squared)
  }

  class(relax) <- "slcmaRelax"
  relax
}

#' Print slcmaRelax
#'
#' Print method for objects of class \code{slcmaRelax}
#'
#' @param x
#'
#' @export
print.slcmaRelax <- function (x) {
  if(is.null(x$totalR2)) {
    cat(sprintf("\nR-squared from relaxed lasso fit: %.3f\n", x$relaxR2))
  }
  else {
    cat(sprintf("\nPartial R-squared from relaxed lasso fit: %.3f\n", x$relaxR2))
    cat(sprintf("  Total R-squared from relaxed lasso fit: %.3f\n", x$totalR2))
  }
  cat("", fill=TRUE)
  if(dim(x$coefficients)[1]==0) {
    cat("No variables selected\n")
  }
  else {
    print(x$coefficients)
  }
  invisible()
}

slcmaBayes <- function(x) {
  y_resid <- x$y_resid
  X_hypos <- x$X_hypos
  p <- dim(X_hypos)[2]
  sigma <- x$multiple$sigma
  numer <- numeric(p)
  for(j in 1:p) {
    residuals <- lm(y_resid ~ X_hypos[,j])$residuals
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

