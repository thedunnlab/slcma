---
title: "Examples demonstrating slcma R functions v0.5.R"
author: Andrew Smith
date: 15 December 2022
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Examples demonstrating slcma R functions v0.5.R}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputnc}
---

## Getting started

Load the library.


```r
library(slcma)
```

```
## Loading required package: lars
```

```
## Warning: package 'lars' was built under R version 4.2.1
```

```
## Loaded lars 1.3
```

```
## Loading required package: selectiveInference
```

```
## Loading required package: glmnet
```

```
## Loading required package: Matrix
```

```
## Loaded glmnet 4.1-6
```

```
## Loading required package: intervals
```

```
## 
## Attaching package: 'intervals'
```

```
## The following object is masked from 'package:Matrix':
## 
##     expand
```

```
## Loading required package: survival
```

```
## Loading required package: adaptMCMC
```

```
## Loading required package: parallel
```

```
## Loading required package: coda
```

```
## Loading required package: MASS
```

```
## Loading required package: mice
```

```
## 
## Attaching package: 'mice'
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```
## The following objects are masked from 'package:base':
## 
##     cbind, rbind
```

```
## Loading required package: mvtnorm
```

## Example 1. SLCMA

Generate some data that matches that seen in Mishra et al (2009).
This is the code I supplied in my 2015 Epidemiology paper.


```r
N <- c(141,20,88,80,40,35,317,367)
BMI <- c(28.7,27.5,29.1,27.8,28.3,27.1,27.6,26.0)
SE <- c(0.5,1.1,0.7,0.7,0.8,0.9,0.3,0.2)
S1 <- c(0,1,0,0,1,1,0,1)
S2 <- c(0,0,1,0,1,0,1,1)
S3 <- c(0,0,0,1,0,1,1,1)
x1 <- rep(S1,N)
x2 <- rep(S2,N)
x3 <- rep(S3,N)
epsilon <- rnorm(sum(N))
e <- lm(epsilon ~ x1 * x2 * x3)$residuals
y <- rep(BMI,N) + rep(SE*sqrt(N),N) * e/sd(e)
```

As in my 2015 Epidemiology paper, perform SLCMA with 3 critical periods (early, middle, late)

* Accumulation

* Early-middle mobility (up and down)

* Middle-late mobility (up and down)

## LARS (stage 1)


```r
myslcma <-
  slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + Mobility(x1,x2) + Mobility(x2,x3))
```

```
##                      Term                             Role
##               (Intercept)       Adjusted for in all models
##                        x1 Available for variable selection
##                        x2 Available for variable selection
##                        x3 Available for variable selection
##  Accumulation(x1, x2, x3) Available for variable selection
##        Mobility(x1, x2)up Available for variable selection
##      Mobility(x1, x2)down Available for variable selection
##        Mobility(x2, x3)up Available for variable selection
##      Mobility(x2, x3)down Available for variable selection
```

Get a summary of the SLMCA, to see what variables are added at which stage.
This is more meaningful that the usual output of lars.


```r
summary(myslcma)
```

```
## 
## Summary of LARS procedure
##  Step        Variable selected Variable removed Variables    R2
##     1 Accumulation(x1, x2, x3)                          1 0.022
##     2                       x1                          2 0.027
##     3     Mobility(x2, x3)down                          3 0.038
##     4                       x3                          4 0.039
##     5       Mobility(x1, x2)up                          5 0.039
##     6                                        x3         4 0.040
##     7     Mobility(x1, x2)down                          5 0.040
```

## Elbow plots

Default plot

```r
plot(myslcma)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

For diagnostic purposes, 
can also show the part of the elbow plot after variables get removed and others added in.


```r
plot(myslcma, show.remove=TRUE, show.labels=FALSE)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

Sometimes it is more helpful to show
an alternative elbow plot with the R2 from the relaxed lasso.


```r
plot(myslcma, relax=TRUE)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

Can add your own labels to the elbow plot

```r
plot(myslcma, labels = c("Early critical period","Middle critical period","Late critical period",
                         "Accumulation","Early-middle mobility (up)","Early-middle mobility (down)",
                         "Middle-late mobility (up)", "Middle-late mobility (down)"))
```

```
##                         Label                     Term
##         Early critical period                       x1
##        Middle critical period                       x2
##          Late critical period                       x3
##                  Accumulation Accumulation(x1, x2, x3)
##    Early-middle mobility (up)       Mobility(x1, x2)up
##  Early-middle mobility (down)     Mobility(x1, x2)down
##     Middle-late mobility (up)       Mobility(x2, x3)up
##   Middle-late mobility (down)     Mobility(x2, x3)down
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

Or if you can't remember all the variable names,
and just want to label the ones that appear on the plot.

```r
plot(myslcma, labels = c("Accumulation","Early critical period","Middle-late mobility (down)",
                         "Late critical period", "Early-middle mobility (up)"), selection.order=TRUE)
```

```
##                        Label                     Term
##                 Accumulation Accumulation(x1, x2, x3)
##        Early critical period                       x1
##  Middle-late mobility (down)     Mobility(x2, x3)down
##         Late critical period                       x3
##   Early-middle mobility (up)       Mobility(x1, x2)up
##                           x3                       x3
##         Mobility(x1, x2)down     Mobility(x1, x2)down
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

The function provides a sanity check so you can compare the default labels with the user-defined labels.

## Inference (stage 2)

Maybe there is an elbow at Step 3, perform selective inference by default

```r
slcmaInfer(myslcma, 3)
```

```
## Warning in lsfit(x, y, intercept = oo): 'X' matrix was collinear
```

```
## 
## Inference for model at Step 3 of LARS procedure
## 
## Number of selected variables: 3
## Lasso R-squared: 0.038
## 
## Results from fixed lasso inference (selective inference):
## 
## Standard deviation of noise (specified or estimated) sigma = 5.255
## 
## Testing results at lambda = 5.553, with alpha = 0.050
## 
##                            Coef P-value  CI.lo  CI.up LoTailArea UpTailArea
## x1                       -0.778   0.158 -1.710  0.862      0.024      0.025
## Accumulation(x1, x2, x3) -0.608   0.011 -1.244 -0.098      0.025      0.024
## Mobility(x2, x3)down      1.196   0.051 -0.296  2.197      0.024      0.024
```

We have the option to report the results of the relaxed lasso instead

```r
slcmaInfer(myslcma, 3, method="relax")
```

```
## 
## Inference for model at Step 3 of LARS procedure
## 
## Number of selected variables: 3
## Lasso R-squared: 0.038
## 
## Results from relaxed lasso:
## 
## Relaxed r-squared: 0.040
## 
##                                Coef        SE      CI.lo      CI.up
## x1                       -0.7781529 0.4660779 -1.6916487  0.1353430
## Accumulation(x1, x2, x3) -0.6082188 0.2327875 -1.0644739 -0.1519637
## Mobility(x2, x3)down      1.1962127 0.5083487  0.1998675  2.1925580
```

We have the option of trying the max-|t| test, but only at Step 1.


```r
slcmaInfer(myslcma, 1, method="maxt")
```

```
## 
## Inference for model at Step 1 of LARS procedure
## 
## Number of selected variables: 1
## Lasso R-squared: 0.022
## 
## Results from max-|t| test:
## 
## Error standard deviation: 5.247 on 1079 degrees of freedom
## 
##                                Coef    P-value
## Accumulation(x1, x2, x3) -0.9641575 6.8507e-09
## 
## Message regarding calculation of P-value: Normal Completion
##  with estimated numeric absolute error: 3.492905e-09
```

We can also get confidence intervals from the max-|t| test, 
but these are time-consuming so they are not produced by default

```r
slcmaInfer(myslcma, 1, method="maxt", do.maxtCI=TRUE)
```

```
## 
## Inference for model at Step 1 of LARS procedure
## 
## Number of selected variables: 1
## Lasso R-squared: 0.022
## 
## Results from max-|t| test:
## 
## Error standard deviation: 5.247 on 1079 degrees of freedom
## Confidence interval coverage 95.0 percent
## 
##                                Coef    P-value     CI.lo     CI.up
## Accumulation(x1, x2, x3) -0.9641575 6.8507e-09 -1.316442 -0.639391
## 
## Message regarding calculation of P-value: Normal Completion
##  with estimated numeric absolute error: 3.492905e-09
```

We have the option of Bayesian inference,
to find the posterior probability of the first-selected variable being the most important

```r
slcmaInfer(myslcma, 1, method="Bayes")
```

```
## 
## Inference for model at Step 1 of LARS procedure
## 
## Number of selected variables: 1
## Lasso R-squared: 0.022
## 
## Posterior probabilities:
##                                    
## x1                       0.01526374
## x2                       0.00000055
## x3                       0.00254404
## Accumulation(x1, x2, x3) 0.98218849
## Mobility(x1, x2)up       0.00000037
## Mobility(x1, x2)down     0.00000001
## Mobility(x2, x3)up       0.00000001
## Mobility(x2, x3)down     0.00000279
```

It is possible to produce results for all inference methods

```r
slcmaInfer(myslcma, 1, method=c("selectiveInference","relax","maxt","Bayes"))
```

```
## Warning in lsfit(x, y, intercept = oo): 'X' matrix was collinear
```

```
## 
## Inference for model at Step 1 of LARS procedure
## 
## Number of selected variables: 1
## Lasso R-squared: 0.022
## 
## Results from fixed lasso inference (selective inference):
## 
## Standard deviation of noise (specified or estimated) sigma = 5.255
## 
## Testing results at lambda = 18.578, with alpha = 0.050
## 
##                            Coef P-value  CI.lo  CI.up LoTailArea UpTailArea
## Accumulation(x1, x2, x3) -0.964       0 -1.277 -0.626      0.024      0.024
## 
## Results from max-|t| test:
## 
## Error standard deviation: 5.247 on 1079 degrees of freedom
## 
##                                Coef    P-value
## Accumulation(x1, x2, x3) -0.9641575 6.8507e-09
## 
## Message regarding calculation of P-value: Normal Completion
##  with estimated numeric absolute error: 3.492905e-09
## 
## Results from relaxed lasso:
## 
## Relaxed r-squared: 0.033
## 
##                                Coef        SE     CI.lo      CI.up
## Accumulation(x1, x2, x3) -0.9641575 0.1575558 -1.272961 -0.6553539
## 
## Posterior probabilities:
##                                    
## x1                       0.01526374
## x2                       0.00000055
## x3                       0.00254404
## Accumulation(x1, x2, x3) 0.98218849
## Mobility(x1, x2)up       0.00000037
## Mobility(x1, x2)down     0.00000001
## Mobility(x2, x3)up       0.00000001
## Mobility(x2, x3)down     0.00000279
```

## Example 2. Adjusting for covariates

Generate some data
(this is the code I supplied in my 2015 IJE paper).

```r
n <- 400
set.seed(1234)
covariate <- rnorm(n)
x1 <- covariate + rnorm(n)
x2 <- x1 + rnorm(n)
x3 <- x2 + 2*rnorm(n)
y <- 2*covariate + x3 + 3*rnorm(n)
```


## Adjust for covariates (stage 1)


Add any covariates into the formula,
then use the option adjust = a numeric vector containing the position(s) of covariate(s) in the formula.
Here there is one covariate (called 'covariate') and it is in the 8th position after the ~.
Hence adjust=8.


```r
slcmaAdjust <-
  slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + Change(x1,x2) + Change(x2,x3) + Change(x1,x3) + covariate,
        adjust=8)
```

```
##                      Term                             Role
##               (Intercept)       Adjusted for in all models
##                        x1 Available for variable selection
##                        x2 Available for variable selection
##                        x3 Available for variable selection
##  Accumulation(x1, x2, x3) Available for variable selection
##            Change(x1, x2) Available for variable selection
##            Change(x2, x3) Available for variable selection
##            Change(x1, x3) Available for variable selection
##                 covariate       Adjusted for in all models
```

The function provides a sanity check so you can check the variables have the desired roles

You can revisit the sanity check

```r
slcmaAdjust$sanity
```

```
##                       Term                             Role
## 1              (Intercept)       Adjusted for in all models
## 2                       x1 Available for variable selection
## 3                       x2 Available for variable selection
## 4                       x3 Available for variable selection
## 5 Accumulation(x1, x2, x3) Available for variable selection
## 6           Change(x1, x2) Available for variable selection
## 7           Change(x2, x3) Available for variable selection
## 8           Change(x1, x3) Available for variable selection
## 9                covariate       Adjusted for in all models
```

Can summarize and produce plots as before

```r
summary(slcmaAdjust) 
```

```
## 
## Summary of LARS procedure
##  Step        Variable selected Variable removed Variables    R2
##     1                       x3                          1 0.316
##     2 Accumulation(x1, x2, x3)                          2 0.359
##     3           Change(x1, x2)                          3 0.360
```

```r
plot(slcmaAdjust)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

Can perform inference (stage 2) as before

```r
slcmaInfer(slcmaAdjust,1, method=c("selectiveInference","relax","maxt","Bayes"))
```

```
## Warning in lsfit(x, y, intercept = oo): 'X' matrix was collinear
```

```
## 
## Inference for model at Step 1 of LARS procedure
## 
## Number of selected variables: 1
## Lasso R-squared: 0.316
## 
## Results from fixed lasso inference (selective inference):
## 
## Standard deviation of noise (specified or estimated) sigma = 3.106
## 
## Testing results at lambda = 15.514, with alpha = 0.050
## 
##     Coef P-value CI.lo CI.up LoTailArea UpTailArea
## x3 0.964       0 0.835 1.093      0.024      0.024
## 
## Results from max-|t| test:
## 
## Error standard deviation: 3.094 on 391 degrees of freedom
## 
##         Coef P-value
## x3 0.9637339       0
## 
## Message regarding calculation of P-value: Normal Completion
##  with estimated numeric absolute error: 0.000000e+00
## 
## Results from relaxed lasso:
## 
## Relaxed r-squared: 0.357
## 
##         Coef         SE     CI.lo    CI.up
## x3 0.9637339 0.06493734 0.8364591 1.091009
## 
## Posterior probabilities:
##                                   
## x1                       0.0000000
## x2                       0.0000000
## x3                       0.9999993
## Accumulation(x1, x2, x3) 0.0000007
## Change(x1, x2)           0.0000000
## Change(x2, x3)           0.0000000
## Change(x1, x3)           0.0000000
```

## Example 3. Covariates and missing values

Some simulated data with missing values
(this is the code I supplied to The Dunn Lab in 2016)


```r
set.seed(123)
covariate <- rnorm(250)
x1 <- covariate + rnorm(250)
x2 <- x1 + 2*rnorm(250)
x3 <- x2 + 2*rnorm(250)
y <- 2*covariate + 2*(x3-x1) + 5*rnorm(250)
covariate[sample(1:250, size=10)] <- NA
x1[sample(1:250, size=10)] <- NA
x2[sample(1:250, size=10)] <- NA
x3[sample(1:250, size=10)] <- NA
y[sample(1:250, size=10)] <- NA
incompleteData <- data.frame(covariate,x1,x2,x3,y)
```

Complete case analysis

```r
CCslcma <- slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + 
                 Change(x1,x3) + Ever(x1,x2,x3) + covariate, incompleteData, adjust=7)
```

```
##                      Term                             Role
##               (Intercept)       Adjusted for in all models
##                        x1 Available for variable selection
##                        x2 Available for variable selection
##                        x3 Available for variable selection
##  Accumulation(x1, x2, x3) Available for variable selection
##            Change(x1, x3) Available for variable selection
##          Ever(x1, x2, x3) Available for variable selection
##                 covariate       Adjusted for in all models
```

Impute missing values (5 imputated datasets to save time)

```r
imputedData <- mice(incompleteData, m=5, maxit=20, seed=321)
```

```
## 
##  iter imp variable
##   1   1  covariate  x1  x2  x3  y
##   1   2  covariate  x1  x2  x3  y
##   1   3  covariate  x1  x2  x3  y
##   1   4  covariate  x1  x2  x3  y
##   1   5  covariate  x1  x2  x3  y
##   2   1  covariate  x1  x2  x3  y
##   2   2  covariate  x1  x2  x3  y
##   2   3  covariate  x1  x2  x3  y
##   2   4  covariate  x1  x2  x3  y
##   2   5  covariate  x1  x2  x3  y
##   3   1  covariate  x1  x2  x3  y
##   3   2  covariate  x1  x2  x3  y
##   3   3  covariate  x1  x2  x3  y
##   3   4  covariate  x1  x2  x3  y
##   3   5  covariate  x1  x2  x3  y
##   4   1  covariate  x1  x2  x3  y
##   4   2  covariate  x1  x2  x3  y
##   4   3  covariate  x1  x2  x3  y
##   4   4  covariate  x1  x2  x3  y
##   4   5  covariate  x1  x2  x3  y
##   5   1  covariate  x1  x2  x3  y
##   5   2  covariate  x1  x2  x3  y
##   5   3  covariate  x1  x2  x3  y
##   5   4  covariate  x1  x2  x3  y
##   5   5  covariate  x1  x2  x3  y
##   6   1  covariate  x1  x2  x3  y
##   6   2  covariate  x1  x2  x3  y
##   6   3  covariate  x1  x2  x3  y
##   6   4  covariate  x1  x2  x3  y
##   6   5  covariate  x1  x2  x3  y
##   7   1  covariate  x1  x2  x3  y
##   7   2  covariate  x1  x2  x3  y
##   7   3  covariate  x1  x2  x3  y
##   7   4  covariate  x1  x2  x3  y
##   7   5  covariate  x1  x2  x3  y
##   8   1  covariate  x1  x2  x3  y
##   8   2  covariate  x1  x2  x3  y
##   8   3  covariate  x1  x2  x3  y
##   8   4  covariate  x1  x2  x3  y
##   8   5  covariate  x1  x2  x3  y
##   9   1  covariate  x1  x2  x3  y
##   9   2  covariate  x1  x2  x3  y
##   9   3  covariate  x1  x2  x3  y
##   9   4  covariate  x1  x2  x3  y
##   9   5  covariate  x1  x2  x3  y
##   10   1  covariate  x1  x2  x3  y
##   10   2  covariate  x1  x2  x3  y
##   10   3  covariate  x1  x2  x3  y
##   10   4  covariate  x1  x2  x3  y
##   10   5  covariate  x1  x2  x3  y
##   11   1  covariate  x1  x2  x3  y
##   11   2  covariate  x1  x2  x3  y
##   11   3  covariate  x1  x2  x3  y
##   11   4  covariate  x1  x2  x3  y
##   11   5  covariate  x1  x2  x3  y
##   12   1  covariate  x1  x2  x3  y
##   12   2  covariate  x1  x2  x3  y
##   12   3  covariate  x1  x2  x3  y
##   12   4  covariate  x1  x2  x3  y
##   12   5  covariate  x1  x2  x3  y
##   13   1  covariate  x1  x2  x3  y
##   13   2  covariate  x1  x2  x3  y
##   13   3  covariate  x1  x2  x3  y
##   13   4  covariate  x1  x2  x3  y
##   13   5  covariate  x1  x2  x3  y
##   14   1  covariate  x1  x2  x3  y
##   14   2  covariate  x1  x2  x3  y
##   14   3  covariate  x1  x2  x3  y
##   14   4  covariate  x1  x2  x3  y
##   14   5  covariate  x1  x2  x3  y
##   15   1  covariate  x1  x2  x3  y
##   15   2  covariate  x1  x2  x3  y
##   15   3  covariate  x1  x2  x3  y
##   15   4  covariate  x1  x2  x3  y
##   15   5  covariate  x1  x2  x3  y
##   16   1  covariate  x1  x2  x3  y
##   16   2  covariate  x1  x2  x3  y
##   16   3  covariate  x1  x2  x3  y
##   16   4  covariate  x1  x2  x3  y
##   16   5  covariate  x1  x2  x3  y
##   17   1  covariate  x1  x2  x3  y
##   17   2  covariate  x1  x2  x3  y
##   17   3  covariate  x1  x2  x3  y
##   17   4  covariate  x1  x2  x3  y
##   17   5  covariate  x1  x2  x3  y
##   18   1  covariate  x1  x2  x3  y
##   18   2  covariate  x1  x2  x3  y
##   18   3  covariate  x1  x2  x3  y
##   18   4  covariate  x1  x2  x3  y
##   18   5  covariate  x1  x2  x3  y
##   19   1  covariate  x1  x2  x3  y
##   19   2  covariate  x1  x2  x3  y
##   19   3  covariate  x1  x2  x3  y
##   19   4  covariate  x1  x2  x3  y
##   19   5  covariate  x1  x2  x3  y
##   20   1  covariate  x1  x2  x3  y
##   20   2  covariate  x1  x2  x3  y
##   20   3  covariate  x1  x2  x3  y
##   20   4  covariate  x1  x2  x3  y
##   20   5  covariate  x1  x2  x3  y
```

Apply LARS to data with same covariance structure as pooled data across imputations 
 (data created via Cholesky decomposition)
 

```r
MIslcma <- slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + 
                 Change(x1,x3) + Ever(x1,x2,x3) + covariate, data=imputedData, adjust=7)
```

```
##                      Term                             Role
##               (Intercept)       Adjusted for in all models
##                        x1 Available for variable selection
##                        x2 Available for variable selection
##                        x3 Available for variable selection
##          Ever(x1, x2, x3) Available for variable selection
##                 covariate       Adjusted for in all models
##  Accumulation(x1, x2, x3) Available for variable selection
##            Change(x1, x3) Available for variable selection
```


Can summarize and produce plots as before

```r
summary(MIslcma) 
```

```
## 
## Summary of LARS procedure
##  Step Variable selected Variable removed Variables    R2
##     1    Change(x1, x3)                          1 0.544
##     2                x2                          2 0.564
##     3                x1                          3 0.565
##     4  Ever(x1, x2, x3)                          4 0.565
```

```r
plot(MIslcma)
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)

Can perform inference (stage 2) as before

```r
slcmaInfer(MIslcma,1, method=c("selectiveInference","relax","maxt","Bayes"))
```

```
## Warning in lsfit(x, y, intercept = oo): 'X' matrix was collinear
```

```
## 
## Inference for model at Step 1 of LARS procedure
## 
## Number of selected variables: 1
## Lasso R-squared: 0.544
## 
## Results from fixed lasso inference (selective inference):
## 
## Standard deviation of noise (specified or estimated) sigma = 5.205
## 
## Testing results at lambda = 15.814, with alpha = 0.050
## 
##                 Coef P-value CI.lo CI.up LoTailArea UpTailArea
## Change(x1, x3) 2.027       0 1.801 2.252      0.024      0.025
## 
## Results from max-|t| test:
## 
## Error standard deviation: 5.195 on 242 degrees of freedom
## 
##                    Coef P-value
## Change(x1, x3) 2.026825       0
## 
## Message regarding calculation of P-value: Normal Completion
##  with estimated numeric absolute error: 0.000000e+00
## 
## Results from relaxed lasso:
## 
## Relaxed r-squared: 0.561
## 
##                    Coef        SE   CI.lo   CI.up
## Change(x1, x3) 2.026825 0.1142905 1.80282 2.25083
## 
## Posterior probabilities:
##                           
## x1                       0
## x2                       0
## x3                       0
## Ever(x1, x2, x3)         0
## Accumulation(x1, x2, x3) 0
## Change(x1, x3)           1
```

Can recreate the dummy dataset used within the above SLCMA

```r
myMIdata <- MI2dummy(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + 
                     Change(x1,x3) + Ever(x1,x2,x3) + covariate, imputedData)
```

Note how sums, sums of squares and cross products are the same for `Accumulation(x1,x2,x3)` and `x1+x2+x3`.


```r
sum(myMIdata$`Accumulation(x1, x2, x3)`)
```

```
## [1] 27.47647
```

```r
sum(myMIdata$x1+myMIdata$x2+myMIdata$x3)
```

```
## [1] 27.47647
```

```r
sum(myMIdata$`Accumulation(x1, x2, x3)`^2)
```

```
## [1] 8938.248
```

```r
sum((myMIdata$x1+myMIdata$x2+myMIdata$x3)^2)
```

```
## [1] 8938.248
```

```r
sum(myMIdata$x1 * myMIdata$`Accumulation(x1, x2, x3)`)
```

```
## [1] 1371.406
```

```r
sum(myMIdata$x1 * (myMIdata$x1+myMIdata$x2+myMIdata$x3))
```

```
## [1] 1371.406
```

```r
sum(myMIdata$y * myMIdata$`Accumulation(x1, x2, x3)`)
```

```
## [1] 7591.831
```

```r
sum(myMIdata$y * (myMIdata$x1+myMIdata$x2+myMIdata$x3))
```

```
## [1] 7591.831
```