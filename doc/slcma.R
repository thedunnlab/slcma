## -----------------------------------------------------------------------------
library(slcma)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
myslcma <-
  slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + Mobility(x1,x2) + Mobility(x2,x3))

## -----------------------------------------------------------------------------
summary(myslcma)

## -----------------------------------------------------------------------------
plot(myslcma)

## -----------------------------------------------------------------------------
plot(myslcma, show.remove=TRUE, show.labels=FALSE)

## -----------------------------------------------------------------------------
plot(myslcma, relax=TRUE)

## -----------------------------------------------------------------------------
plot(myslcma, labels = c("Early critical period","Middle critical period","Late critical period",
                         "Accumulation","Early-middle mobility (up)","Early-middle mobility (down)",
                         "Middle-late mobility (up)", "Middle-late mobility (down)"))

## -----------------------------------------------------------------------------
plot(myslcma, labels = c("Accumulation","Early critical period","Middle-late mobility (down)",
                         "Late critical period", "Early-middle mobility (up)"), selection.order=TRUE)

## -----------------------------------------------------------------------------
slcmaInfer(myslcma, 3)

## -----------------------------------------------------------------------------
slcmaInfer(myslcma, 3, method="relax")

## -----------------------------------------------------------------------------
slcmaInfer(myslcma, 1, method="maxt")

## -----------------------------------------------------------------------------
slcmaInfer(myslcma, 1, method="maxt", do.maxtCI=TRUE)

## -----------------------------------------------------------------------------
slcmaInfer(myslcma, 1, method="Bayes")

## -----------------------------------------------------------------------------
slcmaInfer(myslcma, 1, method=c("selectiveInference","relax","maxt","Bayes"))

## -----------------------------------------------------------------------------
n <- 400
set.seed(1234)
covariate <- rnorm(n)
x1 <- covariate + rnorm(n)
x2 <- x1 + rnorm(n)
x3 <- x2 + 2*rnorm(n)
y <- 2*covariate + x3 + 3*rnorm(n)

## -----------------------------------------------------------------------------
slcmaAdjust <-
  slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + Change(x1,x2) + Change(x2,x3) + Change(x1,x3) + covariate,
        adjust=8)

## -----------------------------------------------------------------------------
slcmaAdjust$sanity

## -----------------------------------------------------------------------------
summary(slcmaAdjust) 
plot(slcmaAdjust)

## -----------------------------------------------------------------------------
slcmaInfer(slcmaAdjust,1, method=c("selectiveInference","relax","maxt","Bayes"))

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
CCslcma <- slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + 
                 Change(x1,x3) + Ever(x1,x2,x3) + covariate, incompleteData, adjust=7)

## -----------------------------------------------------------------------------
imputedData <- mice(incompleteData, m=5, maxit=20, seed=321)

## -----------------------------------------------------------------------------
MIslcma <- slcma(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + 
                 Change(x1,x3) + Ever(x1,x2,x3) + covariate, data=imputedData, adjust=7)

## -----------------------------------------------------------------------------
summary(MIslcma) 
plot(MIslcma)

## -----------------------------------------------------------------------------
slcmaInfer(MIslcma,1, method=c("selectiveInference","relax","maxt","Bayes"))

## -----------------------------------------------------------------------------
myMIdata <- MI2dummy(y ~ x1 + x2 + x3 + Accumulation(x1,x2,x3) + 
                     Change(x1,x3) + Ever(x1,x2,x3) + covariate, imputedData)

## -----------------------------------------------------------------------------
sum(myMIdata$`Accumulation(x1, x2, x3)`)
sum(myMIdata$x1+myMIdata$x2+myMIdata$x3)

sum(myMIdata$`Accumulation(x1, x2, x3)`^2)
sum((myMIdata$x1+myMIdata$x2+myMIdata$x3)^2)

sum(myMIdata$x1 * myMIdata$`Accumulation(x1, x2, x3)`)
sum(myMIdata$x1 * (myMIdata$x1+myMIdata$x2+myMIdata$x3))

sum(myMIdata$y * myMIdata$`Accumulation(x1, x2, x3)`)
sum(myMIdata$y * (myMIdata$x1+myMIdata$x2+myMIdata$x3))

