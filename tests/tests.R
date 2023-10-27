if(!require(testthat)) install.packages("testthat")
library(testthat)

test_that("Check lifecourse hypothesis functions", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  expect_identical(Accumulation(testX1,testX2), testX1+testX2)
  expect_identical(Recency(1:2,testX1,testX2), testX1+2*testX2)
  expect_error(Recency(1:3,testX1,testX2))
  expect_identical(Change(testX1,testX2), I(testX2-testX1))
  expect_identical(Mobility(testX1,testX2)[,1], (1-testX1)*testX2)
  expect_identical(Mobility(testX1,testX2)[,2], (1-testX2)*testX1)
  expect_identical(Always(testX1,testX2), testX1*testX2)
  expect_identical(Ever(testX1,testX2), c(0,0,1,1,1,1,1,1,1))
})

test_that("Check slcma function", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
  expect_error(slcma(~testX1 + testX2))
  expect_error(slcma(cbind(testY,testY) ~ testX1 + testX2))
  expect_error(slcma(testY ~ 1, silent=TRUE))
  expect_equal(test_slcma$base.r.squared, 0)
  expect_identical(test_slcma$adjust.df, 0)
  expect_equal(dim(test_slcma$X_hypos), c(9,3))
  expect_length(test_slcma$y_resid, 9)
  expect_gte(test_slcma$multiple$sigma, 0)
  expect_gte(test_slcma$multiple$r.squared, 0)
  expect_lte(test_slcma$multiple$r.squared, 1)
  expect_gte(test_slcma$multiple$fstatistic, 0)
  expect_identical(test_slcma$multiple$df1, 2)
  expect_identical(test_slcma$multiple$df2, 6)
  expect_gte(test_slcma$multiple$p, 0)
  expect_lte(test_slcma$multiple$p, 1)
  expect_s3_class(test_slcma$fit, "lars")
})

test_that("Check adjusting for covariates", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  testC1 <- c(0,0,0,1,0,0,1,1,1)
  testC2 <- c(0,1,0,1,0,0,1,1,1)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2) + 
    testC1 + testC2, adjust=4:5, silent=TRUE)
  expect_gt(test_slcma$base.r.squared, 0)
  expect_identical(test_slcma$adjust.df, 2)
  expect_equal(dim(test_slcma$X_hypos), c(9,3))
  expect_error(slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2) +
    testC1 + testC2, adjust=1:5, silent=TRUE))
})

test_that("Check summary", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  testC1 <- c(0,0,0,1,0,0,1,1,1)
  testC2 <- c(0,1,0,1,0,0,1,1,1)
  for(check in 1:2) {
    if(check==1) test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
    if(check==2) test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2) + 
                   testC1 + testC2, adjust=4:5, silent=TRUE)
    summary_output <- summary(test_slcma)
    expect_identical(summary_output$sanity, test_slcma$sanity)
    expect_identical(summary_output$base.r.squared, test_slcma$base.r.squared)
    expect_identical(summary_output$adjust.df, test_slcma$adjust.df)
    expect_identical(summary_output$X_hypos, test_slcma$X_hypos)
    expect_identical(summary_output$y_resid, test_slcma$y_resid)
    expect_identical(summary_output$multiple, test_slcma$multiple)
    expect_identical(summary_output$fit, test_slcma$fit)
    expect_vector(summary_output$step)
    expect_identical(summary_output$step[1], 0)
    expect_vector(summary_output$selected)
    expect_identical(summary_output$selected[1],"")
    expect_vector(summary_output$removed)
    expect_identical(summary_output$removed[1:2],c("",""))
    expect_vector(summary_output$vars)
    expect_equal(as.numeric(summary_output$vars[1:2]), 0:1)
  }
})

test_that("Check plot", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
  expect_silent(plot(test_slcma))
  expect_silent(plot(test_slcma, relax=TRUE, show.remove=TRUE, show.labels=FALSE))
  expect_output(plot(test_slcma, labels=c("a","b")))
  expect_output(plot(test_slcma, labels="a", selection.order=TRUE))
})

test_that("Check fixed lasso inference", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
  expect_error(slcmaFLI(test_slcma, -1))
  expect_error(slcmaFLI(test_slcma, 0.5))
  expect_error(slcmaFLI(test_slcma, 3))
  FLI_output <- slcmaFLI(test_slcma, 1)
  expect_length(FLI_output$coef, 1)
  expect_length(FLI_output$ci, 2)
  expect_gte(FLI_output$pv, 0)
  expect_lte(FLI_output$pv, 1)
})
# Gives warning because of issues in selectiveInference package

test_that("Check max-|t| test", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
  Maxt_output <- slcmaMaxt(test_slcma)
  expect_length(Maxt_output$coefficients, 4)
  expect_gte(Maxt_output$error, 0)
  expect_identical(Maxt_output$msg, "Normal Completion")
  Maxt_output <- slcmaMaxt(test_slcma, do.CI = FALSE)
  expect_length(Maxt_output$coefficients, 2)
})

test_that("Check relaxed lasso", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
  expect_error(slcmaRelax(test_slcma, -1))
  expect_error(slcmaRelax(test_slcma, 0.5))
  expect_error(slcmaRelax(test_slcma, 3))
  Relax_output <- slcmaRelax(test_slcma, 1)
  expect_length(Relax_output$coefficients, 4)
  expect_gte(Relax_output$relaxR2, 0)
  Relax_output <- slcmaRelax(test_slcma, 0)
  expect_equal(Relax_output$relaxR2, 0)
})

test_that("Check Bayesian posterior probabilities", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
  Bayes_output <- slcmaBayes(test_slcma)
  expect_length(Bayes_output$probs, 3)
})

test_that("Check wrapper for SLCMA inference", {
  testX1 <- c(0,0,0,0,1,1,1,1,1)
  testX2 <- c(0,0,1,1,0,0,1,1,0)
  testY <- c(3.8,4.9,4.9,5.0,5.0,5.2,5.3,5.9,6.6)
  test_slcma <- slcma(testY ~ testX1 + testX2 + Accumulation(testX1,testX2), silent=TRUE)
  expect_error(slcmaInfer(test_slcma, -1))
  expect_error(slcmaInfer(test_slcma, 0.5))
  expect_error(slcmaInfer(test_slcma, 3))
  Infer_output <- slcmaInfer(test_slcma, 1, method=c("fli","maxt","relax","Bayes"))
  expect_identical(Infer_output$sanity, test_slcma$sanity)
  expect_identical(Infer_output$base.r.squared, test_slcma$base.r.squared)
  expect_identical(Infer_output$adjust.df, test_slcma$adjust.df)
  expect_identical(Infer_output$X_hypos, test_slcma$X_hypos)
  expect_identical(Infer_output$y_resid, test_slcma$y_resid)
  expect_identical(Infer_output$multiple, test_slcma$multiple)
  expect_identical(Infer_output$fit, test_slcma$fit)
  expect_identical(Infer_output$fli, slcmaFLI(test_slcma, 1))
  expect_identical(Infer_output$maxt, slcmaMaxt(test_slcma, do.CI=FALSE))
  expect_identical(Infer_output$relax, slcmaRelax(test_slcma, 1))
  expect_identical(Infer_output$Bayes, slcmaBayes(test_slcma))
})
# Gives warning becuase of issues in selectiveInference package


test_that("Check multiple imputation", {
  test_imputed <- structure(
    list(
      structure(list(c(-0.5,-0.2, 1.5, 0.0,  NA, 1.7, 0.4,-1.2,-0.6,-0.4),
                     c(-0.9,-0.7, 1.2, 0.1, 1.7,  NA, 1.5,-0.6,-0.8,-1.9), 
                     c(-2.1,-2.7,  NA, 1.6,-1.2, 1.4,-0.2,-4.7,-0.5,-2.1), 
                     c(  NA,-2.9, 4.2, 2.0,-1.6, 1.1, 1.7,-5.1,-4.5,-2.5), 
                     c(-2.3, -10, 9.1, 3.3, -19, 7.7, 2.6, 0.4,  NA,-4.2)), 
                names = c("testC1", "testX1", "testX2", "testC2", "testY"), 
                row.names = c(NA, 10L), 
                class = "data.frame"), 
      structure(list(structure(list(-0.5, 0.4), names = c("1","2"), row.names = "5", class = "data.frame"),
                     structure(list(-1.9,-0.7), names = c("1","2"), row.names = "6", class = "data.frame"),
                     structure(list( 3.2, 2.3), names = c("1","2"), row.names = "3", class = "data.frame"), 
                     structure(list( 0.9, 9.0), names = c("1","2"), row.names = "1", class = "data.frame"), 
                     structure(list(-5.4,-4.5), names = c("1","2"), row.names = "9", class = "data.frame")),
                names = c("testC1", "testX1", "testX2", "testC2", "testY")), 
      2, 
      structure(c(F, F, F, F, T, F, F, F, F, F,
                  F, F, F, F, F, T, F, F, F, F, 
                  F, F, T, F, F, F, F, F, F, F, 
                  T, F, F, F, F, F, F, F, F, F,
                  F, F, F, F, F, F, F, F, T, F),
                dim = c(10L, 5L), 
                dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                                c("testC1", "testX1", "testX2", "testC2", "testY")))), 
    names = c("data", "imp", "m", "where"),
    class = "mids")
  test_formula <-   testY ~ testC1 + testX1 + testX2 + Accumulation(testX1,testX2) + 
                              Change(testX1,testX2) + Mobility(testX1,testX2) + testC2
  test_dummy <- MI2dummy(test_formula, test_imputed)
  expect_equal(dim(test_dummy), c(10,9))
  expect_equal(sum(test_dummy$`Accumulation(testX1, testX2)`), sum(test_dummy$testX1+test_dummy$testX2))
  expect_equal(sum(test_dummy$`Accumulation(testX1, testX2)`^2), sum((test_dummy$testX1+test_dummy$testX2)^2))
  expect_equal(sum(test_dummy$`Accumulation(testX1, testX2)`*test_dummy$testX1), 
                    sum((test_dummy$testX1+test_dummy$testX2)*test_dummy$testX1))
  expect_equal(sum(test_dummy$`Accumulation(testX1, testX2)`*test_dummy$testY), 
                    sum((test_dummy$testX1+test_dummy$testX2)*test_dummy$testY))
  expect_error(MI2dummy(test_formula, list(test_imputed)))
  test_slcma <- slcma(test_formula, adjust=c(1,7), data=test_imputed, silent=TRUE)
  test_sanity <- test_slcma$sanity
  expect_identical(as.character(test_sanity$Role[test_sanity$Term=="testC1"]), "Adjusted for in all models")
  expect_identical(as.character(test_sanity$Role[test_sanity$Term=="testC2"]), "Adjusted for in all models")
  expect_equal(sum(test_sanity$Role == "Adjusted for in all models"), 3)
  expect_equal(dim(test_slcma$X_hypos), c(10,6))
})



