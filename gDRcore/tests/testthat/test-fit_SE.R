library(testthat); library(gDRcore);

test_that("fit_SE2 errors as expected", {
  se <- SummarizedExperiment::SummarizedExperiment()
  expect_error(fit_SE2(se))
})
