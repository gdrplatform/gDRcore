library(testthat); library(gDRcore);

test_that("fit_SE errors as expected", {
  se <- SummarizedExperiment::SummarizedExperiment()
  expect_error(fit_SE(se))
})
