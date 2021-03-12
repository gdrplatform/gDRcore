library(testthat); library(gDR)

test_that("get_SE_experiment_metadata works as expected", {
  exp_md <- list("Super" = "Star", "Serena" = "Williams")
  se <- SummarizedExperiment::SummarizedExperiment(metadata = list(experiment_metadata = exp_md))
  oexp_md <- get_SE_experiment_metadata(se)
  expect_equal(oexp_md, exp_md)
})


test_that("get_SE_keys and set_SE_keys works as expected", {
  keys <- list(Keys = list(Day0 = "TEST", Other = "STUFF"))
  se <- SummarizedExperiment::SummarizedExperiment(metadata = keys)
  nkeys <- get_SE_keys(se)
  expect_equal(nkeys$Day0, "TEST")
  expect_equal(nkeys$Other, "STUFF")

  # Test for all keys.
  keys2 <- list("test" = "NEW")
  se <- set_SE_keys(se, keys2)
  nkeys2 <- get_SE_keys(se)
  expect_equal(nkeys2$test, "NEW")
  expect_equal(length(nkeys2), 1)
})
