test_that("average_SE works as expected", {
  # Set up.
  d <- rep(seq(0.1, 0.9, 0.1), each = 4)
  v <- rep(seq(0.1, 0.4, 0.1), 9)
  df <- S4Vectors::DataFrame(Concentration = d,
                             masked = rep(c(TRUE, TRUE, TRUE, FALSE), 9),
                             GRvalue = v,
                             RelativeViability = v)
  normalized <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = df)

  keys <- list(Trt = "Concentration",
               "masked_tag" = "masked")
  assays <- list("Normalized" = normalized)
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays)
  se <- gDRutils::set_SE_keys(se, keys)

  # With masking.
  se1 <-
    average_SE(
      se,
      data_type = "single-agent",
      override_masked = FALSE,
      normalized_assay = "Normalized",
      averaged_assay = "Averaged"
    )
  avg1 <- SummarizedExperiment::assays(se1)[["Averaged"]][1, 1][[1]]
  expect_true(all(avg1$Concentration == seq(0.1, 0.9, 0.1)))
  expect_true(all(avg1$GRvalue == 0.4))
  expect_true(all(avg1$RelativeViability == 0.4))

  # With no masking.
  se2 <-
    average_SE(
      se,
      "single-agent",
      override_masked = TRUE,
      normalized_assay = "Normalized",
      averaged_assay = "Averaged"
    )
  avg2 <- SummarizedExperiment::assays(se2)[["Averaged"]][1, 1][[1]]
  expect_true(all(avg2$Concentration == seq(0.1, 0.9, 0.1)))
  expect_true(all(avg2$GRvalue == 0.25))
  expect_true(all(avg2$RelativeViability == 0.25))
})
