test_that("average_SE works as expected", {
  # Set up.
  d <- rep(seq(0.1, 0.9, 0.1), each = 4)
  v <- rep(seq(0.1, 0.4, 0.1), 9)
  df <- S4Vectors::DataFrame(Concentration = d,
                             normalization_type = rep(c("GR", "RV"),
                                                      length(v) * 2),
                             x = rep(v, 2))
  normalized <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = df)

  keys <- list(Trt = "Concentration")
  assays <- list("Normalized" = normalized)
  se <- SummarizedExperiment::SummarizedExperiment(assays = assays)
  se <- gDRutils::set_SE_keys(se, keys)
  se <- gDRutils::set_SE_identifiers(se, gDRutils::get_env_identifiers())

  # With masking.
  se1 <-
    average_SE(
      se,
      data_type = "single-agent",
      normalized_assay = "Normalized",
      averaged_assay = "Averaged"
    )
  avg1 <- SummarizedExperiment::assays(se1)[["Averaged"]][1, 1][[1]]
  expect_true(all(avg1$Concentration == rep(seq(0.1, 0.9, 0.1), each = 2)))
  expect_true(all(avg1$GRvalue == 0.4))
  expect_true(all(avg1$RelativeViability == 0.4))

  # With no masking.
  se2 <-
    average_SE(
      se,
      "single-agent",
      normalized_assay = "Normalized",
      averaged_assay = "Averaged"
    )
  avg2 <- SummarizedExperiment::assays(se2)[["Averaged"]][1, 1][[1]]
  expect_true(all(avg2$Concentration == rep(seq(0.1, 0.9, 0.1), each = 2)))
  expect_true(all(avg2$GRvalue == 0.25))
  expect_true(all(avg2$RelativeViability == 0.25))
})
