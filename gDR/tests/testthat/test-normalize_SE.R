#library(testthat); library(gDR);

test_that("normalize_SE2 works as expected", {
  # Set up. 
  conc <- rep(seq(0.1, 0.3, 0.1), 2)
  ctrl_df <- S4Vectors::DataFrame(Barcode = c(1, 2),
                                  RefReadout = rep(2, 2),
				  Day0Readout = rep(1, 2), 
				  UntrtReadout = rep(6, 2))

  trt_df <- S4Vectors::DataFrame(CorrectedReadout = rep(seq(1, 3, 1), 2), 
                                 Concentration = conc,
				 Barcode = rep(c(1, 2), each = 3), 
				 masked = rep(FALSE, 6))

  coldata <- S4Vectors::DataFrame(cl_name = "Mickey Mouse", 
                                  ref_time = 1)

  rowdata <- S4Vectors::DataFrame(duration = 2)

  ctrl <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = ctrl_df)
  trted <- BumpyMatrix::splitAsBumpyMatrix(row = 1, column = 1, x = trt_df)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list("RawTreated" = trted, "Controls" = ctrl), 
                                                   colData = coldata, 
                                                   rowData = rowdata)
  keys <- list("nested_keys" = "Barcode", 
               "Trt" = "Concentration", 
               "cellline_name" = "cl_name", 
               "cellline_ref_div_time" = "ref_time", 
               "duration" = "duration", 
               "masked_tag" = "masked")

  se <- set_SE_keys(se, keys) 
  se <- normalize_SE2(se)
  normalized <- SummarizedExperiment::assays(se)[["Normalized"]][1, 1][[1]]

  expect_true(is(normalized, "DataFrame"))
  expect_equal(dim(normalized), c(6, 4))
  expect_true(all(colnames(normalized) %in% c("Concentration", "masked", "GRvalue", "RelativeViability")))
  expect_equal(normalized$Concentration, conc)
})
