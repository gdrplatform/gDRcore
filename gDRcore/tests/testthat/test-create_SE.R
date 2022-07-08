test_that("validate_mapping works as expected", {
  ref_df <- data.frame(Gnumber = paste0("DRUG_", c(rep(10, 3), rep(11, 7))),
                       Gnumber2 = "vehicle", value = runif(10))
  trt_df <- data.frame(Gnumber = paste0("DRUG_", 10),
                       Gnumber2 = paste0("DRUG_", 1:10),
                       value = runif(10))
  merged_df <- validate_mapping(trt_df, ref_df, nested_confounders = NULL)
  expect_equal(dim(merged_df), c(13, 3))  
  
  ref_df$Barcode <- c("A", rep("B", 9))
  trt_df$Barcode <- "A"
  merged_df2 <- validate_mapping(trt_df, ref_df, nested_confounders = "Barcode")
  expect_equal(dim(merged_df2), c(11, 4))
})

test_that("validate_mapping catches reverse single-agent data", {
  trt_df <- data.frame(Gnumber = "DRUG_10",
                       Gnumber_2 = "DRUG_2", value = runif(1),
                       Concentration = runif(1),
                       Concentration_2 = runif(1))
  ref_df <- data.frame(Gnumber = c("DRUG_10", "DRUG_2", "DRUG_3"),
                       Gnumber_2 = rep("vehicle", 3), value = runif(3),
                       Concentration = runif(3), Concentration_2 = rep(0, 3))
  merged_df <- validate_mapping(trt_df, ref_df, nested_confounders = NULL)
  expect_equal(dim(merged_df), c(3, 5))  
})
