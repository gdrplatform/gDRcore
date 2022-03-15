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
  