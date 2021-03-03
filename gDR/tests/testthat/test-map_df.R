library(gDR); library(testthat)
source("setUp.R")

test_that("map_df works as expected", {
  md_df <- unique(md_df)
  ref <- md_df$Gnumber %in% c("vehicle", "untreated")
  ref_df <- md_df[ref, ]
  trt_df <- md_df[!ref, ]
  Keys <- identify_keys(test_df)

  mapping <- map_df(trt_df, 
                    ref_df, 
                    row_endpoint_value_filter = rep(TRUE, nrow(trt_df)), 
                    Keys = Keys, 
                    ref_type = "untrt_Endpoint")

  expect_equal(names(mapping), rownames(trt_df))
  expect_equal(length(mapping), nrow(trt_df))

  out <- lapply(seq(33, 64, 1), function(x) {as.character(c(x %% 16), 16 + (x %% 16))})
  names(out) <- seq(33, 64, 1)

  expect_equal(mapping, out)
  expect_equal(sort(unique(unname(unlist(mapping)))), sort(rownames(ref_df)))
})


test_that("Best match is detected for missing controls", {

})


test_that("NAs are returned for missing controls", {

})
