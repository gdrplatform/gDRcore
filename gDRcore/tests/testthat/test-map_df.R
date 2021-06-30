library(gDRcore)
library(testthat)


test_that("map_df works as expected", {
  md_df <- unique(md_df)
  ref <- md_df$Gnumber %in% c("vehicle", "untreated")
  ref_df <- md_df[ref, ]
  trt_df <- md_df[!ref, ]
  Keys <- identify_keys(test_df)

  ref_type <- "untrt_Endpoint"
  mapping <- map_df(trt_df, 
                    ref_df, 
                    ref_cols = Keys[[ref_type]],
                    ref_type = ref_type)

  expect_equal(names(mapping), rownames(trt_df))
  expect_equal(length(mapping), nrow(trt_df))

  out <- lapply(seq(33, 64, 1), function(x) {
    as.character(c((x - 1) %% 16 + 1, 17 + ((x - 1) %% 16)))})
  names(out) <- seq(33, 64, 1)

  expect_equal(mapping, out)
  expect_equal(sort(unique(unname(unlist(mapping)))), sort(rownames(ref_df)))
})


# TODO: 
# Best match is detected for missing controls


# TODO:
# NAs are returned for missing controls
