test_that("map_df works as expected", {
  md_df <- unique(md_df)
  md_df$rn <- as.character(seq_len(nrow(md_df)))
  ref <- md_df$Gnumber %in% c("vehicle", "untreated")
  
  ref_df <- md_df[ref, ]
  trt_df <- md_df[!ref, ]
  Keys <- identify_keys(test_df)

  ref_type <- "untrt_Endpoint"
  mapping <- map_df(trt_df, 
                    ref_df, 
                    ref_cols = Keys[[ref_type]],
                    ref_type = ref_type)

  expect_equal(names(mapping), trt_df$rn)
  expect_equal(length(mapping), nrow(trt_df))

  out <- lapply(seq(33, 64, 1), function(x) {
    as.character(c((x - 1) %% 16 + 1, 17 + ((x - 1) %% 16)))
    })
  names(out) <- seq(33, 64, 1)

  expect_equal(mapping, out)
  expect_equal(sort(unique(unname(unlist(mapping)))), sort(rownames(ref_df)))
})



# TODO: test_that("Best match is detected for missing controls", {}) # nolint


# TODO: test_that("NAs are returned for missing controls", {}) # nolint


mat_elem <- data.table::data.table(DrugName = rep(c("untreated", "drugA", "drugB", "untreated"), 2),
                                   DrugName_2 = rep(c("untreated", "vehicle", "drugA", "drugB"), 2),
                                   clid = rep(c("C1", "C2"), each = 4))

test_that(".map_references works as expected", {
  # Combination data.
  obs <- .map_references(mat_elem, rowData_colnames = c("DrugName", "DrugName_2"))
  exp <- list("3" = c("2", "4"), "7" = c("6", "8"))
  expect_equal(obs, exp)

  # Single-agent data.
  colname <- c("DrugName", "clid")
  mat_elem2 <- mat_elem[, colname, with = FALSE]
  obs2 <- .map_references(mat_elem2, rowData_colnames = c("DrugName", "DrugName_2"))
  exp2 <- list("2" = NULL, "3" = NULL, "6" = NULL, "7" = NULL)
  expect_equal(obs2, exp2)
})

test_that("map_untreated works as expected", {
  # Combination data.
  obs <- map_untreated(mat_elem)
  expect_equal(sum(obs), 2)
})
