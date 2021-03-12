library(testthat); library(gDR)

test_that(".create_mapping_factors works as expected", {
  rowdata <- data.frame("Billy" = LETTERS[1:3], "Joel" = paste0("A", LETTERS[1:3]))
  coldata <- data.frame("Air" = letters[9:7], "Supply" = paste0("a", letters[9:7]))
  mf <- gDR:::.create_mapping_factors(rowdata = rowdata, coldata = coldata)
  new_len <- nrow(rowdata) * nrow(coldata)
  expect_equal(nrow(mf), new_len)
  expect_equal(rownames(mf), as.character(seq(new_len)))
  expect_equal(ncol(mf), ncol(rowdata) + ncol(coldata) + 2)
})
