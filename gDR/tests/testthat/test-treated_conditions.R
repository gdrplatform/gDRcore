library(testthat); library(gDR)

context("Test that the treated and untreated functions can be appropriately fetched")

test_that(".get_untreated_conditions works", {
  old_drug_name <- gDRutils::get_identifier("drugname")
  on.exit(gDRutils::set_identifier("drugname", old_drug_name))

  gDRutils::set_identifier("drugname", "testdrugname")
  untrt_tags <- gDRutils::get_identifier("untreated_tag")
  n_untrt <- length(untrt_tags)

  data <- data.frame(testdrugname = c(untrt_tags, "treated"), 
                     drug_id = LETTERS[seq_len(n_untrt + 1)])
  data$name_ <- paste0(data$testdrugname, "_", data$drug_id)
  untrt <- gDR:::.get_untreated_conditions(data)
  trt <- gDR:::.get_treated_conditions(data)

  expect_equal(untrt, data$name_[seq_len(n_untrt)])
  expect_equal(trt, data$name_[setdiff(seq_len(length(data)), seq_len(n_untrt))])
})


test_that(".assign_treated_and_untreated_conditions works", {
  drugnames <- c(paste0("G", seq_len(5)), rep(gDRutils::get_identifier("untreated_tag"), 2))
  df <- data.frame(DrugName = drugnames) 
  df <- .assign_treated_and_untreated_conditions(df)
  expect_true(all(dim(df) == c(9, 2)))
  expect_equal(df$treated_untreated, as.factor(c(rep("treated", 5), rep("untreated", 4))))
})
