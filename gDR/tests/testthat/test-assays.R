library(testthat); library(gDR)

context("")

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
