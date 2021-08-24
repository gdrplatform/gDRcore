test_that(".assign_treated_and_untreated_conditions works", {
  drugnames <- c(paste0("G", seq_len(5)), rep(gDRutils::get_identifier("untreated_tag"), 2))
  df <- data.frame(DrugName = drugnames) 
  df <- .assign_treated_and_untreated_conditions(df)
  expect_true(all(dim(df) == c(9, 2)))
  expect_equal(df$treated_untreated, as.factor(c(rep("treated", 5), rep("untreated", 4))))

  # Multiple drugs.
  df_ <- data.frame(DrugName = c("vehicle", "vehicle", "bogus"), DrugName2 = c("vehicle", "bogus", "bogus"))
  df_ <- .assign_treated_and_untreated_conditions(df_, drugname_col = c("DrugName", "DrugName2"))
  expect_true(all(dim(df_) == c(3, 3)))
  expect_equal(df_$treated_untreated, as.factor("untreated", "treated", "treated"))
})


test_that(".assign_treated_and_untreated_conditions errors as expected", {
  drugnames <- c(paste0("G", seq_len(5)), rep("", 2))
  df <- data.frame(DrugName = drugnames) 
  expect_error(.assign_treated_and_untreated_conditions(df), regexp = "no untreated conditions matching pattern")
})
