test_that(".assign_conditions works", {
  a <- 5
  b <- 2
  drugnames <- c(paste0("G", seq_len(a)), rep(gDRutils::get_identifier("untreated_tag"), 2))
  df <- data.frame(DrugName = drugnames, Concentration = c(rep(1, a), rep(0, b * 2)))
  df <- .assign_conditions(df, nested_identifiers = c("Concentration"))
  expect_true(all(dim(df) == c(9, 3)))
  expect_equal(df$treated_untreated, as.factor(c(rep("treated", 5), rep("untreated", 4))))

  # Multiple drugs.
  df_ <- data.frame(DrugName = c("vehicle", "vehicle", "bogus"), DrugName2 = c("vehicle", "bogus", "bogus"),
    Concentration = c(0, 0, 1), Concentration_2 = c(0, 1, 1))
  df_ <- .assign_conditions(df_, nested_identifiers = c("Concentration", "Concentration_2"))
  expect_true(all(dim(df_) == c(3, 5)))
  expect_equal(df_$treated_untreated, as.factor(c("untreated", "reference", "treated")))
})


test_that(".assign_conditions errors as expected", {
  df <- data.frame(Concentration = rep(1, 5), Concentration_2 = rep(0.3333, 5)) 
  expect_error(.assign_conditions(df, "Concentration"), regexp = "no untreated conditions")
  expect_error(.assign_conditions(df, c("Concentration", "Concentration_2")), regexp = "no untreated conditions")
})
