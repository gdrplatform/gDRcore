test_that("grr_matches works as expected", {
  # Combination data.
  mat_elem <- data.table::data.table(DrugName = rep(c("untreated", "drugA", "drugB", "untreated"), 2),
                                     DrugName_2 = rep(c("untreated", "vehicle", "drugA", "drugB"), 2),
                                     clid = rep(c("C1", "C2"), each = 4))
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  ref_idx <- which(mat_elem$DrugName %in% untreated_tag |  mat_elem$DrugName_2 %in% untreated_tag)
  ref <- mat_elem[ref_idx, ]
  treated <- mat_elem[-ref_idx, ]
  valid <- c("DrugName", "DrugName_2")
  clid <- "clid"
  # split data.tables to simple model with clid column and drug column
  trt <- gDRutils::loop(valid, function(x) treated[, c(clid, x), with = FALSE])
  trt <- do.call(paste, do.call(rbind, gDRutils::loop(trt, function(x) setNames(x, names(trt[[1]])))))
  
  ref <- gDRutils::loop(valid, function(x) ref[, c(clid, x), with = FALSE])
  ref <- do.call(paste, do.call(rbind, gDRutils::loop(ref, function(x) setNames(x, names(ref[[1]])))))
  
  matchTrtRef <- grr_matches(trt, ref, list = FALSE, all.y = FALSE)
  expect_equal(dim(matchTrtRef), c(4, 2))
  expect_s3_class(matchTrtRef, "data.table")
  
  matchTrtRefList <- grr_matches(trt, ref, list = TRUE)
  expect_equal(length(matchTrtRefList), 4)
})

test_that("cleanup_metadata works as expected", {
  df <- data.table::data.table(clid = "CELL_LINE",
                               Gnumber = "DRUG_1",
                               Concentration = 3,
                               Duration = 72)
  cleanup_df <- purrr::quietly(cleanup_metadata)(df)
  expect_equal(dim(cleanup_df$result), c(1, 11))
})

test_that("data_model.character works as expected", {
  expect_equal(data_model("co-dilution"), "single-agent")
  expect_equal(data_model("combination"), "combination")
  
  expect_error(data_model("single-a"))
  expect_error(data_model(""))
  expect_error(data_model(NULL))
})

test_that("validate_data_models_availability works as expected", {
  d_type <- c("co-dilution", "combination")

  expect_equal(validate_data_models_availability(d_types = d_type, 
                                                 s_d_models = gDRutils::get_supported_experiments()),
               NULL)
  expect_equal(validate_data_models_availability(d_types = "single-a", 
                                                 s_d_models = gDRutils::get_supported_experiments()),
               NULL)
  expect_error(validate_data_models_availability(d_types = d_type, 
                                                 s_d_models = "single_agent"),
               "'nested_identifiers_l' lacks information for the ")
  expect_error(validate_data_models_availability(d_types = d_type, 
                                                 s_d_models = NULL),
               "'nested_identifiers_l' lacks information for the ")
  
})

