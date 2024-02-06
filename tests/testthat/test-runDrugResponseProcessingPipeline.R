test_that("paste_warnings works as expected", {
  fun1 <- function(x) {
    warning("warning 1")
    warning("warning 2")
    10
  }
  funOutput <- purrr::quietly(fun1)()
  expect_warning(paste_warnings(funOutput$warnings), regexp = "warning 1\nwarning 2")
})


test_that("main pipeline functions works as expected", {
  p_dir <- file.path(tempdir(), "pcheck")
  suppressWarnings(dir.create(p_dir))
  on.exit(unlink(p_dir, TRUE))

  # Define path for data stored in gDR package
  dataDir <- system.file("extdata", "data1", package = "gDRimport")

  # Extract path for example raw_data
  manifest <- list.files(dataDir, pattern = "manifest", full.names = TRUE)
  template <- list.files(dataDir, pattern = "Template", full.names = TRUE)
  raw_data <- list.files(dataDir, pattern = "^RawData", full.names = TRUE)
  l_data <- purrr::quietly(gDRimport::load_data)(manifest, template, raw_data)
  imported_data <-  purrr::quietly(merge_data)(l_data$result$manifest,
                                               l_data$result$treatments,
                                               l_data$result$data)
  
  input_data <- imported_data$result

  input_data <- input_data[input_data$CellLineName %in% unique(input_data$CellLineName)[1] &
                           input_data$DrugName_2 %in% gDRutils::get_env_identifiers("untreated_tag"), ]
  
  ### runDrugResponseProcessingPipeline ###
  expect_true(length(list.files(p_dir)) == 0)

  mae_v1 <- purrr::quietly(runDrugResponseProcessingPipeline)(
    input_data,
    data_dir = p_dir
  )
  expect_true(length(list.files(p_dir)) > 0)
  expect_length(mae_v1$warnings, 2)

  mae_v2 <-
    purrr::quietly(runDrugResponseProcessingPipeline)(
      input_data,
      data_dir = p_dir,
      partial_run = TRUE,
      start_from = "fit_SE"
    )
  expect_length(mae_v2$warnings, 2)

  mae_v3 <-
    purrr::quietly(runDrugResponseProcessingPipeline)(
      input_data,
      data_dir = p_dir,
      partial_run = TRUE,
      start_from = "normalize_SE",
      selected_experiments = c("single-agent")
    )
  expect_length(mae_v3$warnings, 2)

  mae_v4 <-
    purrr::quietly(runDrugResponseProcessingPipeline)(
      mae_v1$result
    )
  expect_length(mae_v4$warnings, 3)

  expect_identical(mae_v1$result, mae_v2$result)
  expect_identical(mae_v2$result, mae_v3$result)
  
  mae_v3$result <- gDRutils::MAEpply(mae_v3$result, function(x) {
    SummarizedExperiment::assay(x, "Controls") <- NULL
    SummarizedExperiment::assay(x, "RawTreated") <- NULL
    x
  })
  
  mae_v4$result <- gDRutils::MAEpply(mae_v4$result, function(x) {
    SummarizedExperiment::assay(x, "Controls") <- NULL
    SummarizedExperiment::assay(x, "RawTreated") <- NULL
    x
  })
  
  # Clear internal metadata (sessionInfo) to not break the tests
  mae_v3$result$`single-agent`@metadata$.internal <- NULL
  mae_v4$result$`single-agent`@metadata$.internal <- NULL
  
  expect_identical(mae_v3$result, mae_v4$result)

  testthat::expect_error(
    runDrugResponseProcessingPipeline(input_data, selected_experiments = "single-agent"),
    "^Selected experiments"
  )

  testthat::expect_error(
    runDrugResponseProcessingPipeline(input_data, partial_run = TRUE),
    "^Path for/to the intermediate data"
  )

  ### prepare_input ###
  nc <- intersect(names(input_data), gDRutils::get_env_identifiers("barcode"))
  inl <- prepare_input(input_data, nc, .get_default_nested_identifiers())
  expect_list(inl)
  inl_names <- c("df_", "df_list", "nested_confounders", "nested_identifiers_l", "exps")
  expect_equal(names(inl), inl_names)

  expect_error(prepare_input(input_data, list(), NULL), "nested_confounders")
  expect_error(prepare_input(input_data, NULL, list()), "nested_identifiers")
})
