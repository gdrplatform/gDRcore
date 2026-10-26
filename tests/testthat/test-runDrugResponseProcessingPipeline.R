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
  imported_data <- purrr::quietly(merge_data)(l_data$result$manifest,
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
  expect_lte(length(mae_v1$warnings), 3)
  expect_true(any(grepl("Switching into 'Barcode' nested confounder.", mae_v1$warnings)))
  expect_true(any(grepl("method L-BFGS-B uses 'factr'.", mae_v1$warnings)))

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
  expect_lte(length(mae_v3$warnings), 3)
  expect_true(any(grepl("Switching into 'Barcode' nested confounder.", mae_v3$warnings)))
  expect_true(any(grepl("method L-BFGS-B uses 'factr'.", mae_v3$warnings)))

  mae_v4 <-
    purrr::quietly(runDrugResponseProcessingPipeline)(
      mae_v1$result
    )
  expect_true(any(grepl("Switching into 'Barcode' nested confounder.", mae_v3$warnings)))
  expect_true(any(grepl("method L-BFGS-B uses 'factr'.", mae_v3$warnings)))

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

test_that("runDrugResponseProcessingPipeline processes time-course data", {
  tc_name <- gDRutils::get_supported_experiments("time-course")
  concs <- c(0.01, 0.1, 1, 10)
  layout <- rbind(
    data.table::data.table(Gnumber = "vehicle", Concentration = 0, replicate = 1:4),
    data.table::CJ(Gnumber = c("DrugA", "DrugB"), Concentration = concs, replicate = 1:2)
  )
  layout[, `:=`(WellRow = LETTERS[seq_len(.N)], WellColumn = "1", DrugName = Gnumber,
                drug_moa = "unknown", clid = "CL1", CellLineName = "CL1", Barcode = "PLATE1")]
  input_dt <- layout[, .(Duration = seq(0, 72, 12)), by = names(layout)]
  potency <- c(vehicle = Inf, DrugA = 0.1, DrugB = 1)
  input_dt[, ReadoutValue := 1000 * 2 ^ (Duration / 24 *
                                           (1 - Concentration / (Concentration + potency[Gnumber])))]
  input_dt[, replicate := NULL]

  expect_message(
    mae <- runDrugResponseProcessingPipeline(data.table::copy(input_dt),
                                             nested_confounders = "Barcode"),
    "Time-course data stops at the LogFoldChange assay"
  )
  expect_equal(names(MultiAssayExperiment::experiments(mae)), tc_name)
  expect_equal(SummarizedExperiment::assayNames(mae[[tc_name]]),
               c("RawTreated", "Controls", "LogFoldChange"))

  # the same SE as the time-course report builds step by step
  .drop_record_id <- function(dt) dt[, setdiff(names(dt), "record_id"), with = FALSE]
  se_ref <- normalize_SE(create_SE(data.table::copy(input_dt), data_type = tc_name),
                         data_type = tc_name)
  for (assay_name in SummarizedExperiment::assayNames(se_ref)) {
    expect_equal(
      .drop_record_id(gDRutils::convert_se_assay_to_dt(mae[[tc_name]], assay_name)),
      gDRutils::convert_se_assay_to_dt(se_ref, assay_name),
      ignore_attr = TRUE
    )
  }
})

test_that("runDrugResponseProcessingPipeline processes co-dilution data", {
  cd_name <- gDRutils::get_supported_experiments("cd")
  set.seed(2)
  cell_lines <- gDRtestData::create_synthetic_cell_lines()
  drugs <- gDRtestData::create_synthetic_drugs()

  mae <- suppressWarnings(generateCodilutionSmall(cell_lines, drugs, save = FALSE))
  expect_equal(names(MultiAssayExperiment::experiments(mae)), cd_name)
  expect_true("Metrics" %in% SummarizedExperiment::assayNames(mae[[cd_name]]))
  metrics <- gDRutils::convert_se_assay_to_dt(mae[[cd_name]], "Metrics")
  expect_gt(NROW(metrics), 0L)
})
