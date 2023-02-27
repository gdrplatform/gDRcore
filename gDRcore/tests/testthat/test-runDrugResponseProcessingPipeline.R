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
dir.create(p_dir) 
on.exit(unlink(p_dir, TRUE))

# Define path for data stored in gDR package
dataDir <- system.file("extdata", "data1", package = "gDRimport")

# Extract path for example raw_data
manifest <- list.files(dataDir, pattern = "manifest", full.names = TRUE)
template <- list.files(dataDir, pattern = "Template", full.names = TRUE)
raw_data <- list.files(dataDir, pattern = "^RawData", full.names = TRUE)
l_data <- gDRimport::load_data(manifest, template, raw_data)
imported_data <- gDRcore::merge_data(l_data$manifest, l_data$treatments, l_data$data)

### runDrugResponseProcessingPipeline ###
expect_true(length(list.files(p_dir)) == 0)

# mae <- purrr::quietly(gDRtestData::generateNoiseRawData)(
#   cell_lines, drugs, e_inf, ec50, hill_coef
# )

mae_v1 <- purrr::quietly(gDRcore:::runDrugResponseProcessingPipeline)(imported_data, data_dir = p_dir,
                                                                      add_raw_data = TRUE)
expect_true(length(list.files(p_dir)) > 0)
expect_length(mae_v1$warnings, 7)

mae_v2 <-
  purrr::quietly(gDRcore:::runDrugResponseProcessingPipeline)(
    imported_data,
    data_dir = p_dir,
    partial_run = TRUE,
    start_from = "fit_SE"
  )
expect_length(mae_v2$warnings, 3)

mae_v3 <-
  purrr::quietly(gDRcore:::runDrugResponseProcessingPipeline)(
    imported_data,
    data_dir = p_dir,
    partial_run = TRUE,
    start_from = "normalize_SE",
    selected_experiments = c("single-agent")
  )
expect_length(mae_v3$warnings, 4)

mae_v4 <-
  purrr::quietly(gDRcore:::runDrugResponseProcessingPipeline)(
    mae_v1$result
  )
expect_length(mae_v4$warnings, 7)

expect_identical(mae_v1$result, mae_v2$result)
expect_identical(mae_v2$result, mae_v3$result)
expect_identical(mae_v3$result, mae_v4$result)


testthat::expect_error(
  gDRcore:::runDrugResponseProcessingPipeline(imported_data, selected_experiments = "single-agent"),
  "^Selected experiments"
)

testthat::expect_error(
  gDRcore:::runDrugResponseProcessingPipeline(imported_data, partial_run = TRUE),
  "^Path for/to the intermediate data"
)

### prepare_input ###
nc <- intersect(names(imported_data), gDRutils::get_env_identifiers("barcode"))
inl <- prepare_input(imported_data, nc, .get_default_nested_identifiers())
expect_list(inl)
inl_names <- c("df_", "df_list", "nested_confounders", "nested_identifiers_l", "exps")
expect_equal(names(inl), inl_names)

expect_error(prepare_input(imported_data, list(), NULL), "nested_confounders")
expect_error(prepare_input(imported_data, NULL, list()), "nested_identifiers")



})
