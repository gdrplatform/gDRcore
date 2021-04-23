skip("Unit tests will be recreated in GDR-853")
testthat::context("Test load_files.R")

library("gDR")

 # Load toy set and reference files
testDataDir <- system.file(package = "gDR", "testdata", "data6")

manifest_file <-
  list.files(testDataDir, pattern = "Manifest_*", full.names = TRUE) %>%
  grep("anonym", ., invert = TRUE, value = TRUE)
template_file = file.path(testDataDir, c('trtmt1.xlsx',
                                         'untreated.xlsx')) %>%
  grep("anonym", ., invert = TRUE, value = TRUE)
results_file <-
  list.files(testDataDir, pattern = "12\\w*", full.names = TRUE) %>%
  grep("anonym", ., invert = TRUE, value = TRUE)

lRef <- read_ref_data(testDataDir)

testthat::test_that("load_data works as expected", {
  
  result <- load_data(manifest_file = manifest_file,
            df_template_files = template_file,
            results_file = results_file,
            instrument = "EnVision")
  
  expect_equal(
    lapply(lRef$lData_data, as.character),
    lapply(result$data, as.character)
  )
  expect_equal(
    lapply(lRef$lData_manifest, as.character),
    lapply(result$manifest, as.character)
  )
  expect_equal(
    lapply(lRef$lData_treatments, as.character),
    lapply(result$treatments, as.character)
  )
})

testthat::test_that("load_data throwing expected errors", {
  # Test assertion:
  expect_error(
    load_data(manifest_file = 1,
              df_template_files = template_file,
              results_file = results_file,
              instrument = "EnVision"),
    "'manifest_file' must be a character vector"
  )
  expect_error(
    load_data(manifest_file = "fake/path",
              df_template_files = template_file,
              results_file = results_file,
              instrument = "EnVision"),
    "'manifest_file' must be a readable path"
  )
  expect_error(
    load_data(manifest_file = manifest_file,
              df_template_files = 1,
              results_file = results_file,
              instrument = "EnVision"),
    "'df_template_files' must be a character vector or data.frame"
  )
  expect_error(
    load_data(manifest_file = manifest_file,
              df_template_files = c("fake/path/1", "fake/path/2"),
              results_file = results_file,
              instrument = "EnVision"),
    "Following path(s) with no read permission found: 'fake/path/1, fake/path/2'",
    fixed = TRUE
  )
  expect_error(
  load_data(manifest_file = manifest_file,
            df_template_files = template_file,
            results_file = 1,
            instrument = "EnVision"),
  "'results_file' must be a character vector or data.frame"
  )
  expect_error(
    load_data(manifest_file = manifest_file,
              df_template_files = template_file,
              results_file = results_file,
              instrument = 1),
    "'instrument' must be a character vector"
  )
  
  # Check 'manifest' incompatible with 'template'

  ## Load template_file from different toy set
  template_file <- 
    list.files(system.file(package = "gDR", "testdata", "data2"), pattern = "Template_*", full.names = TRUE)
  
  expect_error(
    load_data(manifest_file = manifest_file,
              df_template_files = template_file,
              results_file = results_file,
              instrument = "EnVision"),
    "Template does not contains all expected headers for a 'template'. 'Gnumber' is/are required. Please correct your template."
  )
  
})

