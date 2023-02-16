library(testthat)

test_that("merge_data works as expected", {

  manifestFile <- system.file(package = "gDRimport", "extdata", "data1", "manifest.xlsx")
  templateFiles <- list.files(system.file(package = "gDRimport", "extdata", "data1"), pattern = "^Template", full.names = TRUE)
  rawDataFiles <- list.files(system.file(package = "gDRimport", "extdata", "data1"), pattern = "^RawData", full.names = TRUE)

  manifest <- gDRimport::load_manifest(manifestFile)
  template <- gDRimport::load_templates(templateFiles)
  rawData <- gDRimport::load_results(rawDataFiles, manifest$headers, instrument = "EnVision")

  merged <- gDRcore::merge_data(manifest$data, template, rawData)

  expect_equal(
    merged$ReadoutValue[merged$WellRow == "A" & merged$WellColumn == 3],
    rawData$ReadoutValue[rawData$WellRow == "A" & rawData$WellColumn == 3]
  )
})
