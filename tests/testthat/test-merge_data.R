test_that("merge_data works as expected", {
  
  o_cols <- c("WellRow", "WellColumn", "Barcode")
  
  manifestFile <- system.file(package = "gDRimport", "extdata", "data1", "manifest.xlsx")
  templateFiles <- list.files(system.file(package = "gDRimport", "extdata", "data1"),
                              pattern = "^Template", full.names = TRUE)
  rawDataFiles <- list.files(system.file(package = "gDRimport", "extdata", "data1"),
                             pattern = "^RawData", full.names = TRUE)

  manifest <- gDRimport::load_manifest(manifestFile)
  template <- gDRimport::load_templates(templateFiles)
  rawData <- gDRimport::load_results(rawDataFiles, manifest$headers, instrument = "EnVision")
  rawData <- data.table::setorderv(data.table::setDF(rawData), o_cols)

  manifest$data <- data.table::setDT(manifest$data)
  template <- data.table::setDT(template)
  rawData <- data.table::setDT(rawData)
  
  merged_quietly <- purrr::quietly(merge_data)(manifest$data, template, rawData)
  merged <- data.table::setorderv(data.table::setDF(merged_quietly$result), o_cols)

  # checking readout value from rawData
  expect_equal(
    merged$ReadoutValue[merged$WellRow == "A" & merged$WellColumn == 3],
    rawData$ReadoutValue[rawData$WellRow == "A" & rawData$WellColumn == 3]
  )

  # checking templates are loaded
  expect_equal(
    unique(merged$Template),
    unique(manifest$data$Template)
  )

  # checking column names from manifest and rawdata are present
  expect_true(all(names(merged) %in% c(unlist(unname(manifest$headers)), names(rawData))))

  # testing wrong input
  expect_error(
    merge_data("test", template, rawData),
    "Assertion on 'manifest' failed: Must be a data.table",
    fixed = TRUE
  )
})
