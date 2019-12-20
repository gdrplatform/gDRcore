library(magrittr)

standardize_df <- function(df) {
  df %>% dplyr::mutate_all(as.character)
}

read_ref_data <- function(inDir, prefix = "ref") { 
  files <- list.files(inDir, paste0(prefix, "_.+\\.tsv$"), full.names = TRUE)
  lFiles <- lapply(files, function(x) { readr::read_delim(x, delim = "\t")})
  names(lFiles) <- gsub("\\.tsv", "", gsub(paste0("^", prefix, "_"), "", basename(files)))
  refKeys <- yaml::read_yaml(file.path(inDir, paste0(prefix, "_keys.yaml")))
  refRowMaps <- yaml::read_yaml(file.path(inDir, paste0(prefix, "_row_maps.yaml")))
  lFiles$ref_keys <- refKeys
  lFiles$ref_row_maps <- refRowMaps
  lFiles
}

write_ref_data_df <- function(lData, outDir, prefix = "ref") {
  
  myL <- lapply(1:length(lData), function(x) {
    outFile <- file.path(outDir, paste0(prefix, "_lData_", names(lData)[x], ".tsv"))
    readr::write_delim(lData[[x]], outFile, delim = "\t")
  })
  
}

write_ref_data_se <- function(se, outDir, prefix = "ref") {
  
  #assays 
  myL <- lapply(SummarizedExperiment::assayNames(se), function(x) {
    outFile <- file.path(outDir, paste0(prefix, "_assay_", x, ".tsv"))
    readr::write_delim(gDR::assay_to_df(se, x), outFile, delim = "\t")
  })
  
  #df_raw_data from metadata  
  outFile <- file.path(outDir, paste0(prefix, "_df_raw_data.tsv"))
  readr::write_delim(metadata(se)$df_raw_data, outFile, delim = "\t")
    
  keys_yaml <- yaml::as.yaml(metadata(se)$Keys)
  yaml::write_yaml(keys_yaml, file.path(outDir, paste0(prefix,"_keys.yaml")))
  
  row_maps_yaml <- yaml::as.yaml(metadata(se)$row_maps)
  yaml::write_yaml(row_maps_yaml, file.path(outDir, paste0(prefix,"_row_maps.yaml")))
}

test_lData <- function(lData, lRef) {
  testthat::expect_true(all(vapply(lData, function(x) {
    any(class(x) %in% "data.frame")
  }, logical(1L))))
  testthat::expect_equal(names(lData), c("manifest", "treatments", "data"))
  expect_equal(lData$manifest, lRef$lData_manifest)
  expect_equal(standardize_df(lData$data), standardize_df(data.frame(lRef$lData_data)))
  expect_equal(standardize_df(lData$treatments),
               standardize_df(data.frame(lRef$lData_treatments)))
}

test_se <- function(se, lRef) {
  expect_equal(standardize_df(metadata(se)$df_raw_data),
               standardize_df(data.frame(lRef$df_raw_data)))
  expect_equal(yaml::as.yaml(metadata(se)$Keys), paste0(lRef$ref_keys, "\n"))
  expect_equal(yaml::as.yaml(metadata(se)$row_maps),
               paste0(lRef$ref_row_maps, "\n"))
  
  #assays check
  myL <- lapply(SummarizedExperiment::assayNames(se), function(x) {
    print(x)
    xAs <- gDR::assay_to_df(se, x, merge_metrics = FALSE)
    xDf <- lRef[[paste0("assay_", x)]]
    expect_true(nrow(xAs) == nrow(xDf))
    expect_equal(xAs, data.frame(xDf), tolerance = 1e-5)
  })
}

