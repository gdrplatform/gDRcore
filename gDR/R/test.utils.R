library(magrittr)

standardize_df <- function(df) {
  df %>% dplyr::mutate_all(as.character)
}

#' @export
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
    readr::write_delim(gDR::assay_to_df(se, x, merge_metrics = TRUE), outFile, delim = "\t")
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

test_se_normalized <- function(se, lRef) {
  expect_equal(standardize_df(metadata(se)$df_raw_data),
               standardize_df(data.frame(lRef$df_raw_data)))
  expect_equal(yaml::as.yaml(metadata(se)$Keys), paste0(lRef$ref_keys, "\n"))
  expect_equal(yaml::as.yaml(metadata(se)$row_maps),
               paste0(lRef$ref_row_maps, "\n"))
  
    x = "Normalized"
    xAs <- gDR::assay_to_df(se, x, merge_metrics = TRUE)
    xAs$DrugName <- as.character(xAs$DrugName)
    xDf <- lRef[[paste0("assay_", x)]]
    if(x %in% c("Controls", "Avg_Controls")){
      xDf$DivisionTime <- as.numeric(xDf$DivisionTime)
    }
    expect_true(nrow(xAs) == nrow(xDf))
    expect_equal(xAs, data.frame(xDf), tolerance = 1e-5)
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
    xAs <- gDR::assay_to_df(se, x, merge_metrics = TRUE)
    xAs$DrugName <- as.character(xAs$DrugName)
    xDf <- lRef[[paste0("assay_", x)]]
    if(x %in% c("Controls", "Avg_Controls")){
    xDf$DivisionTime <- as.numeric(xDf$DivisionTime)
    }
    expect_true(nrow(xAs) == nrow(xDf))
    expect_equal(xAs, data.frame(xDf), tolerance = 1e-5)
  })
}

test_synthetic_normalization <- function(se, refNormalizedTsv){
  xAs <- gDR::assay_to_df(se, "Normalized", merge_metrics = TRUE)[, c("GRvalue", "RelativeViability")]
  xRef <- refNormalizedTsv[, c("GRvalue", "RelativeViability")]
  expect_true(nrow(xAs) == nrow(xRef))
  expect_equal(xAs, data.frame(xRef), tolerance = 0)
}

#' @export
save_file_type_info <-
  function(v,
           save_dir,
           normKeysFileName = "normKeys.json",
           dfRawDataFileName = "dfRawData.tsv",
           fileTypeInfoName = "fileTypeInfo.csv") {
    checkmate::assert_true(all(c("Manifest", "Template", "RawData") %in% names(v)))
    
    tbl <- tibble::tibble(
      data_type = c(
        rep("manifest", length(v$Manifest$name)),
        rep("template", length(v$Template$name)),
        rep("rawData", length(v$RawData$name)),
        "normKeys",
        "dfRawData"
      ),
      name = c(
        v$Manifest$name,
        v$Template$name,
        v$RawData$name,
        normKeysFileName,
        dfRawDataFileName
      )
    )
    outFile <- file.path(save_dir, fileTypeInfoName)
    write.csv(tbl, outFile, row.names = FALSE)
  }
