library(magrittr)

#' standardize_df
#'
#' Transform all the columns in a dataframe to character type
#'
#' @param df a dataframe for standardization
#'
#' @return a standardized dataframe
#' @export
#'
standardize_df <- function(df) {
  # Assertions:
  stopifnot(inherits(df, "data.frame"))
  df %>% dplyr::mutate_all(as.character)
}

#' read_ref_data
#' 
#' Read reference data 
#'
#' @param inDir a directory path of reference data
#' @param prefix a prefix of reference file names ('ref' by default)
#'
#' @return a list of reference data
#' @export
#'
read_ref_data <- function(inDir, prefix = "ref") {
  # Assertions:
  checkmate::assert_string(inDir)
  checkmate::assert_string(prefix)


  files <- list.files(inDir, paste0(prefix, "_.+\\.tsv$"), full.names = TRUE)
  lFiles <- lapply(files, function(x) { readr::read_delim(x, delim = "\t")})
  names(lFiles) <- gsub("\\.tsv", "", gsub(paste0("^", prefix, "_"), "", basename(files)))
  refKeys <- yaml::read_yaml(file.path(inDir, paste0(prefix, "_keys.yaml")))
  refRowMaps <- yaml::read_yaml(file.path(inDir, paste0(prefix, "_row_maps.yaml")))
  lFiles$ref_keys <- refKeys
  lFiles$ref_row_maps <- refRowMaps
  lFiles
}

#' write_ref_data_df
#' 
#' Write reference dataframe
#'
#' @param lData a list with dataset
#' @param outDir an output directory
#' @param prefix a prefix of reference file names ('ref' by default)
#'
#' @return 
#' @export
#'
write_ref_data_df <- function(lData, outDir, prefix = "ref") {
  # Assertions:
  checkmate::assert_list(lData)
  checkmate::assert_string(outDir)
  checkmate::assert_string(prefix)
  
  myL <- lapply(1:length(lData), function(x) {
    outFile <- file.path(outDir, paste0(prefix, "_lData_", names(lData)[x], ".tsv"))
    readr::write_delim(lData[[x]], outFile, delim = "\t")
  })

}

#' write_ref_data_se
#'
#' @param se a SummarizedExperiment with DR data
#' @param outDir an output directory
#' @param prefix a prefix of reference file names ('ref' by default)
#'
#' @return
#' @export
#'
write_ref_data_se <- function(se, outDir, prefix = "ref") {
  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(outDir)
  checkmate::assert_string(prefix)
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

#' test_lData
#'
#' Test raw data against reference data
#'
#' @param lData a list with dataset
#' @param lRef a list with reference data
#'
#' @return
#' @export
#'
test_lData <- function(lData, lRef) {
  # Assertions:
  testthat::expect_true(all(vapply(lData, function(x) {
    any(class(x) %in% "data.frame")
  }, logical(1L))))
  checkmate::assert_list(lData)
  checkmate::assert_list(lRef)
  
  class(lData$manifest) <- class(lRef$lData_manifest)
  
  x <- standardize_df(lData$data)
  y <- standardize_df(data.frame(lRef$lData_data))
  class(y) <- class(x)
  
  testthat::expect_equal(names(lData), c("manifest", "treatments", "data"))
  expect_equal(lData$manifest, lRef$lData_manifest)
  expect_equal(x, y)
  expect_equal(standardize_df(lData$treatments),
               standardize_df(data.frame(lRef$lData_treatments)))
}

#' test_se_normalized
#'
#' Compare calculated normalization with the reference
#'
#' @param se a SummarizedExperiment object with the Normalized assay
#' @param lRef a list of reference datasets
#'
#' @return
#' @export
#'

test_se_normalized <- function(se, lRef) {
  
  expect_equal(standardize_df(metadata(se)$df_raw_data),
               standardize_df(data.frame(lRef$df_raw_data)))
  expect_equal(yaml::as.yaml(metadata(se)$Keys), paste0(lRef$ref_keys, "\n"))
  expect_equal(yaml::as.yaml(metadata(se)$row_maps),
               paste0(lRef$ref_row_maps, "\n"))
  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_list(lRef)

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

#' test_se
#'
#' Compare all the assays in the SummarizedExperiment object with the reference
#'
#' @param se a SummarizedExperiment object with the calculated assay
#' @param lRef a list of reference datasets
#'
#' @return
#' @export
#'

test_se <- function(se, lRef) {
  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_list(lRef)
  
  x <- standardize_df(metadata(se)$df_raw_data)
  y <- standardize_df(tibble::tibble(lRef$df_raw_data))
  class(y) <- class(x)
  
  expect_equal(x, y)
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

#' test_synthetic_normalization
#'
#' Compare calculated normalization with the reference in synthetic data
#'
#' @param se a SummarizedExperiment object with the calculated assay
#' @param refNormalizedTsv a tsv file with the reference normalized synthetic data
#'
#' @return
#' @export
#'

test_synthetic_normalization <- function(se, refNormalizedTsv) {
  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  stopifnot(inherits(refNormalizedTsv, "data.frame"))
  
  xAs <- gDR::assay_to_df(se, "Normalized")
  xRef <- refNormalizedTsv[refNormalizedTsv$Concentration>0 , ]
  expect_true(nrow(xAs) == nrow(xRef))

  xAs[,grepl('Conc', colnames(xAs))] = round(xAs[,grepl('Conc', colnames(xAs))],5)
  xRef[,grepl('Conc', colnames(xRef))] = round(xRef[,grepl('Conc', colnames(xRef))],5)

  metadata_cols = intersect(colnames(xRef),
    c('clid', 'Barcode',  'Duration', 'WellRow', 'WellColumn',
      colnames(xRef)[grepl('Gnum',colnames(xRef))], colnames(xRef)[grepl('Conc', colnames(xRef))]))
  data_cols = intersect(colnames(xRef), c("GRvalue", "RelativeViability"))

  xMerge = merge(xAs[, c(metadata_cols, data_cols)],
        as.data.frame(xRef[, c(metadata_cols, data_cols)]),
        by = metadata_cols, all = T)

  expect_equal(xMerge[,paste0(data_cols, '.x')], xMerge[,paste0(data_cols, '.y')], tolerance = 1e-4)
}


#' save_file_type_info
#'
#' Save information about file types
#'
#' @param v a list with filenames
#' @param save_dir an output directory
#' @param normKeysFileName a filename with normalization keys ('normKeys.json' by default)
#' @param dfRawDataFileName a filename with rawdata ('dfRawData.tsv' by default)
#' @param fileTypeInfoName a filename with file type information ('fileTypeInfo.csv' by default)
#'
#' @return
#' @export
#'

save_file_type_info <-
  function(v,
           save_dir,
           normKeysFileName = "normKeys.json",
           dfRawDataFileName = "dfRawData.tsv",
           fileTypeInfoName = "fileTypeInfo.csv") {
    # Assertions:
    checkmate::assert_list(v)
    checkmate::assert_string(normKeysFileName)
    checkmate::assert_string(dfRawDataFileName)
    checkmate::assert_string(fileTypeInfoName)
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
