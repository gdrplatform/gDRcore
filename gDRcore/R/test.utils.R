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
  data.frame(lapply(df, as.character))
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
  lFiles <- lapply(files, function(x) { read.table(x, sep = "\t", header = TRUE)})
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
    write.table(lData[[x]], outFile, sep = "\t", quote = FALSE, row.names = FALSE)
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
    write.table(gDRutils::assay_to_dt(se, x, merge_metrics = TRUE), outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  })

  #df_raw_data from metadata
  outFile <- file.path(outDir, paste0(prefix, "_df_raw_data.tsv"))
  write.table(metadata(se)$df_raw_data, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  

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
  
  xAs <- gDRutils::assay_to_dt(se, "Normalized")
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



#' @export
#'
test_se <- function(se, lRef) {
  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_list(lRef)
  
  x <- standardize_df(metadata(se)[["df_"]])
  y <- standardize_df(data.table::data.table(lRef$df_raw_data))
  class(y) <- class(x)
  expect_equal(colnames(y),colnames(x))
  # TODO: 'Medium' column not present
  # now moved to metadata(se)$experiment_metadata, right?
  expect_equal(x, y)
  expect_equal(x, y[colnames(y) != "Medium"])
  
  # TODO: additional cols:
  ## ref_Endpoint (Concentration_2, DrugName_2, Gnumber_2)
  ## untrt_Endpoint (E2)
  expect_equal(metadata(se)$Keys, yaml::yaml.load(paste0(lRef$ref_keys, "\n")))
  
  # TODO: no 'row_maps' metadata
  # any alternative test (or should this chunk be removed?)
  expect_equal(metadata(se)$row_maps,
               yaml::yaml.load(paste0(lRef$ref_row_maps, "\n")))

  # Test metadata keys. 
  cols <- c("Trt", "ref_Endpoint", "untrt_Endpoint", "Day0") 
  a <- lapply(identify_keys(df_)[cols], sort)
  b <- lapply(identify_keys2(df_)[cols], sort)

  pass <- identical(a, b)

  #assays check
  # 1. Assay names
  ## old assay names
  # grep("assay_", names(lRef), value = TRUE)
  # [1] "assay_Averaged"     "assay_Avg_Controls" "assay_Controls"     "assay_Metrics"      "assay_Normalized"  
  ## new assay names 
  # assayNames(se)
  # [1] "RawTreated" "Controls"   "Normalized" "Averaged"   "Metrics"   
  
  # 2. different naming convention for rId/cId
  # ````
  # mets_ref <- lRef$assay_Metrics
  # mets <- gDRutils::assay_to_dt(se = se,assay_name = "Metrics")
  # unique(mets_ref$rId)[2]
  # "G02001876_G02001876.1-4_168_vehicle_0_0.0001_vehicle_9"
  # unique(mets$rId)[2]
  # "0_G02001876.1-4_vehicle_168_0.0001_G02001876_vehicle_9"
  # ```
  # we can convert new rId/cId to old naming scheme with sth like
  # paste(unlist(strsplit(unique(mets$rId)[2], "_"))[c(6, 2, 4, 7, 1, 5, 3, 8)], collapse = "_")
  # but probably it would be safer to stick to the old naming convention
  # (or at least we should check how/if it would affect gDRviz)
 
  # 3. same data in the assays
  # we should probably check that same data is returned for four assays:
  # Controls, Normalized, Averaged, Metrics
  #myL <- lapply(SummarizedExperiment::assayNames(se), function(x) {
  myL <- lapply(c("Controls", "Normalized", "Averaged", "Metrics"), function(x) {
    print(x)
    xAs <- gDRutils::assay_to_dt(se, x, merge_metrics = TRUE)
    xAs$DrugName <- as.character(xAs$DrugName)
    xDf <- lRef[[paste0("assay_", x)]]
    if(x %in% c("Controls", "Avg_Controls")){
    xDf$DivisionTime <- as.numeric(xDf$DivisionTime)
    }
    expect_true(nrow(xAs) == nrow(xDf))
    expect_equivalent(sort(xAs[, order(names(xAs))]), 
                      sort(data.table::as.data.table(xDf)[, order(names(xDf))]), tolerance = 1e-5)
  })
}

#' @export
#'
check_identity_of_dfs <- function(mat1, mat2, i, j, cols = NULL) {
  a <- mat1[i, j][[1]]
  rownames(a) <- NULL

  b <- mat2[i, j][[1]]
  rownames(b) <- NULL

  if (is.null(cols)) {
    cols <- colnames(a)
  }

  a <- a[order(a[, cols]), ]
  b <- b[order(b[, cols]), ]
  testthat::expect_equal(a[, cols, drop = FALSE], b[, cols, drop = FALSE])
}
