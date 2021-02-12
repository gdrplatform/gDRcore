#' createSE
#' TODO: delete me for create_SE function. 
#' @export
createSE <- function(...) {
  .Deprecated("create_SE", package="gDR")
  create_SE(...)
}


#' create_SE
#'
#' Create SummarizedExperiment object from data.frame(s) with dose-reponse data
#'
#' @param df_data a dataframe with DR data
#' @param data_type string type of data to be returned: with untreated conditions only ('untreated'),
#' with treated conditions only ('treated') or all
#' @param readout a character that indiciated readout value in a df ('ReadoutValue' by default)
#' @param discard_keys a vector of keys that should be discarded
#' @param assay_type string assay type for the returned SummarizedExperiment
#'
#' @return SummarizedExperiment object with dose-reponse data
#'
#' @export
create_SE <-
  function(df_data,
           data_type = c("untreated", "treated", "all"),
           readout = 'ReadoutValue',
           discard_keys = NULL,
           assay_type = c("matrix", "BumpyMatrix")) {
    # Assertions:
    stopifnot(any(inherits(df_data, "data.frame"), inherits(df_data, "DataFrame")))
    checkmate::assert_character(data_type)
    checkmate::assert_string(readout)
    checkmate::assert_character(discard_keys, null.ok = TRUE)
    data_type <- match.arg(data_type)
    assay_type <- match.arg(assay_type)

    # stopifnot(all(names(dfList) %in% .assayNames))
    # #dfList must contain first assay (i.e. df_data)
    # stopifnot(.assayNames[1] %in% names(dfList))

    
    mats <- if (assay_type == "matrix") {
      gDRutils::df_to_assay(df_data, data_type = data_type, discard_keys = discard_keys)
    } else if (assay_type == "BumpyMatrix") {
      gDRutils::df_to_bm_assay(df_data, data_type = data_type, discard_keys = discard_keys)
    } else {
      stop(sprintf("bad assay type: '%s'", assay_type))
    }
    
    allMetadata <- getMetaData(df_data, discard_keys = discard_keys)

    seColData <- allMetadata$colData
    rownames(seColData) <- seColData$name_
    seRowData <- allMetadata$rowData
    rownames(seRowData) <- seRowData$name_

    seColData <-
      seColData[colnames(mats), setdiff(colnames(seColData), c('col_id', 'name_'))]
    seRowData <- seRowData[rownames(mats),
                           setdiff(colnames(seRowData), c('row_id', 'name_'))]
    matsL <- list(mats)
    names(matsL) <- readout
    
    se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
                                                     colData = seColData,
                                                     rowData = seRowData)

  }


#' Create a SummarizedExperiment object
#'
#' Create a SummarizedExperiment object from a data.frame where all treatments are on rows,
#' conditions are on columns, treated readouts live in an assay called \code{"treated"},
#' and reference readouts live in an assay called \code{"Controls"}.
#'
#' @param df_ data.frame of raw drug response data containing both treated and untreated values.
#' @param readout string of the name containing the cell viability readout values.
#' @param discard_keys character vector of column names to include in the data.frames in the assays of the resulting \code{SummarizedExperiment} object. 
#' Defaults to \code{NULL}. 
#'
#' @seealso normalize_SE
#' @details 
#' This is most commonly used in preparation for downstream normalization.
#'
#' @export
#'
create_SE2 <- function(df_, readout = "ReadoutValue", discard_keys = NULL) {
  # Assertions:
  stopifnot(any(inherits(df_, "data.frame"), inherits(df_, "DataFrame")))
  checkmate::assert_string(readout)
  checkmate::assert_character(discard_keys, null.ok = TRUE)

  Keys <- identify_keys(df_)
  Keys$discard_keys <- discard_keys

  if (!is.null(discard_keys)) {
    Keys$DoseResp <- setdiff(Keys$DoseResp, discard_keys)
  }

  if (!(gDRutils::get_identifier("masked_tag") %in% colnames(df_))) {
    df_[, gDRutils::get_identifier('masked_tag')] <- FALSE
  }

  # Remove background value from readout (at least 1e-10 to avoid artefactual normalized values).
  df_$CorrectedReadout <- pmax(df_$ReadoutValue - df_$BackgroundValue, 1e-10)

  ## Identify treatments, conditions, and experiment metadata.
  md <- getMetaData(df_, discard_keys = discard_keys)
  coldata <- md$condition_md
  rowdata <- md$treatment_md
  exp_md <- md$experiment_md

  mapping_entries <- .create_mapping_factors(rowdata, coldata)
  mapping_entries$groupings <- rownames(mapping_entries) 

  ## Identify treated and untreated conditions.
  assigned_mapping_entries <- .assign_treated_and_untreated_conditions(mapping_entries)
  split_list <- split(assigned_mapping_entries, f = assigned_mapping_entries$treated_untreated)  
  if (length(split_list) != 2L) {
    stop(sprintf("unexpected conditions found: '%s'", 
      paste(setdiff(levels(assigned_mapping_entries$treated_untreated), c("treated", "untreated")), collapse = ",")))
  }
  treated <- split_list[["treated"]]
  untreated <- split_list[["untreated"]]

  ## Map references.
  # Map untreated references. 
  row_endpoint_value_filter <- rep(TRUE, nrow(untreated))
  Keys$untrt_Endpoint <- setdiff(Keys$untrt_Endpoint, names(key_values))
  untrt_endpoint_map <- map_df(treated, untreated, row_endpoint_value_filter, Keys, ref_type = "untrt_Endpoint")

  # Map day0 references.
  Keys$Day0 <- setdiff(Keys$Day0, names(key_values))
  day0_endpoint_map <- map_df(treated, untreated, row_endpoint_value_filter, Keys, ref_type = "Day0")

  # Map cotreatment references.
  Keys$ref_Endpoint <- setdiff(Keys$ref_Endpoint, names(key_values))
  cotrt_endpoint_map <- map_df(treated, treated, row_endpoint_value_filter, Keys, ref_type = "ref_Endpoint")

  ## Check for failed cotreatment mappings. 

  ## Combine all references with respective treatments.
  
  out <- vector("list", nrow(treated))
  for (trt in rownames(treated)) {
    trt_df <- treated[rownames()]  

    untrt_ref <- untrt_endpoint_map[[trt]]  
    untrt_df <- untreated[rownames(untreated) %in% untrt_ref, ]
    untrt_df <- create_control_df(untrt_df, control_mean_fxn, out_col_name = "UntrtReadout")

    day0_ref <- day0_endpoint_map[[trt]]
    day0_df <- split_list$untreated[rownames(split_list$untreated) %in% day0_ref]
    day0_df <- create_control_df(untrt_df, control_mean_fxn, out_col_name = "day0Readout")

    cotrt_ref <- cotrt_endpoint_map[[trt]]  
    cotrt_df <- split_list$treated[rownames(split_list$treated) %in% cotrt_ref]
    cotrt_df <- create_control_df(untrt_df, control_mean_fxn, out_col_name = "RefReadout")
 
    # Merge all data.frames together.
    ref_df <- merge(untrt)
    out[[i]] <- ref_df

  }
  out_df <- do.call(rbind, out)
  
  ## Create BumpyMatrix.
  df_to_assay(out_df)

  ## create the BumpyMatrix with the proper raw data values.
  # the matrices are named: "RawTreated" and "UntreatedReferences"
  matsL <-
    df_to_assays(
      data = df_,
      meta_data = md,
      endpoint_map = untrt_endpoint_map,
      discard_keys = discard_keys
    )
  
  # Capture important values in experiment metadata.
  experiment_md <- list(exp_md = exp_md, df_ = df_, Keys = Keys)

  # filter out to 'treated' conditions only 
  filtered_rowdata <-
    rowdata[rownames(rowdata) %in% rownames(matsL[[1]]), ] 
  
  # expect at least single treatment condition
  stopifnot(nrow(filtered_rowdata) > 0)
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
    colData = coldata,
    rowData = filtered_rowdata,
    metadata = experiment_md)
  se
}
  # Identify groupings on the original df. 
  df_ <- merge(df_, assigned_mapping_entries, by = c(colnames(rowdata), colnames(coldata)), all.x = TRUE)

  ## Join the metadata mappings back with the original data.
  mapping_entries <- merge(mapping_entries, untrt_endpoint_map) # Check that the other references are filled with NAs. 

