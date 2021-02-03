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
           assay_type = c("matrix", "BumpyMatrix"), 
           aggregate_ref_FXN = function(x) {mean(x, trim = 0.25)}) {
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
#' @param aggregate_ref_FXN function used for averaging references.
#' Defaults to trimmed arithmetic mean with trim = 0.25.
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

  if (!(gDRutils::get_identifier('masked_tag') %in% colnames(df_))) {
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

  # enforced key values for end points (override selected_keys) --> for rows of the SE
  # Keys$untrt_Endpoint <- setdiff(Keys$untrt_Endpoint, names(key_values))
  row_endpoint_value_filter <- rep(TRUE, nrow(split_list[["untreated"]]))

  ## Map the treatments to their references.
  untrt_endpoint_map <- map_df(split_list$treated, split_list$untreated, row_endpoint_value_filter, Keys, ref_type = "untrt_Endpoint")

  # Identify groupings on the original df. 
  df_ <- merge(df_, assigned_mapping_entries, by = c(colnames(rowdata), colnames(coldata)), all.x = TRUE)
  # Split this again to get the references. 
  split_list <- split(df_, f = df_$treated_untreated)
  split_list[["untreated"]]

  # Aggregate where there are multiple references for a single treatment. 
  refs <- unique(untrt_endpoint_map)
  n_refs <- length(refs) # Identify how many unique control groups there are.

  ref_cache <- vector("list", n_refs)
  names(ref_cache) <- vapply(refs, function(x) paste(x, collapse = "_"), character(0))

  for (i in seq_along(untrt_endpoint_map)) {
    trt_refs <- untrt_endpoint_map[[i]]
    if (length(trt_refs > 1L)) {
      agg_readout <- aggregate_ref_FXN(untrt[untrt$groupings %in% trt_refs, readout])
      # TODO: Figure out what values should actually go in here. Looks like we'll need "UntrtReadout". 
      # Will this be just a single value? If so, we don't need a BumpyMatrix, and can just create a matrix and put it in the matrix list.
      ref_cache[[i]] <- aggregate_ref_FXN(untrt[untrt$groupings %in% trt_refs, readout])
    }
  }

  ## Join the metadata mappings back with the original data.
  mapping_entries <- merge(mapping_entries, untrt_endpoint_map) # Check that the other references are filled with NAs. 

  ## TODO: Create the BumpyMatrix with the proper normalized values.
  ###############
  ## FILL ME
  ###############

  matsL <- NULL
  #matsL <- list(mats)
  #names(matsL) <- readout

  # Capture important values in experiment metadata.
  experiment_md <- c(exp_md, list(df_ = df_, Keys = Keys))
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
    colData = coldata,
    rowData = rowdata, 
    metadata = experiment_md)
  se
}
