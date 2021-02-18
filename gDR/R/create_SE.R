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
#' @param control_mean_fxn function indicating how to average controls.
#' Defaults to \code{mean(x, trim = 0.25)}.
#' @param key_values
#' @param discard_keys character vector of column names to include in the data.frames in the assays of the resulting \code{SummarizedExperiment} object. 
#' Defaults to \code{NULL}. 
#'
#' @seealso normalize_SE
#' @details 
#' This is most commonly used in preparation for downstream normalization.
#'
#' @export
#'
create_SE2 <- function(df_, 
                       readout = "ReadoutValue", 
                       control_mean_fxn = function(x) {mean(x, trim = 0.25)}, 
                       key_values = NULL,
                       discard_keys = NULL) {

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
  # First look amongst the untreated.
  cotrt_endpoint_map <- map_df(treated, untreated, row_endpoint_value_filter, Keys, ref_type = "ref_Endpoint")

  if ("Gnumber_2" %in% colnames(treated)) {
    # Remove Gnumber_2, DrugName_2, and Concentration_2.
    Keys$ref_Endpoint <- setdiff(Keys$ref_Endpoint, 
      c(names(key_values), c("Gnumber_2", "DrugName_2", "Concentration_2")))

    # Then look amongst the treated to fill any missing cotrt references.
    missing_cotrt <- vapply(cotrt_endpoint_map, function(x) {length(x) == 0L}, TRUE)
    if (any(missing_cotrt)) {
      missing_mappings <- names(cotrt_endpoint_map)[missing_cotrt]
      missing_trt_mappings <- treated[treated$groupings %in% names(missing_mappings)]

      missing_cotrt_endpoint_map <- map_df(missing_trt_mappings, treated, row_endpoint_value_filter, Keys, ref_type = "ref_Endpoint")

      # Fill found mappings.
      for (grp in names(missing_cotrt_endpoint_map)) {
        if (length(missing_cotrt_endpoint_map[[grp]]) != 0L) {
          cotrt_endpoint_map[grp] <- missing_cotrt_endpoint_map[[grp]]
        }
      }
    }
  }

  ## TODO: Check for failed cotreatment mappings. 

  ## Combine all references with respective treatments.
  # Merge raw data back with groupings.
  dfs <- merge(df_, assigned_mapping_entries, by = c(colnames(rowdata), colnames(coldata)))

  # Remove all rowdata and coldata. 
  dfs <- dfs[!colnames(dfs) %in% c(colnames(rowdata), colnames(coldata), "treated_untreated")]
  df_cols <- colnames(dfs) != "groupings"
  trt_out <- ref_out <- vector("list", nrow(treated))
  for (i in seq_len(nrow(treated))) {
    trt <- rownames(treated)[i]
    trt_df <- dfs[dfs$groupings == trt, df_cols, drop = FALSE]  

    ref_df <- NULL
    if (nrow(trt_df) > 0L) {
      untrt_ref <- untrt_endpoint_map[[trt]]  
      untrt_df <- dfs[dfs$groupings %in% untrt_ref, df_cols, drop = FALSE]
      untrt_df <- create_control_df(
        untrt_df, 
        Keys = Keys, 
        key = "untrt_Endpoint", 
        control_mean_fxn, 
        out_col_name = "UntrtReadout"
      )

      day0_ref <- day0_endpoint_map[[trt]]
      day0_df <- dfs[dfs$groupings %in% day0_ref, df_cols, drop = FALSE]
      day0_df <- create_control_df(
        day0_df, 
        Keys = Keys, 
        key = "Day0", 
        control_mean_fxn, 
        out_col_name = "Day0Readout"
      )

      cotrt_ref <- cotrt_endpoint_map[[trt]]  
      if (length(cotrt_ref) == 0L) {
	# Set the cotrt reference to the untreated reference.
	cotrt_df <- untrt_df 
	colnames(cotrt_df)[grepl("UntrtReadout", colnames(cotrt_df))] <- "RefReadout"
      } else {
	cotrt_df <- dfs[dfs$groupings %in% cotrt_ref, , drop = FALSE]
	cotrt_df <- create_control_df(
	  cotrt_df, 
	  Keys = Keys, 
	  key = "ref_Endpoint", 
	  control_mean_fxn, 
	  out_col_name = "RefReadout"
	)
      }
   
      ## Merge all data.frames together.
      # Try to merge by plate, but otherwise just use mean.. 
      ref_df <- untrt_df
      if (nrow(cotrt_df) > 0L) {
        #ref_conc <- SummarizedExperiment::rowData(normSE)[i, 'Concentration_2']
        #if () # See whether the reference is in treated or untreated.

	merge_cols <- intersect(colnames(cotrt_df), c("Barcode", discard_keys))
	ref_df <- merge(untrt_df, cotrt_df[, c("RefReadout", merge_cols)], by = merge_cols, all = TRUE)
        if (!all(sort(unique(cotrt_df$Barcode)) == sort(unique(untrt_df$Barcode)))) {
          # Merging by barcodes will result in NAs. 
	  ref_df$UntrtReadout[is.na(ref_df$UntrtReadout)] <- mean(ref_df$UntrtReadout, na.rm = TRUE)
	  ref_df$RefReadout[is.na(ref_df$RefReadout)] <- mean(ref_df$RefReadout, na.rm = TRUE)
	  # TODO: Should this also use the control_mean_fxn? 
	}
	  
      } else {
	ref_df$RefReadout <- ref_df$UntrtReadout
      }

      if (nrow(day0_df) > 0L) {
	ref_df <- merge(day0_df[, setdiff(colnames(day0_df), "Barcode"), drop = FALSE], ref_df)
      } else {
        ref_df$Day0Readout <- NA
      }

    } else {
      trt_df <- NULL 
    }

    if (!is.null(ref_df)) {
      row_id <- unique(trt_df$row_id)
      col_id <- unique(trt_df$col_id)
      if (length(row_id) != 1L || length(col_id) != 1L) {
        stop(sprintf("non-unique row_ids: '%s' and col_ids: '%s'", 
          paste0(row_id, collapse = ", "), paste0(col_id, collapse = ", ")))
      }
      ref_df$row_id <- row_id
      ref_df$col_id <- col_id
    }

    ref_out[[i]] <- ref_df
    trt_out[[i]] <- trt_df
  }

  names(ref_out) <- names(trt_out) <- rownames(treated)

  trt_out <- do.call(rbind, trt_out)
  ref_out <- do.call(rbind, ref_out)
  
  keep <- colnames(trt_out)[!colnames(trt_out) %in% c("row_id", "col_id")]
  #treated_mat <- BumpyMatrix::splitAsBumpyMatrix(trt_out[, keep], row = trt_out$row_id, col = trt_out$col_id, sparse = TRUE)
  treated_mat <- BumpyMatrix::splitAsBumpyMatrix(trt_out, row = trt_out$row_id, col = trt_out$col_id)

  keep <- colnames(ref_out)[!colnames(ref_out) %in% c("row_id", "col_id")]
  #reference_mat <- BumpyMatrix::splitAsBumpyMatrix(ref_out[, keep], row = ref_out$row_id, col = ref_out$col_id, sparse = TRUE)
  reference_mat <- BumpyMatrix::splitAsBumpyMatrix(ref_out, row = ref_out$row_id, col = ref_out$col_id)

  matsL <- list(RawTreated = treated_mat, Controls = reference_mat)

  # Capture important values in experiment metadata.
  experiment_md <- list(experiment_metadata = exp_md, df_ = df_, Keys = Keys)

  # filter out to 'treated' conditions only 
  filtered_rowdata <-
    rowdata[rownames(rowdata) %in% rownames(matsL[[1]]), ] 
  
  # expect at least single treatment condition
  stopifnot(nrow(filtered_rowdata) > 0)
 
  se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
    colData = coldata[match(colnames(treated_mat), rownames(coldata)), ],
    rowData = filtered_rowdata[match(rownames(treated_mat), rownames(filtered_rowdata)), ],
    metadata = experiment_md)
  se
}
