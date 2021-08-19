#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
create_SE <- function(df_,
                      readout = "ReadoutValue",
                      control_mean_fxn = function(x) {
                        mean(x, trim = 0.25)
                      },
                      nested_keys = c(gDRutils::get_identifier("concentration"),
                                      gDRutils::get_identifier("barcode")),
                      override_untrt_controls = NULL) {

  # Assertions:
  stopifnot(any(inherits(df_, "data.frame"), inherits(df_, "DataFrame")))
  checkmate::assert_string(readout)
  checkmate::assert_character(nested_keys, null.ok = TRUE)

  if (is(df_, "data.table")) {
    df_ <- S4Vectors::DataFrame(df_)
  }

  identifiers <- gDRutils::get_identifier()
  Keys <- identify_keys(df_, nested_keys, override_untrt_controls, identifiers)

  if (!(identifiers$masked_tag %in% colnames(df_))) {
    df_[, identifiers$masked_tag] <- FALSE
  }

  # Remove background value from readout (at least 1e-10 to avoid artefactual normalized values).
  df_$CorrectedReadout <- pmax(df_$ReadoutValue - df_$BackgroundValue, 1e-10)

  ## Identify treatments, conditions, and experiment metadata.
  md <- gDRutils::split_SE_components(df_, nested_keys = Keys$nested_keys)
  coldata <- md$condition_md
  rowdata <- md$treatment_md
  exp_md <- md$experiment_md

  mapping_entries <- .create_mapping_factors(rowdata, coldata)
  mapping_entries$groupings <- rownames(mapping_entries) 

  ## Identify treated and untreated conditions.
  assigned_mapping_entries <- .assign_treated_and_untreated_conditions(mapping_entries,
    c(identifiers$drugname, paste0(identifiers$drugname, "_2")))
  split_list <- split(mapping_entries, f = assigned_mapping_entries$treated_untreated)  
  if (length(split_list) != 2L) {
    stop(sprintf("unexpected conditions found: '%s'", 
      paste(setdiff(levels(assigned_mapping_entries$treated_untreated), c("treated", "untreated")), collapse = ",")))
  }
  treated <- split_list[["treated"]]
  untreated <- split_list[["untreated"]]

  ## Map references.
  references <- list(untrt_Endpoint = "untrt_Endpoint", Day0 = "Day0")

  ref_maps <- lapply(references, function(ref_type) {
    map_df(treated, untreated, override_untrt_controls = override_untrt_controls,
           ref_cols = Keys[[ref_type]], ref_type = ref_type)
  })

  ## Combine all references with respective treatments.
  # Merge raw data back with groupings.
  dfs <- merge(df_, mapping_entries, by = c(colnames(rowdata), colnames(coldata)))

  # Remove all rowdata and coldata. 
  groupings <- dfs$groupings
  dfs <- dfs[c(md$data_fields, "row_id", "col_id")]

  trt_out <- ref_out <- vector("list", nrow(treated))
  for (i in seq_len(nrow(treated))) {
    trt <- rownames(treated)[i]
    trt_df <- dfs[groupings == trt, , drop = FALSE]  

    if (nrow(trt_df) == 0L) {
      next # do nothing, there is no data to handle
    } else {
      ref_type <- "untrt_Endpoint"
      untrt_ref <- ref_maps[[ref_type]][[trt]]  
      untrt_df <- dfs[groupings %in% untrt_ref, , drop = FALSE]
      untrt_df <- create_control_df(
        untrt_df, 
        control_cols = Keys[[ref_type]], 
        control_mean_fxn, 
        out_col_name = "UntrtReadout"
      )

      ref_type <- "Day0"
      day0_ref <- ref_maps[[ref_type]][[trt]]
      day0_df <- dfs[groupings %in% day0_ref, , drop = FALSE]
      day0_df <- create_control_df(
        day0_df, 
        control_cols = Keys[[ref_type]],
        control_mean_fxn, 
        out_col_name = "Day0Readout"
      )

      ## Merge all data.frames together.
      # Try to merge by plate, but otherwise just use mean. 
      ref_df <- untrt_df 
      merge_cols <- intersect(colnames(day0_df), Keys$nested_keys)
      if (nrow(day0_df) > 0L) {
        ref_df <- merge(untrt_df, day0_df, by = merge_cols, all = TRUE)
      } else {
        ref_df$Day0Readout <- NA
      } 
      
      row_id <- unique(trt_df$row_id)
      col_id <- unique(trt_df$col_id)
      if (length(row_id) != 1L || length(col_id) != 1L) {
        stop(sprintf("non-unique row_ids: '%s' and col_ids: '%s'", 
          paste0(row_id, collapse = ", "), paste0(col_id, collapse = ", ")))
      }
      ref_df$row_id <- row_id
      ref_df$col_id <- col_id
    
      ref_out[[i]] <- ref_df
      trt_out[[i]] <- trt_df
    }

  }

  names(ref_out) <- names(trt_out) <- rownames(treated)
  
  trt_out <- do.call(rbind, trt_out)
  ref_out <- do.call(rbind, ref_out)
  trt_keep <- !colnames(trt_out) %in% c("row_id", "col_id")
  ref_keep <- !colnames(ref_out) %in% c("row_id", "col_id")

  treated_mat <- BumpyMatrix::splitAsBumpyMatrix(trt_out[, trt_keep, drop = FALSE],
                                                 row = trt_out$row_id, col = trt_out$col_id)
  reference_mat <- BumpyMatrix::splitAsBumpyMatrix(ref_out[, ref_keep, drop = FALSE],
                                                   row = ref_out$row_id, col = ref_out$col_id)
  matsL <- list(RawTreated = treated_mat, Controls = reference_mat)

  # Filter out to 'treated' conditions only.
  treated_rowdata <- rowdata[rownames(treated_mat), , drop = FALSE] 
  
  # Assertions.
  stopifnot(nrow(treated_rowdata) > 0)
  stopifnot(nrow(treated_rowdata) == length(unique(trt_out$row_id)))
 
  se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
    colData = coldata[match(colnames(treated_mat), rownames(coldata)), ],
    rowData = treated_rowdata,
    metadata = list(df_ = df_))

  # Capture important values in experiment metadata.
  se <- gDRutils::set_SE_identifiers(se, identifiers)
  se <- gDRutils::set_SE_experiment_metadata(se, exp_md)
  se <- gDRutils::set_SE_keys(se, Keys)

  se
}
