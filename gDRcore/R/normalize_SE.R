#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
normalize_SE <- function(se, 
                         nested_identifiers = gDRutils::get_SE_identifiers(se, "concentration"),
                         nested_confounders = gDRutils::get_SE_identifiers(se, "barcode"),
                         control_assay = "Controls", 
                         raw_treated_assay = "RawTreated", 
                         normalized_assay = "Normalized",
                         ndigit_rounding = 4) {
  
  # Assertions
  checkmate::assert_number(ndigit_rounding)
  gDRutils::validate_se_assay_name(se, control_assay)
  gDRutils::validate_se_assay_name(se, raw_treated_assay)

  # Keys
  nested_keys <- c(nested_identifiers, nested_confounders)
  masked_key <- gDRutils::get_SE_identifiers(se, "masked_tag")
  trt_keys <- gDRutils::get_SE_keys(se, key_type = "Trt")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)

  cl_names <- cdata[, gDRutils::get_SE_identifiers(se, "cellline_name"), drop = FALSE]
  cl_ref_div_times <- cdata[, gDRutils::get_SE_identifiers(se, "cellline_ref_div_time"), drop = FALSE]
  durations <- rdata[, gDRutils::get_SE_identifiers(se, "duration"), drop = FALSE]

  refs <- SummarizedExperiment::assays(se)[[control_assay]]
  trt <- SummarizedExperiment::assays(se)[[raw_treated_assay]]
  
  # Extract common nested_confounders shared by trt_df and ref_df
  nested_confounders <- Reduce(intersect, list(nested_confounders,
                                               names(BumpyMatrix::unsplitAsDataFrame(trt)),
                                               names(BumpyMatrix::unsplitAsDataFrame(refs))))
  
  norm_cols <- c("RelativeViability", "GRvalue")
  out <- vector("list", nrow(se) * ncol(se))
  ref_rel_viability <- ref_GR_value <- div_time <- matrix(NA, nrow = nrow(se), ncol = ncol(se), dimnames = dimnames(se))
  msgs <- NULL
  # Column major order, so go down first.
  for (j in seq_along(colnames(se))) {
    cl_name <- cl_names[j, ]
    ref_div_time <- cl_ref_div_times[j, ]

    for (i in seq_along(rownames(se))) {
      duration <- durations[i, ]

      ref_df <- refs[i, j][[1]]
      trt_df <- trt[i, j][[1]]

      if (length(trt_df) == 0L || nrow(trt_df) == 0L) {
        next # skip if no data
      }

      if (length(ref_df) == 0L || nrow(ref_df) == 0L) {
        msgs <- c(msgs,
          sprintf("Missing control data. Treatment Id: '%s' Cell_line Id: '%s'", 
            rownames(se)[i], colnames(se)[j]))
        next
      }

      # pad the ref_df for missing values based on nested_keys (uses mean across all available values)
      if (!is.null(nested_keys) && length(nested_keys) > 0) {
        ref_df <- fill_NA_ref(ref_df, nested_keys)
      }
      
      # Merge to ensure that the proper discard_key values are mapped.
      all_readouts_df <- merge(trt_df, 
        ref_df, 
        by = nested_confounders,
        all.x = TRUE)

      normalized <- S4Vectors::DataFrame(matrix(NA, nrow = nrow(trt_df), ncol = length(norm_cols)))
      colnames(normalized) <- c(norm_cols)

      # Normalized treated.
      normalized$RelativeViability <- round(all_readouts_df$CorrectedReadout /
                                              all_readouts_df$UntrtReadout, ndigit_rounding)
      normalized$GRvalue <- calculate_GR_value(rel_viability = normalized$RelativeViability, 
        corrected_readout = all_readouts_df$CorrectedReadout, 
        day0_readout = all_readouts_df$Day0Readout, 
        untrt_readout = all_readouts_df$UntrtReadout, 
        ndigit_rounding = ndigit_rounding, 
        duration = duration, 
        ref_div_time = ref_div_time, 
        cl_name = cl_name)

      # Carry over present treated keys.
      keep <- intersect(c(nested_identifiers, trt_keys, masked_key), colnames(all_readouts_df))
      normalized <- cbind(all_readouts_df[keep], normalized) 

      normalized$row_id <- rep(rownames(se)[i], nrow(trt_df))
      normalized$col_id <- rep(colnames(se)[j], nrow(trt_df))

      out[[nrow(se) * (j - 1) + i]] <- normalized
    }
  }

  if (!is.null(msgs)) {
    futile.logger::flog.warn(paste0(msgs, collapse = "\n"))
  }
  out <- do.call("rbind", out)
  
  norm <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(normalized) %in% c("row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assays(se)[[normalized_assay]] <- norm
  se
}


#' @keywords internal
fill_NA_ref <- function(ref_df, nested_keys) {
  data_columns <- setdiff(colnames(ref_df), nested_keys)
  ref_df_mean <- colMeans(data.frame(ref_df[, data_columns, drop = FALSE]), na.rm = TRUE)
  for (col in data_columns) {
    ref_df[is.na(ref_df[, col]), col] <- ref_df_mean[[col]]
  }
  S4Vectors::DataFrame(ref_df)
}
