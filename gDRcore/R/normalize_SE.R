#' normalize_SE
#'
#' Normalize drug response data from treated and untreated pairings.
#'
#' @param se \code{BumpyMatrix} object with assays \code{"RawTreated"} and \code{"Controls"}.
#' @param control_assay string containing the name of the assay representing the controls in the \code{se}.
#' Defaults to \code{"Controls"}.
#' @param raw_treated_assay string containing the name of the assay representing the raw treated data in the \code{se}.
#' Defaults to \code{"RawTreated"}.
#' @param normalized_assay string of the assay name containing the normalized data.
#' Defaults to \code{"Normalized"}.
#' @param ref_GR_assay string of the name of the reference GR assay in the \linkS4class{SummarizedExperiment}.
#' @param ref_RV_assay string of the name of the reference Relative Viability assay in the \linkS4class{SummarizedExperiment}.
#' @param ndigit_rounding integer specifying number of digits of rounding during calculations.
#' Defaults to \code{4}.
#'
#' @return \code{BumpyMatrix} object with a new assay named \code{"Normalized"} containing \code{DataFrame}s 
#' holding \code{RelativeViability} and \code{GRvalue}, as well as new assays named \code{RefRelativeViability}, \code{RefGRvalue}, and \code{DivisionTime} values.
#'
#' @family runDrugResponseProcessingPipelineFxns
#' @export
#'
normalize_SE <- function(se, 
                          control_assay = "Controls", 
                          raw_treated_assay = "RawTreated", 
                          normalized_assay = "Normalized", 
                          ref_GR_assay = "RefGRvalue", 
                          ref_RV_assay = "RefRelativeViability", 
                          ndigit_rounding = 4) {

  # Assertions
  checkmate::assert_number(ndigit_rounding)
  if (!all(c(raw_treated_assay, control_assay) %in% SummarizedExperiment::assayNames(se))) {
    stop(sprintf("missing expected assays: '%s'", 
      setdiff(c(raw_treated_assay, control_assay), SummarizedExperiment::assayNames(se))))
  }

  refs <- SummarizedExperiment::assays(se)[[control_assay]]
  trt <- SummarizedExperiment::assays(se)[[raw_treated_assay]]

  nested_keys <- gDRutils::get_SE_keys(se, key_type = "nested_keys")
  trt_keys <- gDRutils::get_SE_keys(se, key_type = "Trt")
  cl_name_key <- gDRutils::get_SE_keys(se, "cellline_name")
  cl_ref_div_time_key <- gDRutils::get_SE_keys(se, "cellline_ref_div_time")
  duration_key <- gDRutils::get_SE_keys(se, "duration")
  masked_key <- gDRutils::get_SE_keys(se, "masked_tag")

  norm_cols <- c("RelativeViability", "GRvalue")
  out <- vector("list", nrow(se) * ncol(se))

  ref_rel_viability <- ref_GR_value <- div_time <- matrix(NA, nrow = nrow(se), ncol = ncol(se), dimnames = dimnames(se))

  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)

  # Column major order, so go down first.
  for (j in seq_along(colnames(se))) {
    cl_name <- cdata[j, cl_name_key]
    ref_div_time <- cdata[j, cl_ref_div_time_key]

    for (i in seq_along(rownames(se))) {
      duration <- rdata[i, duration_key]

      ref_df <- refs[i, j][[1]]
      trt_df <- trt[i, j][[1]]

      if (length(trt_df) == 0L || nrow(trt_df) == 0L) {
        next # skip if no data
        # TODO: Does this still need to initialize an empty DFrame with appropriate colnames?
      }

      if (length(ref_df) == 0L || nrow(ref_df) == 0L) {
        futile.logger::flog.warn(
          sprintf("Missing control data. Treatment Id: '%s' Cell_line Id: '%s'", 
            rownames(se)[i], colnames(se)[j])
        )
        next
      }

      # pad the ref_df for missing values based on nested_keys (uses mean across all available values)
      if (!is.null(nested_keys) && length(nested_keys) > 0) {
        ref_df_complete <- unique(trt_df[,nested_keys,drop=FALSE])
        ref_df_complete <- merge(ref_df_complete, ref_df, by = nested_keys)
        data_columns <- setdiff(colnames(ref_df), nested_keys)
        ref_df_mean <- lapply(ref_df[, data_columns, drop=FALSE], function(x) mean(x, na.rm = TRUE))
        for (col in data_columns) {
          ref_df_complete[is.na(ref_df_complete[,col]), col] <- ref_df_mean[[col]]
        }
      } else {
        ref_df_complete <- ref_df
      }

      # Merge to ensure that the proper discard_key values are mapped.
      all_readouts_df <- merge(trt_df, 
        ref_df_complete, 
        by = nested_keys,
        all.x = TRUE)

      normalized <- S4Vectors::DataFrame(matrix(NA, nrow = nrow(trt_df), ncol = length(norm_cols)))
      colnames(normalized) <- c(norm_cols)

      # Normalized treated.
      normalized$RelativeViability <- round(all_readouts_df$CorrectedReadout/all_readouts_df$UntrtReadout, ndigit_rounding)
      normalized$GRvalue <- calculate_GR_value(rel_viability = normalized$RelativeViability, 
        corrected_readout = all_readouts_df$CorrectedReadout, 
        day0_readout = all_readouts_df$Day0Readout, 
        untrt_readout = all_readouts_df$UntrtReadout, 
        ndigit_rounding = ndigit_rounding, 
        duration = duration, 
        ref_div_time = ref_div_time, 
        cl_name = cl_name)

      # Carry over present treated keys.
      normalized <- cbind(all_readouts_df[intersect(c(trt_keys, masked_key), colnames(all_readouts_df))], normalized) 

      normalized$row_id <- rep(rownames(se)[i], nrow(trt_df))
      normalized$col_id <- rep(colnames(se)[j], nrow(trt_df))

      out[[nrow(se) * (j - 1) + i]] <- normalized

      ## This assignment can be used to make comparisons between old and refactored methods (GDR-621).
      ## However, the aggregation in this manner technically biases the references to their 
      ## corresponding treated entries, so is incorrect.
      #ref_df <- all_readouts_df 

      ## Perform the calculations on all references.
      ## Then, take the mean to report the final reference normalized value.
      RV_vec <- ref_df_complete$RefReadout/ref_df_complete$UntrtReadout
      GR_vec <- calculate_GR_value(rel_viability = RV_vec, 
        corrected_readout = ref_df_complete$RefReadout, 
        day0_readout = ref_df_complete$Day0Readout, 
        untrt_readout = ref_df_complete$UntrtReadout, 
        ndigit_rounding = ndigit_rounding, 
        duration = duration, 
        ref_div_time = ref_div_time, 
        cl_name = cl_name)

      ref_rel_viability[i, j] <- round(mean(RV_vec, na.rm=TRUE), ndigit_rounding)
      ref_GR_value[i, j] <- round(mean(GR_vec, na.rm = TRUE), ndigit_rounding)
      div_time[i, j] <- round(duration / log2(mean(ref_df$UntrtReadout/ref_df$Day0Readout, na.rm = TRUE)), ndigit_rounding)
    }
  }

  out <- do.call("rbind", out)
  norm <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(normalized) %in% c("row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assays(se)[[normalized_assay]] <- norm

  SummarizedExperiment::assays(se)[[ref_GR_assay]] <- ref_GR_value
  SummarizedExperiment::assays(se)[[ref_RV_assay]] <- ref_rel_viability
  SummarizedExperiment::assays(se)[["DivisionTime"]] <- div_time

  return(se)
}
