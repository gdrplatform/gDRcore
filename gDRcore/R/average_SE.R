#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
average_SE <- function(se,
                       nested_identifiers = gDRutils::get_SE_identifiers(se, "concentration"),
                       override_masked = FALSE,
                       normalized_assay = "Normalized", 
                       averaged_assay = "Averaged") {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  gDRutils::validate_se_assay_name(se, normalized_assay)

  trt_keys <- gDRutils::get_SE_keys(se, "Trt")
  masked_tag_key <- gDRutils::get_SE_keys(se, "masked_tag")

  if (!(length(trt_keys) > 0L && trt_keys != "")) {
    stop("unexpected treated keys on 'se' object")
  }

  if (!(length(masked_tag_key) > 0L && masked_tag_key != "")) {
    stop("unexpected masked_tag on 'se' object")
  }

  normalized <- SummarizedExperiment::assay(se, normalized_assay)

  std_cols <- c("GRvalue", "RelativeViability")
  out <- vector("list", nrow(se) * ncol(se))
  count <- 0
  for (i in seq_len(nrow(se))) {
    for (j in seq_len(ncol(se))) {
      count <- count + 1
      norm_df <- normalized[i, j][[1]]
      if (nrow(norm_df) == 0L) {
        next
        }

      # bypass 'masked' filter
      masked <- norm_df[[masked_tag_key]] & !override_masked

      p_trt_keys <- intersect(trt_keys, colnames(norm_df))

      if (sum(!masked) > 0) {
        if (length(missing <- setdiff(std_cols, colnames(norm_df))) > 0L) {
          stop(sprintf("missing expected columns in nested normalized dataframe: '%s'", 
            paste0(missing, collapse = ", ")))
        }

        avg_df <- stats::aggregate(norm_df[!masked, std_cols],
          by = as.list(norm_df[!masked, nested_identifiers, drop = FALSE]), 
          function(x) mean(x, na.rm = TRUE))

        std_df <- stats::aggregate(norm_df[!masked, std_cols],
          by = as.list(norm_df[!masked, nested_identifiers, drop = FALSE]), 
          function(x) stats::sd(x, na.rm = TRUE))
        colnames(std_df)[colnames(std_df) %in% std_cols] <-
          paste0("std_", colnames(std_df)[colnames(std_df) %in% std_cols])

        agg_df <- S4Vectors::DataFrame(merge(avg_df, std_df, by = nested_identifiers))
      } else {
        # only one or no unmasked value --> create degenerated dataframe
        all_cols <- c(std_cols, p_trt_keys, paste0("std_", std_cols), "row_id", "col_id")
        agg_df <- S4Vectors::DataFrame(matrix(NA, 1, length(all_cols)))
        colnames(agg_df) <- all_cols
      }

      if (nrow(agg_df) != 0L) {
        agg_df$row_id <- rep(rownames(se)[i], nrow(agg_df))
        agg_df$col_id <- rep(colnames(se)[j], nrow(agg_df))
      }
      out[[count]] <- agg_df
    }
  }
  
  out <- S4Vectors::DataFrame(do.call("rbind", out))
  averaged <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c(masked_tag_key, "row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assays(se)[[averaged_assay]] <- averaged
  se
}
