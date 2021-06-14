#' average_SE
#'
#' Average the normalized data within the nested \code{DataFrame}s. 
#'
#' @param se a \linkS4class{SummarizedExperiment} with drug response data that has a 
#' \code{normalized_assay}.
#' @param override_masked boolean indicating whether or not to override the masked wells
#' in the averaging and include all wells. 
#' Defaults to \code{FALSE}.
#' @param normalized_assay string of the assay name containing the normalized data.
#' Defaults to \code{"Normalized"}.
#' @param averaged_assay string of the assay to output averaged values in.
#' Defaults to \code{"Averaged"}.
#'
#' @return a \linkS4class{SummarizedExperiment} with an additional assay 
#' specified by \code{averaged_assay} containing normalized data with 
#' mean and standard deviation calculations for each unique treatment in the nested
#' \code{DataFrame}s. 
#'
#' @details Expects that \code{gDRutils::get_SE_keys(se)} will have values for both 
#' \code{"Trt"} and \code{"masked_tag"}.
#'
#' @family runDrugResponseProcessingPipelineFxns
#'
#' @export
#'
average_SE <- function(se, 
                        override_masked = FALSE, 
                        normalized_assay = "Normalized", 
                        averaged_assay = "Averaged") {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  if (!(normalized_assay %in% SummarizedExperiment::assayNames(se))) {
    stop(sprintf("missing expected assays: '%s'", normalized_assay))
  }

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
      if (nrow(norm_df) == 0L) {next}

      # bypass 'masked' filter
      masked <- norm_df[[masked_tag_key]] & !override_masked

      p_trt_keys <- intersect(trt_keys, colnames(norm_df))

      if (sum(!masked) > 0) {        
        if (length(missing <- setdiff(std_cols, colnames(norm_df))) > 0L) {
          stop(sprintf("missing expected columns in nested normalized dataframe: '%s'", 
            paste0(missing, collapse = ", ")))
        }

        avg_df <- stats::aggregate(norm_df[!masked, std_cols],
          by = as.list(norm_df[!masked, p_trt_keys, drop = FALSE]), 
          function(x) mean(x, na.rm = TRUE))

        std_df <- stats::aggregate(norm_df[!masked, std_cols],
          by = as.list(norm_df[!masked, p_trt_keys, drop = FALSE]), 
          function(x) stats::sd(x, na.rm = TRUE))
        colnames(std_df)[colnames(std_df) %in% std_cols] <-
          paste0("std_", colnames(std_df)[colnames(std_df) %in% std_cols])

        agg_df <- S4Vectors::DataFrame(merge(avg_df, std_df, by = p_trt_keys))
      } else {
        # only one or no unmasked value --> create degenerated dataframe
        all_cols <- c(std_cols, p_trt_keys, paste0("std_", std_cols), 'row_id', 'col_id')
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
  return(se)
}
