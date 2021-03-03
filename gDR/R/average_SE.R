#' average_SE
#'
#' Avereage normalized SummarizedExperiment of DR data
#'
#' @param normSE a SummarizedExperiment with normalized DR data
#' @param TrtKeys a vector of keys used for averaging (NULL by default)
#'
#' @return a SummarizedExperiment with additional assay with averaged DR data
#'
#' @export
#'
average_SE <- function(normSE, TrtKeys = NULL, include_masked = F) {
  .Deprecated(msg = "see average_SE2 for similar, but not identical functionality")

  # Assertions:
  checkmate::assert_class(normSE, "SummarizedExperiment")
  checkmate::assert_vector(TrtKeys, null.ok = TRUE)

    avgSE <- normSE
    if (is.null(TrtKeys)) {
        if ("Keys" %in% names(metadata(normSE))) {
          TrtKeys <- metadata(normSE)$Keys$Trt
          TrtKeys <- setdiff(TrtKeys, metadata(normSE)$Keys$discard_keys)
        } else {
          TrtKeys <- identify_keys(normSE)$Trt
        }
    }
    metadata(normSE)$Keys$Trt <- TrtKeys

    SummarizedExperiment::assay(avgSE, "Averaged") <- SummarizedExperiment::assay(avgSE, "Normalized")
    avgSE <- aapply(avgSE, function(x) {
        # bypass 'masked' filter
        x$masked <- x$masked & !include_masked

        subKeys <- intersect(TrtKeys, colnames(x))
        if (sum(!x$masked) >= 1) {
            df_av <- aggregate(x[ !x$masked ,
                                  c("GRvalue", "RelativeViability","CorrectedReadout")],
                            by = as.list(x[ !x$masked , subKeys, drop = FALSE]),
                            FUN = function(y) mean(y, na.rm = TRUE))
            df_std <- aggregate(x[!x$masked, c("GRvalue", "RelativeViability")],
                                by = as.list(x[ !x$masked, subKeys, drop = FALSE]),
                                FUN = function(x) sd(x, na.rm = TRUE))
            colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")] =
                paste0("std_",
                    colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")])
            return( merge(df_av, df_std, by = subKeys) )
        } else { # case: (nrow(x) == 0 || all(x$masked))
            df_ = as.data.frame(matrix(0,0,length(subKeys)+5))
            colnames(df_) = c(subKeys,
                  c("GRvalue", "RelativeViability","CorrectedReadout"),
                  paste0("std_", c("GRvalue", "RelativeViability")))
            return(df_)
        } 
    }, "Averaged")

    SummarizedExperiment::assay(avgSE, "Avg_Controls") <- SummarizedExperiment::assay(avgSE, "Controls")
    avgSE <- aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys <- intersect(TrtKeys, colnames(x))
            df_av <- DataFrame(lapply(x[, c("Day0Readout", "UntrtReadout",
                    "RefGRvalue", "RefRelativeViability",
                    "RefReadout", "DivisionTime")], FUN = function(y) mean(y, na.rm = TRUE)))
            return( df_av )
        } else return(x)
    }, "Avg_Controls")

    return(avgSE)
}


# TODO: Mention something about the std (standard deviation). 
#' average_SE2
#'
#' Average the normalized data per
#'
#' @param se a \linkS4class{SummarizedExperiment} with drug response data.
#' @param include_masked boolean indicating whether or not to include masked wells
#' in the averaging. 
#' This is used as an override to whatever wells have been masked in the original data.
#' @param normalized_assay string of the assay name containing the normalized data.
#' Defaults to \code{Normalized}.
#'
#' @return a SummarizedExperiment with additional assay with averaged DR data
#' @seealso runDrugResponseProcessingPipeline2
#'
#' @export
#'
average_SE2 <- function(se, 
                        include_masked = FALSE, 
                        normalized_assay = "Normalized") {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")

  normalized <- SummarizedExperiment::assay(se, normalized_assay)

  trt_keys <- get_SE_keys(se, "Trt")

  out <- vector("list", nrow(se) * ncol(se))
  count <- 1
  for (i in seq_len(nrow(se))) {
    for (j in seq_len(ncol(se))) {
      norm_df <- normalized[i, j][[1]]

      # bypass 'masked' filter
      masked <- norm_df[[gDRutils::get_identifier("masked_tag")]] & !include_masked

      if (nrow(norm_df[!masked, ]) > 1L) {
        p_trt_keys <- intersect(trt_keys, colnames(norm_df))
	std_cols <- c("GRvalue", "RelativeViability")

	avg_df <- aggregate(norm_df[!masked, std_cols],
	  by = as.list(norm_df[!masked, p_trt_keys, drop = FALSE]), 
	  function(x) mean(x, na.rm = TRUE))

	std_df <- aggregate(norm_df[!masked, std_cols],
	  by = as.list(norm_df[!masked, p_trt_keys, drop = FALSE]), 
	  function(x) {sd(x, na.rm = TRUE)})
        colnames(std_df)[colnames(std_df) %in% std_cols] <-
          paste0("std_", colnames(std_df)[colnames(std_df) %in% std_cols])

	agg_df <- merge(avg_df, std_df, by = p_trt_keys) 
      } else {
        agg_df <- norm_df
      }

      if (nrow(agg_df) != 0L) {
	agg_df$row_id <- rep(rownames(se)[i], nrow(agg_df))
	agg_df$col_id <- rep(colnames(se)[j], nrow(agg_df))
      }
      out[[count]] <- agg_df
      count <- count + 1
    }
  }

  out <- DataFrame(do.call("rbind", out))
  averaged <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c("masked", "row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assays(se)[["Averaged"]] <- averaged
  return(se)
}
