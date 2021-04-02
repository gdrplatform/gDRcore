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
        if ("Keys" %in% names(S4Vectors::metadata(normSE))) {
          TrtKeys <- S4Vectors::metadata(normSE)$Keys$Trt
          TrtKeys <- setdiff(TrtKeys, S4Vectors::metadata(normSE)$Keys$discard_keys)
        } else {
          TrtKeys <- identify_keys(normSE)$Trt
        }
    }
    S4Vectors::metadata(normSE)$Keys$Trt <- TrtKeys

    SummarizedExperiment::assay(avgSE, "Averaged") <- SummarizedExperiment::assay(avgSE, "Normalized")
    avgSE_assay <- assay(avgSE, "Averaged")
    
    for (i in seq_len(nrow(avgSE))) {
        for (j in seq_len(ncol(avgSE))) {
        x <- avgSE_assay[[i,j]]
        # bypass 'masked' filter
        x$masked <- x$masked & !include_masked

        subKeys <- intersect(TrtKeys, colnames(x))
        if (sum(!x$masked) >= 1) {
            df_av <- stats::aggregate(x[ !x$masked ,
                                  c("GRvalue", "RelativeViability","CorrectedReadout")],
                            by = as.list(x[ !x$masked , subKeys, drop = FALSE]),
                            FUN = function(y) mean(y, na.rm = TRUE))
            df_std <- stats::aggregate(x[!x$masked, c("GRvalue", "RelativeViability")],
                                by = as.list(x[ !x$masked, subKeys, drop = FALSE]),
                                FUN = function(x) stats::sd(x, na.rm = TRUE))
            colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")] =
                paste0("std_",
                    colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")])
            df_ <- merge(df_av, df_std, by = subKeys) 
        } else { # case: (nrow(x) == 0 || all(x$masked))
            df_ <- as.data.frame(matrix(0,0,length(subKeys)+5))
            colnames(df_) = c(subKeys,
                  c("GRvalue", "RelativeViability","CorrectedReadout"),
                  paste0("std_", c("GRvalue", "RelativeViability")))
        } 
        avgSE_assay[[i,j]] <- df_
    }}
    assay(avgSE, "Averaged") <- avgSE_assay

    SummarizedExperiment::assay(avgSE, "Avg_Controls") <- SummarizedExperiment::assay(avgSE, "Controls")
    avgSE <- aapply(avgSE, function(x) {
        if (nrow(x) > 1) {
            subKeys <- intersect(TrtKeys, colnames(x))
            df_av <- S4Vectors::DataFrame(lapply(x[, c("Day0Readout", "UntrtReadout",
                    "RefGRvalue", "RefRelativeViability",
                    "RefReadout", "DivisionTime")], FUN = function(y) mean(y, na.rm = TRUE)))
            return( df_av )
        } else return(x)
    }, "Avg_Controls")

    return(avgSE)
}


#' average_SE2
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
#' @details Expects that \code{get_SE_keys(se)} will have values for both 
#' \code{"Trt"} and \code{"masked_tag"}.
#'
#' @family runDrugResponseProcessingPipelineFxns
#'
#' @export
#'
average_SE2 <- function(se, 
                        override_masked = FALSE, 
                        normalized_assay = "Normalized", 
                        averaged_assay = "Averaged") {

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  if (!(normalized_assay %in% SummarizedExperiment::assayNames(se))) {
    stop(sprintf("missing expected assays: '%s'", normalized_assay))
  }

  trt_keys <- get_SE_keys(se, "Trt")
  masked_tag_key <- get_SE_keys(se, "masked_tag")

  if (!(length(trt_keys) > 0L && trt_keys != "")) {
    stop("unexpected treated keys on 'se' object")
  }

  if (!(length(masked_tag_key) > 0L && masked_tag_key != "")) {
    stop("unexpected masked_tag on 'se' object")
  }

  normalized <- SummarizedExperiment::assay(se, normalized_assay)

  std_cols <- c("GRvalue", "RelativeViability")
  out <- vector("list", nrow(se) * ncol(se))
  count <- 1
  for (i in seq_len(nrow(se))) {
    for (j in seq_len(ncol(se))) {
      norm_df <- normalized[i, j][[1]]

      # bypass 'masked' filter
      masked <- norm_df[[masked_tag_key]] & !override_masked

      if (sum(!masked) > 0) {
	p_trt_keys <- intersect(trt_keys, colnames(norm_df))
	
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
	agg_df <- S4Vectors::DataFrame(matrix(NA, 1, length(std_cols)*2 + length(p_trt_keys) + 2))
	colnames(agg_df) <- c(std_cols, p_trt_keys, paste0("std_", std_cols), 'row_id', 'col_id')
      }

      if (nrow(agg_df) != 0L) {
	agg_df$row_id <- rep(rownames(se)[i], nrow(agg_df))
	agg_df$col_id <- rep(colnames(se)[j], nrow(agg_df))
      }
      out[[count]] <- agg_df
      count <- count + 1
    }
  }
  # account for cases in which all values for a given condition are masked
  out <- S4Vectors::DataFrame(do.call("rbind", out))
  averaged <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(out) %in% c(masked_tag_key, "row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  # needs to specify rownames/colnames(averaged) in case that the matrix is not complete (e.g all values for a given cell line are masked)
  SummarizedExperiment::assays(se)[[averaged_assay]] <- averaged
  return(se)
}
