#' average_SE
#'
#' Average the assays of a SummarizedExperiment of drug response data
#'
#' @param se a \linkS4class{SummarizedExperiment} with drug response data.
#' @param TrtKeys a vector of keys used for averaging (NULL by default)
#' @param aggregate_FXN function used for averaging data that should be considered replicates.
#' Defaults to trimmed arithmetic mean with trim = 0.25.
#'
#' @return a SummarizedExperiment with additional assay with averaged DR data
#'
#' @export
#'
average_SE <- function(se, TrtKeys = NULL, include_masked = FALSE, 
  aggregate_FXN = function(x) {mean(x, na.rm = TRUE, trim = 0.25)}, treated_assay = "RawTreated", reference_assay = "UntreatedReferences") {
#
#  # Aggregate where there are multiple references for a single treatment. 
#  refs <- unique(untrt_endpoint_map)
#  n_refs <- length(refs) # Identify how many unique control groups there are.
#
#  ref_cache <- vector("list", n_refs)
#  names(ref_cache) <- vapply(refs, function(x) paste(x, collapse = "_"), character(0))
#
#  for (i in seq_along(untrt_endpoint_map)) {
#    trt_refs <- untrt_endpoint_map[[i]]
#    if (length(trt_refs > 1L)) {
#      # Note that the metadata no longer needs to be carried, as the only relevant information at this point is the mapping
#      # which will be captured through the matrix indices.  
#      agg_readout <- aggregate_FXN(untrt[untrt$groupings %in% trt_refs, readout])
#      ref_cache[[i]] <- agg_readout
#    }
#  }

  # Assertions:
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_vector(TrtKeys, null.ok = TRUE)

  if (is.null(TrtKeys)) {
    if ("Keys" %in% names(metadata(se))) {
      TrtKeys <- metadata(se)$Keys$Trt
      TrtKeys <- setdiff(TrtKeys, metadata(se)$Keys$discard_keys)
    } else {
      TrtKeys <- identify_keys(se)$Trt
    }
  }
  metadata(se)$Keys$Trt <- TrtKeys

  trt_fields <- c("GRvalue", "RelativeViability")
  ref_fields <- c("Day0Readout", "UntrtReadout", "RefGRvalue", "RefRelativeViability", "RefReadout", "DivisionTime")

  trt <- SummarizedExperiment::assay(se, treated_assay)
  avg_trt <- trt
  ref <- SummarizedExperiment::assay(se, reference_assay)
  avg_ref <- ref

  for (i in rownames(se)) {
    for (j in colnames(se)) {
      trt_df <- trt[i, j]
      ref_df <- ref[i, j]

      # bypass 'masked' filter
      masked <- trt_df$masked & !include_masked

      subKeys <- intersect(TrtKeys, colnames(trt_df))
      if (sum(!masked) >= 1) {
	df_av <- aggregate(trt_df[!masked, trt_fields],
	  by = as.list(trt_df[!masked, subKeys, drop = FALSE]), aggregate_FXN)
	df_std <- aggregate(trt_df[!masked, trt_fields],
	  by = as.list(trt_df[!masked, subKeys, drop = FALSE]), aggregate_FXN)
	colnames(df_std) <- paste0("std_", colnames(df_std))
	avg_trt_df <- merge(df_av, df_std, by = subKeys)
      } else { # case: (nrow(trt_df) == 0 || all(masked))
	avg_trt_df <- as.data.frame(matrix(0, 0, length(subKeys) + 5)) # TODO: remove this hardcoded number.
	colnames(avg_trt_df) <- c(subKeys, trt_fields, paste0("std_", trt_fields))
      } 
 
      if (nrow(ref_df) > 1) {
	subKeys <- intersect(TrtKeys, colnames(ref_df))
	avg_ref_df <- DataFrame(lapply(ref_df[, ref_fields], aggregate_FXN))
      } else {
        avg_ref_df <- ref_df
      }
      
      avg_trt[i, j] <- avg_trt_df
      avg_ref[i, j] <- avg_ref_df
    }
  }
  # TODO: Add the two BumpyMatrices back into the SE.
  # paste0("Averaged", )
  return(se)
}
