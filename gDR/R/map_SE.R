#' map_SE
#'
#' Perfmorm mapping for normalization
#'
#' @param normSE a SummarizedExperiment with normalization assaay
#' @param ctrlSE a SummarizedExperiment object with information for controls
#' @param row_endpoint_value_filter an array with key values for end points
#' @param Keys a list of all identified keys
#' @param T0 a logical indicating if the mapping should be performer for Time=0 (FALSE by default)
#'
#' @return a list of mapping
#' @export
#'
map_SE <- function(normSE, ctrlSE, row_endpoint_value_filter, Keys, T0 = FALSE) {
  .Deprecated(msg = "see map_df for similar, but not identical functionality")

  # Assertions:
  checkmate::assert_class(normSE, "SummarizedExperiment")
  checkmate::assert_class(ctrlSE, "SummarizedExperiment")
  checkmate::assert_logical(row_endpoint_value_filter)
  checkmate::assert_list(Keys)
  checkmate::assert_logical(T0)

  duration_col <- gDRutils::get_identifier("duration")
  if (T0) {
    ctrl <- "Day0"
    keyValuesList <- list(T0 = SummarizedExperiment::rowData(ctrlSE)[, duration_col] == 0)
    matchFactor <- "T0"
    ref_type <- "T=0"
  } else {
    ctrl <- "untrt_Endpoint"
    keyValuesList <- list(key_values = row_endpoint_value_filter)
    matchFactor <- duration_col 
    ref_type <- "endpoint"
  }

  ctrl_rdata <- SummarizedExperiment::rowData(ctrlSE)
  norm_rdata <- SummarizedExperiment::rowData(normSE)
  norm_rnames <- rownames(normSE)

  out <- vector("list", length(norm_rnames))
  for (i in seq_len(length(norm_rnames))) {
    treatment <- norm_rnames[i]

    # define matrix with matching metadata
    ctrl_md <- intersect(Keys[[ctrl]], names(ctrl_rdata))
    names(ctrl_md) <- ctrl_md

    conc <- apply(cbind(array(0, nrow(ctrlSE)), # padding to avoid empty df;
			ctrl_rdata[, agrep("Concentration", colnames(ctrl_rdata)), drop = FALSE]), 1,
		  function(z) {all(z == 0)})

    match_mx <- IRanges::LogicalList(c(
      lapply(ctrl_md, function(y) ctrl_rdata[, y] == norm_rdata[treatment, y]),
      c(keyValuesList, list(conc = conc))))

    match_idx <- which(apply(as.matrix(match_mx), 2, all)) # test matching conditions
    if (length(match_idx) == 0) {
      # If no exact match, try to find best match (as many metadata fields as possible)
      futile.logger::flog.warn("Missing untreated controls '%s' for: ('%s')", ref_type, treatment)
      idx <- apply(as.matrix(match_mx), 2, function(y) sum(y, na.rm = TRUE)) 
      idx <- idx * match_mx[[matchFactor]]

      if (any(idx > 0)) {
	match_idx <- which.max(idx)
	futile.logger::flog.warn("Found partial match:", rownames(ctrlSE)[match_idx])
      } else { # failed to find any potential match
	futile.logger::flog.warn("No partial match found")
      }
    }

    out[[i]] <- rownames(ctrlSE)[match_idx]
  }
  names(out) <- norm_rnames
  out
}
