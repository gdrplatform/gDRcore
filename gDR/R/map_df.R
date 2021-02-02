#' Map treated conditions to their respective references.
#'
#' Map treated conditions to their respective Day0, untreated, or single-agent references using condition metadata.
#'
#' @param trt_rdata data.frame of treated metadata. 
#' @param ref_rdata data.frame of untreated metadata.
#' @param row_endpoint_value_filter
#' @param Keys named list of keys and values.
#' Likely obtained from \code{identify_keys()}.
#' @param ref_type string of the reference type to map to.
#' Should be one of \code(c("Day0", "untrt_Endpoint", "ref_Endpoint")).
#'
#' @return named list where the names are the treated conditions and the values are the corresponding references for each condition.
#'
#' @export
#'
map_df <- function(trt_rdata, ref_rdata, row_endpoint_value_filter, Keys, ref_type = c("Day0", "untrt_Endpoint", "ref_Endpoint")) {
  # Assertions:
  checkmate::assert_class(normSE, "data.frame")
  checkmate::assert_class(ctrlSE, "data.frame")
  checkmate::assert_logical(row_endpoint_value_filter)
  checkmate::assert_list(Keys)

  duration_col <- gDRutils::get_identifier("duration")

  conc <- cbind(array(0, nrow(ref_rdata)), # padding to avoid empty df;
    ref_rdata[, agrep("Concentration", colnames(ref_rdata)), drop = FALSE])
  is_ref_conc <- apply(conc, 1, function(z) {all(z == 0)})

  if (ref_type == "Day0") {
    matching_list <- list(T0 = ref_rdata[, duration_col] == 0, conc = is_ref_conc) # Identifying which of the durations have a value of 0.
    matchFactor <- "T0"
  } else if (ref_type == "untrt_Endpoint"){
    matching_list <- list(key_values = row_endpoint_value_filter, conc = is_ref_conc)
    matchFactor <- duration_col 
  }

  trt_rnames <- rownames(trt_rdata)

  # define matrix with matching metadata
  ref_md <- intersect(Keys[[ref_type]], names(ref_rdata))
  names(ref_md) <- ref_md

  out <- vector("list", length(trt_rnames))
  for (i in seq_len(length(trt_rnames))) {
    treatment <- trt_rnames[i]
    refs <- lapply(ref_md, function(y) {ref_rdata[, y] == trt_rdata[treatment, y]}),

    match_mx <- as.matrix(c(refs, matching_list))
    match_idx <- which(apply(match_mx), 2, all) # test matching conditions
    if (length(match_idx) == 0) {
      # No exact match, try to find best match (as many metadata fields as possible).
      futile.logger::flog.warn("Missing untreated controls '%s' for: ('%s')", ref_type, treatment)
      idx <- apply(match_mx, 2, function(y) sum(y, na.rm = TRUE)) 
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
  names(out) <- trt_rnames
  out
}
