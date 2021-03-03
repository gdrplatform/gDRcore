#' Map treated conditions to their respective references.
#'
#' Map treated conditions to their respective Day0, untreated, or single-agent references using condition metadata.
#'
#' @param trt_md data.frame of treated metadata. 
#' @param ref_md data.frame of untreated metadata.
#' @param row_endpoint_value_filter boolean 
#' @param Keys named list of keys and values.
#' Likely obtained from \code{identify_keys()}.
#' @param ref_type string of the reference type to map to.
#' Should be one of \code{c("Day0", "untrt_Endpoint", "ref_Endpoint")}.
#'
#' @return named list mapping treated metadata to untreated metadata.
#'
#' @seealso identify_keys2
#' @export
#'
map_df <- function(trt_md, ref_md, row_endpoint_value_filter, Keys, ref_type = c("Day0", "untrt_Endpoint", "ref_Endpoint")) {
  # Assertions:
  checkmate::assert_class(trt_md, "data.frame")
  checkmate::assert_class(ref_md, "data.frame")
  checkmate::assert_logical(row_endpoint_value_filter)
  checkmate::assert_list(Keys)
  ref_type <- match.arg(ref_type)

  duration_col <- gDRutils::get_identifier("duration")

  conc <- cbind(array(0, nrow(ref_md)), # padding to avoid empty df;
    ref_md[, agrep("Concentration", colnames(ref_md)), drop = FALSE])
  is_ref_conc <- apply(conc, 1, function(z) {all(z == 0)})

  if (ref_type == "Day0") {
    matching_list <- list(T0 = ref_md[, duration_col] == 0, conc = is_ref_conc) # Identifying which of the durations have a value of 0.
    matchFactor <- "T0"
  } else if (ref_type == "untrt_Endpoint") {
    matching_list <- list(key_values = row_endpoint_value_filter, conc = is_ref_conc)
    matchFactor <- duration_col 
  } else if (ref_type == "ref_Endpoint") {
    matching_list <- NULL
    matchFactor <- NULL
  }

  trt_rnames <- rownames(trt_md)

  # define matrix with matching metadata
  # TODO: Is this necessary? I think all the keys should now only include keys that are present.
  present_ref_cols <- intersect(Keys[[ref_type]], names(ref_md))
  names(present_ref_cols) <- present_ref_cols

  out <- list("vector", length(trt_rnames))
  for (i in seq_len(length(trt_rnames))) {
    treatment <- trt_rnames[i]
    refs <- lapply(present_ref_cols, function(y) {ref_md[, y] == trt_md[treatment, y]})

    all_checks <- c(refs, matching_list)
    match_mx <- do.call("rbind", all_checks)
    rownames(match_mx) <- names(all_checks)
    match_idx <- which(apply(match_mx, 2, all)) # test matching conditions
    if (length(match_idx) == 0) {
      # No exact match, try to find best match (as many metadata fields as possible).
      futile.logger::flog.warn("Missing reference controls '%s' for: ('%s')", ref_type, treatment)
      idx <- apply(match_mx, 2, function(y) sum(y, na.rm = TRUE)) 
      # TODO: Sort this out so that it also takes the average in case multiple are found.
      idx <- idx * match_mx[matchFactor, ]

      if (any(idx > 0)) {
	match_idx <- which.max(idx)
	futile.logger::flog.warn("Found partial match:", rownames(ref_md)[match_idx])
      } else { # failed to find any potential match
	futile.logger::flog.warn("No partial match found")
      }
    }
    out[[i]] <- rownames(ref_md)[match_idx] # TODO: Check that this properly handles NAs. 
  }
  names(out) <- trt_rnames
  out
}
