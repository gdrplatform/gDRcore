#' identify_keys
#'
#' Group columns from a data.frame that correspond to different 
#'
#' @param df_ a data.frame to identify keys for.
#' @param nested_keys character vector of keys to exclude from the returned list. 
#' The keys discarded should be identical to the keys in the third
#' dimension of the SummarizedExperiment.
#' Defaults to the \code{"Barcode"} and the \code{masked} identifier.
#' @param override_untrt_controls named list containing defining factors in the treatments.
#' Defaults to \code{NULL}.
#' @param identifiers named list containing all identifiers to use during processing.
#' By default, this value will be obtained by the environment.
#'
#' @return named list of key types and their corresponding key values. 
#'
#' @details This is most likely to be used for provenance tracking and will 
#' be placed on the SummarizedExperiment metadata for downstream analyses
#' to reference. 
#'
#' @seealso map_df, create_SE
#' @export
#'
identify_keys <- function(df_,  
                          nested_keys = NULL,
                          override_untrt_controls = NULL,
                          identifiers = gDRutils::get_env_identifiers()) {
  # Assertions:
  stopifnot(inherits(df_, c("data.frame", "DataFrame")))

  all_keys <- colnames(df_)
  dropped_nested_keys <- setdiff(nested_keys, all_keys)
  if (length(dropped_nested_keys) != 0L) {
    warning(sprintf("ignoring nested_keys input: '%s' which are not present in data.frame",
      paste0(dropped_nested_keys, collapse = ", ")))
    nested_keys <- intersect(nested_keys, all_keys)
  }

  dropped_override_untrt_controls <- setdiff(names(override_untrt_controls), all_keys)
  if (length(dropped_override_untrt_controls) != 0L) {
    warning(sprintf("ignoring override_untrt_controls input: '%s' which are not present in data.frame",
      paste0(dropped_override_untrt_controls, collapse = ", ")))
    override_untrt_controls <- intersect(override_untrt_controls, all_keys)
  }

  x <- c("Concentration", 
    identifiers$drug, 
    identifiers$drug_name,
    identifiers$drug_moa)
  pattern <- sprintf("%s|%s|%s|%s", x[1], x[2], x[3], x[4])
  pattern_keys <- grepl(pattern, all_keys)
  
  duration_col <- identifiers$duration

  keys <- list(Trt = setdiff(all_keys, c(nested_keys, override_untrt_controls)),
    ref_Endpoint = setdiff(all_keys, c(x, override_untrt_controls)),
    untrt_Endpoint = setdiff(all_keys[!pattern_keys], override_untrt_controls),
    Day0 = setdiff(all_keys[!pattern_keys], duration_col),
    nested_keys = nested_keys
  )

  keys <- lapply(keys, function(x) setdiff(x, c(gDRutils::get_header("raw_data"),
    gDRutils::get_header("normalized_results"), 
    "Template", 
    identifiers$well_position, 
    gDRutils::get_header("averaged_results"),
    gDRutils::get_header("metrics_results"), 
    identifiers$cellline_ref_div_time)))

  keys$masked_tag <- identifiers$masked_tag
  keys$cellline_name <- identifiers$cellline_name
  keys$cellline_ref_div_time <- identifiers$cellline_ref_div_time
  keys$duration <- duration_col 
  keys$untreated_tag <- identifiers$untreated_tag

  t0 <- df_[, duration_col] == 0
  # Remove keys where all values are NA.
  # TODO: Improve this.
  for (k in keys[["untrt_Endpoint"]]) {
    if (all(is.na(df_[, k]))) {
      keys <- lapply(keys, function(x) setdiff(x, k))
    }
    if (all(is.na(df_[t0, k]))) {
      keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
    }
  }

  return(keys)
}
