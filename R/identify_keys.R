#' identify_keys
#'
#' Group columns from a data.table that correspond to different 
#'
#' @param df_ a data.table to identify keys for.
#' @param nested_keys character vector of keys to exclude from the returned 
#' list. The keys discarded should be identical to the keys in the third
#' dimension of the SummarizedExperiment.
#' Defaults to the \code{"Barcode"} and the \code{masked} identifier.
#' @param override_untrt_controls named list containing defining factors in the 
#' treatments. Defaults to \code{NULL}.
#' @param identifiers named list containing all identifiers to use during 
#' processing. By default, this value will be obtained by the environment.
#' 
#' @examples 
#' n <- 64
#' md_df <- data.table::data.table(
#'   Gnumber = rep(c("vehicle", "untreated", paste0("G", seq(2))), each = 16), 
#'   DrugName = rep(c("vehicle", "untreated", paste0("GN", seq(2))), each = 16), 
#'   clid = paste0("C", rep_len(seq(4), n)),
#'   CellLineName = paste0("N", rep_len(seq(4), n)),
#'   replicates = rep_len(paste0("R", rep(seq(4), each = 4)), 64),
#'   drug_moa = "inhibitor",
#'   ReferenceDivisionTime = rep_len(c(120, 60), n),
#'   Tissue = "Lung",
#'   parental_identifier = "CL12345",
#'   Duration = 160
#' )
#' md_df <- unique(md_df)
#' ref <- md_df$Gnumber %in% c("vehicle", "untreated")
#' trt_df <- md_df[!ref, ]
#' identify_keys(trt_df)
#'
#' @return named list of key types and their corresponding key values. 
#'
#' @details This is most likely to be used for provenance tracking and will 
#' be placed on the SummarizedExperiment metadata for downstream analyses
#' to reference. 
#'
#' @seealso map_df, create_SE
#' @keywords identify_keys
#' @export
#'
identify_keys <- function(df_,  
                          nested_keys = NULL,
                          override_untrt_controls = NULL,
                          identifiers = gDRutils::get_env_identifiers()) {
  # Assertions:
  stopifnot(inherits(df_, c("data.table", "DataFrame")))
  
  all_keys <- colnames(df_)
  dropped_nested_keys <- setdiff(nested_keys, all_keys)
  if (length(dropped_nested_keys) != 0L) {
    warning(
      sprintf(
        "ignoring nested_keys input: '%s' which are not present in data.table",
        toString(dQuote(dropped_nested_keys, q = FALSE))
      )
    )
    nested_keys <- intersect(nested_keys, all_keys)
  }

  dropped_override_untrt_controls <- 
    setdiff(names(override_untrt_controls), all_keys)
  if (length(dropped_override_untrt_controls) != 0L) {
    warning(
      sprintf(
        "ignoring override_untrt_controls input: '%s' which are 
        not present in data.table",
        toString(dQuote(dropped_override_untrt_controls, q = FALSE))
      )
    )
    override_untrt_controls <- intersect(override_untrt_controls, all_keys)
  }
  
  x <- c("concentration", "drug", "drug_name", "drug_moa")
  pattern_keys <- all_keys %in% 
    identifiers[grep(paste(x, collapse = "|"), names(identifiers))]
  
  duration_col <- identifiers$duration

  keys <- list(Trt = setdiff(all_keys, c(nested_keys, override_untrt_controls)),
    ref_Endpoint = setdiff(
      all_keys, 
      c(identifiers[x], override_untrt_controls)
    ),
    untrt_Endpoint = setdiff(all_keys[!pattern_keys], override_untrt_controls),
    Day0 = setdiff(all_keys[!pattern_keys], duration_col),
    nested_keys = nested_keys
  )

  keys <- gDRutils::loop(
    keys, 
    function(x) {
      setdiff(x, c(gDRutils::get_header("raw_data"),
      gDRutils::get_header("normalized_results"), 
      identifiers$template, 
      identifiers$well_position, 
      gDRutils::get_header("averaged_results"),
      gDRutils::get_header("metrics_results"), 
      identifiers$cellline_ref_div_time))
    })

  keys$masked_tag <- identifiers$masked_tag
  keys$cellline_name <- identifiers$cellline_name
  keys$cellline_ref_div_time <- identifiers$cellline_ref_div_time
  keys$duration <- duration_col 
  keys$untreated_tag <- identifiers$untreated_tag

  t0 <- df_[, duration_col, with = FALSE] == 0
  # Remove keys where all values are NA.
  # TODO: Improve this.
  for (k in keys[["untrt_Endpoint"]]) {
    if (all(is.na(df_[, k, with = FALSE]))) {
      keys <- gDRutils::loop(keys, function(x) setdiff(x, k))
    }
    if (all(is.na(df_[which(t0), k, with = FALSE]))) {
      keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
    }
  }

  return(keys)
}
