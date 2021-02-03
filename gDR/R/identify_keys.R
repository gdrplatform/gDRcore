#' identify_keys
#'
#' Identify keys in the DR data represented by dataframe or SummarizedExperiment objects
#'
#' @param obj a dataframe or SummarizedExperiment with keys
#'
#' @return a list of keys
#' @export
#'
identify_keys <- function(obj) {
  # Assertions:
  stopifnot(inherits(obj, c("data.frame", "DataFrame", "SummarizedExperiment")))

  if (inherits(obj, "SummarizedExperiment")) {
    rdata <- SummarizedExperiment::rowData(obj)
    cdata <- SummarizedExperiment::colData(obj)

    se_untrt <- NULL
    all_keys <- unique(c(colnames(rdata), colnames(cdata),
      unlist(lapply(SummarizedExperiment::assay(obj), colnames))))
  } else { # case of a data frame
    all_keys <- colnames(obj)
  }

  keys <- list(Trt = setdiff(all_keys, "Barcode"),
    DoseResp = setdiff(all_keys,  "Barcode"),
    ref_Endpoint = setdiff(all_keys, c("Concentration",
      gDRutils::get_identifier("drug"),
      gDRutils::get_identifier("drugname"))),
    untrt_Endpoint = all_keys[ c(-agrep("Concentration", all_keys),
      -agrep(gDRutils::get_identifier("drug"), all_keys),
      -agrep(gDRutils::get_identifier("drugname"), all_keys))])
  keys[["Day0"]] <- setdiff(keys[["untrt_Endpoint"]], gDRutils::get_identifier("duration"))
  keys <- lapply(keys, function(x) setdiff(x, c(gDRutils::get_header("raw_data"),
    gDRutils::get_header("normalized_results"), "Template", 
    gDRutils::get_identifier("well_position"), gDRutils::get_header("averaged_results"),
	  gDRutils::get_header("metrics_results"), "ReferenceDivisionTime")))
  keys <- lapply(keys, sort)

  # check if all values of a key is NA
  for (k in keys[["untrt_Endpoint"]]) {
    if (inherits(obj, "SummarizedExperiment")) {
      # check the metadata fields for NA
      if (k %in% colnames(rdata)) {
	df_ <- rdata
      } else if (k %in% colnames(cdata)) {
	df_ <- cdata
      } else {
	next # not a metadata
      }

      if (all(is.na(df_[,k]))) {
        keys <- lapply(keys, function(x) setdiff(x, k))
      }

      if (!is.null(se_untrt) && k %in% colnames(SummarizedExperiment::rowData(se_untrt))) {
	df_ <- SummarizedExperiment::rowData(se_untrt)
	if (all(is.na(df_[df_[, gDRutils::get_identifier("duration")] == 0, k]))) {
	    keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
	}
      }
    } else { # case of a data frame
      if (all(is.na(obj[, k]))) {
	keys <- lapply(keys, function(x) setdiff(x, k))
      }
      if (all(is.na(obj[obj[,gDRutils::get_identifier("duration")] == 0, k]))) {
	keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
      }
    }
  }
  return(keys)
}
