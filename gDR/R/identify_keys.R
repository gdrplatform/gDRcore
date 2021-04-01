#' identify_keys
#'
#' Identify keys in the DR data represented by dataframe or SummarizedExperiment or MultiAssayExperiment objects
#'
#' @param df_se_mae a dataframe or SummarizedExperiment or MultiassayExperiment with keys
#'
#' @return a list of keys
#' @export
#'
identify_keys <- function(df_se_mae) {
  .Deprecated(msg = "see identify_keys2 for similar, but not identical functionality")

  # Assertions:
  stopifnot(inherits(df_se_mae, c("data.frame", "MultiAssayExperiment", "SummarizedExperiment")))


    if (any(class(df_se_mae) %in% c("MultiAssayExperiment", "SummarizedExperiment"))) {
        if ("MultiAssayExperiment" %in% class(df_se_mae)) {
            # if MAE, convert to SE based on the treated SE (could be optimized)
            df_se_mae <- df_se_mae[["treated"]]
            se_untrt <-  df_se_mae[["untreated"]]
        } else se_untrt <- NULL
        all_keys <- unique(c(
            colnames(SummarizedExperiment::rowData(df_se_mae)),
            colnames(SummarizedExperiment::colData(df_se_mae)),
            unlist(lapply(SummarizedExperiment::assay(df_se_mae), colnames))))
    } else { # case of a data frame
        all_keys <- colnames(df_se_mae)
    }

    keys <- list(Trt = setdiff(all_keys, "Barcode"),
            DoseResp = setdiff(all_keys,  "Barcode"),
            ref_Endpoint = setdiff(all_keys, c("Concentration",
                                            gDRutils::get_identifier("drug"),
                                            gDRutils::get_identifier("drugname"),
                                            gDRutils::get_identifier("drug_moa"))),
            untrt_Endpoint = all_keys[ c(-agrep("Concentration", all_keys),
                                            -agrep(gDRutils::get_identifier("drug"), all_keys),
                                            -agrep(gDRutils::get_identifier("drugname"), all_keys),
                                            -agrep(gDRutils::get_identifier("drug_moa"), all_keys))])
    keys[["Day0"]] <- setdiff(keys[["untrt_Endpoint"]], gDRutils::get_identifier("duration"))
    keys <- lapply(keys, function(x) setdiff(x, c(gDRutils::get_header("raw_data"),
        gDRutils::get_header("normalized_results"), "Template", gDRutils::get_identifier("WellPosition"), gDRutils::get_header("averaged_results"),
            gDRutils::get_header("metrics_results"), "ReferenceDivisionTime"
    )))
    keys <- lapply(keys, sort)

    # check if all values of a key is NA
    for (k in keys[["untrt_Endpoint"]]) {

        if ("SummarizedExperiment" %in% class(df_se_mae)) {
            # check the metadata fields for NA
            if (k %in% colnames(SummarizedExperiment::rowData(df_se_mae))) df_ <- SummarizedExperiment::rowData(df_se_mae)
            else if (k %in% colnames(SummarizedExperiment::colData(df_se_mae))) df_ <- SummarizedExperiment::colData(df_se_mae)
            else next # not a metadata

            if (all(is.na(df_[,k]))) keys <- lapply(keys, function(x) setdiff(x, k))

            if (!is.null(se_untrt) && k %in% colnames(SummarizedExperiment::rowData(se_untrt))) {
                df_ <- SummarizedExperiment::rowData(se_untrt)
                if (all(is.na(df_[df_[,gDRutils::get_identifier("duration")] == 0, k]))) {
                    keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
                }
            }
        } else { # case of a data frame
            if (all(is.na(df_se_mae[, k]))) {
                keys <- lapply(keys, function(x) setdiff(x, k))
            }
            if (all(is.na(df_se_mae[df_se_mae[,gDRutils::get_identifier("duration")] == 0, k]))) {
                keys[["Day0"]] <- setdiff(keys[["Day0"]], k)
            }
        }
    }
  return(keys)
}


#' identify_keys2
#'
#' Group columns from a data.frame that correspond to different 
#'
#' @param df_ a data.frame to identify keys for.
#' @param nested_keys character vector of keys to exclude from the returned list. 
#' The keys discarded should be identical to the keys in the third
#' dimension of the SummarizedExperiment.
#' Defaults to the \code{"Barcode"} and the \code{masked} identifier.
#' @param override_controls named list containing defining factors in the treatments.
#' Defaults to \code{NULL}.
#'
#' @return named list of key types and their corresponding key values. 
#'
#' @details This is most likely to be used for provenance tracking and will 
#' be placed on the SummarizedExperiment metadata for downstream analyses
#' to reference. 
#'
#' @seealso map_df, create_SE2
#' @export
#'
identify_keys2 <- function(df_,  
                           nested_keys = c("Barcode", gDRutils::get_identifier("masked_tag")),
                           override_controls = NULL) {
  # Assertions:
  stopifnot(inherits(df_, c("data.frame", "DataFrame")))

  all_keys <- colnames(df_)

  x <- c("Concentration", 
    gDRutils::get_identifier("drug"), 
    gDRutils::get_identifier("drugname"),
    gDRutils::get_identifier("drug_moa"))
  pattern <- sprintf("%s|%s|%s|%s", x[1], x[2], x[3], x[4])
  pattern_keys <- grepl(pattern, all_keys)
  
  duration_col <- gDRutils::get_identifier("duration")

  keys <- list(Trt = setdiff(all_keys, c(nested_keys, override_controls)),
    ref_Endpoint = setdiff(all_keys, c(x, override_controls)),
    untrt_Endpoint = setdiff(all_keys[!pattern_keys], override_controls),
    Day0 = setdiff(all_keys[!pattern_keys], duration_col),
    nested_keys = intersect(nested_keys, colnames(df_))
  )

  keys <- lapply(keys, function(x) setdiff(x, c(gDRutils::get_header("raw_data"),
    gDRutils::get_header("normalized_results"), 
    "Template", 
    gDRutils::get_identifier("well_position"), 
    gDRutils::get_header("averaged_results"),
    gDRutils::get_header("metrics_results"), 
    gDRutils::get_identifier("cellline_ref_div_time"))))

  keys$masked_tag <- gDRutils::get_identifier("masked_tag")
  keys$cellline_name <- gDRutils::get_identifier("cellline_name")
  keys$cellline_ref_div_time <- gDRutils::get_identifier("cellline_ref_div_time")
  keys$duration <- duration_col 
  keys$untreated_tag <- gDRutils::get_identifier("untreated_tag")

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
