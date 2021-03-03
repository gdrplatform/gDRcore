#' @import S4Vectors

#### AUXILIARY FUNCTIONS ####
.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$", 
                          gDRutils::get_identifier("drugname"), 
                          gDRutils::get_identifier("drugname"))

.untreated_tag_patterns <- vapply(gDRutils::get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
.untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse="|")

.assayNames <-
  c("df_raw_data",
    "df_normalized",
    "df_averaged",
    "df_metrics")


#' aapply
#'
#' Works like sapply but on each nested dataframe of the assay of an SE
#'
#' @param SE a SummarizedExperiment object
#' @param fx any function
#' @param assay_type a name of an assay
#'
#' @return the same SE object with updated nested dataframe
#' @export
aapply <-
  function(SE, fx, assay_type = 1) {
    # Assertions:
    checkmate::assert_class(SE, "SummarizedExperiment")
    checkmate::assert_function(fx)
    checkmate::assert_scalar(assay_type)

    SummarizedExperiment::assay(SE, assay_type, withDimnames=FALSE) = matrix(sapply(SummarizedExperiment::assay(SE, assay_type), fx), nrow = nrow(SE), ncol = ncol(SE))
    return(SE)
  }
