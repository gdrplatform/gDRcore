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


#' .assign_treated_and_untreated_conditions
#'
#' Assign the treated and untreated conditions to a data.frame
#'
#' @param data.frame containing the drug name field.
#'
#' @return data.frame containing an additional factor column called \code{treated_untreated} specifying
#' whether the condition is a treated or untreated entry.
#'
#' @export
#'
.assign_treated_and_untreated_conditions <- function(df_) {
  drugnames <- tolower(as.data.frame(df_)[, gDRutils::get_identifier("drugname")])
  untreated <- grepl(.untreatedDrugNameRegex, drugnames)
  df_$treated_untreated <- "treated"
  df_$treated_untreated[untreated] <- "untreated"
  df_$treated_untreated <- as.factor(df_$treated_untreated)
  df_
}


#' .get_untreated_conditions
#'
#' Get untreated conditions
#'
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with untreated conditions
#'
#' @export
#'
.get_untreated_conditions <-
  function(drug_data) {
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    drugnames <- tolower(as.data.frame(drug_data)[, gDRutils::get_identifier("drugname")])
    drug_data[grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }

#' .get_treated_conditions
#'
#' Get treated conditions
#'
#' @param drug_data data.frame or DataFrame with treatment information
#'
#' @return character vector with treated conditions
#'
#' @export
#'
.get_treated_conditions <-
  function(drug_data) {
    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    drugnames <- tolower(as.data.frame(drug_data)[, gDRutils::get_identifier("drugname")])
    drug_data[!grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }


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
