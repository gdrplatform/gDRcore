.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$", 
                          gDRutils::get_identifier("drugname"), 
                          gDRutils::get_identifier("drugname"))

.untreated_tag_patterns <- vapply(gDRutils::get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
.untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse="|")


#' .assign_treated_and_untreated_conditions
#'
#' Assign the treated and untreated conditions to a data.frame
#'
#' @param df_ containing the drug name field as specified by \code{gDRutils::get_identifier("drugname")}..
#'
#' @return data.frame containing an additional factor column called \code{treated_untreated} specifying
#' whether the condition is a treated or untreated entry.
#'
#' @export
#'
.assign_treated_and_untreated_conditions <- function(df_) {
  drugnames <- tolower(as.data.frame(df_)[, gDRutils::get_identifier("drugname")])
  untreated <- grepl(.untreatedDrugNameRegex, drugnames)
  if (!any(untreated)) {
    stop(sprintf("no untreated conditions matching pattern: '%s'", .untreatedDrugNameRegex))
  }

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
    .Deprecated(msg = "use '.assign_treated_and_untreated_conditions' instead")

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
    .Deprecated(msg = "use '.assign_treated_and_untreated_conditions' instead")

    # Assertions:
    stopifnot(any(inherits(drug_data, "data.frame"), inherits(drug_data, "DataFrame")))
    drugnames <- tolower(as.data.frame(drug_data)[, gDRutils::get_identifier("drugname")])
    drug_data[!grepl(.untreatedDrugNameRegex, drugnames), "name_"]
  }
