.drugNameRegex <- sprintf("^%s$|^%s_[[:digit:]]+$", 
                          gDRutils::get_identifier("drugname"), 
                          gDRutils::get_identifier("drugname"))

.untreated_tag_patterns <- vapply(gDRutils::get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
.untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")


#' .assign_treated_and_untreated_conditions
#'
#' Assign the treated and untreated conditions to a data.frame
#'
#' @param df_ data.frame containing the drug name field.
#' @param drugname_col string of the column in \code{df_} corresponding to the drug name. 
#' Defaults to \code{gDRutils::get_identifier("drugname")}..
#'
#' @return data.frame containing an additional factor column called \code{treated_untreated} specifying
#' whether the condition is a treated or untreated entry.
#'
#' @export
#'
.assign_treated_and_untreated_conditions <- function(df_, drugname_col = gDRutils::get_identifier("drugname")) {
  if (!drugname_col %in% colnames(df_)) {
    stop(sprintf("missing drug name column: %s", drugname_col))
  }
  drugnames <- tolower(as.data.frame(df_)[, drugname_col])
  untreated <- grepl(.untreatedDrugNameRegex, drugnames)
  if (!any(untreated)) {
    stop(sprintf("no untreated conditions matching pattern: '%s'", .untreatedDrugNameRegex))
  }

  df_$treated_untreated <- "treated"
  df_$treated_untreated[untreated] <- "untreated"
  df_$treated_untreated <- as.factor(df_$treated_untreated)
  df_
}
