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
.assign_treated_and_untreated_conditions <- function(df_, drugname_col) {
#.assign_treated_and_untreated_conditions <- function(df_, drugname_col = gDRutils::get_identifier(c("drugname", "drugname2"))) {
# TODO: Replace once GDR-911 has been merged.
  valid <- intersect(drugname_col, colnames(df_))
  if (length(valid) == 0L) {
    stop(sprintf("missing drug name column(s): %s", paste0(valid, ", ")))
  }

  .untreated_tag_patterns <- vapply(gDRutils::get_identifier("untreated_tag"), sprintf, fmt = "^%s$", character(1))
  .untreatedDrugNameRegex <- paste(.untreated_tag_patterns, collapse = "|")

  validate <- as.data.frame(df_)[, valid, drop = FALSE]
  is_untreated <- data.frame(lapply(validate, function(x) grepl(.untreatedDrugNameRegex, tolower(x))))
  untreated <- rowSums(is_untreated) == length(valid)
  if (!any(untreated)) {
    stop(sprintf("no untreated conditions matching pattern: '%s' in columns: '%s'",
      .untreatedDrugNameRegex, paste0(valid, ", ")))
  }

  df_$treated_untreated <- "treated"
  df_$treated_untreated[untreated] <- "untreated"
  df_$treated_untreated <- as.factor(df_$treated_untreated)
  df_
}
