#' .assign_conditions
#'
#' Assign the treated, untreated, and reference conditions to a data.frame
#'
#' @param df_ data.frame containing the drug name field.
#' @param nested_identifiers Character vector of column names in \code{df_}
#' that 
#'
#' @return data.frame containing an additional factor column called \code{condition} specifying
#' whether the condition is a treated, untreated, or reference entry.
#' 
#' @details Untreated rows are identified by having both nested_identifiers set to 0.
#' When nested_identifiers are drugs, this represent entries where all drugs are at 0 concentration.
#'
#' @export
#'
.assign_conditions <- function(df_, nested_identifiers) {
  valid <- intersect(nested_identifiers, colnames(df_))
  if (length(valid) == 0L) {
    stop(sprintf("input 'df_' does not contain nested_identifiers: '%s'", paste0(valid, ", ")))
  }

  nzero <- rowSums(data.frame(lapply(df_[valid], function(x) x == 0L)))
  untreated <- nzero == length(valid)
  treated <- nzero == 0L

  if (!any(untreated)) {
    stop(sprintf("no untreated conditions [conc = 0] found in columns: '%s'",
      paste0(valid, ", ")))
  }

  df_$conditions <- "reference"
  df_$conditions[untreated] <- "untreated"
  df_$conditions[treated] <- "treated"
  df_$conditions <- as.factor(df_$conditions)
  df_
}
