#' .assign_conditions
#'
#' Assign every row in a code{data.frame} to be a \code{"treated"}, 
#' \code{"untreated"}, or \code{"reference"} condition.
#'
#' @param df_ \code{data.frame} containing the \code{nested_identifiers}.
#' @param nested_identifiers Character vector of column names in \code{df_}
#'                           that will determine a condition type.
#'
#' @return data.frame containing an additional factor column called 
#' \code{"conditions"} specifying whether the condition is a \code{"treated"} 
#' (drug treated), \code{"untreated"}, or \code{"reference"} entry.
#' 
#' @details 
#' Conditions are considered \code{"untreated"} when both 
#' \code{nested_identifiers} have a value of \code{0}. 
#' When nested_identifiers are drugs, this represent entries where all drugs 
#' are at \code{0} concentration. Conditions are considered \code{"reference"} 
#' when only one of the \code{nested_identifiers} is \code{0}. 
#' In the context of drugs, this represents a single-agent response.
#'
#' @export
#'
.assign_conditions <- function(df_, nested_identifiers) {
  checkmate::assert_data_frame(df_)
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  valid <- intersect(nested_identifiers, colnames(df_))
  if (length(valid) == 0L) {
    stop(
      sprintf(
        "input 'df_' does not contain nested_identifiers: '%s'", 
        toString(valid)
      )
    )
  }

  nzero <- rowSums(data.frame(lapply(df_[valid], function(x) x == 0L)))
  untreated <- nzero == length(valid)
  treated <- nzero == 0L

  if (!any(untreated)) {
    stop(
      sprintf(
        "no untreated conditions [conc = 0] found in columns: '%s'", 
        toString(valid)
      )
    )
  }

  df_$conditions <- "reference"
  df_$conditions[untreated] <- "untreated"
  df_$conditions[treated] <- "treated"
  df_$conditions <- as.factor(df_$conditions)
  df_
}
