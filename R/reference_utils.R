#' Create a control dataframe for a treatment-cell line combination.
#'
#' Create an aggregated control data.frame.
#'
#' @param df_ data.frame
#' @param control_cols character vector of columns to include in the resulting 
#' control data.frame.
#' @param control_mean_fxn function indicating how to average controls.
#' Defaults to \code{mean(x, trim = 0.25)}.
#' @param out_col_name string of the output readout that will replace 
#' \code{CorrectedReadout}.
#'
#' @return data.frame of values aggregated by
#'
create_control_df <- function(df_,
                              control_cols,
                              control_mean_fxn = function(x) mean(x, trim = 0.25),
                              out_col_name) {

  checkmate::assert_data_frame(df_)
  checkmate::assert_character(control_cols)
  checkmate::assert_function(control_mean_fxn)
  checkmate::assert_string(out_col_name)

  if (nrow(df_) != 0L) {
    # Rename CorrectedReadout.
    col_intersect <- intersect(control_cols, colnames(df_))
    df_ <- df_[, c("CorrectedReadout", col_intersect), drop = FALSE]
    colnames(df_)[grepl("CorrectedReadout", colnames(df_))] <- out_col_name

    # Aggregate by all non-readout data (the metadata).
    if (ncol(df_) > 1) {
      df_ <- stats::aggregate(df_[, out_col_name, drop = FALSE],
              by = as.list(df_[, colnames(df_) != out_col_name, drop = FALSE]),
              function(x) control_mean_fxn(x))
    } else if (ncol(df_) == 1) {
      # only ReadoutValue column exists (i.e. no 'Barcode')
      df_ <- data.frame(control_mean_fxn(df_[, out_col_name]))
      colnames(df_) <- out_col_name
    } else {
      stop(sprintf("unexpected columns in data.frame: '%s'",
        paste0(colnames(df_), collapse = ", ")))
    }
  }
  df_
}
