#' Get predicted values for a given fit and input.
#'
#' Map fittings to identifiers and compute the predicted values for corresponding fits.
#'
#' @param pred numeric vector for which you want predictions.
#' @param match_col vector to match on \code{fittings} to get the correct fit.
#' @param fittings data.frame of fit metrics.
#' @param fitting_id_col string of the column name in \code{fittings} that should be
#' used to match with \code{match_col} .
#'
#' @return numeric vector of predicted values given \code{pred} inputs and \code{fittings} values.
#'
#' @export
map_ids_to_fits <- function(pred, match_col, fittings, fitting_id_col) {
  ridx <- S4Vectors::match(round(log10(match_col), 2), round(log10(fittings[[fitting_id_col]]), 2))
  metrics <- fittings[ridx, c(fitting_id_col, "x_inf", "x_0", "ec50", "h")]
  # Extrapolate fitted values.
  out <- gDRutils::logistic_4parameters(pred,
    metrics$x_inf,
    metrics$x_0,
    metrics$ec50,
    metrics$h)
  out
}


#' Calculate the difference between values in two data.frames
#'
#' Calculate the difference between values, likely representing the same metric, from two data.frames.
#'
#' @param metric data.frame often representing readouts derived by calculating some metric.
#' Examples of this could include hsa or bliss calculations from single-agent data. 
#' @param measured data.frame often representing measured data from an experiment.
#' @param series_identifiers character vector of identifiers in \code{measured} or \code{metric}
#' which define a unique data point.
#' @param metric_col string of the column in \code{metric} to use in the excess calculation.
#' @param measured_col string of the column in \code{measured} to use in the excess calculation.
#'
#' @return DataFrame of \code{measured}, now with an additional column named \code{excess} (positive values for synergy/benefit).
#' @export
calculate_excess <- function(metric, measured, series_identifiers, metric_col, measured_col) {
  # TODO: Ensure same dims metric, measured
  # TODO: Ensure there is a unique entry for series_id1, series_id2 and no repeats
  # TODO: Check that there are no NAs
  idx <- S4Vectors::match(DataFrame(measured[, series_identifiers]), DataFrame(metric[, series_identifiers]))
  
  out <- measured[, series_identifiers]
  excess <- metric[idx, metric_col] - measured[, measured_col]
  excess[ apply(as.matrix(out[, series_identifiers])==0, 1, any) ] <- 0 # single agent should be 0
  out$excess <- excess
  as.data.frame(out)
}


#' @keywords internal
#' @noRd
convertDFtoBumpyMatrixUsingIds <- function(df, row_id = "row_id", col_id = "col_id") {
  BumpyMatrix::splitAsBumpyMatrix(
    df[!colnames(df) %in% c(row_id, col_id)],
    row = df[[row_id]], col = df[[col_id]])
}
