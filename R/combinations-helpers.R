#' Get predicted values for a given fit and input.
#'
#' Map fittings to identifiers and compute the predicted values 
#' for corresponding fits.
#'
#' @param pred numeric vector for which you want predictions.
#' @param match_col vector to match on \code{fittings} to get the correct fit.
#' @param fittings data.table of fit metrics.
#' @param fitting_id_col string of the column name in \code{fittings} that 
#' should be used to match with \code{match_col} .
#'
#' @return 
#' Numeric vector of predicted values given \code{pred} inputs 
#' and \code{fittings} values.
#' 
#' @examples
#' pred <- c(1, 5, 5)
#' match_col <- c(1, 1, 2)
#' fitting_id_col <- "match_on_me"
#' 
#' fit1 <- data.table::data.table(h = 2.09, x_inf = 0.68, x_0 = 1, ec50 = 0.003)
#' fit2 <- data.table::data.table(h = 0.906, x_inf = 0.46, x_0 = 1, ec50 = 0.001)
#' fittings <- do.call(rbind, list(fit1, fit2))
#' fittings[[fitting_id_col]] <- c(1, 2)
#' 
#' map_ids_to_fits(pred, match_col, fittings, fitting_id_col)
#'
#' @keywords map_df
#' @export
map_ids_to_fits <- function(pred, match_col, fittings, fitting_id_col) {
  
  ridx <- S4Vectors::match(
    round(log10(match_col), 2), round(log10(fittings[[fitting_id_col]]), 2)
  )
  colnames <- c(fitting_id_col, "x_inf", "x_0", "ec50", "h")
  metrics <- fittings[ridx, colnames, with = FALSE]
  # Extrapolate fitted values.
  out <- gDRutils::predict_efficacy_from_conc(
    pred,
    metrics$x_inf,
    metrics$x_0,
    metrics$ec50,
    metrics$h
  )
  out
}


#' Calculate the difference between values in two data.tables
#'
#' Calculate the difference between values, likely representing 
#' the same metric, from two data.tables.
#'
#' @param metric data.table often representing 
#'               readouts derived by calculating some metric. 
#'               Examples of this could include hsa or bliss calculations 
#'               from single-agent data. 
#' @param measured data.table often representing 
#'                 measured data from an experiment.
#' @param series_identifiers character vector of identifiers in 
#'                           \code{measured} or \code{metric} 
#'                           which define a unique data point.
#' @param metric_col string of the column in \code{metric} 
#'                   to use in excess calculation.
#' @param measured_col string of the column in \code{measured} 
#'                     to use in excess calculation.
#'
#' @examples
#' metric <- data.table::data.table(
#'   Concentration = c(1, 2, 3, 1, 2, 3),
#'   Concentration_2 = c(1, 1, 1, 2, 2, 2),
#'   GRvalue = c(100, 200, 300, 400, 500, 600)
#' )
#' measured <- data.table::data.table(
#'   Concentration = c(3, 1, 2, 2, 1, 3),
#'   Concentration_2 = c(1, 1, 1, 2, 2, 2),
#'   testvalue = c(200, 0, 100, 400, 300, 500)
#' )
#' series_identifiers <- c("Concentration", "Concentration_2")
#' metric_col <- "GRvalue"
#' measured_col <- "testvalue"
#' calculate_excess(
#'   metric, 
#'   measured, 
#'   series_identifiers, 
#'   metric_col, 
#'   measured_col
#' )
#'
#' @return data.table of \code{measured}, now with an additional column named
#' \code{excess} (positive values for synergy/benefit).
#' @keywords combinations
#' @export
calculate_excess <- function(metric, 
                             measured, 
                             series_identifiers, 
                             metric_col, 
                             measured_col) {
  # TODO: Ensure same dims metric, measured
  # TODO: Ensure there is a unique entry for series_id1, series_id2 and 
  # no repeats
  # TODO: Check that there are no NAs
  idx <- S4Vectors::match(
    DataFrame(measured[, series_identifiers, with = FALSE]), 
    DataFrame(metric[, series_identifiers, with = FALSE])
  )
  
  out <- measured[, series_identifiers, with = FALSE]
  excess <- metric[idx, metric_col, with = FALSE] - measured[, measured_col, with = FALSE]
  excess[apply(as.matrix(out[, series_identifiers, with = FALSE]) == 0, 1, any)] <- NA
  out$x <- excess
  out
}

#' Calculate score for HSA and Bliss
#'
#' @param excess numeric vector with excess
#'
#' @examples
#' metric <- data.table::data.table(
#'   Concentration = c(1, 2, 3, 1, 2, 3),
#'   Concentration_2 = c(1, 1, 1, 2, 2, 2),
#'   GRvalue = c(100, 200, 300, 400, 500, 600)
#' )
#' measured <- data.table::data.table(
#'   Concentration = c(3, 1, 2, 2, 1, 3),
#'   Concentration_2 = c(1, 1, 1, 2, 2, 2),
#'   testvalue = c(200, 0, 100, 400, 300, 500)
#' )
#' series_identifiers <- c("Concentration", "Concentration_2")
#' metric_col <- "GRvalue"
#' measured_col <- "testvalue"
#' x <- calculate_excess(
#'   metric, 
#'   measured, 
#'   series_identifiers, 
#'   metric_col, 
#'   measured_col
#' )
#' calculate_score(x$x)
#'
#' @return numeric vector with calculated score
#' @keywords combinations
#' @export
calculate_score <- function(excess) {
  checkmate::assert_numeric(excess)
  # average the top 10-percentile excess to get a single value  for the excess
  mean(excess[excess >= stats::quantile(excess, 0.9, na.rm = TRUE)], na.rm = TRUE)
}


#' @keywords internal
#' @noRd
convertDFtoBumpyMatrixUsingIds <- function(df, 
                                           row_id = "row_id", 
                                           col_id = "col_id") {
  BumpyMatrix::splitAsBumpyMatrix(
    df[!colnames(df) %in% c(row_id, col_id)],
    row = df[[row_id]], col = df[[col_id]])
}


#' Standardize concentrations.
#'
#' Utilize a map to standardize concentrations.
#'
#' @param original_concs numeric vector of concentrations to replace 
#'                       using \code{conc_map}.
#' @param conc_map data.table of two columns named \code{original_conc_col} 
#'                 and \code{standardized_conc_col}.
#' @param original_conc_col string of the name of the column in \code{conc_map}
#'                          containing the original concentrations to replace.
#' @param standardized_conc_col string of the name of the column 
#'                              in \code{conc_map} containing the standardized 
#'                              concentrations to use for replacement.
#'
#' @examples
#' conc_map <- data.table::data.table(
#'   orig = c(0.99, 0.6, 0.456, 0.4), 
#'   std = c(1, 0.6, 0.46, 0.4)
#' )
#' original_concs <- c(0.456, 0.456, 0.4, 0.99)
#' exp <- c(0.46, 0.46, 0.4, 1)
#' obs <- replace_conc_with_standardized_conc(
#'   original_concs, 
#'   conc_map,
#'   original_conc_col = "orig", 
#'   standardized_conc_col = "std"
#' )
#' @return numeric vector of standardized concentrations.
#'
#' @seealso map_conc_to_standardized_conc
#' @keywords utils
#' @export
replace_conc_with_standardized_conc <- function(original_concs, 
                                                conc_map, 
                                                original_conc_col, 
                                                standardized_conc_col) {
  
  out <- unlist(conc_map[
    match(
      original_concs, 
      conc_map[[original_conc_col]]
    ), 
    standardized_conc_col, with = FALSE])
  if (length(out) != length(original_concs)) {
    stop("standardized output is not the same length as the input")
  }
  out
}


#' Create a mapping of concentrations to standardized concentrations.
#'
#' Create a mapping of concentrations to standardized concentrations.
#'
#' @param conc1 numeric vector of the concentrations for drug 1.
#' @param conc2 numeric vector of the concentrations for drug 2.
#'
#' @examples
#' 
#' ratio <- 0.5
#' conc1 <- c(0, 10 ^ (seq(-3, 1, ratio)))
#' 
#' shorter_range <- conc1[-1]
#' noise <- runif(length(shorter_range), 1e-12, 1e-11)
#' conc2 <- shorter_range + noise
#' 
#' map_conc_to_standardized_conc(conc1, conc2)
#'
#' @return data.table of 2 columns named \code{"concs"} and \code{"rconcs"}
#' containing the original concentrations and their closest matched 
#' standardized concentrations respectively. and their new standardized 
#' concentrations.
#'
#' @seealso replace_conc_w_standardized_conc
#' @details The concentrations are standardized in that they will contain 
#' regularly spaced dilutions and close values will be rounded.
#' @keywords utils
#' @export
map_conc_to_standardized_conc <- function(conc1, conc2) {
  # Remove single-agent.
  
  conc_1 <- setdiff(conc1, 0)
  conc_2 <- setdiff(conc2, 0)
  
  conc_1 <- sort(unique(conc_1))
  rconc1 <- .standardize_conc(conc_1)
  
  conc_2 <- sort(unique(conc_2))
  rconc2 <- .standardize_conc(conc_2)
  
  rconc <- c(0, unique(c(rconc1, rconc2)))
  .find_closest_match <- function(x) {
    rconc[abs(rconc - x) == min(abs(rconc - x))]
  }
  concs <- unique(c(conc1, conc2))
  mapped_rconcs <- vapply(concs, .find_closest_match, numeric(1))
  map <- unique(data.table::data.table(concs = concs, rconcs = mapped_rconcs))
  
  tol <- 1
  
  # Check if standardized values are within 5% of the original values
  round_diff <- which(abs(map$concs - map$rconcs) / map$concs > 0.05)
  
  map$rconcs[round_diff] <- map$concs[round_diff]
  mismatched <- which(
    gDRutils::round_concentration(map$conc, tol) != gDRutils::round_concentration(map$rconc, tol)
  )
  for (i in mismatched) {
    warning(sprintf("mapping original concentration '%s' to '%s'",
                    map[i, "concs"], map[i, "rconcs"]))
  }
  map
}


#' Standardize concentration values.
#'
#' Standardize concentration values.
#'
#' @param conc numeric vector of the concentrations
#'
#' @examples
#' concs <- 10 ^ (seq(-1, 1, 0.9))
#' .standardize_conc(concs)
#'
#' @return vector of standardized concentrations
#' @details If no \code{conc} are passed, \code{NULL} is returned.
#'
#' @keywords utils
#' @export
.standardize_conc <- function(conc) {
  rconc <- if (S4Vectors::isEmpty(conc)) {
    NULL
  } else if (length(unique(gDRutils::round_concentration(conc, 3))) > 4) {
    # 4 is determined by the fewest number of concentrations required to be 
    # considered a "matrix".
    log10_step <- .calculate_dilution_ratio(conc)
    num_steps <- round((log10(max(conc)) - log10(min(conc)) / log10_step), 0)
    seqs <- log10(max(conc)) - (log10_step * c(0:num_steps))
    sort(gDRutils::round_concentration(10 ^ seqs, 3))
  } else {
    # Few enough concentrations, don't need to calculate a series.
    gDRutils::round_concentration(conc, 3)
  }
  rconc
}



#' Calculate a dilution ratio.
#'
#' Calculation a dilution ratio from a vector of concentrations.
#'
#' @param concs numeric vector of concentrations.
#'
#' @return numeric value of the dilution ratio for a given set of 
#' concentrations.
#' @keywords internal
#' @noRd
.calculate_dilution_ratio <- function(concs) {
  first_removed <- concs[-1]
  first_two_removed <- first_removed[-1]
  last_removed <- concs[-length(concs)]
  last_two_removed <- last_removed[-length(last_removed)]
  
  dil_ratios <- c(
    log10(first_removed / last_removed), 
    log10(first_two_removed / last_two_removed)
  )
  rounded_dil_ratios <- gDRutils::round_concentration(dil_ratios, 2)
  
  # Get most frequent dilution ratio.
  highest_freq_ratio <- names(
    sort(table(rounded_dil_ratios), decreasing = TRUE)
  )[[1]]
  as.numeric(highest_freq_ratio)
}
