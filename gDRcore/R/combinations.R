#' @keywords internal
.create_combo_control <- function(nested_identifiers) {
  out <- data.frame()
  for (i in nested_identifiers) {
    out[[i]] <- 0
  }
  out$GRvalue <- 1
  out$std_GRvalue <- 0
  out$RelativeViability <- 1
  out$std_RelativeViability <- 0
  out
}

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


#' @details HSA takes the minimum of the two single agents readouts.
calculate_HSA <- function(sa1, series_id1, sa2, series_id2, metric) {
  .calculate_matrix_metric(sa1, series_id1, sa2, series_id2, metric, FXN = pmin)
}


#' @details Bliss is the mulitplication of the single agent readouts.
calculate_Bliss <- function(sa1, series_id1, sa2, series_id2, metric) {
  .calculate_matrix_metric(sa1, series_id1, sa2, series_id2, metric, FXN = "*")
}


#' @keywords internal
.calculate_matrix_metric <- function(sa1, series_id1, sa2, series_id2, metric, FXN) {
  colnames(sa1)[colnames(sa1) == metric] <- "metric1"
  colnames(sa2)[colnames(sa2) == metric] <- "metric2"

  # TODO: ensure they're unique?
  u <- expand.grid(sa1[[series_id1]], sa2[[series_id2]])
  colnames(u) <- c(series_id1, series_id2)

  idx <- match(u[[series_id1]], sa1[[series_id1]])
  u <- base::merge(u, sa1[, c(series_id1, "metric1")], by = series_id1)
  u <- base::merge(u, sa2[, c(series_id2, "metric2")], by = series_id2)

  metric <- do.call(FXN, list(u$metric1, u$metric2))
  cbind(u, metric)
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
#' @return DataFrame of \code{measured}, now with an additional column named \code{excess}.
#' @export
calculate_excess <- function(metric, measured, series_identifiers, metric_col, measured_col) {
  # TODO: Ensure same dims metric, measured
  # TODO: Ensure there is a unique entry for series_id1, series_id2 and no repeats
  # TODO: Check that there are no NAs
  idx <- S4Vectors::match(DataFrame(measured[, series_identifiers]), DataFrame(metric[, series_identifiers]))
  
  out <- measured[, series_identifiers]
  excess <- measured[, measured_col] - metric[idx, metric_col]
  out$excess <- excess
  out
}


define_matrix_position <- function(mean_matrix,
                      conc_margin = 10 ^ 0.5,
                      log2_pos_offset = log10(3) / 2
              ) {

  checkmate::assert_number(conc_margin)
  checkmate::assert_number(log2_pos_offset)

  axis_1 <- data.frame(conc_1 = round(as.numeric(rownames(mean_matrix)), 5),
          log10conc_1 = 0, pos_y = 0, marks_y = 0)
  axis_1$log10conc_1 <- log10(axis_1$conc)
  axis_1$pos_y <- axis_1$log10conc_1
  axis_1$pos_y[1] <- 2 * axis_1$pos_y[2] - axis_1$pos_y[3] - log10(1.5)
  axis_1$marks_y <- sprintf("%.2g", axis_1$conc_1)

  # drug_2 is diluted along the columns and will be the x-axis of the matrix plots
  axis_2 <- data.frame(conc_2 = round(as.numeric(colnames(mean_matrix)), 5),
                log10conc_2 = 0, pos_x = 0, marks_x = 0)
  axis_2$log10conc_2 <- log10(axis_2$conc_2)
  axis_2$pos_x <- axis_2$log10conc_2
  axis_2$pos_x[1] <- 2 * axis_2$pos_x[2] - axis_2$pos_x[3] - log10(1.5)
  axis_2$marks_x <- sprintf("%.2g", axis_2$conc_2)

  list(axis_1 = axis_1, axis_2 = axis_2)
}



calculate_Loewe <- function(mean_matrix, 
                      row_fittings, 
                      col_fittings, 
                      codilution_fittings, 
                      series_identifiers,
                      normalization_type,
                      conc_margin = 10 ^ 0.5,
                      log2_pos_offset = log10(3) / 2
                    ) {

  checkmate::assert_number(conc_margin)
  checkmate::assert_number(log2_pos_offset)
  checkmate::assert_character(normalization_type) 

  axes <- define_matrix_position(mean_matrix, conc_margin = conc_margin, log2_pos_offset = log2_pos_offset)
  axis_1 <- axes$axis_1
  axis_2 <- axes$axis_2

  if (min(mean_matrix, na.rm = TRUE) > 0.7) {
    iso_cutoff <- NULL
  } else {
    if (normalization_type == "GR") {
      iso_cutoff <- seq(max(-0.25, ceiling(20 * min(mean_matrix + 0.08, na.rm = TRUE)) / 20), 0.8,  0.05)
    } else {
      iso_cutoff <- seq(max(0.2, ceiling(20 * min(mean_matrix + 0.08, na.rm = TRUE)) / 20), 0.8,  0.05)
    }
    names(iso_cutoff) <- as.character(iso_cutoff)
  }

  # create the variable for the different isobolograms
  all_iso <- vector("list", length(iso_cutoff))
  
  # TOOD: by_col and by_row should become encoded by the IRanges object length

  # get the IC50/GR50 on the marginal (single agent) based on the fit functions and capping
  ref_x50 <- c(conc_1 = min(col_fittings[col_fittings$cotrt_value == 0, "xc50"], max(axis_1$conc_1) * conc_margin),
              conc_2 = min(row_fittings[1, "xc50"], max(axis_2$conc_2) * conc_margin))

  names(all_iso) <- iso_cutoff
  for (isobol_value in iso_cutoff) { # run through the different isobolograms
    # cutoff point by row
    row_fittings <- row_fittings[order(row_fittings$cotrt_value, decreasing = TRUE), ]
    df_iso <- cbind(conc_1 = row_fittings[, "cotrt_value"], data.frame(conc_2 =
      ifelse(row_fittings$x_0 < isobol_value, 0, ifelse(row_fittings$x_inf > isobol_value,
        Inf,
        row_fittings$ec50 * ((((row_fittings$x_0 - row_fittings$x_inf) / (isobol_value - row_fittings$x_inf)) - 1) ^
            (1 / pmax(row_fittings$h, 0.01))))
          ),
          fit_type = "by_row"))

    # cutoff point by columns
    col_fittings <- col_fittings[order(row_fittings$cotrt_value, decreasing = FALSE), ]
    df_iso <- rbind(df_iso, cbind(data.frame(conc_1 =
      ifelse(col_fittings$x_0 < isobol_value, 0, ifelse(col_fittings$x_inf > isobol_value,
        Inf,
        col_fittings$ec50 * ((((col_fittings$x_0 - col_fittings$x_inf) / (isobol_value - col_fittings$x_inf)) - 1) ^
            (1 / pmax(col_fittings$h, 0.01))))
          ),
      conc_2 = col_fittings[, "cotrt_value"],
        fit_type = "by_col")))

    # capping
    df_iso$conc_1 <- pmin(df_iso$conc_1, max(axis_1$conc_1) * conc_margin)
    df_iso$conc_2 <- pmin(df_iso$conc_2, max(axis_2$conc_2) * conc_margin)

    ref_conc_1 <- pmin(df_iso$conc_1[df_iso$conc_2 == 0 & df_iso$fit_type == "by_col"],
                    max(axis_1$conc_1) * conc_margin)
    ref_conc_2 <- pmin(df_iso$conc_2[df_iso$conc_1 == 0 & df_iso$fit_type == "by_row"],
                    max(axis_2$conc_2) * conc_margin)

    # cutoff point by diagonal (co-dilution)
    # co-dil is given as concentration of drug 1
    if (nrow(codilution_fittings) > 1) {
      codilution_fittings <- codilution_fittings[order(codilution_fittings$ratio, decreasing = TRUE), ]
      codilution_fittings <- codilution_fittings[codilution_fittings$fit_type %in% "DRC3pHillFitModelFixS0", ]
      conc_mix <- ifelse (codilution_fittings$x_0 < isobol_value, 
          0,
          ifelse(codilution_fittings$x_inf > isobol_value,
            Inf,
            codilution_fittings$ec50 * ((((codilution_fittings$x_0 - codilution_fittings$x_inf) /
                                        (isobol_value - codilution_fittings$x_inf)) - 1) ^
                                      (1 / codilution_fittings$h))
          )
        )
        
      df_iso_codil <- data.frame(conc_1 = conc_mix / (1 + 1 / codilution_fittings$ratio),
                        conc_2 = conc_mix / (1 + codilution_fittings$ratio), fit_type = "by_codil")
      # avoid extrapolation
      capped_idx <- df_iso_codil$conc_1 > (max(axis_1$conc_1) * conc_margin) |
            df_iso_codil$conc_2 > (max(axis_2$conc_2) * conc_margin)
      df_iso_codil$conc_1[capped_idx] <- pmin(max(axis_1$conc_1) * conc_margin,
                  codilution_fittings$conc_ratio * (max(axis_2$conc_2) * conc_margin))[capped_idx]
      df_iso_codil$conc_2[capped_idx] <- pmin(max(axis_2$conc_2) * conc_margin,
                  (max(axis_1$conc_1) * conc_margin) / codilution_fittings$conc_ratio)[capped_idx]

      df_iso <- rbind(df_iso, df_iso_codil)
    }
    # remove low concentration values
    df_iso <- df_iso[!is.na(df_iso$conc_1) & !is.na(df_iso$conc_2), ]
    df_iso <- df_iso[(df_iso$conc_1 > axis_1$conc_1[2] / 2 | df_iso$fit_type == "by_row") &
                       (df_iso$conc_2 > axis_2$conc_2[2] / 2 | df_iso$fit_type == "by_col"), ]

    if (nrow(df_iso) < 5) {
      next
    }

    # cap the values for plotting and smoothing (pos_ is in the log10 space)
    df_iso$pos_x <- pmax(log10(df_iso$conc_2), min(axis_2$pos_x))
    df_iso$pos_y <- pmax(log10(df_iso$conc_1), min(axis_1$pos_y))

    df_iso$pos_x <- pmin(df_iso$pos_x, max(axis_2$pos_x) + log10(conc_margin))
    df_iso$pos_y <- pmin(df_iso$pos_y, max(axis_1$pos_y) + log10(conc_margin))

    # rotate 45 degree to calculate smooth curve (x_ is in the rotated log10 space):
    df_iso$x1 <- (df_iso$pos_x - min(axis_2$pos_x) -
                              (df_iso$pos_y - min(axis_1$pos_y))) / sqrt(2) # conc_ratio
    df_iso$x2 <- (df_iso$pos_x - min(axis_2$pos_x) +
                              (df_iso$pos_y - min(axis_1$pos_y))) / sqrt(2) # new response value
    # offset helps with smoothing of the edges (x2_off is in the rotated log10 space with offset):
    x2_extra_offset <- 1 / 4 
    df_iso$x2_off <- df_iso$x2 + abs(df_iso$x1) * x2_extra_offset
    isobol_x1 <- seq(min(df_iso$x1) - log2_pos_offset, max(df_iso$x1) + log2_pos_offset, .1)

    # perform the smoothing
    df_iso_curve <- data.frame(x1 = isobol_x1, x2_off = zoo::rollmean(
      rowMeans(do.call(cbind, lapply(names(which(table(df_iso$fit_type) > 1)), function(x)
        stats::approx(x = df_iso$x1[df_iso$fit_type == x], 
               y = df_iso$x2_off[df_iso$fit_type == x], xout = isobol_x1)$y)),
        na.rm = TRUE), 5, fill = NA))
    df_iso_curve <- df_iso_curve[!is.na(df_iso_curve$x2_off), ]
    df_iso_curve$x2 <- df_iso_curve$x2_off - abs(df_iso_curve$x1) * x2_extra_offset # remove offset

    # rotate back the position
    df_iso_curve$pos_x <- (df_iso_curve$x1 + df_iso_curve$x2) / sqrt(2) + min(axis_2$pos_x)
    df_iso_curve$pos_y <- (-df_iso_curve$x1 + df_iso_curve$x2) / sqrt(2) + min(axis_1$pos_y)

    # calculate the reference (additive model in the rotated space)
    c2 <- ref_conc_2 / (1 + (ref_conc_2 / ref_conc_1) * (10 ^ (
                -(sqrt(2) * df_iso_curve$x1 + min(axis_2$pos_x) - min(axis_1$pos_y)))))
    df_iso_curve$x2_ref <- (log10(c2) - min(axis_2$pos_x) +
              (log10(ref_conc_1 * (1 - c2 / ref_conc_2)) - min(axis_1$pos_y))) / sqrt(2)

    # cap the concentrations for the reference
    over_edge <- pmax(0, (-df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
                  min(axis_1$pos_y) -  (max(axis_1$pos_y) + conc_margin)) +
                pmax(0, (df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
                  min(axis_2$pos_x) -  (max(axis_2$pos_x) + conc_margin))
    df_iso_curve$x2_ref_cap <- df_iso_curve$x2_ref - over_edge * sqrt(2)

    # rotate back the reference
    df_iso_curve$pos_x_ref <- (df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
      min(axis_2$pos_x)
    df_iso_curve$pos_y_ref <- (-df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
      min(axis_1$pos_y)
    df_iso_curve <- rbind(NA, df_iso_curve, NA)
    df_iso_curve[1, c("pos_x", "pos_x_ref")] <- min(axis_2$pos_x)
    df_iso_curve[1, c("pos_y", "pos_y_ref")] <- log10(ref_conc_1)
    df_iso_curve[1, "x1"] <-
                (df_iso_curve$pos_x[1] - min(axis_2$pos_x) -
                    (df_iso_curve$pos_y[1] - min(axis_1$pos_y))) / sqrt(2)
    df_iso_curve[1, c("x2", "x2_ref", "x2_ref_cap")] <-
                (df_iso_curve$pos_x[1] - min(axis_2$pos_x) +
                    (df_iso_curve$pos_y[1] - min(axis_1$pos_y))) / sqrt(2)

    df_iso_curve[nrow(df_iso_curve), c("pos_x", "pos_x_ref")] <- log10(ref_conc_2)
    df_iso_curve[nrow(df_iso_curve), c("pos_y", "pos_y_ref")] <- min(axis_1$pos_y)
    df_iso_curve[nrow(df_iso_curve), c("x2", "x2_ref", "x2_ref_cap")] <-
                (df_iso_curve$pos_x[nrow(df_iso_curve)] - min(axis_2$pos_x) +
                    (df_iso_curve$pos_y[nrow(df_iso_curve)] -
                    min(axis_1$pos_y))) / sqrt(2)
    df_iso_curve[nrow(df_iso_curve), "x1"] <-
                (df_iso_curve$pos_x[nrow(df_iso_curve)] - min(axis_2$pos_x) -
                    (df_iso_curve$pos_y[nrow(df_iso_curve)] -
                    min(axis_1$pos_y))) / sqrt(2)

    # calculate CI across range to concentration ratios (in the rotated log10 space)
    df_iso_curve$log10_ratio_conc <- df_iso_curve$x1
    # delta in the log10 space is the ratio of 
    df_iso_curve$log2_CI <- zoo::rollmean( 
        log2(10) * (df_iso_curve$x2 - df_iso_curve$x2_ref) / sqrt(2), 4,
        fill = c(0, 0, 0))

    # cap position for plotting the isobolograms
    df_iso$pos_y <- pmin(df_iso$pos_y, max(axis_1$pos_y) + log2_pos_offset)
    df_iso$pos_x <- pmin(df_iso$pos_x, max(axis_2$pos_x) + log2_pos_offset)
    df_iso_curve$pos_y <- pmin(df_iso_curve$pos_y, max(axis_1$pos_y) + log2_pos_offset)
    df_iso_curve$pos_x <- pmin(df_iso_curve$pos_x, max(axis_2$pos_x) + log2_pos_offset)
    df_iso_curve$pos_y_ref <- pmin(df_iso_curve$pos_y_ref,
        max(axis_1$pos_y) + log2_pos_offset)
    df_iso_curve$pos_x_ref <- pmin(df_iso_curve$pos_x_ref,
        max(axis_2$pos_x) + log2_pos_offset)

    range <- 2 # in log10 domain --> 100-fold range for calculating averaged CI 
    ratio_idx <- which(
        (df_iso_curve$log10_ratio_conc > (min(df_iso_curve$log10_ratio_conc) + range / 2)) &
        (df_iso_curve$log10_ratio_conc < (max(df_iso_curve$log10_ratio_conc) - range / 2)))

    if (length(ratio_idx) == 0) {
      ratio_idx <- round(nrow(df_iso_curve) / 2)
    }

    df_100x_AUC <- data.frame(log10_ratio_conc = df_iso_curve$log10_ratio_conc[ratio_idx],
      AUC_log2CI <- vapply(ratio_idx, function(x) mean(df_iso_curve$log2_CI[
        (df_iso_curve$log10_ratio_conc > (df_iso_curve$log10_ratio_conc[x] - range / 2)) &
        (df_iso_curve$log10_ratio_conc <= (df_iso_curve$log10_ratio_conc[x] + range / 2))]),
        FUN.VALUE = double(1)
      )
    )

    all_iso[[as.character(isobol_value)]] <- list(
      df_iso = df_iso,
      df_iso_curve = df_iso_curve,
      ref_conc_1 = ref_conc_1,
      ref_conc_2 = ref_conc_2,
      # max_log2CI = max(df_iso_curve$log2_CI),
      # min_log2CI = min(df_iso_curve$log2_CI),
      AUC_log2CI = min(df_100x_AUC$AUC_log2CI)
      )
    }
    all_iso <- all_iso[!vapply(all_iso, is.null, FUN.VALUE = logical(1))]

    df_all_iso_points <- do.call(rbind, lapply(names(all_iso), 
            function(x) cbind(iso_level = x, all_iso[[x]]$df_iso[, c('conc_1', 'conc_2', 'pos_x', 'pos_y', 'fit_type')])))
    df_all_iso_curves <- do.call(rbind, lapply(names(all_iso), 
            function(x) cbind(iso_level = x, all_iso[[x]]$df_iso_curve[, c('pos_x', 'pos_y', 'pos_x_ref', 'pos_y_ref', 'log10_ratio_conc', 'log2_CI')])))
    df_all_AUC_log2CI <- do.call(rbind, lapply(names(all_iso), 
            function(x) data.frame(iso_level = x, CI_100x = all_iso[[x]]$AUC_log2CI, ref_conc_1 = all_iso[[x]]$ref_conc_1, ref_conc_2 = all_iso[[x]]$ref_conc_2)))
    
    isobologram_out <- list(
      df_all_iso_points = df_all_iso_points,
      df_all_iso_curves = df_all_iso_curves,
      df_all_AUC_log2CI = df_all_AUC_log2CI
    )
}
