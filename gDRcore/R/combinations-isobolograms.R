#' Define matrix grid positions
#'
#' @param conc1 drug_1 concentration
#' @param conc2 drug_2 concentration
#'
#' @details
#' drug_1 is diluted along the rows as the y-axis and
#' drug_2 is diluted along the columns and will be the x-axis.
define_matrix_grid_positions <- function(conc1, conc2) {
  .generate_gap_for_single_agent <- function(x) {
    2 * x[2] - x[3] - log10(1.5)
  } 

  conc_1 <- sort(unique(round_concentration(conc1)))
  pos_y <- log10conc_1 <- log10(conc_1)
  pos_y[1] <- .generate_gap_for_single_agent(log10conc_1)
  axis_1 <- data.frame(conc_1 = conc_1,
    log10conc_1 = log10conc_1,
    pos_y = pos_y,
    marks_y = sprintf("%.2g", conc_1)
  )

  conc_2 <- sort(unique(round_concentration(conc2)))
  pos_x <- log10conc_2 <- log10(conc_2)
  pos_x[1] <- .generate_gap_for_single_agent(log10conc_2)
  axis_2 <- data.frame(conc_2 = conc_2,
    log10conc_2 = log10conc_2,
    pos_x = pos_x,
    marks_x = sprintf("%.2g", conc_2)
  )

  list(axis_1 = axis_1, axis_2 = axis_2)
}


#' @importFrom zoo rollmean
calculate_Loewe <- function(df_mean, 
                      row_fittings, 
                      col_fittings, 
                      codilution_fittings, 
                      series_identifiers,
                      normalization_type,
                      conc_margin = 10 ^ 0.5,
                      log2_pos_offset = log10(3) / 2,
                      min_n_conc = 6,
                      min_r2 = 0.8
                    ) {

  checkmate::assert_number(conc_margin)
  checkmate::assert_number(log2_pos_offset)
  checkmate::assert_character(normalization_type) 
  if (length(series_identifiers) != 2L) {
    stop("only series_identifiers of length 2 are currently supported")
  }
 
  iso_cutoffs <- get_isocutoffs(df_mean, normalization_type)

  all_iso <- vector("list", length(iso_cutoffs))
  names(all_iso) <- iso_cutoffs

  axes <- define_matrix_grid_positions(df_mean[[series_identifiers[1]]], df_mean[[series_identifiers[2]]])
  axis_1 <- axes$axis_1
  axis_2 <- axes$axis_2
  
  max1_cap <- round_concentration(max(axis_1$conc_1) * conc_margin)
  max2_cap <- round_concentration(max(axis_2$conc_2) * conc_margin)

  row_fittings <- row_fittings[order(row_fittings$cotrt_value, decreasing = TRUE), ]
  col_fittings <- col_fittings[order(col_fittings$cotrt_value, decreasing = FALSE), ]
  if (!is.null(codilution_fittings)) {
    codilution_fittings <- codilution_fittings[order(codilution_fittings$ratio, decreasing = TRUE), ]
    codilution_fittings <- codilution_fittings[codilution_fittings$fit_type %in% "DRC3pHillFitModelFixS0", ]
    # more restrictive to be sure that it adds fit quality
    codilution_fittings <-
      codilution_fittings[codilution_fittings$N_conc >= min_n_conc &
                            codilution_fittings$r2 > min_r2, ]
  }
  futile.logger::flog.debug("codilution fittings after transformation/filtering:")
  futile.logger::flog.debug(knitr::kable(codilution_fittings))
  
  for (isobol_value in iso_cutoffs) { # run through the different isobolograms
    
    # avoid getting a fit result that is close to the asymptotic value (and thus less precise)
    # filtered_cf can be empty (i.e. no rows) - this edge case is handled in calculate_isobolograms
    filtered_cf <- codilution_fittings[codilution_fittings$x_inf < (isobol_value - .1), ]
    
    df_iso <- calculate_isobolograms(row_fittings, col_fittings, filtered_cf, isobol_value, max1_cap, max2_cap)
    ref_conc_1 <- pmin(df_iso$conc_1[df_iso$conc_2 == 0 & df_iso$fit_type == "by_col"], max1_cap)
    ref_conc_2 <- pmin(df_iso$conc_2[df_iso$conc_1 == 0 & df_iso$fit_type == "by_row"], max2_cap)

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
                - (sqrt(2) * df_iso_curve$x1 + min(axis_2$pos_x) - min(axis_1$pos_y)))))
    if (length(c2) == 0) {
      return(NULL)
    }
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
    df_iso_curve$log2_CI <- zoo::rollmean(log2(10) * (df_iso_curve$x2 - df_iso_curve$x2_ref) / sqrt(2),
                                          4,
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
      AUC_log2CI = min(df_100x_AUC$AUC_log2CI)
      )
    }
    all_iso <- all_iso[!vapply(all_iso, is.null, FUN.VALUE = logical(1))]

    df_all_iso_points <- do.call(rbind, lapply(names(all_iso), 
            function(x) cbind(iso_level = x,
                              all_iso[[x]]$df_iso[, c("conc_1", "conc_2", "pos_x", "pos_y", "fit_type")])))
    df_all_iso_curves <- do.call(rbind, lapply(names(all_iso), 
            function(x) cbind(iso_level = x,
                              all_iso[[x]]$df_iso_curve[, c("pos_x", "pos_y", "pos_x_ref", "pos_y_ref",
                                                            "log10_ratio_conc", "log2_CI")])))
    df_all_AUC_log2CI <- do.call(rbind, lapply(names(all_iso), 
            function(x) data.frame(iso_level = x, CI_100x = all_iso[[x]]$AUC_log2CI,
                                   ref_conc_1 = all_iso[[x]]$ref_conc_1, ref_conc_2 = all_iso[[x]]$ref_conc_2)))
    
    isobologram_out <- list(
      df_all_iso_points = df_all_iso_points,
      df_all_iso_curves = df_all_iso_curves,
      df_all_AUC_log2CI = df_all_AUC_log2CI
    )
}


get_isocutoffs <- function(df_mean, normalization_type) {
  if (min(df_mean[[normalization_type]], na.rm = TRUE) > 0.7) {
    iso_cutoffs <- NULL
  } else {
    if (normalization_type == "GRvalue") {
      max_val <- -0.25
    } else {
      max_val <- 0.2
    }
    iso_cutoffs <- seq(max(max_val, ceiling(20 * min(df_mean[[normalization_type]] + 0.08,
                                                    na.rm = TRUE)) / 20), 0.8,  0.05)
    names(iso_cutoffs) <- as.character(iso_cutoffs)
  }
  iso_cutoffs
}


calculate_isobolograms <- function(row_fittings, col_fittings, codilution_fittings, isobol_value,
  max1_cap, max2_cap) {
    # cutoff point by row
    conc2 <- gDRutils::predict_conc_from_efficacy(
      efficacy = isobol_value,
      x_inf = row_fittings$x_inf,
      x_0 = row_fittings$x_0,
      ec50 = row_fittings$ec50,
      h = row_fittings$h
    )
    row_iso <- data.frame(conc_1 = row_fittings[, "cotrt_value"],
      conc_2 = conc2,
      fit_type = "by_row")

    # cutoff point by columns
    conc1 <- gDRutils::predict_conc_from_efficacy(
      efficacy = isobol_value,
      x_inf = col_fittings$x_inf,
      x_0 = col_fittings$x_0,
      ec50 = col_fittings$ec50,
      h = col_fittings$h
    )
    col_iso <- data.frame(conc_1 = conc1,
      conc_2 = col_fittings[, "cotrt_value"],
      fit_type = "by_col")

    xy_iso <- do.call("rbind", list(row_iso, col_iso))
    xy_iso$conc_1 <- pmin(xy_iso$conc_1, max1_cap)
    xy_iso$conc_2 <- pmin(xy_iso$conc_2, max2_cap)

    # cutoff point by diagonal (co-dilution)
    # co-dil is given as concentration of drug 1
    if (!is.null(codilution_fittings) && NROW(codilution_fittings) > 1) {
      conc_mix <- gDRutils::predict_conc_from_efficacy(
        efficacy = isobol_value,
        x_inf = codilution_fittings$x_inf,
        x_0 = codilution_fittings$x_0,
        ec50 = codilution_fittings$ec50,
        h = codilution_fittings$h
      ) 
      codil_iso <- data.frame(conc_1 = conc_mix / (1 + codilution_fittings$ratio),
        conc_2 = conc_mix / (1 + (1 / codilution_fittings$ratio)),
        fit_type = "by_codil")

      # Keep isobologram values within matrix values.
      capped_idx <- codil_iso$conc_1 > max1_cap | codil_iso$conc_2 > max2_cap
      codil_iso$conc_1[capped_idx] <- pmin(max1_cap, codilution_fittings$conc_ratio * max2_cap)[capped_idx]
      codil_iso$conc_2[capped_idx] <- pmin(max2_cap, max1_cap / codilution_fittings$conc_ratio)[capped_idx]

      isobologram_values <- do.call("rbind", list(xy_iso, codil_iso))
    } else {
      isobologram_values <- xy_iso
    }
    isobologram_values
  }
