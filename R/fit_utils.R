####
# apply_combo_sa_fits — Step 1+2 of fit_SE.combinations: SA fits → Metrics assay
####

#' apply_combo_sa_fits
#'
#' Fit single-agent dose-response curves for each co-treatment concentration
#' in a combination SE and store the results in the \code{Metrics} assay.
#' This is \strong{step 1+2} of the \code{\link{fit_SE.combinations}} pipeline:
#'
#' \enumerate{
#'   \item Standardise concentrations and build the complete dose-response matrix.
#'   \item Fit SA curves per co-treatment concentration via
#'     \code{fit_combo_cotreatments()} and \code{fit_combo_codilutions()}.
#'   \item Compute \strong{smooth} SA predictions at every combo point via
#'     \code{map_ids_to_fits()} and store them together with the SA fit parameters.
#' }
#'
#' The resulting \code{Metrics} assay contains \code{dilution_drug},
#' \code{cotrt_value}, \code{ec50}, \code{h}, \code{x_inf}, \code{x_0},
#' and \code{normalization_type} — the same schema as produced by
#' \code{fit_SE.combinations()}.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   an \code{Averaged} assay (combination experiment).
#' @param series_identifiers character vector of length 2: column names for the
#'   two concentration axes (e.g. \code{c("Concentration", "Concentration_2")}).
#'   Defaults to \code{get_default_nested_identifiers(se, data_model("combination"))}.
#' @param normalization_types character vector of normalization types to process.
#'   Default \code{c("GR", "RV")}.
#' @param averaged_assay string; name of the input assay. Default \code{"Averaged"}.
#' @param metrics_assay string; name of the output assay. Default \code{"Metrics"}.
#' @param fit_source string recorded in the \code{fit_source} column. Default \code{"gDR"}.
#'
#' @return Updated \code{SummarizedExperiment} with a \code{metrics_assay}
#'   assay containing SA fit parameters per co-treatment concentration.
#'
#' @seealso \code{\link{apply_combo_excess}}, \code{\link{apply_combo_isobolograms}},
#'   \code{\link{apply_combo_scores}}, \code{\link{fit_SE.combinations}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' SummarizedExperiment::assays(combo_se) <-
#'   SummarizedExperiment::assays(combo_se)["Averaged"]
#' combo_se_fitted <- apply_combo_sa_fits(combo_se[1, 1])
#' "Metrics" %in% SummarizedExperiment::assayNames(combo_se_fitted)
#'
#' @importFrom checkmate assert_class assert_string assert_character
#' @export
apply_combo_sa_fits <- function(se,
                                series_identifiers = NULL,
                                normalization_types = c("GR", "RV"),
                                averaged_assay = "Averaged",
                                metrics_assay = "Metrics",
                                fit_source = "gDR") {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE, len = 2L)
  checkmate::assert_character(normalization_types, min.len = 1L)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_string(fit_source)
  gDRutils::validate_se_assay_name(se, averaged_assay)

  if (is.null(series_identifiers)) {
    series_identifiers <- get_default_nested_identifiers(
      se, data_model(gDRutils::get_supported_experiments("combo"))
    )
  }

  id <- series_identifiers[1]
  id2 <- series_identifiers[2]

  avg <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, averaged_assay)
  )
  iterator <- unique(avg[, c("column", "row")])

  out_list <- gDRutils::loop(seq_len(NROW(iterator)), function(row_idx) {
    # Each list element must be named "metrics" for rbindParallelList compatibility
    x <- iterator[row_idx, ]
    i <- x[["row"]]
    j <- x[["column"]]

    avg_combo <- data.table::as.data.table(
      avg[avg[["row"]] == i & avg[["column"]] == j, ]
    )

    if (all(is.na(avg_combo[, .SD, .SDcols = !c("row", "column", "normalization_type")]))) {
      return(list(metrics = data.table::data.table(row_id = i, col_id = j)))
    }

    # Standardise concentrations (handles numeric precision mismatches)
    conc_map <- gDRutils::map_conc_to_standardized_conc(avg_combo[[id]], avg_combo[[id2]])
    avg_combo[[id]] <- replace_conc_with_standardized_conc(avg_combo[[id]], conc_map, "concs", "rconcs")
    avg_combo[[id2]] <- replace_conc_with_standardized_conc(avg_combo[[id2]], conc_map, "concs", "rconcs")

    mean_avg_combo <- avg_combo[, lapply(.SD, mean),
                                by = c(id, id2, "normalization_type"),
                                .SDcols = c("x", "x_std")]

    # Build complete concentration matrix
    conc1 <- table(avg_combo[[id]])
    conc1 <- sort(as.numeric(names(conc1)[conc1 > (max(conc1[names(conc1) != "0"]) / 2)]))
    conc2 <- table(avg_combo[[id2]])
    conc2 <- sort(as.numeric(names(conc2)[conc2 > (max(conc2[names(conc2) != "0"]) / 2)]))
    complete_cj <- data.table::CJ(unique(c(0, conc1)), unique(c(0, conc2)))
    data.table::setnames(complete_cj, c(id, id2))
    complete <- mean_avg_combo[complete_cj, on = c(id, id2)]

    avg_combo_dt <- data.table::as.data.table(mean_avg_combo)
    metrics_list <- vector("list", length(normalization_types))

    for (nt_idx in seq_along(normalization_types)) {
      norm_type <- normalization_types[[nt_idx]]
      avg_subset <- avg_combo_dt[normalization_type == norm_type]
      complete_subset <- complete[normalization_type == norm_type | is.na(normalization_type)]
      complete_subset[is.na(normalization_type), normalization_type := norm_type]

      # SA fits per co-treatment concentration (drug_1: series=conc1, cotrt=conc2)
      drug_1 <- fit_combo_cotreatments(avg_subset, series_id = id, cotrt_id = id2, norm_type)
      drug_1 <- if (all(is.na(drug_1$fit_type))) {
        x_tmp <- drug_1[1, ]
        x_tmp[, "cotrt_value"] <- NA
        x_tmp
      } else {
        na.omit(drug_1, col = "fit_type")
      }

      # SA fits per co-treatment concentration (drug_2: series=conc2, cotrt=conc1)
      drug_2 <- na.omit(
        fit_combo_cotreatments(avg_subset, series_id = id2, cotrt_id = id, norm_type),
        col = "fit_type"
      )

      # Codilution fits (diagonal)
      codilution <- na.omit(
        fit_combo_codilutions(avg_subset, series_identifiers, norm_type),
        col = "fit_type"
      )

      metrics_names <- c("drug_1", "drug_2")

      # Smooth SA predictions at every combo concentration
      complete_subset$col_values <- map_ids_to_fits(
        pred = complete_subset[[id]], match_col = complete_subset[[id2]],
        drug_1, "cotrt_value"
      )
      complete_subset$row_values <- map_ids_to_fits(
        pred = complete_subset[[id2]], match_col = complete_subset[[id]],
        drug_2, "cotrt_value"
      )
      if (!is.null(codilution)) {
        complete_subset$codil_values <- map_ids_to_fits(
          pred = complete_subset[[id2]] + complete_subset[[id]],
          match_col = gDRutils::round_concentration(complete_subset[[id2]] / complete_subset[[id]], 1),
          codilution, "ratio"
        )
        metrics_names <- c(metrics_names, "codilution")
      }

      metrics_merged <- data.table::rbindlist(mget(metrics_names), fill = TRUE)
      metrics_merged$dilution_drug <- rep(
        metrics_names,
        vapply(mget(metrics_names), NROW, numeric(1))
      )
      metrics_merged <- metrics_merged[
        !(metrics_merged$fit_type %in% c("DRCInvalidFitResult", "DRCTooFewPointsToFit")), ]
      if (NROW(metrics_merged) == 0L) {
        metrics_merged <- metrics_merged[, lapply(.SD, function(x) NA)]
        metrics_merged[, `:=`(normalization_type = norm_type, fit_source = fit_source)]
      } else {
        metrics_merged[, fit_source := fit_source]
      }

      metrics_merged$row_id <- i
      metrics_merged$col_id <- j
      metrics_list[[nt_idx]] <- metrics_merged
    }

    list(metrics = data.table::rbindlist(metrics_list, fill = TRUE))
  })

  all_metrics <- rbindParallelList(out_list, "metrics")
  SummarizedExperiment::assays(se)[[metrics_assay]] <-
    convertDFtoBumpyMatrixUsingIds(all_metrics)
  se
}


####
# .compute_smooth_matrix — shared helper: rebuild smooth values from Metrics SA fits
####

# Rebuild the smooth-values matrix for one (drug-combo × cell-line × norm_type)
# cell from the SA fit parameters stored in the Metrics assay.
# Returns a data.table with columns: id, id2, normalization_type, smooth.
#' @keywords internal
.compute_smooth_matrix <- function(avg_combo, met_dt, id, id2, norm_type) {
  avg_sub <- avg_combo[avg_combo$normalization_type == norm_type, ]
  met_sub <- met_dt[met_dt$normalization_type == norm_type, ]
  if (NROW(avg_sub) == 0L || NROW(met_sub) == 0L) return(NULL)

  # Standardize concentrations
  conc_map <- gDRutils::map_conc_to_standardized_conc(avg_sub[[id]], avg_sub[[id2]])
  avg_sub[[id]] <- replace_conc_with_standardized_conc(avg_sub[[id]], conc_map, "concs", "rconcs")
  avg_sub[[id2]] <- replace_conc_with_standardized_conc(avg_sub[[id2]], conc_map, "concs", "rconcs")
  mean_avg <- avg_sub[, lapply(.SD, mean),
                      by = c(id, id2, "normalization_type"), .SDcols = c("x", "x_std")]

  conc1 <- table(avg_sub[[id]])
  conc1 <- sort(as.numeric(names(conc1)[conc1 > (max(conc1[names(conc1) != "0"]) / 2)]))
  conc2 <- table(avg_sub[[id2]])
  conc2 <- sort(as.numeric(names(conc2)[conc2 > (max(conc2[names(conc2) != "0"]) / 2)]))
  complete_cj <- data.table::CJ(unique(c(0, conc1)), unique(c(0, conc2)))
  data.table::setnames(complete_cj, c(id, id2))
  complete_sub <- mean_avg[complete_cj, on = c(id, id2)]
  complete_sub[is.na(normalization_type), normalization_type := norm_type]

  drug1_p <- met_sub[met_sub$dilution_drug == "drug_1", ]
  drug2_p <- met_sub[met_sub$dilution_drug == "drug_2", ]
  codil_p <- met_sub[met_sub$dilution_drug == "codilution", ]

  complete_sub$col_values <- map_ids_to_fits(complete_sub[[id]], complete_sub[[id2]], drug1_p, "cotrt_value")
  complete_sub$row_values <- map_ids_to_fits(complete_sub[[id2]], complete_sub[[id]], drug2_p, "cotrt_value")
  if (NROW(codil_p) > 0L) {
    complete_sub$codil_values <- map_ids_to_fits(
      complete_sub[[id2]] + complete_sub[[id]],
      gDRutils::round_concentration(complete_sub[[id2]] / complete_sub[[id]], 1),
      codil_p, "ratio"
    )
  }

  keep <- intersect(colnames(complete_sub), c("x", "row_values", "col_values", "codil_values"))
  mat <- as.matrix(complete_sub[, keep, with = FALSE])
  complete_sub$smooth <- rowMeans(mat, na.rm = TRUE)
  av_matrix <- complete_sub[, c(id, id2, "smooth"), with = FALSE]
  av_matrix[, normalization_type := norm_type]
  av_matrix[get(id) == 0 & get(id2) == 0, smooth := 1]
  av_matrix
}


####
# apply_combo_excess — Step 3 of fit_SE.combinations: smooth → excess assay
####

#' apply_combo_excess
#'
#' Compute per-point HSA and Bliss excess values for a combination SE,
#' writing the result to the \code{excess} assay.
#' This is \strong{step 3} of the \code{\link{fit_SE.combinations}} pipeline.
#'
#' Requires the \code{Averaged} and \code{Metrics} assays to be present
#' (the latter as produced by \code{\link{apply_combo_sa_fits}} or
#' \code{fit_SE.combinations}).
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   \code{Averaged} and \code{Metrics} assays.
#' @param series_identifiers character vector of length 2: concentration column
#'   names. Defaults to \code{get_default_nested_identifiers(se, ...)}.
#' @param normalization_types character vector. Default \code{c("GR", "RV")}.
#' @param averaged_assay string. Default \code{"Averaged"}.
#' @param metrics_assay string. Default \code{"Metrics"}.
#' @param excess_assay string; name of the output assay. Default \code{"excess"}.
#'
#' @return Updated SE with an \code{excess_assay} assay containing per-point
#'   \code{smooth}, \code{hsa_excess}, and \code{bliss_excess} columns.
#'
#' @seealso \code{\link{apply_combo_sa_fits}}, \code{\link{apply_combo_isobolograms}},
#'   \code{\link{apply_combo_scores}}, \code{\link{fit_SE.combinations}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_excess <- apply_combo_excess(combo_se[1, 1])
#' "excess" %in% SummarizedExperiment::assayNames(combo_se_excess)
#'
#' @importFrom checkmate assert_class assert_string assert_character
#' @export
apply_combo_excess <- function(se,
                               series_identifiers = NULL,
                               normalization_types = c("GR", "RV"),
                               averaged_assay = "Averaged",
                               metrics_assay = "Metrics",
                               excess_assay = "excess") {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE, len = 2L)
  checkmate::assert_character(normalization_types, min.len = 1L)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_string(excess_assay)
  gDRutils::validate_se_assay_name(se, averaged_assay)
  gDRutils::validate_se_assay_name(se, metrics_assay)

  if (is.null(series_identifiers)) {
    series_identifiers <- get_default_nested_identifiers(
      se, data_model(gDRutils::get_supported_experiments("combo"))
    )
  }
  id <- series_identifiers[1]
  id2 <- series_identifiers[2]

  avg <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, averaged_assay)
  ))
  met <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, metrics_assay),
    row.field = "row", column.field = "column"
  ))
  iterator <- unique(avg[, c("column", "row")])

  out_list <- gDRutils::loop(seq_len(NROW(iterator)), function(row_idx) {
    i <- iterator[row_idx, "row"][[1L]]
    j <- iterator[row_idx, "column"][[1L]]
    avg_combo <- avg[avg$row == i & avg$column == j, ]
    met_cell <- met[met$row == i & met$column == j, ]

    if (all(is.na(avg_combo[, .SD, .SDcols = !c("row", "column", "normalization_type")]))) {
      return(list(excess = data.table::data.table(row_id = i, col_id = j)))
    }

    excess_list <- vector("list", length(normalization_types))
    for (nt_idx in seq_along(normalization_types)) {
      norm_type <- normalization_types[[nt_idx]]

      av_matrix <- .compute_smooth_matrix(avg_combo, met_cell, id, id2, norm_type)
      if (is.null(av_matrix) || NROW(av_matrix) == 0L) next

      cols <- c(id, id2, "smooth")
      sa1 <- av_matrix[av_matrix[[id2]] == 0, cols, with = FALSE]
      sa2 <- av_matrix[av_matrix[[id]] == 0, cols, with = FALSE]

      hsa <- calculate_HSA(sa1, id, sa2, id2, norm_type)
      h_excess <- calculate_excess(hsa, av_matrix,
                                   series_identifiers = c(id, id2),
                                   metric_col = "metric", measured_col = "smooth")
      data.table::setnames(h_excess, "x", "hsa_excess")

      bliss <- calculate_Bliss(sa1, id, sa2, id2, norm_type)
      bliss_excess <- calculate_excess(bliss, av_matrix,
                                       series_identifiers = c(id, id2),
                                       metric_col = "metric", measured_col = "smooth")
      data.table::setnames(bliss_excess, "x", "bliss_excess")

      excess <- data.table::merge.data.table(
        data.table::merge.data.table(av_matrix, h_excess, all = TRUE),
        bliss_excess, all = TRUE
      )
      excess$row_id <- i
      excess$col_id <- j
      excess_list[[nt_idx]] <- excess
    }

    list(excess = data.table::rbindlist(excess_list, fill = TRUE))
  })

  all_excess <- rbindParallelList(out_list, "excess")
  SummarizedExperiment::assays(se)[[excess_assay]] <-
    convertDFtoBumpyMatrixUsingIds(all_excess)
  se
}


####
# apply_combo_isobolograms — Step 4 of fit_SE.combinations: Loewe CI → isobolograms assay
####

#' apply_combo_isobolograms
#'
#' Compute Loewe combination index isobolograms for a combination SE,
#' writing results to the \code{isobolograms} and \code{all_iso_points} assays.
#' This is \strong{step 4} of the \code{\link{fit_SE.combinations}} pipeline.
#'
#' Requires the \code{Averaged} and \code{Metrics} assays (as produced by
#' \code{\link{apply_combo_sa_fits}} or \code{fit_SE.combinations}).
#' Isobologram analysis is only performed when at least 9 combo concentration
#' points are available — otherwise \code{NA} is written, matching
#' \code{fit_SE.combinations} behaviour.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   \code{Averaged} and \code{Metrics} assays.
#' @param series_identifiers character vector of length 2. Defaults to
#'   \code{get_default_nested_identifiers(se, ...)}.
#' @param normalization_types character vector. Default \code{c("GR", "RV")}.
#' @param averaged_assay string. Default \code{"Averaged"}.
#' @param metrics_assay string. Default \code{"Metrics"}.
#' @param isobolograms_assay string; name of the Loewe curves output assay.
#'   Default \code{"isobolograms"}.
#' @param iso_points_assay string; name of the iso points output assay.
#'   Default \code{"all_iso_points"}.
#'
#' @return Updated SE with \code{isobolograms_assay} and \code{iso_points_assay}
#'   assays.
#'
#' @seealso \code{\link{apply_combo_sa_fits}}, \code{\link{apply_combo_excess}},
#'   \code{\link{apply_combo_scores}}, \code{\link{fit_SE.combinations}}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_iso <- apply_combo_isobolograms(combo_se[1, 1])
#' "isobolograms" %in% SummarizedExperiment::assayNames(combo_se_iso)
#'
#' @importFrom checkmate assert_class assert_string assert_character
#' @export
apply_combo_isobolograms <- function(se,
                                     series_identifiers = NULL,
                                     normalization_types = c("GR", "RV"),
                                     averaged_assay = "Averaged",
                                     metrics_assay = "Metrics",
                                     isobolograms_assay = "isobolograms",
                                     iso_points_assay = "all_iso_points") {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE, len = 2L)
  checkmate::assert_character(normalization_types, min.len = 1L)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_string(isobolograms_assay)
  checkmate::assert_string(iso_points_assay)
  gDRutils::validate_se_assay_name(se, averaged_assay)
  gDRutils::validate_se_assay_name(se, metrics_assay)

  if (is.null(series_identifiers)) {
    series_identifiers <- get_default_nested_identifiers(
      se, data_model(gDRutils::get_supported_experiments("combo"))
    )
  }
  id <- series_identifiers[1]
  id2 <- series_identifiers[2]

  avg <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, averaged_assay)
  ))
  met <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, metrics_assay),
    row.field = "row", column.field = "column"
  ))
  iterator <- unique(avg[, c("column", "row")])

  out_list <- gDRutils::loop(seq_len(NROW(iterator)), function(row_idx) {
    i <- iterator[row_idx, "row"][[1L]]
    j <- iterator[row_idx, "column"][[1L]]
    avg_combo <- avg[avg$row == i & avg$column == j, ]
    met_cell <- met[met$row == i & met$column == j, ]

    if (all(is.na(avg_combo[, .SD, .SDcols = !c("row", "column", "normalization_type")]))) {
      return(list(
        isobolograms = data.table::data.table(row_id = i, col_id = j),
        iso_points = data.table::data.table(row_id = i, col_id = j)
      ))
    }

    iso_curves_list <- vector("list", length(normalization_types))
    iso_points_list <- vector("list", length(normalization_types))

    for (nt_idx in seq_along(normalization_types)) {
      norm_type <- normalization_types[[nt_idx]]

      av_matrix <- .compute_smooth_matrix(avg_combo, met_cell, id, id2, norm_type)
      if (is.null(av_matrix) || NROW(av_matrix) == 0L) next

      drug1_p <- met_cell[met_cell$normalization_type == norm_type &
                          met_cell$dilution_drug == "drug_1", ]
      drug2_p <- met_cell[met_cell$normalization_type == norm_type &
                          met_cell$dilution_drug == "drug_2", ]
      codil_p <- met_cell[met_cell$normalization_type == norm_type &
                          met_cell$dilution_drug == "codilution", ]
      if (NROW(codil_p) == 0L) codil_p <- NULL

      # Determine whether enough dense combo points exist for Loewe
      discard1 <- names(which(
        table(av_matrix[!is.na(av_matrix$smooth), id, with = FALSE]) < 3
      ))
      discard2 <- names(which(
        table(av_matrix[!is.na(av_matrix$smooth), id2, with = FALSE]) < 3
      ))
      av_dense <- av_matrix[
        !(av_matrix[[id]] %in% discard1) & !(av_matrix[[id2]] %in% discard2), ]

      isobologram_out <- if (
        sum((av_dense[[id]] * av_dense[[id2]]) > 0) > 9 &&
        NROW(drug2_p) > 0L && any(drug2_p$cotrt_value == 0, na.rm = TRUE)
      ) {
        calculate_Loewe(
          av_matrix, drug2_p, drug1_p, codil_p,
          series_identifiers = c(id, id2),
          normalization_type = norm_type
        )
      } else {
        list(
          df_all_iso_points = data.table::data.table(row_id = NA, col_id = NA),
          df_all_iso_curves = data.table::data.table(row_id = NA, col_id = NA)
        )
      }

      isobologram_out$df_all_iso_points$row_id <- i
      isobologram_out$df_all_iso_points$col_id <- j
      isobologram_out$df_all_iso_curves$row_id <- i
      isobologram_out$df_all_iso_curves$col_id <- j
      isobologram_out$df_all_iso_points$normalization_type <- norm_type
      isobologram_out$df_all_iso_curves$normalization_type <- norm_type

      iso_points_list[[nt_idx]] <- isobologram_out$df_all_iso_points
      iso_curves_list[[nt_idx]] <- isobologram_out$df_all_iso_curves
    }

    list(
      isobolograms = data.table::rbindlist(iso_curves_list, fill = TRUE),
      iso_points = data.table::rbindlist(iso_points_list, fill = TRUE)
    )
  })

  all_isobolograms <- rbindParallelList(out_list, "isobolograms")
  all_iso_points <- rbindParallelList(out_list, "iso_points")
  SummarizedExperiment::assays(se)[[isobolograms_assay]] <-
    convertDFtoBumpyMatrixUsingIds(all_isobolograms)
  SummarizedExperiment::assays(se)[[iso_points_assay]] <-
    convertDFtoBumpyMatrixUsingIds(all_iso_points)
  se
}


####
# apply_combo_scores — high-level wrapper replicating fit_SE.combinations scoring
####

#' apply_combo_scores
#'
#' Compute Bliss and HSA synergy scores for a combination SE, replicating the
#' exact scoring logic of \code{\link{fit_SE.combinations}}.
#'
#' Unlike the low-level \code{\link{bliss_fit_fn}} and \code{\link{hss_fit_fn}}
#' (which work on raw Averaged data per triplet), this function uses the SA fit
#' parameters from the \code{Metrics} assay to generate curve-smoothed single-agent
#' responses via \code{\link[gDRutils]{predict_efficacy_from_conc}} before computing
#' excess.  This produces results numerically identical to \code{fit_SE.combinations}.
#'
#' @details
#' The function requires a \code{Metrics} assay containing columns
#' \code{dilution_drug}, \code{cotrt_value}, \code{ec50}, \code{h},
#' \code{x_inf}, \code{x_0}, and \code{normalization_type} — as produced by
#' \code{fit_SE.combinations} or \code{\link{apply_fit_to_se}} with
#' \code{fit_drug_response_metrics}.
#'
#' Scoring steps (per drug-combo × cell-line × normalization_type):
#' \enumerate{
#'   \item Predict smooth SA responses at every combo concentration using
#'     \code{predict_efficacy_from_conc} on \code{drug_1} and \code{drug_2}
#'     parameter sets from \code{Metrics}.
#'   \item Average \code{col_values} (drug-1-along-conc1) and
#'     \code{row_values} (drug-2-along-conc2) → \code{smooth}.
#'   \item Compute Bliss-expected using \code{\link{calculate_Bliss}} and
#'     HSA-expected using \code{\link{calculate_HSA}} on the smooth SA edges.
#'   \item Compute excess via \code{\link{calculate_excess}}.
#'   \item Score = mean of top 10-percentile excess via \code{\link{calculate_score}}.
#' }
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} with
#'   \code{"Averaged"} and \code{"Metrics"} assays (combination experiment).
#' @param scores_assay string; name of the output assay to write scores into.
#'   Default \code{"scores"}.
#' @param averaged_assay string; name of the input assay. Default \code{"Averaged"}.
#' @param metrics_assay string; name of the assay containing SA fit parameters.
#'   Default \code{"Metrics"}.
#' @param excess_assay string or \code{NULL}; if the SE already contains a
#'   pre-computed excess assay (e.g. from \code{\link{apply_combo_excess}}),
#'   pass its name here to skip excess re-computation and score directly from it.
#'   Default \code{NULL} (excess is computed internally).
#' @param normalization_types character vector of normalization types to process.
#'   Default \code{c("GR", "RV")}.
#' @param fit_source string recorded in the \code{fit_source} column. Default
#'   \code{"gDR"}.
#' @param score_FUN function to reduce per-point excess values to a scalar
#'   score. Default \code{\link{calculate_score}} (mean of top 10-percentile).
#'
#' @return Updated \code{SummarizedExperiment} with a \code{scores_assay}
#'   assay containing \code{bliss_score} and \code{hsa_score} per triplet.
#'
#' @seealso \code{\link{bliss_fit_fn}}, \code{\link{hss_fit_fn}} for the
#'   simplified raw-data variants.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_out <- apply_combo_scores(combo_se)
#' "scores" %in% SummarizedExperiment::assayNames(combo_se_out)
#'
#' @importFrom checkmate assert_class assert_string assert_character assert_function
#' @export
apply_combo_scores <- function(se,
                               scores_assay = "scores",
                               averaged_assay = "Averaged",
                               metrics_assay = "Metrics",
                               excess_assay = NULL,
                               normalization_types = c("GR", "RV"),
                               fit_source = "gDR",
                               score_FUN = calculate_score) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(scores_assay)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_string(excess_assay, null.ok = TRUE)
  checkmate::assert_character(normalization_types, min.len = 1L)
  checkmate::assert_function(score_FUN)
  gDRutils::validate_se_assay_name(se, averaged_assay)
  gDRutils::validate_se_assay_name(se, metrics_assay)

  # Fast path: if excess assay already exists, score directly from it.
  # This gives results numerically identical to fit_SE.combinations()
  # because both use the same excess values computed in apply_combo_excess().
  if (!is.null(excess_assay) && excess_assay %in% SummarizedExperiment::assayNames(se)) {
    exc <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
      SummarizedExperiment::assay(se, excess_assay),
      row.field = "row", column.field = "column"
    ))
    score_rows <- lapply(
      split(exc, by = c("row", "column", "normalization_type"), sorted = TRUE),
      function(dt) {
        data.table::data.table(
          row = dt$row[1L],
          column = dt$column[1L],
          normalization_type = dt$normalization_type[1L],
          fit_source = fit_source,
          bliss_score = if (all(is.na(dt$bliss_excess))) NA_real_ else
            score_FUN(dt$bliss_excess),
          hsa_score = if (all(is.na(dt$hsa_excess))) NA_real_ else
            score_FUN(dt$hsa_excess)
        )
      }
    )
    scores_dt <- data.table::rbindlist(score_rows)
    se <- gDRutils::persist_fit_assay(se, scores_dt, "replace", scores_assay, "row", "column",
                         upsert_key_cols = c("fit_source", "normalization_type"))
    return(se)
  }

  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")
  series_identifiers <- c(conc1_col, conc2_col)

  avg_bm <- SummarizedExperiment::assay(se, averaged_assay)
  met_bm <- SummarizedExperiment::assay(se, metrics_assay)

  # Iterate over each (drug-combo × cell-line) BumpyMatrix cell
  score_rows <- vector("list", prod(dim(se)) * length(normalization_types))
  idx <- 0L

  for (ri in seq_len(NROW(se))) {
    for (ci in seq_len(NCOL(se))) {
      avg_dt <- data.table::as.data.table(avg_bm[ri, ci][[1L]])
      met_dt <- data.table::as.data.table(met_bm[ri, ci][[1L]])

      if (NROW(avg_dt) == 0L || NROW(met_dt) == 0L) next

      row_nm <- rownames(se)[ri]
      col_nm <- colnames(se)[ci]

      for (norm_type in normalization_types) {
        avg_sub <- avg_dt[avg_dt$normalization_type == norm_type, ]
        met_sub <- met_dt[met_dt$normalization_type == norm_type, ]
        if (NROW(avg_sub) == 0L || NROW(met_sub) == 0L) next

        # Smooth SA predictions for every point in the averaged combo data
        # drug_1 params: SA along conc1 for each cotrt conc2 value
        drug1_params <- met_sub[met_sub$dilution_drug == "drug_1", ]
        drug2_params <- met_sub[met_sub$dilution_drug == "drug_2", ]

        if (NROW(drug1_params) == 0L || NROW(drug2_params) == 0L) next

        # col_values: predict drug-1 response at conc1 for each cotrt conc2
        avg_sub$col_values <- vapply(seq_len(NROW(avg_sub)), function(i) {
          cotrt <- avg_sub[[conc2_col]][i]
          conc <- avg_sub[[conc1_col]][i]
          p <- drug1_params[abs(drug1_params$cotrt_value - cotrt) ==
                              min(abs(drug1_params$cotrt_value - cotrt)), ][1L, ]
          if (NROW(p) == 0L || any(is.na(c(p$ec50, p$h, p$x_inf, p$x_0)))) {
            NA_real_
          } else {
            gDRutils::predict_efficacy_from_conc(conc, p$x_inf, p$x_0, p$ec50, p$h)
          }
        }, numeric(1))

        # row_values: predict drug-2 response at conc2 for each cotrt conc1
        avg_sub$row_values <- vapply(seq_len(NROW(avg_sub)), function(i) {
          cotrt <- avg_sub[[conc1_col]][i]
          conc <- avg_sub[[conc2_col]][i]
          p <- drug2_params[abs(drug2_params$cotrt_value - cotrt) ==
                              min(abs(drug2_params$cotrt_value - cotrt)), ][1L, ]
          if (NROW(p) == 0L || any(is.na(c(p$ec50, p$h, p$x_inf, p$x_0)))) {
            NA_real_
          } else {
            gDRutils::predict_efficacy_from_conc(conc, p$x_inf, p$x_0, p$ec50, p$h)
          }
        }, numeric(1))

        avg_sub$smooth <- rowMeans(
          avg_sub[, c("col_values", "row_values"), with = FALSE],
          na.rm = TRUE
        )
        avg_sub[avg_sub[[conc1_col]] == 0 & avg_sub[[conc2_col]] == 0, smooth := 1]

        # SA edges for calculate_HSA / calculate_Bliss
        sa1 <- avg_sub[avg_sub[[conc2_col]] == 0, c(series_identifiers, "smooth"), with = FALSE]
        sa2 <- avg_sub[avg_sub[[conc1_col]] == 0, c(series_identifiers, "smooth"), with = FALSE]

        if (NROW(sa1) == 0L || NROW(sa2) == 0L) next

        # HSA score
        hsa_expected <- calculate_HSA(sa1, conc1_col, sa2, conc2_col, norm_type)
        h_excess <- calculate_excess(
          hsa_expected, avg_sub,
          series_identifiers = series_identifiers,
          metric_col = "metric", measured_col = "smooth"
        )
        if (all(is.na(h_excess$x))) {
          hsa_score <- NA_real_
        } else {
          hsa_score <- score_FUN(h_excess$x)
        }

        # Bliss score
        bliss_expected <- calculate_Bliss(sa1, conc1_col, sa2, conc2_col, norm_type)
        bliss_excess <- calculate_excess(
          bliss_expected, avg_sub,
          series_identifiers = series_identifiers,
          metric_col = "metric", measured_col = "smooth"
        )
        if (all(is.na(bliss_excess$x))) {
          bliss_score <- NA_real_
        } else {
          bliss_score <- score_FUN(bliss_excess$x)
        }

        idx <- idx + 1L
        score_rows[[idx]] <- data.table::data.table(
          row = row_nm,
          column = col_nm,
          normalization_type = norm_type,
          fit_source = fit_source,
          bliss_score = bliss_score,
          hsa_score = hsa_score
        )
      }
    }
  }

  if (idx > 0L) {
    scores_dt <- data.table::rbindlist(score_rows[seq_len(idx)])
    se <- gDRutils::persist_fit_assay(se, scores_dt, "replace", scores_assay, "row", "column",
                         upsert_key_cols = c("fit_source", "normalization_type"))
  }
  se
}


#' bliss_fit_fn
#'
#' Reference fit function for Bliss independence synergy scoring on combination
#' data.  Computes Bliss-expected response from single-agent edges and derives
#' excess and score from raw averaged data.
#'
#' Intended for use with \code{\link{apply_fit}} on combination SEs:
#'
#' \preformatted{
#'   gDRutils::apply_fit(
#'     combo_se, bliss_fit_fn, "combination",
#'     output_assay = "custom_bliss", fit_source = "bliss"
#'   )
#' }
#'
#' @param avg_dt \code{data.table} for one (drug1 x drug2 x cell line x
#'   normalization_type) combination.  Must contain columns \code{Concentration},
#'   \code{Concentration_2}, \code{x}, and \code{normalization_type}.
#'
#' @return Named list with \code{bliss_score}, \code{bliss_excess_mean},
#'   \code{n_combo_points}, and \code{normalization_type}.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_out <- gDRutils::apply_fit(
#'   combo_se, bliss_fit_fn, "combination",
#'   output_assay = "custom_bliss", fit_source = "bliss"
#' )
#'
#' @export
bliss_fit_fn <- function(avg_dt) {
  norm_type <- as.character(avg_dt$normalization_type[1])
  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")

  sa1 <- avg_dt[avg_dt[[conc2_col]] == 0 & avg_dt[[conc1_col]] > 0, ]
  sa2 <- avg_dt[avg_dt[[conc1_col]] == 0 & avg_dt[[conc2_col]] > 0, ]
  combo <- avg_dt[avg_dt[[conc1_col]] > 0 & avg_dt[[conc2_col]] > 0, ]

  n_combo <- NROW(combo)
  if (n_combo == 0L || NROW(sa1) == 0L || NROW(sa2) == 0L) {
    list(
      normalization_type = norm_type,
      bliss_score = NA_real_,
      bliss_excess_mean = NA_real_,
      n_combo_points = n_combo
    )
  } else {
    # Match SA responses to each combo point by concentration (per-point comparison)
    combo_ordered <- combo[order(combo[[conc1_col]], combo[[conc2_col]]), ]
    sa1_match <- sa1[["x"]][match(combo_ordered[[conc1_col]], sa1[[conc1_col]])]
    sa2_match <- sa2[["x"]][match(combo_ordered[[conc2_col]], sa2[[conc2_col]])]

    if (norm_type == "RV") {
      expected <- sa1_match * sa2_match
    } else {
      # GR adaptation: Holbeck et al., Cancer Res 2017
      # Guard for cytotoxic drugs where GR < 0: negative log2 inputs produce a
      # positive product → false synergy.  Use pmin (more toxic agent) instead.
      expected <- ifelse(
        sa1_match < 0 | sa2_match < 0,
        pmin(sa1_match, sa2_match),
        2^(log2(sa1_match + 1) * log2(sa2_match + 1)) - 1
      )
    }

    excess <- expected - combo_ordered[["x"]]
    q90 <- stats::quantile(excess, 0.9, na.rm = TRUE)

    list(
      normalization_type = norm_type,
      bliss_score = mean(excess[excess >= q90], na.rm = TRUE),
      bliss_excess_mean = mean(excess, na.rm = TRUE),
      n_combo_points = n_combo
    )
  }
}


#' hss_fit_fn
#'
#' Reference fit function for Highest Single Agent (HSA) synergy scoring on
#' combination data.  Uses raw averaged single-agent edges.
#'
#' Intended for use with \code{\link{apply_fit}} on combination SEs:
#'
#' \preformatted{
#'   gDRutils::apply_fit(
#'     combo_se, hss_fit_fn, "combination",
#'     output_assay = "custom_hss", fit_source = "hss"
#'   )
#' }
#'
#' @param avg_dt \code{data.table} for one (drug1 x drug2 x cell line x
#'   normalization_type) combination.  Must contain columns \code{Concentration},
#'   \code{Concentration_2}, \code{x}, and \code{normalization_type}.
#'
#' @return Named list with \code{hss_score}, \code{hss_excess_mean},
#'   \code{n_combo_points}, and \code{normalization_type}.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_out <- gDRutils::apply_fit(
#'   combo_se, hss_fit_fn, "combination",
#'   output_assay = "custom_hss", fit_source = "hss"
#' )
#'
#' @export
hss_fit_fn <- function(avg_dt) {
  norm_type <- as.character(avg_dt$normalization_type[1])
  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")

  sa1 <- avg_dt[avg_dt[[conc2_col]] == 0 & avg_dt[[conc1_col]] > 0, ]
  sa2 <- avg_dt[avg_dt[[conc1_col]] == 0 & avg_dt[[conc2_col]] > 0, ]
  combo <- avg_dt[avg_dt[[conc1_col]] > 0 & avg_dt[[conc2_col]] > 0, ]

  n_combo <- NROW(combo)
  if (n_combo == 0L || NROW(sa1) == 0L || NROW(sa2) == 0L) {
    list(
      normalization_type = norm_type,
      hss_score = NA_real_,
      hss_excess_mean = NA_real_,
      n_combo_points = n_combo
    )
  } else {
    # Match SA responses to each combo point by concentration (per-point comparison)
    combo_ordered <- combo[order(combo[[conc1_col]], combo[[conc2_col]]), ]
    sa1_match <- sa1[["x"]][match(combo_ordered[[conc1_col]], sa1[[conc1_col]])]
    sa2_match <- sa2[["x"]][match(combo_ordered[[conc2_col]], sa2[[conc2_col]])]

    # HSA expected: the more potent single agent at each combo grid point
    expected <- pmin(sa1_match, sa2_match)
    excess <- expected - combo_ordered[["x"]]
    q90 <- stats::quantile(excess, 0.9, na.rm = TRUE)

    list(
      normalization_type = norm_type,
      hss_score = mean(excess[excess >= q90], na.rm = TRUE),
      hss_excess_mean = mean(excess, na.rm = TRUE),
      n_combo_points = n_combo
    )
  }
}



####
# Deprecated aliases
####

#' Deprecated alias for apply_fit
#'
#' @param ... arguments passed to \code{\link[gDRutils]{apply_fit}}
#'
#' @return updated \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'
#' @examples
#' # superseded — call gDRutils::apply_fit() instead
#' @keywords internal
#' @export
apply_custom_fit <- function(...) {
  .Deprecated("gDRutils::apply_fit")
  gDRutils::apply_fit(...)
}

#' Deprecated alias for apply_fits
#'
#' @param ... arguments passed to \code{\link[gDRutils]{apply_fits}}
#'
#' @return updated \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'
#' @examples
#' # superseded — call gDRutils::apply_fits() instead
#' @keywords internal
#' @export
apply_custom_fits <- function(...) {
  .Deprecated("gDRutils::apply_fits")
  gDRutils::apply_fits(...)
}
