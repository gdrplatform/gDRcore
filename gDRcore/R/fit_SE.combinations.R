#' fit_SE for combination screens
#'
#' Perform fittings for combination screens.
#'
#' @param se \code{SummarizedExperiment} object with a BumpyMatrix assay containing averaged data.
#' @param series_identifiers character vector of the column names in the nested \code{DFrame}
#' corresponding to nested identifiers.
#' @param metrics_assay string of the name of the metrics assay to output
#' in the returned \linkS4class{SummarizedExperiment}
#' whose combination represents a unique series for which to fit curves. 
#' @param conc_margin margin for calculation and plots as fold-change over highest test conc for calculation
#' @param log2_pos_offset max offset for conc
#' @param norm_types character vector of normalization types used for calculating combo matrix
#'
#' @details
#' This function assumes that the combination is set up with both concentrations nested in the assay.
#'
#' @return A code{SummarizedExperiment} object with an additional assay containing the combination metrics.
#' @export
#'
fit_SE.combinations <- function(se,
                                series_identifiers = NULL,
                                normalization_types = c("GRvalue", "RelativeViability"),
                                averaged_assay = "Averaged",
                                metrics_assay = "Metrics"
                                ) {

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE)
  checkmate::assert_character(normalization_types)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  
  if (is.null(series_identifiers)) {
    series_identifiers <- get_nested_default_identifiers(se, averaged_assay)
  }
  
  if (length(series_identifiers) != 2L) {
    stop("gDR only supports 'series_identifiers' arguments with length '2' for the 'fit_SE.combinations' function")
  }
  
  avg <- BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se, averaged_assay))

  rdata <- SummarizedExperiment::rowData(se)
  cdata <- SummarizedExperiment::colData(se)
  cl_name <- gDRutils::get_SE_identifiers(se, "cellline_name")

  id <- series_identifiers[1]
  id2 <- series_identifiers[2]

  iterator <- unique(avg[, c("column", "row")])
  bliss_excess <- hsa_excess <- metrics <- isobolograms <- smooth_mx <- vector("list", nrow(iterator))
  bliss_score <- hsa_score <- CIScore_50 <- CIScore_80 <- S4Vectors::DataFrame(matrix(NA, nrow(iterator), 0))
  for (row in seq_len(nrow(iterator))) {
    x <- iterator[row, ]
    i <- x[["row"]]
    j <- x[["column"]]

    avg_combo <- avg[avg$row == i & avg$column == j, ]
    conc_map <- map_conc_to_standardized_conc(avg_combo[[id]], avg_combo[[id2]])

    avg_combo[[id]] <- replace_conc_with_standardized_conc(avg_combo[[id]], conc_map, "concs", "rconcs")
    avg_combo[[id2]] <- replace_conc_with_standardized_conc(avg_combo[[id2]], conc_map, "concs", "rconcs")
    for (metric in normalization_types) {
      metric_name <- switch(metric,
                            "GRvalue" = "GR",
                            "RelativeViability" = "RV",
                            "unknown")

      # fit by column: the series in the primary identifier, the cotrt is the secondary one
      col_fittings <- gDRcore:::fit_combo_cotreatments(avg_combo, series_id = id, cotrt_id = id2, metric)
      col_fittings <- col_fittings[!is.na(col_fittings$fit_type), ]

      # fit by row (flipped): the series in the secondary identifier, the cotrt is the primary one
      row_fittings <- gDRcore:::fit_combo_cotreatments(avg_combo, series_id = id2, cotrt_id = id, metric)
      row_fittings <- row_fittings[!is.na(row_fittings$fit_type), ]

      # fit by codilution (diagonal)
      codilution_fittings <- gDRcore:::fit_combo_codilutions(avg_combo, series_identifiers, metric)
      codilution_fittings <- codilution_fittings[!is.na(codilution_fittings$fit_type), ]

      # apply the fit to get smoothed data: results per column
      # (along primary identifier for each value of the secondary identifier)
      avg_combo$col_values <- map_ids_to_fits(pred = avg_combo[[id]],
                                             match_col = avg_combo[[id2]], col_fittings, "cotrt_value")
      # apply the fit the get smoothed data: results per row
      # (along secondary identifier for each value of the primary identifier)
      avg_combo$row_values <- map_ids_to_fits(pred = avg_combo[[id2]],
                                             match_col = avg_combo[[id]], row_fittings, "cotrt_value")
      metrics_names <- c("col_fittings", "row_fittings")
      if (!is.null(codilution_fittings)) {
        # apply the fit the get smoothed data: codilution results (along sum of identifiers for each ratio)
        avg_combo$codil_values <- map_ids_to_fits(pred = avg_combo[[id2]] + avg_combo[[id]],
                                                 match_col = round_concentration(avg_combo[[id2]] / avg_combo[[id]], 1),
                                                 codilution_fittings, "ratio")
        metrics_names <- c(metrics_names, "codilution_fittings")
      } 
      metrics_merged <- do.call(plyr::rbind.fill, mget(metrics_names))
      metrics_merged$fit_type <- sub("(.*)(\\..*)", "\\1", rownames(metrics_merged))

      keep <- intersect(colnames(avg_combo), c(metric, "row_values", "col_values", "codil_values"))
      mat <- as.matrix(avg_combo[, keep])
      avg_combo$average <- rowMeans(mat, na.rm = TRUE)
      
      # remove rows/columns with less than 2 values
      discard_conc_1 <- names(which(table(avg_combo[[id]][!is.na(avg_combo$average)]) < 3))
      discard_conc_2 <- names(which(table(avg_combo[[id2]][!is.na(avg_combo$average)]) < 3))
      avg_combo <- avg_combo[!(avg_combo[[id]] %in% discard_conc_1) & !(avg_combo[[id2]] %in% discard_conc_2), ]

      # just keep the relevant columns and change to the metric name
      avg <- as.data.frame(avg_combo[, c(id, id2, "average")])
      avg[avg[[id]] == 0 & avg[[id2]] == 0] <- 0
      colnames(avg)[3] <- metric

      if (NROW(avg) == 0) {
        avg <- h_excess <- b_excess <- NULL
      } else {
        sa1 <- avg[avg[[id2]] == 0, c(id, id2, metric)]
        sa2 <- avg[avg[[id]] == 0, c(id, id2, metric)]
        
        hsa <- calculate_HSA(sa1, id, sa2, id2, metric)
        h_excess <- calculate_excess(hsa, avg, series_identifiers = c(id, id2), metric_col = "metric",
                                    measured_col = metric)

        bliss <- calculate_Bliss(sa1, id, sa2, id2, metric)
        b_excess <- calculate_excess(bliss, avg, series_identifiers = c(id, id2), metric_col = "metric",
                                    measured_col = metric)
      }
      
      # call calculate_Loewe and calculate_isobolograms: 
      isobologram_out <- if (sum((avg[[id]] * avg[[id2]]) > 0) > 9 && min(row_fittings$cotrt_value) == 0) {
        calculate_Loewe(avg, row_fittings, col_fittings, codilution_fittings,
                        series_identifiers = c(id, id2), normalization_type = metric
                        ) 
      } else {
        NULL
      }

      # average the top 10-percentile excess to get a single value for the excess
      hsa_score[row, metric_name] <- ifelse(is.null(h_excess), NA, mean(
        h_excess$excess[h_excess$excess >= quantile(h_excess$excess, 0.9, na.rm = TRUE)], na.rm = TRUE)
      )
      bliss_score[row, metric_name] <- ifelse(is.null(b_excess), NA, mean(
        b_excess$excess[b_excess$excess >= quantile(b_excess$excess, 0.9, na.rm = TRUE)], na.rm = TRUE)
      )

      if (all(vapply(isobologram_out, is.null, logical(1)))) {
        CIScore_50[row, metric_name] <- CIScore_80[row, metric_name] <- NA
      } else {
        CIScore_50[row, metric_name] <- isobologram_out$df_all_AUC_log2CI$CI_100x[
          isobologram_out$df_all_AUC_log2CI$iso_level ==
            min(isobologram_out$df_all_AUC_log2CI$iso_level[isobologram_out$df_all_AUC_log2CI$iso_level >= 0.5])]
        CIScore_80[row, metric_name] <- isobologram_out$df_all_AUC_log2CI$CI_100x[
          isobologram_out$df_all_AUC_log2CI$iso_level == 
            min(isobologram_out$df_all_AUC_log2CI$iso_level[isobologram_out$df_all_AUC_log2CI$iso_level >= 0.2])]
      }
      
      b_excess$row_id <- avg$row_id <- h_excess$row_id <- isobologram_out$df_all_iso_curves$row_id <-
        col_fittings$row_id <- hsa_score[row, "row_id"] <- bliss_score[row, "row_id"] <-
        CIScore_50[row, "row_id"] <- CIScore_80[row, "row_id"] <- metrics_merged$row_id <- i
      b_excess$col_id <- avg$col_id <- h_excess$col_id <- isobologram_out$df_all_iso_curves$col_id <- 
        col_fittings$col_id <- hsa_score[row, "col_id"] <- bliss_score[row, "col_id"] <-
        CIScore_50[row, "col_id"] <- CIScore_80[row, "col_id"] <- metrics_merged$col_id <- j
      b_excess$normalization_type <- h_excess$normalization_type <- 
        isobologram_out$df_all_iso_curves$normalization_type <- metric_name
      
      hsa_excess[[row]] <- plyr::rbind.fill(hsa_excess[[row]], as.data.frame(h_excess))
      bliss_excess[[row]] <- plyr::rbind.fill(bliss_excess[[row]], as.data.frame(b_excess))
      if (!is.null(smooth_mx[[row]]) && length(smooth_mx[[row]]) != 2) { # check if it does not contain only ids
        smooth_mx[[row]] <- as.data.frame(merge(smooth_mx[[row]], avg, all = TRUE,
                                                by = c(id, id2, "row_id", "col_id")))
      } else {
        smooth_mx[[row]] <- as.data.frame(avg)
      }
      isobolograms[[row]] <- plyr::rbind.fill(isobolograms[[row]],
                                              as.data.frame(isobologram_out$df_all_iso_curves))
      metrics[[row]] <- rbind(metrics[[row]], metrics_merged)
    }
  }

  all_smooth_mx <- rbindList(smooth_mx)
  all_hsa_excess <- rbindList(hsa_excess)
  all_b_excess <- rbindList(bliss_excess)
  all_isobolograms <- rbindList(isobolograms)
  all_metrics <- rbindList(metrics)
  
  SummarizedExperiment::assays(se)[["SmoothMatrix"]] <- convertDFtoBumpyMatrixUsingIds(all_smooth_mx)
  SummarizedExperiment::assays(se)[["BlissExcess"]] <- convertDFtoBumpyMatrixUsingIds(all_b_excess)
  SummarizedExperiment::assays(se)[["HSAExcess"]] <- convertDFtoBumpyMatrixUsingIds(all_hsa_excess)
  SummarizedExperiment::assays(se)[["isobolograms"]] <- convertDFtoBumpyMatrixUsingIds(all_isobolograms)

  SummarizedExperiment::assays(se)[["BlissScore"]] <- convertDFtoBumpyMatrixUsingIds(bliss_score)
  SummarizedExperiment::assays(se)[["HSAScore"]] <- convertDFtoBumpyMatrixUsingIds(hsa_score)
  SummarizedExperiment::assays(se)[["CIScore_50"]] <- convertDFtoBumpyMatrixUsingIds(CIScore_50)
  SummarizedExperiment::assays(se)[["CIScore_80"]] <- convertDFtoBumpyMatrixUsingIds(CIScore_80)

  SummarizedExperiment::assays(se)[[metrics_assay]] <- convertDFtoBumpyMatrixUsingIds(all_metrics)
  se
}
