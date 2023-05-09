#' fit_SE for combination screens
#'
#' Perform fittings for combination screens.
#' @param se \code{SummarizedExperiment} object with a BumpyMatrix assay 
#' containing averaged data.
#' @param data_type single-agent vs combination
#' @param series_identifiers character vector of the column names in the 
#' nested \code{DFrame} corresponding to nested identifiers.
#' @param normalization_types character vector of normalization types used for 
#' calculating combo matrix.
#' @param averaged_assay string of the name of the averaged assay to use as 
#' input. in the \code{se}.
#' @param metrics_assay string of the name of the metrics assay to output
#' in the returned \linkS4class{SummarizedExperiment}.
#' whose combination represents a unique series for which to fit curves. 
#'
#' @details
#' This function assumes that the combination is set up with both 
#' concentrations nested in the assay.
#' 
#' @examples 
#' fmae_cms_path <- system.file(
#'   package = "gDRtestData", "testdata", "finalMAE_combo_matrix_small.RDS"
#' )
#' fmae_cms <- readRDS(fmae_cms_path)
#' se1 <- fmae_cms[["matrix"]]
#' SummarizedExperiment::assays(se1) <- 
#'   SummarizedExperiment::assays(se1)["Averaged"]
#' fit_SE.combinations(se1[1, 1])
#'
#' @return A code{SummarizedExperiment} object with an additional assay 
#' containing the combination metrics.
#' 
#' @export
#'
fit_SE.combinations <- function(se,
                                data_type = "matrix",
                                series_identifiers = NULL,
                                normalization_types = c("GR", "RV"),
                                averaged_assay = "Averaged",
                                metrics_assay = "Metrics"
                                ) {

  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_character(series_identifiers, null.ok = TRUE)
  checkmate::assert_character(normalization_types)
  checkmate::assert_string(averaged_assay)
  checkmate::assert_string(metrics_assay)
  
  if (is.null(series_identifiers)) {
    series_identifiers <- get_default_nested_identifiers(
      se, 
      data_model(data_type)
    )
  }
  
  if (length(series_identifiers) != 2L) {
    stop(
      "gDR only supports 'series_identifiers' arguments with length '2' 
      for the 'fit_SE.combinations' function"
    )
  }
  
  avg <- BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, averaged_assay)
  )

  rdata <- SummarizedExperiment::rowData(se)
  cdata <- SummarizedExperiment::colData(se)
  cl_name <- gDRutils::get_SE_identifiers(se, "cellline_name")

  id <- series_identifiers[1]
  id2 <- series_identifiers[2]

  iterator <- unique(avg[, c("column", "row")])

  out <- gDRutils::loop(seq_len(nrow(iterator)), function(row) {
  
    bliss_excess <- hsa_excess <- metrics <- all_iso_points <- 
      isobolograms <- smooth_mx <- NULL
    bliss_score <- hsa_score <- CIScore_50 <- CIScore_80 <- 
      S4Vectors::DataFrame(matrix(NA, 1, 0))
    x <- iterator[row, ]
    i <- x[["row"]]
    j <- x[["column"]]

    avg_combo <- avg[avg[["row"]] == i & avg[["column"]] == j, ]
   
    omit_c_idx <- which(colnames(avg_combo) %in% c("row", "column"))
    if (all(is.na(avg_combo[, -omit_c_idx]))) { # omit masked values (all NAs)
      smooth_mx <- hsa_excess <- bliss_excess <- isobolograms <- metrics <- 
        bliss_score[, c("row_id", "col_id")] <- 
        hsa_score[, c("row_id", "col_id")] <-
        CIScore_50[, c("row_id", "col_id")] <- 
        CIScore_80[, c("row_id", "col_id")] <-
        all_iso_points <-
        data.table::data.table(row_id = i, col_id = j)
      return(list(bliss_excess = bliss_excess,
           hsa_excess = hsa_excess,
           metrics = metrics,
           all_iso_points = all_iso_points,
           isobolograms = isobolograms,
           smooth_mx = smooth_mx,
           bliss_score = bliss_score,
           hsa_score = hsa_score,
           CIScore_50 = CIScore_50,
           CIScore_80 = CIScore_80))
    }
    
    conc_map <- map_conc_to_standardized_conc(avg_combo[[id]], avg_combo[[id2]])

    avg_combo[[id]] <- replace_conc_with_standardized_conc(
      avg_combo[[id]], 
      conc_map, 
      "concs", 
      "rconcs"
    )
    avg_combo[[id2]] <- replace_conc_with_standardized_conc(
      avg_combo[[id2]], 
      conc_map, 
      "concs", 
      "rconcs"
    )
    
    mean_avg_combo <- stats::aggregate(stats::as.formula(
        paste(
          "x ~ ",
          paste(c(id, id2, "normalization_type"), collapse = " + ")
        )
      ),
      function(x) mean(x, na.rm = TRUE),
      data = avg_combo, 
      na.action = stats::na.pass
    )
    # deal with cases of multiple concentrations mapped to the same value 
    # when rounded create a complete matrix with the most frequence combo 
    # concentrations
    conc1 <- table(avg_combo[[id]])
    conc1 <- sort(as.numeric(
      names(conc1)[conc1 > (max(conc1[names(conc1) != "0"]) / 2)]
    ))
    conc2 <- table(avg_combo[[id2]])
    conc2 <- sort(as.numeric(
      names(conc2)[conc2 > (max(conc2[names(conc2) != "0"]) / 2)]
    ))
    
    # create matrix with single agent
    complete <- merge(unique(c(0, conc1)), unique(c(0, conc2)), by = NULL)
    colnames(complete) <- c(id, id2)
    complete <- data.table::as.data.table(merge(complete, mean_avg_combo, all.x = TRUE, by = c(id, id2)))

    for (metric in normalization_types) {
      
      metric_name <- gDRutils::extend_normalization_type_name(metric)
      avg_combo <- data.table::as.data.table(avg_combo)
      avg_subset <- avg_combo[normalization_type == metric]
      
      # fit by column: the series in the primary identifier, the cotrt is the 
      # secondary one
      
      col_fittings <- fit_combo_cotreatments(
        avg_subset,
        series_id = id, 
        cotrt_id = id2, 
        metric
      )
      col_fittings <- if (all(is.na(col_fittings$fit_type))) {
        x <- col_fittings[1, ]
        x[, "cotrt_value"] <- NA
        x
      } else {
        col_fittings[!is.na(col_fittings$fit_type), ]
      }
    
      # fit by row (flipped): the series in the secondary identifier, the 
      # cotrt is the primary one
      row_fittings <- fit_combo_cotreatments(
        avg_subset,
        series_id = id2, 
        cotrt_id = id, 
        metric
      )
      row_fittings <- row_fittings[!is.na(row_fittings$fit_type), ]

      # fit by codilution (diagonal)
      codilution_fittings <- fit_combo_codilutions(
        avg_subset,
        series_identifiers, 
        metric
      )
      codilution_fittings <- 
        codilution_fittings[!is.na(codilution_fittings$fit_type), ]

      # apply the fit to get smoothed data: results per column
      # (along primary identifier for each value of the secondary identifier)
      complete$col_values <- map_ids_to_fits(
        pred = complete[[id]],
        match_col = complete[[id2]], 
        col_fittings, 
        "cotrt_value"
      )
      # apply the fit the get smoothed data: results per row
      # (along secondary identifier for each value of the primary identifier)
      complete$row_values <- map_ids_to_fits(
        pred = complete[[id2]],
        match_col = complete[[id]], 
        row_fittings, "cotrt_value"
      )
      metrics_names <- c("col_fittings", "row_fittings")
      if (!is.null(codilution_fittings)) {
        # apply the fit the get smoothed data: codilution results (along sum 
        # of identifiers for each ratio)
        complete$codil_values <- map_ids_to_fits(
          pred = complete[[id2]] + complete[[id]],
          match_col = round_concentration(complete[[id2]] / complete[[id]], 1),
          codilution_fittings, 
          "ratio"
        )
        metrics_names <- c(metrics_names, "codilution_fittings")
      } 
      metrics_merged <- do.call(plyr::rbind.fill, mget(metrics_names))
      # we need it to distinguish which rows are col_fittings/row_fittings 
      # and codilution_fitting
      metrics_merged$source <- rep(
        metrics_names, 
        vapply(mget(metrics_names), nrow, numeric(1))
      )
      # remove degenerated fits
      metrics_merged <-
        metrics_merged[
          !(metrics_merged$fit_type %in% 
              c("DRCInvalidFitResult", "DRCTooFewPointsToFit")), 
        ]
      if (nrow(metrics_merged) == 0) {
        metrics_merged[1, ] <- NA
      }
      keep <- intersect(
        colnames(complete), 
        c(metric, "row_values", "col_values", "codil_values")
      )
      mat <- as.matrix(complete[, ..keep])
      complete$average <- rowMeans(mat, na.rm = TRUE)
      
      # just keep the relevant columns and change to the metric name
      cols <- c(id, id2, "average")
      av_matrix <- complete[, ..cols]
      colnames(av_matrix)[3] <- metric
      av_matrix[(av_matrix[[id]] == 0) & (av_matrix[[id2]] == 0), metric] <- 1

      if (NROW(av_matrix) == 0) {
        av_matrix <- h_excess <- b_excess <- NULL
      } else {
        
        cols <- c(id, id2, metric)
        sa1 <- av_matrix[av_matrix[[id2]] == 0, ..cols]
        sa2 <- av_matrix[av_matrix[[id]] == 0, ..cols]

        hsa <- calculate_HSA(sa1, id, sa2, id2, metric)
        h_excess <- calculate_excess(
          hsa, 
          av_matrix, 
          series_identifiers = c(id, id2), 
          metric_col = "metric",
          measured_col = metric
        )

        bliss <- calculate_Bliss(sa1, id, sa2, id2, metric)
        b_excess <- calculate_excess(
          bliss, 
          av_matrix, 
          series_identifiers = c(id, id2), 
          metric_col = "metric",
          measured_col = metric
        )
      }
      
      # call calculate_Loewe and calculate_isobolograms: 
      # remove rows/columns with less than 2 values
      discard_conc_1 <- names(which(
        table(av_matrix[[id]][!is.na(av_matrix[[metric]])]) < 3
      ))
      discard_conc_2 <- names(which(
        table(av_matrix[[id2]][!is.na(av_matrix[[metric]])]) < 3
      ))
      av_matrix_dense <- av_matrix[
        !(av_matrix[[id]] %in% discard_conc_1) & 
          !(av_matrix[[id2]] %in% discard_conc_2), 
      ]
      isobologram_out <- if (
          sum((av_matrix_dense[[id]] * av_matrix_dense[[id2]]) > 0) > 9 &&
          min(row_fittings$cotrt_value) == 0
        ) {
          calculate_Loewe(
            av_matrix, 
            row_fittings, 
            col_fittings, 
            codilution_fittings,
            series_identifiers = c(id, id2), 
            normalization_type = metric
          )
      } else {
        list(df_all_iso_points = data.table::data.table(row_id = NA, col_id = NA),
             df_all_iso_curves = data.table::data.table(row_id = NA, col_id = NA))
      }

      # average the top 10-percentile excess to get a single value 
      # for the excess
      hsa_score[, metric_name] <- ifelse(
        is.null(h_excess), 
        NA, 
        mean(
          h_excess$excess[
            h_excess$excess >= 
              stats::quantile(h_excess$excess, 0.9, na.rm = TRUE)
          ], 
          na.rm = TRUE
        )
      )
      bliss_score[, metric_name] <- ifelse(
        is.null(b_excess), 
        NA, 
        mean(
          b_excess$excess[
            b_excess$excess >= 
              stats::quantile(b_excess$excess, 0.9, na.rm = TRUE)
          ], 
          na.rm = TRUE
        )
      )

      if (all(vapply(isobologram_out, is.null, logical(1)))) {
        CIScore_50[, metric_name] <- CIScore_80[, metric_name] <- NA
      } else {
        CIScore_50[, metric_name] <- isobologram_out$df_all_AUC_log2CI$CI_100x[
          isobologram_out$df_all_AUC_log2CI$iso_level ==
            min(isobologram_out$df_all_AUC_log2CI$iso_level[
              isobologram_out$df_all_AUC_log2CI$iso_level >= 0.5
            ])]
        CIScore_80[, metric_name] <- isobologram_out$df_all_AUC_log2CI$CI_100x[
          isobologram_out$df_all_AUC_log2CI$iso_level == 
            min(isobologram_out$df_all_AUC_log2CI$iso_level[
              isobologram_out$df_all_AUC_log2CI$iso_level >= 0.2
            ])]
      }
      
      b_excess$row_id <- av_matrix$row_id <- h_excess$row_id <- 
        isobologram_out$df_all_iso_points$row_id <-
        isobologram_out$df_all_iso_curves$row_id <-
        col_fittings$row_id <-
        hsa_score[, "row_id"] <- bliss_score[, "row_id"] <-
        CIScore_50[, "row_id"] <-
        CIScore_80[, "row_id"] <- metrics_merged$row_id <- i
      b_excess$col_id <-
        av_matrix$col_id <-
        h_excess$col_id <- isobologram_out$df_all_iso_points$col_id <-
        isobologram_out$df_all_iso_curves$col_id <-
        col_fittings$col_id <-
        hsa_score[, "col_id"] <- bliss_score[, "col_id"] <-
        CIScore_50[, "col_id"] <-
        CIScore_80[, "col_id"] <- metrics_merged$col_id <- j
      
      b_excess$normalization_type <- h_excess$normalization_type <- 
        isobologram_out$df_all_iso_points$normalization_type <-
        isobologram_out$df_all_iso_curves$normalization_type <- metric_name
      
      hsa_excess <- plyr::rbind.fill(hsa_excess, as.data.frame(h_excess))
      bliss_excess <- plyr::rbind.fill(bliss_excess, as.data.frame(b_excess))
      # check if it does not contain only ids
      
      if (!is.null(smooth_mx) && ncol(smooth_mx) != 2) {
        smooth_mx <-
          merge(smooth_mx,
                av_matrix,
                all = TRUE,
                by = c(id, id2, "row_id", "col_id"))
      } else {
        smooth_mx <- av_matrix
      }
    
      if (is.null(all_iso_points)) {
        all_iso_points <- data.table::as.data.table(isobologram_out$df_all_iso_points)
      } else {
        all_iso_points <- plyr::rbind.fill(
          all_iso_points,
          data.table::as.data.table(isobologram_out$df_all_iso_points)
        )
      }
      
      if (is.null(isobolograms)) {
        isobolograms <- data.table::as.data.table(isobologram_out$df_all_iso_curves)
      } else {
        isobolograms <- plyr::rbind.fill(
          isobolograms,
          data.table::as.data.table(isobologram_out$df_all_iso_curves)
        )
      }
      
      metrics <- rbind(metrics, metrics_merged)
    }
    list(bliss_excess = bliss_excess,
         hsa_excess = hsa_excess,
         metrics = metrics,
         all_iso_points = all_iso_points,
         isobolograms = isobolograms,
         smooth_mx = smooth_mx,
         bliss_score = bliss_score,
         hsa_score = hsa_score,
         CIScore_50 = CIScore_50,
         CIScore_80 = CIScore_80)
  })

  all_smooth_mx <- rbindParallelList(out, "smooth_mx")
  all_hsa_excess <- rbindParallelList(out, "hsa_excess")
  all_b_excess <- rbindParallelList(out, "bliss_excess")
  all_iso_points <- rbindParallelList(out, "all_iso_points")
  all_isobolograms <- rbindParallelList(out, "isobolograms")
  all_metrics <- rbindParallelList(out, "metrics")
  
  bliss_score <- rbindParallelList(out, "bliss_score")
  hsa_score <- rbindParallelList(out, "hsa_score")
  CIScore_50 <- rbindParallelList(out, "CIScore_50")
  CIScore_80 <- rbindParallelList(out, "CIScore_80")
  
  SummarizedExperiment::assays(se)[["SmoothMatrix"]] <- 
    convertDFtoBumpyMatrixUsingIds(all_smooth_mx)
  SummarizedExperiment::assays(se)[["BlissExcess"]] <- 
    convertDFtoBumpyMatrixUsingIds(all_b_excess)
  SummarizedExperiment::assays(se)[["HSAExcess"]] <- 
    convertDFtoBumpyMatrixUsingIds(all_hsa_excess)
  SummarizedExperiment::assays(se)[["all_iso_points"]] <- 
    convertDFtoBumpyMatrixUsingIds(all_iso_points)
  SummarizedExperiment::assays(se)[["isobolograms"]] <- 
    convertDFtoBumpyMatrixUsingIds(all_isobolograms)

  SummarizedExperiment::assays(se)[["BlissScore"]] <- 
    convertDFtoBumpyMatrixUsingIds(bliss_score)
  SummarizedExperiment::assays(se)[["HSAScore"]] <- 
    convertDFtoBumpyMatrixUsingIds(hsa_score)
  SummarizedExperiment::assays(se)[["CIScore_50"]] <- 
    convertDFtoBumpyMatrixUsingIds(CIScore_50)
  SummarizedExperiment::assays(se)[["CIScore_80"]] <- 
    convertDFtoBumpyMatrixUsingIds(CIScore_80)

  SummarizedExperiment::assays(se)[[metrics_assay]] <- 
    convertDFtoBumpyMatrixUsingIds(all_metrics)
  se
}
