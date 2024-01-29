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
#' fmae_cms <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' 
#' se1 <- fmae_cms[[gDRutils::get_experiment_groups("combination")]]
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
                                data_type = gDRutils::get_experiment_groups("combination"),
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
    print(row)
  
    metrics <- all_iso_points <- isobolograms <- excess_full <- NULL
    excess <- scores <- 
      data.table::data.table(normalization_type = normalization_types)
    x <- iterator[row, ]
    i <- x[["row"]]
    j <- x[["column"]]

    avg_combo <- data.table::as.data.table(avg[avg[["row"]] == i & avg[["column"]] == j, ])
   
    if (all(is.na(avg_combo[, -c("row", "column", "normalization_type")]))) { # omit masked values (all NAs)
      excess <- isobolograms <- metrics <- scores <- all_iso_points <-
        data.table::data.table(row_id = i, col_id = j)
      return(list(excess = excess,
           metrics = metrics,
           all_iso_points = all_iso_points,
           isobolograms = isobolograms,
           scores = scores))
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
    
    mean_avg_combo <-  avg_combo[, lapply(.SD, mean), by = c(id, id2, "normalization_type"), .SDcols = "x"]
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
    complete <- data.table::CJ(unique(c(0, conc1)), unique(c(0, conc2)))
    colnames(complete) <- c(id, id2)
    
    ### because the function below expect an exact match, some values are not passed
    ### from mean_avg_combo to complete if there is a numerical error. This leads to unnecessary NAs in compete
    ### This should never be the case, but I found it in one example with Conc = 0.7407403 -->
    ### mean_avg_combo[,Concentration] - complete_subset[,Concentration] =
    ###   0.00000000000000361 --> FALSE for the match --> NA
    
    complete <- mean_avg_combo[complete, on = c(id, id2)]
  
  
    for (norm_type in normalization_types) {
      print(norm_type)

      avg_combo <- data.table::as.data.table(avg_combo)
      avg_subset <- avg_combo[normalization_type == norm_type]
      complete_subset <- complete[normalization_type == norm_type | is.na(normalization_type)]
      complete_subset[is.na(normalization_type), normalization_type := norm_type]
      
      # fit by column: the series in the primary identifier, the cotrt is the 
      # secondary one
      
      col_fittings <- fit_combo_cotreatments(
        avg_subset,
        series_id = id, 
        cotrt_id = id2, 
        norm_type
      )
      col_fittings <- if (all(is.na(col_fittings$fit_type))) {
        x <- col_fittings[1, ]
        x[, "cotrt_value"] <- NA
        x
      } else {
        na.omit(col_fittings, col = "fit_type")
      }
    
      # fit by row (flipped): the series in the secondary identifier, the 
      # cotrt is the primary one
      row_fittings <- fit_combo_cotreatments(
        avg_subset,
        series_id = id2, 
        cotrt_id = id, 
        norm_type
      )
      
      row_fittings <- na.omit(row_fittings, col = "fit_type")

      # fit by codilution (diagonal)
      codilution_fittings <- fit_combo_codilutions(
        avg_subset,
        series_identifiers, 
        norm_type
      )
      codilution_fittings <- na.omit(codilution_fittings, col = "fit_type")

      # apply the fit to get smoothed data: results per column
      # (along primary identifier for each value of the secondary identifier)
      complete_subset$col_values <- map_ids_to_fits(
        pred = complete_subset[[id]],
        match_col = complete_subset[[id2]], 
        col_fittings, 
        "cotrt_value"
      )
      # apply the fit the get smoothed data: results per row
      # (along secondary identifier for each value of the primary identifier)
      complete_subset$row_values <- map_ids_to_fits(
        pred = complete_subset[[id2]],
        match_col = complete_subset[[id]], 
        row_fittings, "cotrt_value"
      )
      metrics_names <- c("col_fittings", "row_fittings")
      if (!is.null(codilution_fittings)) {
        # apply the fit the get smoothed data: codilution results (along sum 
        # of identifiers for each ratio)
        complete_subset$codil_values <- map_ids_to_fits(
          pred = complete_subset[[id2]] + complete_subset[[id]],
          match_col = round_concentration(complete_subset[[id2]] / complete_subset[[id]], 1),
          codilution_fittings, 
          "ratio"
        )
        metrics_names <- c(metrics_names, "codilution_fittings")
      } 

      metrics_merged <- data.table::rbindlist(mget(metrics_names), fill = TRUE)
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
        metrics_merged <- metrics_merged[, lapply(.SD, function(x) NA)]
        metrics_merged[, `:=`(normalization_type = norm_type, fit_source = "gDR")]
        metrics_merged
      }
      keep <- intersect(
        colnames(complete_subset), 
        c("x", "row_values", "col_values", "codil_values")
      )
      mat <- as.matrix(complete_subset[, keep, with = FALSE])
      # store the values from the Averaged assay for reference
      complete_subset$av_values <- complete_subset$x
      # and replace by the Smoothed values
      complete_subset$smooth <- rowMeans(mat, na.rm = TRUE)
      
      # just keep the relevant columns and change to the metric name
      cols <- c(id, id2, "smooth")
      av_matrix <- complete_subset[, cols, with = FALSE]
      av_matrix[, "normalization_type" := norm_type]
      
      av_matrix[get(id) == 0 & get(id2) == 0, smooth := 1]
      if (NROW(av_matrix) == 0) {
        excess[excess$normalization_type == norm_type,
               c("smooth", "hsa_excess", "bliss_excess")] <- NA
      } else {
        sa1 <- av_matrix[get(id2) == 0, cols, with = FALSE]
        sa2 <- av_matrix[get(id) == 0, cols, with = FALSE]

        hsa <- calculate_HSA(sa1, id, sa2, id2, norm_type)
        h_excess <- calculate_excess(
          hsa, 
          av_matrix, 
          series_identifiers = c(id, id2), 
          metric_col = "metric",
          measured_col = "smooth"
        )
        data.table::setnames(h_excess, "x", "hsa_excess")
        bliss <- calculate_Bliss(sa1, id, sa2, id2, norm_type)
        bliss_excess <- calculate_excess(
          bliss, 
          av_matrix, 
          series_identifiers = c(id, id2), 
          metric_col = "metric",
          measured_col = "smooth"
        )
        data.table::setnames(bliss_excess, "x", "bliss_excess")
        excess <- Reduce(function(x, y) merge(x, y, all = TRUE), list(av_matrix, h_excess, bliss_excess))
      }
      
      # call calculate_Loewe and calculate_isobolograms: 
      # remove rows/columns with less than 2 values
      discard_conc_1 <- names(which(
        table(av_matrix[!is.na(smooth) & normalization_type == norm_type, id, with = FALSE]) < 3
      ))
      discard_conc_2 <- names(which(
        table(av_matrix[!is.na(smooth) & normalization_type == norm_type, id2, with = FALSE]) < 3
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
            normalization_type = norm_type
          )
      } else {
        list(df_all_iso_points = data.table::data.table(row_id = NA, col_id = NA),
             df_all_iso_curves = data.table::data.table(row_id = NA, col_id = NA))
      }

      # average the top 10-percentile excess to get a single value 
      # for the excess
      scores[scores$normalization_type == norm_type, "hsa_score"] <- ifelse(
        is.null(h_excess), 
        NA, 
        mean(
          h_excess$hsa_excess[
            h_excess$hsa_excess >= 
              stats::quantile(h_excess$hsa_excess, 0.9, na.rm = TRUE)
          ], 
          na.rm = TRUE
        )
      )
      scores[scores$normalization_type == norm_type, "bliss_score"] <- ifelse(
        is.null(bliss_excess), 
        NA, 
        mean(
          bliss_excess$bliss_excess[
            bliss_excess$bliss_excess >= 
              stats::quantile(bliss_excess$bliss_excess, 0.9, na.rm = TRUE)
          ], 
          na.rm = TRUE
        )
      )

      if (all(vapply(isobologram_out, function(x) is.null(x) || all(is.na(x)), logical(1)))) {
        scores[scores$normalization_type == norm_type, "CIScore_50"] <-
          scores[scores$normalization_type == norm_type, "CIScore_80"] <- NA
      } else {
        scores[scores$normalization_type == norm_type, "CIScore_50"] <-
          unique(isobologram_out$df_all_AUC_log2CI$CI_100x[
          isobologram_out$df_all_AUC_log2CI$iso_level ==
            min(isobologram_out$df_all_AUC_log2CI$iso_level[
              isobologram_out$df_all_AUC_log2CI$iso_level >= 0.5
            ])])
        scores[scores$normalization_type == norm_type, "CIScore_80"] <-
          unique(isobologram_out$df_all_AUC_log2CI$CI_100x[
          isobologram_out$df_all_AUC_log2CI$iso_level == 
            min(isobologram_out$df_all_AUC_log2CI$iso_level[
              isobologram_out$df_all_AUC_log2CI$iso_level >= 0.2
            ])])
      }
      
      excess$row_id <-
        isobologram_out$df_all_iso_points$row_id <-
        isobologram_out$df_all_iso_curves$row_id <-
        col_fittings$row_id <-
        scores$row_id <-
        metrics_merged$row_id <- i
      excess$col_id <-
        isobologram_out$df_all_iso_points$col_id <-
        isobologram_out$df_all_iso_curves$col_id <-
        col_fittings$col_id <-
        scores$col_id <-
        metrics_merged$col_id <- j
      
      isobologram_out$df_all_iso_points$normalization_type <-
      isobologram_out$df_all_iso_curves$normalization_type <- norm_type
      
      # check if it does not contain only ids
      if (is.null(all_iso_points)) {
        all_iso_points <- isobologram_out$df_all_iso_points
      } else {
        all_iso_points <- data.table::rbindlist(list(
          all_iso_points,
          isobologram_out$df_all_iso_points
        ), fill = TRUE)
      }
      
      if (is.null(isobolograms)) {
        isobolograms <- isobologram_out$df_all_iso_curves
      } else {
        isobolograms <- data.table::rbindlist(list(
          isobolograms,
          isobologram_out$df_all_iso_curves
        ), fill = TRUE)
      }
      
      metrics <- data.table::rbindlist(list(metrics, metrics_merged), fill = TRUE)
      excess_full <- data.table::rbindlist(list(excess_full, excess), fill = TRUE)
    }
    list(metrics = metrics,
         all_iso_points = all_iso_points,
         isobolograms = isobolograms,
         excess = excess_full,
         scores = scores)
  })

  excess <- rbindParallelList(out, "excess")
  all_iso_points <- rbindParallelList(out, "all_iso_points")
  all_isobolograms <- rbindParallelList(out, "isobolograms")
  all_metrics <- rbindParallelList(out, "metrics")
  scores <- rbindParallelList(out, "scores")

  
  SummarizedExperiment::assays(se)[["excess"]] <- 
    convertDFtoBumpyMatrixUsingIds(excess)
  SummarizedExperiment::assays(se)[["all_iso_points"]] <- 
    convertDFtoBumpyMatrixUsingIds(all_iso_points)
  SummarizedExperiment::assays(se)[["isobolograms"]] <- 
    convertDFtoBumpyMatrixUsingIds(all_isobolograms)

  SummarizedExperiment::assays(se)[["scores"]] <- 
    convertDFtoBumpyMatrixUsingIds(scores)

  SummarizedExperiment::assays(se)[[metrics_assay]] <- 
    convertDFtoBumpyMatrixUsingIds(all_metrics)
  se
}
