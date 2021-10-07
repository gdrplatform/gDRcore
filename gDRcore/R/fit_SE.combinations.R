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
  checkmate::assert_character(normalization_types)
  
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

    sa <- avg_combo[[id]] == 0 | avg_combo[[id2]] == 0
    single_agent <- avg_combo[sa, ]
    measured <- avg_combo#[!sa, ] # we can keep the single agent to have a full matrix

    for (metric in normalization_types) {
      sa1 <- single_agent[single_agent[[id]] == 0, c(id, id2, metric)]
      sa2 <- single_agent[single_agent[[id2]] == 0, c(id, id2, metric)]

      # fit by column: the series in the primary identifier, the cotrt is the secondary one
      col_fittings <- gDRcore:::fit_combo_cotreatments(avg_combo, series_id = id, cotrt_id = id2, metric)
      # fit by row (flipped): the series in the secondary identifier, the cotrt is the primary one
      row_fittings <- gDRcore:::fit_combo_cotreatments(avg_combo, series_id = id2, cotrt_id = id, metric)
      codilution_fittings <- gDRcore:::fit_combo_codilutions(measured, series_identifiers, metric)
      
      # apply the fit the get smoothed data: results per column
      # (along primary identifier for each value of the secondary identifier)
      measured$col_values <- map_ids_to_fits(pred = measured[[id]],
                                             match_col = measured[[id2]], col_fittings, "cotrt_value")
      # apply the fit the get smoothed data: results per row
      # (along secondary identifier for each value of the primary identifier)
      measured$row_values <- map_ids_to_fits(pred = measured[[id2]],
                                             match_col = measured[[id]], row_fittings, "cotrt_value")
      if (!is.null(codilution_fittings)) {
        # apply the fit the get smoothed data: codilution results (along sum of identifiers for each ratio)
        measured$codil_values <- map_ids_to_fits(pred = measured[[id2]] + measured[[id]],
                                                 match_col = measured[[id2]] / measured[[id]],
                                                 codilution_fittings, "ratio")
        data.table::setnames(codilution_fittings, "ratio", "cotrt_value")
        metrics_names <- c("col_fittings", "row_fittings", "codilution_fittings")
        metrics_merged <- do.call(rbind, mget(metrics_names))
        metrics_merged$fit_type <- sub("(.*)(\\..*)", "\\1", rownames(metrics_merged))
        data.table::setnames(codilution_fittings, "cotrt_value", "ratio")
      } else {
        metrics_names <- c("col_fittings", "row_fittings")
        metrics_merged <- do.call(rbind, mget(metrics_names))
        metrics_merged$fit_type <- sub("(.*)(\\..*)", "\\1", rownames(metrics_merged))
        }
    
      keep <- intersect(colnames(measured), c(metric, "row_values", "col_values", "codil_values"))
      mat <- as.matrix(measured[, keep])
      measured$average <- rowMeans(mat, na.rm = TRUE)
      
      hsa <- calculate_HSA(sa1, id2, sa2, id, metric)
      h_excess <- calculate_excess(hsa, measured, series_identifiers = c(id, id2), metric_col = "metric",
                                  measured_col = "average")

      bliss <- calculate_Bliss(sa1, id2, sa2, id, metric)
      b_excess <- calculate_excess(bliss, measured, series_identifiers = c(id, id2), metric_col = "metric",
                                  measured_col = "average")

      # TODO: call calculate_Loewe and calculate_isoobolograms
      ## TODO: Create another assay in here with each spot in the matrix as the 2 series_identifier concentrations
      ## and then each new metric that should go for each spot is another column in the nested DataFrame
      
      metric_name <- switch(metric,
                            "GRvalue" = "GR",
                            "RelativeViability" = "RV",
                            "unknown")
      hsa_score[row, metric_name] <- mean(
        h_excess$excess[h_excess$excess <= quantile(h_excess$excess, 0.1, na.rm = TRUE)])
      bliss_score[row, metric_name] <- mean(
        b_excess$excess[b_excess$excess <= quantile(b_excess$excess, 0.1, na.rm = TRUE)])
      

      # contruct full matrix with single agent
      mean_matrix <- reshape2::acast(as.data.frame(measured[, c("average", id, id2)]),
                                     formula = sprintf("%s ~ %s", id, id2), value.var = "average")
      mean_matrix["0", "0"] <- 1 # add the corner of the matrix
      # remove empty rows/columns
      mean_matrix <- mean_matrix[rowSums(!is.na(mean_matrix)) > 2, colSums(!is.na(mean_matrix)) > 2]

      # TO DO : measured$average or mean_matrix should be saved for further plots in some manner
      mean_df <- reshape2::melt(mean_matrix, varnames = c(id, id2), value.name = metric)
      if (nrow(mean_df) == 0) {
        mean_df <- NULL
      }
      # call calculate_Loewe and calculate_isobolograms: 
      isobologram_out <- if (ncol(mean_matrix) > 3) {
        calculate_Loewe(mean_matrix, row_fittings, col_fittings,
                        codilution_fittings, normalization_type = metric) 
      } else {
          NULL
        }
      ## TODO: Create another assay in here with each spot in the matrix as the 2 series_identifier concentrations
      ## and then each new metric that should go for each spot is another column in the nested DataFrame
      hsa_score[row, metric_name] <- mean(
        h_excess$excess[h_excess$excess <= quantile(h_excess$excess, 0.1, na.rm = TRUE)], na.rm = TRUE)
      bliss_score[row, metric_name] <- mean(
        b_excess$excess[b_excess$excess <= quantile(b_excess$excess, 0.1, na.rm = TRUE)], na.rm = TRUE)
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
      
      b_excess$row_id <- mean_df$row_id <- h_excess$row_id <- isobologram_out$df_all_iso_curves$row_id <- 
        col_fittings$row_id <- hsa_score[row, "row_id"] <- bliss_score[row, "row_id"] <-
        CIScore_50[row, "row_id"] <- CIScore_80[row, "row_id"] <- metrics_merged$row_id <- i
      b_excess$col_id <- mean_df$col_id <- h_excess$col_id <- isobologram_out$df_all_iso_curves$col_id <- 
        col_fittings$col_id <- hsa_score[row, "col_id"] <- bliss_score[row, "col_id"] <-
        CIScore_50[row, "col_id"] <- CIScore_80[row, "col_id"] <- metrics_merged$col_id <- j
      b_excess$normalization_type <- h_excess$normalization_type <-
        isobologram_out$df_all_iso_curves$normalization_type <- metric_name
      
      hsa_excess[[row]] <- rbind(hsa_excess[[row]], h_excess)
      bliss_excess[[row]] <- rbind(bliss_excess[[row]], b_excess)
      if (!is.null(smooth_mx[[row]])) {
        smooth_mx[[row]] <- merge(smooth_mx[[row]], mean_df, all = TRUE)
      } else {
        smooth_mx[[row]] <- mean_df
      }
      isobolograms[[row]] <- rbind(isobolograms[[row]], as.data.frame(isobologram_out$df_all_iso_curves))
      metrics[[row]] <- rbind(metrics[[row]], metrics_merged)
    }
  }


  all_hsa_excess <- S4Vectors::DataFrame(do.call("rbind", hsa_excess))
  all_b_excess <- S4Vectors::DataFrame(do.call("rbind", bliss_excess))
  all_smooth_mx <- S4Vectors::DataFrame(do.call("rbind", smooth_mx))
  all_isobolograms <- S4Vectors::DataFrame(do.call("rbind", isobolograms))
  all_CIScore_50 <- S4Vectors::DataFrame(do.call("rbind", CIScore_50))
  all_CIScore_80 <- S4Vectors::DataFrame(do.call("rbind", CIScore_80))
  all_metrics <- S4Vectors::DataFrame(do.call("rbind", metrics))

  SummarizedExperiment::assays(se)[["SmoothMatrix"]] <- BumpyMatrix::splitAsBumpyMatrix(
    all_smooth_mx[!colnames(all_smooth_mx) %in% c("row_id", "col_id")],
    row = all_smooth_mx$row_id, col = all_smooth_mx$col_id)
  SummarizedExperiment::assays(se)[["BlissExcess"]] <- BumpyMatrix::splitAsBumpyMatrix(
    all_b_excess[!colnames(all_b_excess) %in% c("row_id", "col_id")],
    row = all_b_excess$row_id, col = all_b_excess$col_id)
  SummarizedExperiment::assays(se)[["HSAExcess"]] <- BumpyMatrix::splitAsBumpyMatrix(
    all_hsa_excess[!colnames(all_hsa_excess) %in% c("row_id", "col_id")],
    row = all_hsa_excess$row_id, col = all_hsa_excess$col_id)
  SummarizedExperiment::assays(se)[["isobolograms"]] <- BumpyMatrix::splitAsBumpyMatrix(
    all_isobolograms[!colnames(all_isobolograms) %in% c("row_id", "col_id")],
    row = all_isobolograms$row_id, col = all_isobolograms$col_id)
  SummarizedExperiment::assays(se)[["BlissScore"]] <- BumpyMatrix::splitAsBumpyMatrix(
    bliss_score[!colnames(bliss_score) %in% c("row_id", "col_id")],
    row = bliss_score$row_id, col = bliss_score$col_id)
  
  SummarizedExperiment::assays(se)[["HSAScore"]] <- BumpyMatrix::splitAsBumpyMatrix(
    hsa_score[!colnames(hsa_score) %in% c("row_id", "col_id")],
    row = hsa_score$row_id, col = hsa_score$col_id)
  SummarizedExperiment::assays(se)[["CIScore_50"]] <- BumpyMatrix::splitAsBumpyMatrix(
    CIScore_50[!colnames(CIScore_50) %in% c("row_id", "col_id")],
    row = CIScore_50$row_id, col = CIScore_50$col_id)
  SummarizedExperiment::assays(se)[["CIScore_80"]] <- BumpyMatrix::splitAsBumpyMatrix(
    CIScore_80[!colnames(CIScore_80) %in% c("row_id", "col_id")],
    row = CIScore_80$row_id, col = CIScore_80$col_id)
  
  SummarizedExperiment::assays(se)[[metrics_assay]] <- BumpyMatrix::splitAsBumpyMatrix(
    all_metrics[!colnames(all_metrics) %in% c("row_id", "col_id")],
    row = all_metrics$row_id, col = all_metrics$col_id)

  se
}
  

#' @export
calculate_combo_matrix <- function(se,
                                   series_identifiers,
                                   conc_margin = 10 ^ 0.5,
                                   log2_pos_offset = log10(3) / 2,
                                   norm_types = c("RelativeViability", "GRvalue"),
                                   averaged_assay = "Averaged") {
  .Deprecated(new = "fit_SE.combinations")
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_number(conc_margin)
  checkmate::assert_number(log2_pos_offset)
  checkmate::assert_character(norm_types)
  
  agg_results_norm <- list()
  all_combo_variables_norm <- list()
  # run through both metrics
  for (norm_method in norm_types) {
    
    # create empty matrices to store the metrics across all drug combo and cell lines in the SE
    agg_results <- list(
      bliss_q10 = matrix(NA, length(S4Vectors::metadata(se)$drug_combinations), 
      ncol(se),
      dimnames = list(vapply(S4Vectors::metadata(se)$drug_combinations,
                function(x) gsub("Duration", "T", x$name), character(1)),
                paste(SummarizedExperiment::colData(se)$CellLineName, SummarizedExperiment::colData(se)$clid)))
    )
    agg_results$hsa_q10 <- agg_results$CI_100x_50 <- agg_results$CI_100x_80 <- agg_results$bliss_q10

    # create empty matrix of lists for storing all results necessary for plotting
    all_combo_variables <- matrix(list(), length(S4Vectors::metadata(se)$drug_combinations), ncol(se),
      dimnames = list(vapply(S4Vectors::metadata(se)$drug_combinations,
                function(x) gsub("Duration", "T", x$name), FUN.VALUE = character(1)),
                paste(SummarizedExperiment::colData(se)$CellLineName, SummarizedExperiment::colData(se)$clid))
    )
    # run through all drug combinations
    for (idc in seq_len(length(S4Vectors::metadata(se)$drug_combinations))) {
      combo <- S4Vectors::metadata(se)$drug_combinations[[idc]]

      # run through all cell lines
      for (iCL in seq_len(ncol(se))) {
        # get all data as flat data.frame
        
        flat_data <- gDRutils::convert_se_assay_to_dt(se[combo$rows, iCL], averaged_assay)


        # skip if not matrix combo data or too few measured values
        if (nrow(flat_data) == 0L ||
          !all(combo$condition[c("DrugName", "DrugName_2")] %in%
            unique(c(flat_data$DrugName, flat_data$DrugName_2))) ||
            sum(!is.na(flat_data[, ..norm_method])) < 10) {
          print(sprintf("Skipping %s for %s x %s", SummarizedExperiment::colData(se)$CellLineName[iCL],
              combo$condition["DrugName"], combo$condition["DrugName_2"]))
          next
        }

        # avoid mismatch because of numerical rounding
        flat_data$Concentration <- 10 ^ (.25 * round(4 * log10(flat_data$Concentration), 1))
        flat_data$Concentration_2 <- 10 ^ (.25 * round(4 * log10(flat_data$Concentration_2), 1))

        print(sprintf("Calculation for cell line %s treated with %s x %s", 
                      SummarizedExperiment::colData(se)$CellLineName[iCL], 
                      combo$condition["DrugName"], combo$condition["DrugName_2"]))

        # Secondary drug for some combinations becomes primary drug when viewed as single-agent (conc1=0), so swap. 
        # swap the data into the Gnumber and Gnumber_2 such that it can form a matrix
        swap_idx <- flat_data[, gDRutils::get_env_identifiers("drugname")] == 
          combo$condition[[paste0(gDRutils::get_env_identifiers("drugname"), "_2")]] & 
          flat_data$Gnumber_2 %in% gDRutils::get_env_identifiers("untreated_tag")
        temp_df <- flat_data[swap_idx, ]
        temp_df[, c("Gnumber", "DrugName", "drug_moa", "Concentration",
                    "Gnumber_2", "DrugName_2", "drug_moa_2", "Concentration_2")] <- 
          temp_df[, c("Gnumber_2", "DrugName_2", "drug_moa_2", "Concentration_2", 
                      "Gnumber", "DrugName", "drug_moa", "Concentration")]

        # Append the bottom left entry for full matrix creation. 
        temp_df2 <- flat_data[flat_data$Concentration == 0 & flat_data$Concentration_2 == 0, ][1, ]
        temp_df2$GRvalue <- 1
        temp_df2$std_GRvalue <- 0
        temp_df2$RelativeViability <- 1
        temp_df2$std_RelativeViability <- 0
        temp_df2$Concentration <- 0
        temp_df2$Concentration_2 <- 0

        # Then, do things normally downstream. 

        flat_data <- rbind(temp_df2, flat_data[!swap_idx, ], temp_df)
  
        # check the there is only one drug for DrugName and one for DrugName_2
        stopifnot(all(unique(flat_data$DrugName) %in% c(combo$condition["DrugName"],
              gDRutils::get_env_identifiers()$untreated_tag)))
        stopifnot(all(unique(flat_data$DrugName_2) %in% c(combo$condition["DrugName_2"],
              gDRutils::get_env_identifiers()$untreated_tag)))

        # check what are the combo concentrations
        conc_combo <- unique(flat_data[flat_data$Concentration > 0 & flat_data$Concentration_2 > 0, 
                            c("Concentration", "Concentration_2")])

        
        # Get the single agent fit for drug 1
        df_ <- flat_data[flat_data$Concentration_2 %in% 0 & flat_data$Concentration > 0, ]
        fit_drug1 <- gDRutils::fit_curves(
          df_ = df_[!is.na(df_[, norm_method]), ],
          series_identifiers = series_identifiers,
          force_fit = TRUE,
          cap = 0.2,
          normalization_type = ifelse(norm_method == "RelativeViability", "RV", "GR")
        )

        # check if missing concentration for the single agent for drug 1 and infer if needed
        if (!all(conc_combo$Concentration %in% flat_data$Concentration[flat_data$Concentration_2 %in% 0])) {
          df_inf <- data.frame(
            Concentration = sort(setdiff(conc_combo$Concentration, 
                                         flat_data$Concentration[flat_data$Concentration_2 %in% 0])), 
            val = 0)
          df_inf$val <- gDRutils::logistic_4parameters(
            df_inf[, 1], 
            fit_drug1$x_inf, 
            fit_drug1$x_0, 
            fit_drug1$ec50, fit_drug1$h)
          colnames(df_inf)[2] <- norm_method
          
          df_out <- as.data.frame(df_[df_$Concentration_2 == 0, ])[seq_len(nrow(df_inf)), ] # get the metadata
          df_out[, colnames(df_inf)] <- df_inf
          df_out[, grepl("std_", colnames(df_out)) | 
            colnames(df_out) %in% c("CorrectedReadout", setdiff(c("RelativeViability", "GRvalue"), norm_method))] <- 0
          flat_data <- rbind(flat_data, df_out)
        }
        
        # Get the single agent fit for drug 2
        df_ <- flat_data[flat_data$Concentration %in% 0 & flat_data$Concentration_2 > 0, ]
        df_$Concentration <- df_$Concentration_2 # necessary for the fit
        fit_drug2 <- gDRutils::fit_curves(
          df_ = df_[!is.na(df_[, norm_method]), ],
          series_identifiers = series_identifiers,
          normalization_type = ifelse(norm_method == "RelativeViability", "RV", "GR"),
          force_fit = TRUE,
          cap = 0.2
        )

        # check if missing concentration for the single agent for drug 2 and infer if needed
        if (!all(conc_combo$Concentration_2 %in% flat_data$Concentration_2[flat_data$Concentration %in% 0])) {
          df_inf <- data.frame(
            Concentration_2 = sort(setdiff(conc_combo$Concentration_2, 
                                           flat_data$Concentration_2[flat_data$Concentration %in% 0])),
            val = 0)
          df_inf$val <- gDRutils::logistic_4parameters(
            df_inf[, 1], 
            fit_drug2$x_inf, 
            fit_drug2$x_0, 
            fit_drug2$ec50, 
            fit_drug2$h)
          colnames(df_inf)[2] <- norm_method
          
          df_out <- as.data.frame(df_[df_$Concentration == 0, ])[seq_len(nrow(df_inf)), ]
          df_out[, colnames(df_inf)] <- df_inf
          df_out[, grepl("std_", colnames(df_out)) | 
            colnames(df_out) %in% c("CorrectedReadout", setdiff(c("RelativeViability", "GRvalue"), norm_method))] <- 0
          flat_data <- rbind(flat_data, df_out)
        }

        # create the combo matrix for the two drugs with measured data (rows as concentration of drug_1 
        # and columns as drug_2)
        mx_response <- reshape2::acast(flat_data, 
          Concentration ~ Concentration_2,
          value.var = norm_method, 
          fun.aggregate = function(x) mean(x, na.rm = TRUE))
        mx_response[1, 1] <- 1
        # remove empty rows/columns
        mx_response <- mx_response[rowSums(!is.na(mx_response)) > 2, colSums(!is.na(mx_response)) > 2]

        # drug_1 is diluted along the rows and will be the y-axis of the matrix plots
        drug1_axis <- data.frame(conc_1 = round(as.numeric(rownames(mx_response)), 5),
          log10conc1 = 0, pos_y = 0, marks_y = 0)
        drug1_axis$log10conc1 <- log10(drug1_axis$conc)
        drug1_axis$pos_y <- drug1_axis$log10conc1
        drug1_axis$pos_y[1] <- 2 * drug1_axis$pos_y[2] - drug1_axis$pos_y[3] - log10(1.5)
        drug1_axis$marks_y <- sprintf("%.2g", drug1_axis$conc_1)

        # drug_2 is diluted along the columns and will be the x-axis of the matrix plots
        drug2_axis <- data.frame(conc_2 = round(as.numeric(colnames(mx_response)), 5),
                      log10conc2 = 0, pos_x = 0, marks_x = 0)
        drug2_axis$log10conc2 <- log10(drug2_axis$conc_2)
        drug2_axis$pos_x <- drug2_axis$log10conc2
        drug2_axis$pos_x[1] <- 2 * drug2_axis$pos_x[2] - drug2_axis$pos_x[3] - log10(1.5)
        drug2_axis$marks_x <- sprintf("%.2g", drug2_axis$conc_2)


        mx_full <- mx_response # mx_full is the matrix with imputed data for NAs in mx_response
        # impute missing values
        idx_na <- which(is.na(mx_response), arr.ind = TRUE)
        if (nrow(idx_na) > 0) {
          for (i in seq_len(nrow(idx_na))) {
            r_idx <- idx_na[i, "row"] + (-1:1)
            r_idx <- pmin(pmax(r_idx, 1), nrow(mx_response))
            c_idx <- idx_na[i, "col"] + (-1:1)
            c_idx <- pmin(pmax(c_idx, 1), ncol(mx_response))
            close_mx <- mx_response[r_idx, c_idx]
            mx_full[idx_na[i, "row"], idx_na[i, "col"]] <- mean(close_mx, na.rm = TRUE) # use a smoothing strategy
          }
        }

        # function to fit the data by rows or columns and replace measured values by fitted values
        replace_by_fit <- function(mx, idx, by_row = TRUE) {
          if (!by_row) {
            mx <- t(mx)
          }
          df_ <- data.frame(as.numeric(colnames(mx))[-1], mx[idx, -1])
          colnames(df_) <- c("Concentration", norm_method)


          fit_res <- gDRutils::fit_curves(
            df_,
            series_identifiers = series_identifiers,
            e_0 = mx[idx, 1], # use single agent fit
            GR_0 = mx[idx, 1], # use single agent fit
            normalization_type = ifelse(norm_method == "RelativeViability", "RV", "GR"),
            force_fit = TRUE,
            cap = 0.2
          )

          # if failed and values are increasing, fit a reverse sigmoidal
          if (is.na(fit_res$x_0) && stats::median(diff(mx[idx, -1])) > 0 && mean(diff(mx[idx, -1])) > 0) {
              df_[, norm_method] <- 1 - df_[, norm_method]
            fit_res <- gDRutils::fit_curves(
              df_,
              series_identifiers = series_identifiers,
              e_0 = 1 - mx[idx, 1], # use single agent fit
              GR_0 = 1 - mx[idx, 1], # use single agent fit
              normalization_type = ifelse(norm_method == "RelativeViability", "RV", "GR"),
              force_fit = TRUE,
              cap = 0.2
            )
            fit_res[grepl("x_", names(fit_res))] <- 1 - as.matrix(fit_res[grepl("x_", names(fit_res))])
          }

          mx[idx, -1] <- gDRutils::logistic_4parameters(as.numeric(colnames(mx))[-1],
            fit_res$x_inf, fit_res$x_0, fit_res$ec50, fit_res$h)
          if (!by_row) {
            mx <- t(mx)
          }
          return(list(mx = mx, fit_res = fit_res))
        }

        mx_fit <- mx_full # mx_fit are three matrices with fitted data (row, columns, diagonals) based on mx_full
        
        # get the fits for the first row (single agent) and create empty matrices
        mx_fit[1, -1] <- gDRutils::logistic_4parameters(as.numeric(colnames(mx_fit))[-1],
            fit_drug2$x_inf, fit_drug2$x_0, fit_drug2$ec50, fit_drug2$h)
        all_fits <- list(by_row = cbind(data.frame(conc_1 = 0), fit_drug2))
        
        # get the fits for the first column (single agent)
        mx_fit[-1, 1] <- gDRutils::logistic_4parameters(as.numeric(rownames(mx_fit))[-1],
            fit_drug1$x_inf, fit_drug1$x_0, fit_drug1$ec50, fit_drug1$h)
        all_fits[["by_col"]] <- cbind(data.frame(conc_2 = 0), fit_drug1)
        # matrices with first row/column populated
        mx_fit <- lapply(1:3, function(x) mx_fit)
        names(mx_fit) <- c("by_row", "by_col", "by_codil")

        # replace with fits by rows
        for (i in 2:nrow(mx_response)) {  # fit by rows
          out <- replace_by_fit(mx_fit[["by_row"]], i)
          mx_fit[["by_row"]] <- out$mx
          # add fitted data 
          all_fits$by_row <- rbind(all_fits$by_row,
                cbind(data.frame(conc_1 = as.numeric(rownames(mx_fit[["by_row"]])[i])),
            out$fit_res))
        }
        # replace with fits by columns
        for (i in 2:ncol(mx_response)) {  # fit by columns
          out <- replace_by_fit(mx_fit[["by_col"]], i, FALSE)
          mx_fit[["by_col"]] <- out$mx
          # add fitted data 
          all_fits$by_col <- rbind(all_fits$by_col,
                cbind(data.frame(conc_2 = as.numeric(colnames(mx_fit[["by_col"]])[i])),
            out$fit_res))
        }

        # replace with fits by diagonals (co-dilutions with different ratios)
        all_fits[["by_codil"]] <- data.frame()

        # if not nrow != ncol; start from the highest concentration
        n_conc_codil <- min(ncol(mx_response), nrow(mx_response)) - 1
        conc_ratio <- 2 ^ unique(round(log2(
          as.numeric(colnames(mx_response))[ncol(mx_response):(ncol(mx_response) - n_conc_codil + 1)] / 
            as.numeric(rownames(mx_response))[nrow(mx_response):(nrow(mx_response) - n_conc_codil + 1)]), 1))

        if (length(conc_ratio) != 1 || n_conc_codil < 4) { # check that ratios are respected
          mx_fit[["by_codil"]] <- NA
        } else {
          for (i in -(nrow(mx_response) - 1):(ncol(mx_response) - 1)) {
              # get the data for a given diagonal (offset between the concentration and concentration_2)
              idx <- data.frame(col_idx = seq_len(nrow(mx_response)) - i,
                            row_idx = seq_len(nrow(mx_response)))
              idx <- idx[idx$col_idx > 1 & idx$col_idx <= ncol(mx_response) &
                        idx$row_idx > 1 & idx$row_idx <= nrow(mx_response), ]

              if (nrow(idx) < 4) {
                for (j in seq_len(nrow(idx))) {
                  mx_fit[["by_codil"]][idx$row_idx[j], idx$col_idx[j]] <- NA
                }
              } else {
                # fit if enough data points
                resp_value <- array(0, nrow(idx))
                for (j in seq_len(nrow(idx))) {
                  resp_value[j] <- mx_fit[["by_codil"]][idx$row_idx[j], idx$col_idx[j]]
                }

                df_ <- data.frame(as.numeric(rownames(mx_fit[["by_codil"]]))[idx$row_idx] +
                    as.numeric(colnames(mx_fit[["by_codil"]]))[idx$col_idx],
                    resp_value)
                colnames(df_) <- c("Concentration", norm_method)

                fit_res <- gDRutils::fit_curves(
                  df_,
                  series_identifiers = series_identifiers,
                  e_0 = 1,
                  GR_0 = 1,
                  normalization_type = ifelse(norm_method == "RelativeViability", "RV", "GR"),
                  force_fit = TRUE
                )

                fit_resp <- gDRutils::logistic_4parameters(
                  as.numeric(rownames(mx_fit[["by_codil"]]))[idx$row_idx] +
                      as.numeric(colnames(mx_fit[["by_codil"]]))[idx$col_idx],
                  fit_res$x_inf, fit_res$x_0, fit_res$ec50, fit_res$h)
                for (j in seq_len(nrow(idx))) {
                  mx_fit[["by_codil"]][idx$row_idx[j], idx$col_idx[j]] <- fit_resp[j]
                }
                # add fitted data 
                all_fits$by_codil <- rbind(all_fits$by_codil, cbind(data.frame(
                  conc_ratio = as.numeric(rownames(mx_fit[["by_codil"]])[idx$row_idx[1]]) / 
                        as.numeric(colnames(mx_fit[["by_codil"]])[idx$col_idx[1]])),
                      fit_res))
            }
          }
        }
        # store all matrices
        all_mx <- c(list(mx_response = mx_response, mx_full = mx_full), mx_fit)
        # calculate the mean matrix (measured, imputed, and row/column/diagonal fitted data)
        mx_mean <- matrix(sapply(seq_len(length(mx_response)), function(i) {
          mean(sapply(seq_len(length(all_mx)), function(j) all_mx[[j]][i]), na.rm = TRUE)}),
          nrow(mx_response), ncol(mx_response),
          dimnames = list(rownames(mx_response), colnames(mx_response)))
        all_mx <- c(all_mx, list(mx_mean = mx_mean))

        # calculated excess over single agent (HSA model)
        hsa_mx <- all_mx$mx_mean
        hsa_mx[-1, -1] <- pmin(t(matrix(t(hsa_mx[1, -1]), ncol(hsa_mx) - 1, nrow(hsa_mx) - 1)),
                            matrix(hsa_mx[-1, 1], nrow(hsa_mx) - 1, ncol(hsa_mx) - 1))
        hsa_excess <- mx_mean - hsa_mx

        # calculated excess over independence (Bliss model)
        bliss_mx <- all_mx$mx_mean
        bliss_mx[-1, -1] <- t(matrix(t(bliss_mx[1, -1]), ncol(bliss_mx) - 1, nrow(bliss_mx) - 1)) *
                            matrix(bliss_mx[-1, 1], nrow(bliss_mx) - 1, ncol(bliss_mx) - 1)
        bliss_excess <- mx_mean - bliss_mx

        # add combo models and excess for HSA and bliss
        all_mx <- c(all_mx, list(hsa_mx = hsa_mx, hsa_excess = hsa_excess,
                    bliss_mx = bliss_mx, bliss_excess = bliss_excess
            ))

        # decide on the isobologram cutoff levels for the Loewe model 
        if (min(mx_mean, na.rm = TRUE) > 0.7) {
          iso_cutoff <- NULL
        } else {
          if (norm_method == "GRvalue") {
            iso_cutoff <- seq(max(-0.25, ceiling(20 * min(mx_mean + .08, na.rm = TRUE)) / 20), 0.8,  0.05)
          } else {
            iso_cutoff <- seq(max(0.2, ceiling(20 * min(mx_mean + .08, na.rm = TRUE)) / 20), 0.8,  0.05)
          }
          names(iso_cutoff) <- as.character(iso_cutoff)
        }

        # create the variable for the different isobolograms
        all_iso <- vector("list", length(iso_cutoff))
        ref_x50 <- c(conc_1 = min(all_fits[["by_col"]][1, "xc50"], max(drug1_axis$conc_1) * conc_margin),
                    conc_2 = min(all_fits[["by_row"]][1, "xc50"], max(drug2_axis$conc_2) * conc_margin))
        names(all_iso) <- iso_cutoff
        for (isobol_value in iso_cutoff) { # run through the different isobolograms
          # cutoff point by row
          df_fit <- all_fits[["by_row"]]
          df_fit <- df_fit[nrow(df_fit):1, ]
          df_iso <- cbind(df_fit[, "conc_1", drop = FALSE], data.frame(conc_2 =
            ifelse(df_fit$x_0 < isobol_value, 0, ifelse(df_fit$x_inf > isobol_value,
              Inf,
              df_fit$ec50 * ((((df_fit$x_0 - df_fit$x_inf) / (isobol_value - df_fit$x_inf)) - 1) ^
                  (1 / pmax(df_fit$h, 0.01))))
                ),
                fit_type = "by_row"))

          # cutoff point by columns
          df_fit <- all_fits[["by_col"]]
          df_iso <- rbind(df_iso, cbind(data.frame(conc_1 =
            ifelse(df_fit$x_0 < isobol_value, 0, ifelse(df_fit$x_inf > isobol_value,
              Inf,
              df_fit$ec50 * ((((df_fit$x_0 - df_fit$x_inf) / (isobol_value - df_fit$x_inf)) - 1) ^
                  (1 / pmax(df_fit$h, 0.01))))
                ),
            df_fit[, "conc_2", drop = FALSE],
              fit_type = "by_col")))


          df_iso$conc_1 <- pmin(df_iso$conc_1, max(drug1_axis$conc_1) * conc_margin)
          df_iso$conc_2 <- pmin(df_iso$conc_2, max(drug2_axis$conc_2) * conc_margin)

          ref_conc_1 <- pmin(df_iso$conc_1[df_iso$conc_2 == 0 & df_iso$fit_type == "by_col"],
                          max(drug1_axis$conc_1) * conc_margin)
          ref_conc_2 <- pmin(df_iso$conc_2[df_iso$conc_1 == 0 & df_iso$fit_type == "by_row"],
                          max(drug2_axis$conc_2) * conc_margin)

          # cutoff point by diagonal (co-dilution)
          # co-dil is given as concentration of drug 1
          df_fit <- all_fits[["by_codil"]]
          if (nrow(df_fit) > 1) {
            df_fit <- df_fit[nrow(df_fit):1, ]
            df_fit <- df_fit[df_fit$fit_type %in% "DRC3pHillFitModelFixS0", ]
            conc_mix <- ifelse(df_fit$x_0 < isobol_value, 0, ifelse(df_fit$x_inf > isobol_value,
              Inf,
              df_fit$ec50 * ((((df_fit$x_0 - df_fit$x_inf) / (isobol_value - df_fit$x_inf)) - 1) ^
                  (1 / df_fit$h)))
                )
            df_iso_codil <- data.frame(conc_1 = conc_mix / (1 + 1 / df_fit$conc_ratio),
                              conc_2 = conc_mix / (1 + df_fit$conc_ratio), fit_type = "by_codil")
            # avoid extrapolation
            capped_idx <- df_iso_codil$conc_1 > (max(drug1_axis$conc_1) * conc_margin) |
                  df_iso_codil$conc_2 > (max(drug2_axis$conc_2) * conc_margin)
            df_iso_codil$conc_1[capped_idx] <- pmin(max(drug1_axis$conc_1) * conc_margin,
                        df_fit$conc_ratio * (max(drug2_axis$conc_2) * conc_margin))[capped_idx]
            df_iso_codil$conc_2[capped_idx] <- pmin(max(drug2_axis$conc_2) * conc_margin,
                        (max(drug1_axis$conc_1) * conc_margin) / df_fit$conc_ratio)[capped_idx]

            df_iso <- rbind(df_iso, df_iso_codil)
          }
          # remove low concentration values
          df_iso <- df_iso[!is.na(df_iso$conc_1) & !is.na(df_iso$conc_2), ]
          df_iso <- df_iso[(df_iso$conc_1 > drug1_axis$conc_1[2] / 2 | df_iso$fit_type == "by_row") &
                             (df_iso$conc_2 > drug2_axis$conc_2[2] / 2 | df_iso$fit_type == "by_col"), ]

          if (nrow(df_iso) < 5) {
            next
          }
  
          # cap the values for plotting and smoothing
          df_iso$pos_x <- pmax(log10(df_iso$conc_2), min(drug2_axis$pos_x))
          df_iso$pos_y <- pmax(log10(df_iso$conc_1), min(drug1_axis$pos_y))

          df_iso$pos_x <- pmin(df_iso$pos_x, max(drug2_axis$pos_x) + log10(conc_margin))
          df_iso$pos_y <- pmin(df_iso$pos_y, max(drug1_axis$pos_y) + log10(conc_margin))

          # rotate 45 degree to calculate smooth curve:
          df_iso$x1 <- (df_iso$pos_x - min(drug2_axis$pos_x) -
                                    (df_iso$pos_y - min(drug1_axis$pos_y))) / sqrt(2) # conc_ratio
          df_iso$x2 <- (df_iso$pos_x - min(drug2_axis$pos_x) +
                                    (df_iso$pos_y - min(drug1_axis$pos_y))) / sqrt(2) # new response value
          x2_extra_offset <- 1 / 4 # offset helps with smoothing of the edges
          df_iso$x2_off <- df_iso$x2 + abs(df_iso$x1) * x2_extra_offset
          isobol_x1 <- seq(min(df_iso$x1) - log2_pos_offset, max(df_iso$x1) + log2_pos_offset, .1)

          # perform the smoothing
          df_iso_curve <- data.frame(x1 = isobol_x1, x2_off = zoo::rollmean(
            rowMeans(do.call(cbind, lapply(names(which(table(df_iso$fit_type) > 1)), function(x)
              stats::approx(x = df_iso$x1[df_iso$fit_type == x], 
                     y = df_iso$x2_off[df_iso$fit_type == x], xout = isobol_x1)$y)),
              na.rm = TRUE), 5, fill = NA))
          df_iso_curve <- df_iso_curve[!is.na(df_iso_curve$x2_off), ]
          df_iso_curve$x2 <- df_iso_curve$x2_off - abs(df_iso_curve$x1) * x2_extra_offset

          # rotate back the position
          df_iso_curve$pos_x <- (df_iso_curve$x1 + df_iso_curve$x2) / sqrt(2) + min(drug2_axis$pos_x)
          df_iso_curve$pos_y <- (-df_iso_curve$x1 + df_iso_curve$x2) / sqrt(2) + min(drug1_axis$pos_y)

          # calculate the reference (additive model in the rotated space)
          c2 <- ref_conc_2 / (1 + (ref_conc_2 / ref_conc_1) * (10 ^ (
                      -(sqrt(2) * df_iso_curve$x1 + min(drug2_axis$pos_x) - min(drug1_axis$pos_y)))))
          df_iso_curve$x2_ref <- (log10(c2) - min(drug2_axis$pos_x) +
                    (log10(ref_conc_1 * (1 - c2 / ref_conc_2)) - min(drug1_axis$pos_y))) / sqrt(2)

          # cap the concentrations for the reference
          over_edge <- pmax(0, (-df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
                        min(drug1_axis$pos_y) -  (max(drug1_axis$pos_y) + conc_margin)) +
                      pmax(0, (df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
                        min(drug2_axis$pos_x) -  (max(drug2_axis$pos_x) + conc_margin))
          df_iso_curve$x2_ref_cap <- df_iso_curve$x2_ref - over_edge * sqrt(2)

          # rotate back the reference
          df_iso_curve$pos_x_ref <- (df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
            min(drug2_axis$pos_x)
          df_iso_curve$pos_y_ref <- (-df_iso_curve$x1 + df_iso_curve$x2_ref) / sqrt(2) +
            min(drug1_axis$pos_y)
          df_iso_curve <- rbind(NA, df_iso_curve, NA)
          df_iso_curve[1, c("pos_x", "pos_x_ref")] <- min(drug2_axis$pos_x)
          df_iso_curve[1, c("pos_y", "pos_y_ref")] <- log10(ref_conc_1)
          df_iso_curve[1, "x1"] <-
                      (df_iso_curve$pos_x[1] - min(drug2_axis$pos_x) -
                          (df_iso_curve$pos_y[1] - min(drug1_axis$pos_y))) / sqrt(2)
          df_iso_curve[1, c("x2", "x2_ref", "x2_ref_cap")] <-
                      (df_iso_curve$pos_x[1] - min(drug2_axis$pos_x) +
                          (df_iso_curve$pos_y[1] - min(drug1_axis$pos_y))) / sqrt(2)

          df_iso_curve[nrow(df_iso_curve), c("pos_x", "pos_x_ref")] <- log10(ref_conc_2)
          df_iso_curve[nrow(df_iso_curve), c("pos_y", "pos_y_ref")] <- min(drug1_axis$pos_y)
          df_iso_curve[nrow(df_iso_curve), c("x2", "x2_ref", "x2_ref_cap")] <-
                      (df_iso_curve$pos_x[nrow(df_iso_curve)] - min(drug2_axis$pos_x) +
                          (df_iso_curve$pos_y[nrow(df_iso_curve)] -
                          min(drug1_axis$pos_y))) / sqrt(2)
          df_iso_curve[nrow(df_iso_curve), "x1"] <-
                      (df_iso_curve$pos_x[nrow(df_iso_curve)] - min(drug2_axis$pos_x) -
                          (df_iso_curve$pos_y[nrow(df_iso_curve)] -
                          min(drug1_axis$pos_y))) / sqrt(2)

          # calculate CI across range to concentration ratios (in the rotated space)
          df_iso_curve$log10_ratio_conc <- df_iso_curve$x1
          df_iso_curve$log2_CI <- zoo::rollmean(
              log2(10) * (df_iso_curve$x2 - df_iso_curve$x2_ref) / sqrt(2), 4,
              fill = c(0, 0, 0))

          # cap position for plotting the isobolograms
          df_iso$pos_y <- pmin(df_iso$pos_y, max(drug1_axis$pos_y) + log2_pos_offset)
          df_iso$pos_x <- pmin(df_iso$pos_x, max(drug2_axis$pos_x) + log2_pos_offset)
          df_iso_curve$pos_y <- pmin(df_iso_curve$pos_y, max(drug1_axis$pos_y) + log2_pos_offset)
          df_iso_curve$pos_x <- pmin(df_iso_curve$pos_x, max(drug2_axis$pos_x) + log2_pos_offset)
          df_iso_curve$pos_y_ref <- pmin(df_iso_curve$pos_y_ref,
              max(drug1_axis$pos_y) + log2_pos_offset)
          df_iso_curve$pos_x_ref <- pmin(df_iso_curve$pos_x_ref,
              max(drug2_axis$pos_x) + log2_pos_offset)

          range <- 2 # in log10 domain --> 100-fold range for calculating averaged CI 
          ratio_idx <- which(
              (df_iso_curve$log10_ratio_conc > (min(df_iso_curve$log10_ratio_conc) + range / 2)) &
              (df_iso_curve$log10_ratio_conc < (max(df_iso_curve$log10_ratio_conc) - range / 2)))

          if (length(ratio_idx) == 0) {
            ratio_idx <- round(nrow(df_iso_curve) / 2)
          }

          df_100x_AUC <- data.frame(log10_ratio_conc = df_iso_curve$log10_ratio_conc[ratio_idx],
            AUC_CI = vapply(ratio_idx, function(x) mean(df_iso_curve$log2_CI[
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
            max_log2CI = max(df_iso_curve$log2_CI),
            min_log2CI = min(df_iso_curve$log2_CI),
            CI_100x = min(df_100x_AUC$AUC_CI)
            )
        }
        all_iso <- all_iso[!vapply(all_iso, is.null, FUN.VALUE = logical(1))]

        df_CI_100x <- unlist(lapply(all_iso, "[[", "CI_100x"))
        df_CI_100x <- data.frame(level = as.numeric(names(df_CI_100x)), log2_CI = df_CI_100x)

        # individual results
        all_combo_variables[[idc, iCL]] <- list(all_mx = all_mx, all_iso = all_iso, df_CI_100x = df_CI_100x,
          drug1_axis = drug1_axis, drug2_axis = drug2_axis, norm_method = norm_method, ref_x50 = ref_x50,
          condition = c(combo$condition, CellLineName = SummarizedExperiment::colData(se)$CellLineName[iCL], 
                        CLID = SummarizedExperiment::colData(se)$clid[iCL]))

        # store the aggregated metrics for each drug pair/cell line
        agg_results$CI_100x_50[idc, iCL] <- ifelse(length(df_CI_100x) == 0, 1,
            ifelse("0.5" %in% df_CI_100x$level, 2 ^ df_CI_100x$log2_CI[df_CI_100x$level == 0.5], 1))
        agg_results$CI_100x_80[idc, iCL] <- ifelse(length(df_CI_100x) == 0, 1,
            ifelse("0.2" %in% df_CI_100x$level, 2 ^ df_CI_100x$log2_CI[df_CI_100x$level == 0.2], 1))
        agg_results$hsa_q10[idc, iCL] <- mean(
          all_mx$hsa_excess[all_mx$hsa_excess <= quantile(all_mx$hsa_excess, .1, na.rm = TRUE)])
        agg_results$bliss_q10[idc, iCL] <- mean(
          all_mx$bliss_excess[all_mx$bliss_excess <= quantile(all_mx$bliss_excess, .1, na.rm = TRUE)])

      }
    }
    agg_results_norm[[norm_method]] <- agg_results
    all_combo_variables_norm[[norm_method]] <- all_combo_variables
  }

  return(list(all_combo_variables_norm <- all_combo_variables_norm, agg_results_norm = agg_results_norm))

}
