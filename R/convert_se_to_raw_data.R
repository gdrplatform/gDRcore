#' Transform mae into raw_data
#'
#' @param mae MultiAssayExperiment object with SummarizedExperiments containing "RawTreated" and "Controls" assays
#'
#' @return data.frame with raw data
#' @export
convert_mae_to_raw_data <- function(mae) {
  
  checkmate::assert_class(mae, "MultiAssayExperiment")
  lapply(MAEpply(mae, SummarizedExperiment::assayNames), function(x) all(c("RawTreated", "Controls") %in% x))
  
  # Get the data from se
  
  data <- gDRutils::MAEpply(mae, convert_se_to_raw_data)
  data_df <- data.table::rbindlist(data, fill = TRUE)
  
  
  # Identifiers
  conc_cols <- intersect(unlist(get_SE_identifiers(mae[[1]], c("concentration", "concentration2"),
                                                   simplify = FALSE)), names(data_df))
  conc_cols_rev <- intersect(unlist(get_SE_identifiers(mae[[1]], c("concentration2", "concentration"),
                                                   simplify = FALSE)), names(data_df))
  
  drug_cols <- intersect(unlist(get_SE_identifiers(mae[[1]], c("drug", "drug_name", "drug_moa",
                                                         "drug2", "drug_name2", "drug_moa2"),
                                                   simplify = FALSE)), names(data_df))
  
  drug_cols_rev <- intersect(unlist(get_SE_identifiers(mae[[1]], c("drug2", "drug_name2", "drug_moa2",
                                                               "drug", "drug_name", "drug_moa"),
                                                   simplify = FALSE)), names(data_df))
  
  untreated_tag <- get_SE_identifiers(mae[[1]], "untreated_tag")[1]
  masked_tag <- get_SE_identifiers(mae[[1]], "masked_tag")[1]
  
  data_df[, (conc_cols) := lapply(.SD, function(x) replace(x, is.na(x), 0)), .SDcols = conc_cols]
  data_df[, (drug_cols) := lapply(.SD, function(x) replace(x, is.na(x), untreated_tag)), .SDcols = drug_cols]
  
  data_df <- data_df[!duplicated(data_df$record_id), -("record_id")]
  
  
  if (gDRutils::get_env_identifiers("concentration2") %in% colnames(data_df)) {
    single_agent_idx <- 
      data_df[[gDRutils::get_env_identifiers("concentration")]] %in% 0
    data_sa <- data_df[single_agent_idx, ]
    data_df <- data_df[!single_agent_idx, ]
    data.table::setnames(data_sa, c(drug_cols, conc_cols), c(drug_cols_rev, conc_cols_rev))
    data_df <- data.table::rbindlist(list(data_sa, data_df), use.names = TRUE)
  }
  as.data.frame(data_df)
}


#' Transform se into raw_data
#'
#' @param se SummarizedExperiment object with "RawTreated" and "Controls" assays
#'
#' @return data.frame with raw data
#' @export
convert_se_to_raw_data <- function(se) {
  
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_true(all(c("RawTreated", "Controls") %in% SummarizedExperiment::assayNames(se)))
  
  # Get the data from se
  trt <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
  ctrl <- gDRutils::convert_se_assay_to_dt(se, "Controls")
  
  # Identifiers
  conc_cols1 <- intersect(unlist(get_SE_identifiers(se, "concentration",
                                                   simplify = FALSE)), names(trt))
  conc_cols2 <- intersect(unlist(get_SE_identifiers(se, "concentration2",
                                                    simplify = FALSE)), names(trt))
  
  conc_cols <- c(conc_cols1, conc_cols2)
  
  drug_cols1 <- intersect(unlist(get_SE_identifiers(se, c("drug", "drug_name", "drug_moa"),
                                                   simplify = FALSE)), names(trt))
  drug_cols2 <- intersect(unlist(get_SE_identifiers(se, c("drug2", "drug_name2", "drug_moa2"),
                                                   simplify = FALSE)), names(trt))
  drug_cols <- c(drug_cols1, drug_cols2)
  
  
  untreated_tag <- get_SE_identifiers(se, "untreated_tag")[1]
  masked_tag <- get_SE_identifiers(se, "masked_tag")[1]
  
  # Add required cols and correct the data
  ctrl[, (drug_cols) := NULL]
  
  ctrl$Duration <- ifelse(!is.na(ctrl$isDay0) & ctrl$isDay0, 0, ctrl$Duration)
  ctrl[, c("control_type", "isDay0") := NULL]
  ctrl[, c(eval(conc_cols), "BackgroundValue") := 0]
  ctrl[, eval(drug_cols) := untreated_tag]
  ctrl[, eval(masked_tag) := FALSE]
  data.table::setnames(ctrl, "CorrectedReadout", "ReadoutValue")
  
  # Merge and adjust the data
  merged_df <- data.table::rbindlist(list(trt, ctrl), fill = TRUE)
  merged_df$rId <- merged_df$cId <-  merged_df$CorrectedReadout <- NULL
  merged_df <- merged_df[!duplicated(merged_df$record_id), ]

  as.data.frame(merged_df[!is.na(merged_df$ReadoutValue), ])
}
