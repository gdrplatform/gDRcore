#' Transform mae into raw data
#'
#' @param mae MultiAssayExperiment object with SummarizedExperiments containing "RawTreated" and "Controls" assays
#'
#'
#' @examples 
#' mae <- gDRutils::get_synthetic_data("finalMAE_small")
#' convert_mae_to_raw_data(mae)
#' 
#' @return data.table with raw data
#' @export
convert_mae_to_raw_data <- function(mae) {
  
  checkmate::assert_class(mae, "MultiAssayExperiment")
  lapply(gDRutils::MAEpply(mae, SummarizedExperiment::assayNames),
         function(x) all(c("RawTreated", "Controls") %in% x))
  
  # Get the data from se
  
  data <- gDRutils::MAEpply(mae, convert_se_to_raw_data)
  
  data_df <- data.table::rbindlist(data, fill = TRUE)
  
  data_df <- replace_NA_in_raw_data(data_df, mae)
  
  data.table::setorder(data_df)
  
  data_df <- data_df[!duplicated(data_df$record_id), ]
  selected_columns <- !names(data_df) %in% c("record_id", "BackgroundValue",
                                          "WellColumn", "WellRow", "Template", "swap_sa")
  data.table::setorder(data_df, record_id)
  data_df[, selected_columns, with = FALSE]
}


#' Transform se into raw_data
#'
#' @param se SummarizedExperiment object with "RawTreated" and "Controls" assays
#'
#' @return data.table with raw data
#' 
#' @examples 
#' mae <- gDRutils::get_synthetic_data("finalMAE_small")
#' se <- mae[[1]]
#' convert_se_to_raw_data(se)
#' 
#' @export
convert_se_to_raw_data <- function(se) {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_true(all(c("RawTreated", "Controls") %in% SummarizedExperiment::assayNames(se)))
  
  # Get the data from se
  trt <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
  ctrl <- gDRutils::convert_se_assay_to_dt(se, "Controls")
  
  # Identifiers
  conc_cols1 <- intersect(unlist(gDRutils::get_SE_identifiers(se, "concentration",
                                                   simplify = FALSE)), names(trt))
  conc_cols2 <- intersect(unlist(gDRutils::get_SE_identifiers(se, "concentration2",
                                                    simplify = FALSE)), names(trt))
  conc_cols <- c(conc_cols1, conc_cols2)
  drug_cols1 <- intersect(unlist(gDRutils::get_SE_identifiers(se, c("drug", "drug_name", "drug_moa"),
                                                   simplify = FALSE)), names(trt))
  drug_cols2 <- intersect(unlist(gDRutils::get_SE_identifiers(se, c("drug2", "drug_name2", "drug_moa2"),
                                                   simplify = FALSE)), names(trt))
  drug_cols <- c(drug_cols1, drug_cols2)
  untreated_tag <- gDRutils::get_SE_identifiers(se, "untreated_tag")[1]
  masked_tag <- gDRutils::get_SE_identifiers(se, "masked_tag")[1]
  ref_div_time <- gDRutils::get_SE_identifiers(se, "cellline_ref_div_time")
  readout_val <- "ReadoutValue"
  background_val <- "BackgroundValue"
  duration <- gDRutils::get_SE_identifiers(se, "duration")
  
  # Add required cols and correct the data
  ctrl[, (drug_cols) := NULL]
  ctrl[, Duration := ifelse(!is.na(isDay0) & isDay0, 0, Duration)]
  ctrl[, c("control_type", "isDay0") := NULL]
  ctrl[, c(eval(conc_cols), "BackgroundValue") := 0]
  ctrl[, eval(drug_cols) := untreated_tag]
  ctrl[, eval(masked_tag) := FALSE]
  
  trt[get(conc_cols1) == 0, (drug_cols1) := untreated_tag]
  
  
  if (length(conc_cols2) > 0 && length(drug_cols2) > 0) {
    trt[get(conc_cols2) == 0, (drug_cols2) := untreated_tag]
  }
  
  # Merge and adjust the data
  merged_df <- data.table::rbindlist(list(trt, ctrl), fill = TRUE)
  merged_df$rId <- merged_df$cId <-  merged_df$ReadoutValue <- NULL
  data.table::setnames(merged_df, "CorrectedReadout", "ReadoutValue")
  
  exp_metadata <- gDRutils::get_SE_experiment_metadata(se)

  cbind(na.omit(merged_df, col = "ReadoutValue"), data.table::as.data.table(exp_metadata))
}


#' @keywords internal
replace_NA_in_raw_data <- function(df, mae) {
  df_ <- data.table::copy(df)
  # Identifiers
  conc_cols <- intersect(unlist(gDRutils::get_SE_identifiers(mae[[1]], c("concentration", "concentration2"),
                                                             simplify = FALSE)), names(df_))
  
  drug_cols <- intersect(unlist(gDRutils::get_SE_identifiers(mae[[1]], c("drug", "drug_name", "drug_moa",
                                                                         "drug2", "drug_name2", "drug_moa2"),
                                                             simplify = FALSE)), names(df_))
  
  untreated_tag <- gDRutils::get_SE_identifiers(mae[[1]], "untreated_tag")
  
  df_[, (conc_cols) := lapply(.SD, function(x) replace(x, is.na(x), 0)), .SDcols = conc_cols]
  df_[, (drug_cols) := lapply(.SD, function(x) replace(x, is.na(x), untreated_tag[[1]])), .SDcols = drug_cols]
  df_
}
