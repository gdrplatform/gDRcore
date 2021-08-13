#' Calculate combo co-dilution
#'
#' @param se a BumpyMatrix SE with drug response data
#'
#' @return
#' @export
calculate_combo_codilution <- function(se) {
  .Deprecated("fit_SE.combinations")

  checkmate::assert_class(se, "SummarizedExperiment")

  # create new SE for the codilution
  co_dilution_metrics <- S4Vectors::metadata(se)$drug_combinations

  codil_SE <- se[seq_along(co_dilution_metrics), ]
  rownames(codil_SE) <- vapply(co_dilution_metrics, function(x)
      paste0(x$condition[["DrugName"]], "_", x$condition[["DrugName_2"]]),
      character(1))
  SummarizedExperiment::rowData(codil_SE) <- 
      S4Vectors::DataFrame(t(sapply(co_dilution_metrics, "[[", "condition")))
  S4Vectors::metadata(codil_SE) <- list()
  SummarizedExperiment::assays(codil_SE) <-
      SummarizedExperiment::assays(codil_SE)[c("Averaged", "Metrics")]

  mSE_m <- SummarizedExperiment::assay(codil_SE, "Metrics") # not being used currently
  a_SE <- SummarizedExperiment::assay(codil_SE, "Averaged")
  flat_data <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  flat_data$rId <- as.factor(flat_data$rId)
  flat_data$cId <- as.factor(flat_data$cId)

  for (ic in seq_along(co_dilution_metrics)) {
      for (iCL in colnames(se)) {
          combo <- co_dilution_metrics[[ic]]
          df_ <- flat_data[flat_data$rId %in% combo$rows & flat_data$cId == iCL, ]
          df_ <- df_[df_$Concentration_2 != 0, ]
          df_$Concentration_1 <- df_$Concentration
          df_$Concentration <- df_$Concentration_1 + df_$Concentration_2
          a_SE[[ic, iCL]] <- df_
      }
  }
  SummarizedExperiment::assay(codil_SE, "Averaged") <- a_SE

  codil_SE <- gDR::metrics_SE(codil_SE)
  single_SE <- se[SummarizedExperiment::rowData(se)$Concentration_2 == 0, ]
  SummarizedExperiment::assays(single_SE) <- 
    SummarizedExperiment::assays(single_SE)[names(SummarizedExperiment::assays(codil_SE))]
  SummarizedExperiment::rowData(single_SE)$Concentration_2 <- NULL
  SummarizedExperiment::rowData(single_SE)$conc_log10_ratio <- NA
  
  codil_SE <- rbind(codil_SE, single_SE)

}
