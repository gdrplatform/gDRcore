

calculate_combo_codilution <- function(SE){

  # create new SE for the codilution
  co_dilution_metrics = metadata(SE)$drug_combinations

  codil_SE = SE[1:length(co_dilution_metrics),]
  rownames(codil_SE) = sapply(co_dilution_metrics, function(x)
      paste0(x$condition[['DrugName']],'_',x$condition[['DrugName_2']]))
  SummarizedExperiment::rowData(codil_SE) = 
      DataFrame(t(sapply(co_dilution_metrics, '[[', 'condition')))
  metadata(codil_SE) = list()
  SummarizedExperiment::assays(codil_SE) = 
      SummarizedExperiment::assays(codil_SE)[c('Averaged', 'Avg_Controls', 'Metrics')]

  mSE_m <- SummarizedExperiment::assay(codil_SE, "Metrics")
  a_SE = SummarizedExperiment::assay(codil_SE, "Averaged")
  a_ctrl_SE = SummarizedExperiment::assay(codil_SE, "Avg_Controls")
  flat_data = gDRutils::assay_to_dt(SE, 'Averaged')
  flat_data$rId = as.factor(flat_data$rId)
  flat_data$cId = as.factor(flat_data$cId)

  for (ic in 1:length(co_dilution_metrics)) {
      for (iCL in colnames(SE)) {

          combo = co_dilution_metrics[[ic]]
          df_ = flat_data[flat_data$rId %in% combo$rows & flat_data$cId == iCL,]
          
          df_ = df_[df_$Concentration_2 != 0 ,]
          df_$Concentration_1 = df_$Concentration
          df_$Concentration =  df_$Concentration_1 + df_$Concentration_2
          a_SE[[ic, iCL]] = df_
          a_ctrl_SE[[ic, iCL]]$RefReadout = a_ctrl_SE[[ic, iCL]]$UntrtReadout
          a_ctrl_SE[[ic, iCL]]$RefRelativeViability = 1
          a_ctrl_SE[[ic, iCL]]$RefGRvalue = 1
      }
  }
  SummarizedExperiment::assay(codil_SE, "Averaged") <- a_SE
  SummarizedExperiment::assay(codil_SE, "Avg_Controls") <- a_ctrl_SE

  codil_SE = gDR::metrics_SE(codil_SE)
  single_SE = SE[rowData(SE)$Concentration_2 == 0, ]
  assays(single_SE) = assays(single_SE)[names(assays(codil_SE))]
  rowData(single_SE)$Concentration_2 = NULL
  rowData(single_SE)$conc_log10_ratio = NA
  
  codil_SE = rbind(codil_SE, single_SE)

}