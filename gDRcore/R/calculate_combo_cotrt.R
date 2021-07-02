#' Calculate combo co-treatment
#'
#' @param se a BumpyMatrix SE with drug response data
#'
#' @return
#' @export
calculate_combo_cotrt <- function(se) {
    
    checkmate::assert_class(se, "SummarizedExperiment")
    
    drug_combos <- S4Vectors::metadata(se)$drug_combinations
    df_all_combo <- do.call(rbind, 
                           lapply(drug_combos, function(x) as.data.frame(x$condition)))
    rownames(df_all_combo) <- paste0(df_all_combo$Gnumber, 
                                    "_", 
                                    df_all_combo$Gnumber_2, 
                                    "_", 
                                    df_all_combo$Concentration_2)
    
    fixed_metrics <- SummarizedExperiment::SummarizedExperiment( 
        assays = list(delta_RV_AOC_range = matrix(NA, nrow(df_all_combo), ncol(se)),
                      RV_AOC_range = matrix(NA, nrow(df_all_combo), ncol(se)),
                      RV_AOC_range_ref = matrix(NA, nrow(df_all_combo), ncol(se)),
                      
                      delta_RV_mean = matrix(NA, nrow(df_all_combo), ncol(se)),
                      RV_mean = matrix(NA, nrow(df_all_combo), ncol(se)),
                      RV_mean_ref = matrix(NA, nrow(df_all_combo), ncol(se)),
                      
                      delta_GR_AOC_range = matrix(NA, nrow(df_all_combo), ncol(se)),
                      GR_AOC_range = matrix(NA, nrow(df_all_combo), ncol(se)),
                      GR_AOC_range_ref = matrix(NA, nrow(df_all_combo), ncol(se)),
                      
                      log10_ratio_IC50 = matrix(NA, nrow(df_all_combo), ncol(se)),
                      IC50 = matrix(NA, nrow(df_all_combo), ncol(se)),
                      IC50_ref = matrix(NA, nrow(df_all_combo), ncol(se)),
                      
                      log10_ratio_GR50 = matrix(NA, nrow(df_all_combo), ncol(se)),
                      GR50 = matrix(NA, nrow(df_all_combo), ncol(se)),
                      GR50_ref = matrix(NA, nrow(df_all_combo), ncol(se)),
                      
                      delta_E_max = matrix(NA, nrow(df_all_combo), ncol(se)),
                      E_max = matrix(NA, nrow(df_all_combo), ncol(se)),
                      E_max_ref = matrix(NA, nrow(df_all_combo), ncol(se)),
                      
                      delta_GR_max = matrix(NA, nrow(df_all_combo), ncol(se)),
                      GR_max = matrix(NA, nrow(df_all_combo), ncol(se)),
                      GR_max_ref = matrix(NA, nrow(df_all_combo), ncol(se))
        ),
        rowData = df_all_combo,
        colData = SummarizedExperiment::colData(se))
    
    
    # run through all drug combinations
    for (idc in seq_along(S4Vectors::metadata(se)$drug_combinations)) {
        combo <- S4Vectors::metadata(se)$drug_combinations[[idc]]
        print(sprintf("%s x %s (%i/%i)", combo$condition$DrugName, combo$condition$DrugName_2,
                      idc, length(S4Vectors::metadata(se)$drug_combinations)))
        combo_idxes <- df_all_combo$Gnumber == combo$condition$Gnumber & 
            df_all_combo$Gnumber_2 == combo$condition$Gnumber_2        
        
        # run through all cell lines
        for (iCL in seq_len(ncol(se))) {
            # get all data as flat data.frame
            
            flat_data <- gDRutils::convert_se_assay_to_dt(se[combo$rows, iCL], "Metrics")
            if (nrow(flat_data) < 2) next
            
            ref_data <- flat_data[flat_data$Gnumber %in% combo$condition$Gnumber & 
                                      flat_data$Concentration_2 == 0, ]
            
            for (c2 in combo$condition$Concentration_2) {
                combo_idx <- which(combo_idxes & (df_all_combo$Concentration_2 == c2))
                
                cotrt_data <- flat_data[flat_data$Gnumber %in% combo$condition$Gnumber & 
                                            flat_data$Gnumber_2 %in% combo$condition$Gnumber_2 & 
                                            flat_data$Concentration_2 == c2, ]
                
                if (nrow(ref_data) != 1 || nrow(cotrt_data) != 1) next
                
                for (i in c("RV_AOC_range", "RV_mean", "GR_AOC_range", "IC50", "GR50", "E_max", "GR_max")) {
                    SummarizedExperiment::assay(fixed_metrics, i)[combo_idx, iCL] <- as.numeric(cotrt_data[1, ..i])
                    SummarizedExperiment::assay(fixed_metrics, paste0(i, "_ref"))[combo_idx, iCL] <- as.numeric(ref_data[1, ..i])
                }
                for (i in c("RV_AOC_range", "RV_mean", "GR_AOC_range", "E_max", "GR_max")) { # delta metrics
                    SummarizedExperiment::assay(fixed_metrics, paste0("delta_", i))[combo_idx, iCL] <- 
                        as.numeric(cotrt_data[, ..i] - ref_data[, ..i])
                }
                
                SummarizedExperiment::assay(fixed_metrics, "log10_ratio_IC50")[combo_idx, iCL] <- log10(
                    min(2 * 10 ^ ref_data$maxlog10Concentration, ref_data$IC50) / 
                        min(2 * 10 ^ cotrt_data$maxlog10Concentration, cotrt_data$IC50)
                )
                SummarizedExperiment::assay(fixed_metrics, "log10_ratio_GR50")[combo_idx, iCL] <- log10(
                    min(2 * 10 ^ ref_data$maxlog10Concentration, ref_data$GR50) / 
                        min(2 * 10 ^ cotrt_data$maxlog10Concentration, cotrt_data$GR50)
                )
            }
        }
    }
}
