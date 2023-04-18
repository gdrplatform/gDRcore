#' generateNoNoiseRawData
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateNoNoiseRawData <- function(cell_lines, drugs, save = TRUE) {
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # generate the data for the 1st test set: no noise
  #   only for testing purpuses not displayed as example
  df_merged <- prepareMergedData(cell_lines[2:11, ], drugs[2:11, ], 0)
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_small_no_noise.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateNoiseRawData
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateNoiseRawData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 1st test set with noise
  df_merged <- prepareMergedData(cell_lines[2:11, ], drugs[2:11, ])

  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_small.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateLigandData
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateLigandData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 1st test set with ligand as reference
  df_merged <- prepareMergedData(cell_lines[2:6, ], drugs[2:5, ], 0)
  df_merged$Ligand <- 0.1
  df_merged2 <- df_merged[df_merged$Gnumber %in% c("vehicle", "G00002", "G00003"), ]
  df_merged2$Ligand <- 0
  
  idx1 <- df_merged2$clid %in% paste0("CL000", 11:12)
  idx2 <- df_merged2$clid %in% paste0("CL000", 13:14)
  df_merged2$ReadoutValue <- 105 - pmax(0, pmin(104, (105 - df_merged2$ReadoutValue) ^ 1.1))
  df_merged2$ReadoutValue[idx1] <- 0.8 * df_merged2$ReadoutValue[idx1]
  df_merged2$ReadoutValue[idx2] <- 0.5 * df_merged2$ReadoutValue[idx2]
  df_merged2$ReadoutValue <- round(df_merged2$ReadoutValue, 1)
  
  df_merged2$Barcode <- paste0(df_merged2$Barcode, "1")
  df_merged <- rbind(df_merged, df_merged2)
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged, 
      override_untrt_controls = c(Ligand = 0.1),
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_wLigand.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateMediumData
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateMediumData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 2nd (medium size) test set with single agent
  df_merged <- prepareMergedData(cell_lines[seq_len(15), ], drugs[seq_len(40), ])
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_medium.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}


#' generateManyLinesData
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateManyLinesData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 2nd (medium size) test set with single agent
  df_merged <- prepareMergedData(cell_lines, drugs[seq_len(40), ])
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_many_lines.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateManyDrugsData
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateManyDrugsData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with single agent (many drugs)
  df_merged <- prepareMergedData(cell_lines[seq_len(10), ], drugs[seq_len(40), ])
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_many_drugs.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateComboNoNoiseData
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateComboNoNoiseData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with combo (two single dose)
  #   co-treatment drug is only as DrugName_2
  df_merged <- prepareComboMergedData(cell_lines[2:4, ], drugs, noise = 0)
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_2dose_nonoise.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateComboNoNoiseData2
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateComboNoNoiseData2 <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with combo (two single dose)
  #   co-treatment drug is also as single agent as DrugName
  df_merged <- prepareComboMergedData(cell_lines[2:4, ], drugs, drugsIdx1 = c(2:4, 26), noise = 0)
  df_merged <- df_merged[!(df_merged$Gnumber %in% c("vehicle", drugs$Gnumber[26]) &
    df_merged$Gnumber_2 == drugs$Gnumber[26]), ]
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_2dose_nonoise2.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateComboNoNoiseData3
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateComboNoNoiseData3 <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 3rd test set with combo (two single dose)
  #   co-treatment drug does NOT have single agent response
  df_merged <- prepareComboMergedData(
    cell_lines = cell_lines[2:4, ], 
    drugs = drugs, 
    drugsIdx1 = 2:4, 
    noise = 0, 
    modifyDf2 = TRUE
  )
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_2dose_nonoise3.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateComboManyDrugs
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateComboManyDrugs <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with combo (unique dose; many drug)
  df_merged <- prepareComboMergedData(
    cell_lines = cell_lines[2:4, ], 
    drugs = drugs,
    drugsIdx1 = -1,
    drugsIdx2 = c(1, 1),
    concentration = c(0, 2)
  )
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_1dose_many_drugs.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateComboMatrixSmall
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateComboMatrixSmall <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with combo matrix (small, no noise)
  concentration <- 10^ (seq(-3, .5, .5))
  df_layout <- prepareData(cell_lines[7:8, ], drugs[c(4:6), ], concentration)
  df_2 <- prepareData(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26), ], concentration)
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- merge(df_layout, df_2)
  
  df_merged <- gDRtestData::generate_response_data(df_layout_2, 0)
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_matrix_small.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}
  
#' generateComboMatrix
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateComboMatrix <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with combo matrix (mid-size)
  df_layout <- prepareData(cell_lines[seq(1, 30, 4), ], drugs[c(1, 2, 11), ])
  df_2 <- prepareData(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26, 31), ])
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- merge(df_layout, df_2)
  
  browser()
  
  df_merged <- gDRtestData::generate_response_data(df_layout_2)
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_matrix.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateTripleComboMatrix
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateTripleComboMatrix <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with triple combo  (no noise)
  concentration <- 10^ (seq(-3, .5, .5))
  df_layout <- prepareData(cell_lines[7:8, ], drugs[c(4:6), ], concentration)
  
  df_2 <- prepareData(
    cell_lines[cell_lines$clid %in% df_layout$clid, ], 
    drugs[c(21, 26), ],
    c(0, concentration)
  )
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- merge(df_layout, df_2)
  
  df_3 <- prepareData(
    cell_lines[cell_lines$clid %in% df_layout$clid, ], 
    drugs[10, ],
    c(0, .1, 1)
  )
  df_3 <- changeColNames(df_3, drugs, "_3")
  df_layout_3 <- merge(merge(df_layout, df_2), df_3)
  
  df_merged <- gDRtestData::generate_response_data(df_layout_3, 0)
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_triple.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateCodilutionSmall
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateCodilutionSmall <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with combo co-dilution (small)
  df_layout <- prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])

  df_2 <- cbind(drugs[1, , drop = FALSE], df_layout[, "Concentration", drop = FALSE])
  df_layout_2 <- prepareCodilutionData(df_2, df_layout)
  
  df_merged <- gDRtestData::generate_response_data(df_layout_2, 0)
  
  if (requireNamespace("gDRcore", quietly = TRUE)) {
    
    mae <- gDRcore::runDrugResponseProcessingPipeline(
      df_merged,
      nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
    )
    
    if (save) {
      save_rds(
        rdsObj = mae,
        rdsName = "finalMAE_combo_codilution_small.RDS"
      )
    }
    invisible(mae)
  } else {
    df_merged
  }
}

#' generateCodilution
#' 
#' @keywords internal
#' @return data.frame with raw input data or MAE with processed data
generateCodilution <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with combo (co-dilution)
  df_layout <- prepareData(cell_lines[seq(1, 15, 2), ], drugs[seq_len(12), ])
  
  df_2 <- cbind(drugs[c(1, 1), ], df_layout[, "Concentration", drop = FALSE])
  df_layout_2 <- prepareCodilutionData(df_2, df_layout)

  df_merged <- gDRtestData::generate_response_data(df_layout_2)
  mae <- gDRcore::runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_rds(
      rdsObj = mae,
      rdsName = "finalMAE_combo_codilution.RDS"
    )
  }
  
  invisible(mae)
}
