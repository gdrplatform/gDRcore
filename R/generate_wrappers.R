#' generateNoNoiseRawData
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateNoNoiseRawData <- function(cell_lines, drugs, save = TRUE) {
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # generate the data for the 1st test set: no noise
  #   only for testing purpuses not displayed as example
  df_merged <- gDRtestData::prepareMergedData(cell_lines[2:11, ], drugs[2:11, ], 0)
  
  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_small_no_noise.qs"
    )
  }
  invisible(mae)
}

#' generateNoiseRawData
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateNoiseRawData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 1st test set with noise
  
  df_merged <- gDRtestData::prepareMergedData(cell_lines[2:11, ], drugs[2:11, ])

  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_small.qs"
    )
  }
  invisible(mae)
}

#' generateLigandData
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateLigandData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 1st test set with ligand as reference
  df_merged <- gDRtestData::prepareMergedData(cell_lines[2:6, ], drugs[2:5, ], 0)
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
  
  mae <- runDrugResponseProcessingPipeline(
    df_merged, 
    override_untrt_controls = c(Ligand = 0.1),
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_wLigand.qs"
    )
  }
  invisible(mae)
}

#' generateMediumData
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateMediumData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 2nd (medium size) test set with single agent
  df_merged <- gDRtestData::prepareMergedData(cell_lines[seq_len(15), ], drugs[seq_len(40), ])
    
  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_medium.qs"
    )
  }
  invisible(mae)
}


#' generateComboNoNoiseData
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateComboNoNoiseData <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with combo (two single dose)
  #   co-treatment drug is only as DrugName_2
  df_merged <- gDRtestData::prepareComboMergedData(cell_lines[2:4, ], drugs, noise = 0)
    
  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_2dose_nonoise.qs"
    )
  }
  invisible(mae)
}

#' generateComboNoNoiseData2
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateComboNoNoiseData2 <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with combo (two single dose)
  #   co-treatment drug is also as single agent as DrugName
  df_merged <- gDRtestData::prepareComboMergedData(cell_lines[2:4, ], drugs, drugsIdx1 = c(2:4, 26), noise = 0)
  df_merged <- df_merged[!(df_merged$Gnumber %in% c("vehicle", drugs$Gnumber[26]) &
    df_merged$Gnumber_2 == drugs$Gnumber[26]), ]
    
  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_2dose_nonoise2.qs"
    )
  }
  invisible(mae)
}

#' generateComboNoNoiseData3
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateComboNoNoiseData3 <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the 3rd test set with combo (two single dose)
  #   co-treatment drug does NOT have single agent response
  df_merged <- gDRtestData::prepareComboMergedData(
    cell_lines = cell_lines[2:4, ], 
    drugs = drugs, 
    drugsIdx1 = 2:4, 
    noise = 0, 
    modifyDf2 = TRUE
  )

  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_2dose_nonoise3.qs"
    )
  }
  invisible(mae)
}

#' generateComboMatrixSmall
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateComboMatrixSmall <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with combo matrix (small, no noise)
  concentration <- 10^ (seq(-3, .5, .5))
  df_layout <- gDRtestData::prepareData(cell_lines[7:8, ], drugs[c(4:6), ], concentration)
  df_2 <- gDRtestData::prepareData(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26), ], concentration)
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- df_layout[df_2, on = intersect(names(df_layout), names(df_2)),
                           allow.cartesian = TRUE]
  
  df_merged <- gDRtestData::generate_response_data(df_layout_2, 0)

  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_matrix_small.qs"
    )
  }
  invisible(mae)
}
  
#' generateComboMatrix
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateComboMatrix <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with combo matrix (mid-size)
  df_layout <- gDRtestData::prepareData(cell_lines[seq(1, 30, 4), ], drugs[c(1, 2, 11), ])
  df_2 <- gDRtestData::prepareData(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26, 31), ])
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- df_layout[df_2, on = intersect(names(df_layout), names(df_2)),
                          allow.cartesian = TRUE]
  
  df_merged <- gDRtestData::generate_response_data(df_layout_2)

  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_matrix.qs"
    )
  }
  invisible(mae)
}

#' generateTripleComboMatrix
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateTripleComboMatrix <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with triple combo  (no noise)
  concentration <- 10^ (seq(-3, .5, .5))
  df_layout <- gDRtestData::prepareData(cell_lines[7:8, ], drugs[c(4:6), ], concentration)
  
  df_2 <- gDRtestData::prepareData(
    cell_lines[cell_lines$clid %in% df_layout$clid, ], 
    drugs[c(21, 26), ],
    c(0, concentration)
  )
  df_2 <- changeColNames(df_2, drugs, "_2")
  df_layout_2 <- df_layout[df_2, on = intersect(names(df_layout), names(df_2)),
                           allow.cartesian = TRUE]
  
  df_3 <- gDRtestData::prepareData(
    cell_lines[cell_lines$clid %in% df_layout$clid, ], 
    drugs[10, ],
    c(0, .1, 1)
  )
  df_3 <- changeColNames(df_3, drugs, "_3")
df_layout_3 <- df_layout[df_2, on = intersect(names(df_layout), names(df_2)), 
                         allow.cartesian = TRUE][
                           df_3, on = intersect(names(df_layout), names(df_3)), 
                           allow.cartesian = TRUE]
  
  df_merged <- gDRtestData::generate_response_data(df_layout_3, 0)

  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_triple.qs"
    )
  }
  invisible(mae)
}

#' generateCodilutionSmall
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateCodilutionSmall <- function(cell_lines, drugs, save = TRUE) {
  # generate the data with combo co-dilution (small)
  df_layout <- gDRtestData::prepareData(cell_lines[seq_len(2), ], drugs[seq_len(4), ])

  df_2 <- cbind(drugs[1, , drop = FALSE], df_layout[, "Concentration", drop = FALSE])
  df_layout_2 <- gDRtestData::prepareCodilutionData(df_2, df_layout)
  
  df_merged <- gDRtestData::generate_response_data(df_layout_2, 0)
  
  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_codilution_small.qs"
    )
  }
  invisible(mae)
}

#' generateCodilution
#' 
#' @keywords internal
#' @return data.table with raw input data or MAE with processed data
generateCodilution <- function(cell_lines, drugs, save = TRUE) {
  # generate the data for the test set with combo (co-dilution)
  df_layout <- gDRtestData::prepareData(cell_lines[seq(1, 15, 2), ], drugs[seq_len(12), ])
  
  df_2 <- cbind(drugs[c(1, 1), ], df_layout[, "Concentration", drop = FALSE])
  df_layout_2 <- gDRtestData::prepareCodilutionData(df_2, df_layout)

  df_merged <- gDRtestData::generate_response_data(df_layout_2)
  mae <- runDrugResponseProcessingPipeline(
    df_merged,
    nested_confounders = gDRutils::get_env_identifiers("barcode")[1]
  )
  
  if (save) {
    save_qs(
      qsObj = mae,
      qsName = "finalMAE_combo_codilution.qs"
    )
  }
  
  invisible(mae)
}

#' @keywords internal
save_qs <- function(qsObj, qsName) {
  qs::qsave(
    qsObj,
    file.path(system.file("testdata", package = "gDRtestData"), qsName)
  )
}

#' @keywords internal
changeColNames <- function(df, drugs, suffix) {
  cols <- colnames(df) %in% c(colnames(drugs), "Concentration")
  colnames(df)[cols] <- paste0(colnames(df)[cols], suffix)
  
  df
}
