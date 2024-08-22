test_that("identify_data_type and split_raw_data works as expected", {
  df_layout <-
    data.table::as.data.table(merge.data.frame(cell_lines[7:8, ], drugs[c(4:6), ], by = NULL))
  df_layout <- gDRtestData::add_data_replicates(df_layout)
  df_layout <- gDRtestData::add_concentration(df_layout, concentrations = 10 ^ (seq(-3, .5, .5)))
  
  df_2 <-
    data.table::as.data.table(merge.data.frame(
      cell_lines[cell_lines$clid %in% df_layout$clid, ], 
      drugs[c(21, 26), ], 
      by = NULL
    ))
  df_2 <- gDRtestData::add_data_replicates(df_2)
  df_2 <- gDRtestData::add_concentration(df_2, concentrations = 10 ^ (seq(-3, .5, .5)))
  colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")] <- 
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")], "_2")
  
  df_layout_2 <- df_layout[df_2, on = intersect(names(df_layout), names(df_2)),
                           allow.cartesian = TRUE]
  
  df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)
  df <- identify_data_type(df_merged_data)
  expect_equal(ncol(df_merged_data), ncol(df))
  expect_true("type" %in% names(df))
  
  df_list <- split_raw_data(df)
  expect_true(inherits(df_list, "list"))
  expect_true(all(names(df_list) %in% c(gDRutils::get_supported_experiments("combo"),
                                        gDRutils::get_supported_experiments("sa"))))
  
  
  df2 <- data.table::data.table(Gnumber = c(rep("DrugA", 9), "DrugB"),
                                DrugName = c(rep("DrugA", 9), "DrugB"),
                                drug_moa = "unknown",
                                drug_moa_2 = "unknown",
                                Gnumber_2 = c(rep("DrugB", 9), "DrugA"),
                                DrugName_2 = c(rep("DrugB", 9), "DrugA"),
                                Concentration = runif(10) + 0.01,
                                Concentration_2 = runif(10) + 0.01,
                                clid = "CL1")
  df2 <- identify_data_type(df2)
  split_df2 <- split_raw_data(df2)
  expect_true(all(split_df2[["Gnumber"]] == "DrugA"))
})


test_that("collapse drugs works as expected", {
  idfs <- gDRutils::get_env_identifiers()
  cols <- idfs[c("drug", "drug2", "drug3",
                 "drug_name", "drug_name2", "drug_name3",
                 "drug_moa", "drug_moa2", "drug_moa3",
                 "concentration", "concentration2", "concentration3")]
  dt <- data.table::data.table("vehicle",
                               "drug2",
                               "drug3",
                               "vehicle",
                               "drug2",
                               "drug3",
                               "vehicle",
                               "moa1",
                               "moa1",
                               0,
                               0.35,
                               0.2)
  data.table::setnames(dt, unlist(cols))
  dt_collapsed <- collapse_drugs(dt)
  expect_equal(dt_collapsed[[cols$drug]], "drug2")
  expect_equal(dt_collapsed[[cols$drug2]], "drug3")
  expect_equal(dt_collapsed[[cols$drug3]], "vehicle")
})


test_that("process_perturbations works as expected", {
  dt <- data.table::data.table(
    DrugName = c("vehicle", "drugA", "drugA"),
    Concentration = c(0, 10, 0),
    DrugName_2 = c("vehicle", "drugB", "drugB"),
    Concentration_2 = c(0, 20, 0)
  )
  
  drugs_cotrt_ids <- c("DrugName", "DrugName_2")
  conc_cotrt_ids <- c("Concentration", "Concentration_2")
  
  result <- process_perturbations(dt, drugs_cotrt_ids, conc_cotrt_ids)
  
  expected <- data.table::data.table(
    drugA = c(0, 10, 0),
    drugB = c(0, 20, 0)
  )
  
  expect_equal(result, expected)
  
  
  dt2 <- data.table::data.table(
    drug1 = c("vehicle", "drugA", "drugA"),
    conc1 = c(0, 10, 0),
    drug2 = c("vehicle", "drugB", "drugB"),
    conc2 = c(0, 20, 0),
    drug3 = c("vehicle", "drugC", "drugC"),
    conc3 = c(0, 30, 0)
  )
  
  drugs_cotrt_ids <- c("drug1", "drug2", "drug3")
  conc_cotrt_ids <- c("conc1", "conc2", "conc3")
  
  result <- process_perturbations(dt2, drugs_cotrt_ids, conc_cotrt_ids)
  
  expected <- data.table::data.table(
    drugA = c(0, 10, 0),
    drugB = c(0, 20, 0),
    drugC = c(0, 30, 0)
  )
  expect_equal(result, expected)
  
  
  dt3 <- data.table::data.table(
    drug1 = c("vehicle", "drugA", "drugB"),
    conc1 = c(0, 10, 2),
    drug2 = c("vehicle", "drugB", "drugB"),
    conc2 = c(0, 20, 0),
    drug3 = c("vehicle", "drugC", "drugC"),
    conc3 = c(0, 30, 0)
  )
  
  drugs_cotrt_ids <- c("drug1", "drug2", "drug3")
  conc_cotrt_ids <- c("conc1", "conc2", "conc3")
  
  result <- process_perturbations(dt3, drugs_cotrt_ids, conc_cotrt_ids)
  
  expected <- data.table::data.table(
    drug1 = c("vehicle", "drugA", "drugB"),
    conc1 = c(0, 10, 2),
    drugB = c(0, 20, 0),
    drugC = c(0, 30, 0)
  )
  expect_equal(result, expected)
  
  
  dt4 <- data.table::data.table(
    Gnumber = c("vehicle", "drugA", "drugB"),
    Concentration = c(0, 10, 2),
    Gnumber_2 = c("vehicle", "drugB", "drugB"),
    Concentration_2 = c(0, 20, 0),
    drug_moa_2 = c("vehicle", "moa_A", "moa_A"),
    DrugName_2 = c("vehicle", "drugB", "drugB")
  )
  
  drugs_cotrt_ids <- "Gnumber_2"
  conc_cotrt_ids <- "Concentration_2"
  
  result <- process_perturbations(dt4, drugs_cotrt_ids, conc_cotrt_ids)
  
  expected <- data.table::data.table(
    Gnumber = c("vehicle", "drugA", "drugB"),
    Concentration = c(0, 10, 2),
    drugB = c(0, 20, 0)
  )
  expect_equal(result, expected)
})

