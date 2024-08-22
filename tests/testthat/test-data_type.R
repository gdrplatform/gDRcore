test_that("identify_data_type and split_raw_data works as expected", {
  dt_layout <-
    data.table::as.data.table(merge.data.frame(cell_lines[7:8, ], drugs[c(4:6), ], by = NULL))
  dt_layout <- gDRtestData::add_data_replicates(dt_layout)
  dt_layout <- gDRtestData::add_concentration(dt_layout, concentrations = 10 ^ (seq(-3, .5, .5)))
  
  dt_2 <-
    data.table::as.data.table(merge.data.frame(
      cell_lines[cell_lines$clid %in% dt_layout$clid, ], 
      drugs[c(21, 26), ], 
      by = NULL
    ))
  dt_2 <- gDRtestData::add_data_replicates(dt_2)
  dt_2 <- gDRtestData::add_concentration(dt_2, concentrations = 10 ^ (seq(-3, .5, .5)))
  colnames(dt_2)[colnames(dt_2) %in% c(colnames(drugs), "Concentration")] <- 
    paste0(colnames(dt_2)[colnames(dt_2) %in% c(colnames(drugs), "Concentration")], "_2")
  
  dt_layout_2 <- dt_layout[dt_2, on = intersect(names(dt_layout), names(dt_2)),
                           allow.cartesian = TRUE]
  
  dt_merged_data <- gDRtestData::generate_response_data(dt_layout_2, 0)
  dt <- identify_data_type(dt_merged_data)
  expect_equal(ncol(dt_merged_data), ncol(dt))
  expect_true("type" %in% names(dt))
  
  dt_list <- split_raw_data(dt)
  expect_true(inherits(dt_list, "list"))
  expect_true(all(names(dt_list) %in% c(gDRutils::get_supported_experiments("combo"),
                                        gDRutils::get_supported_experiments("sa"))))
  
  
  dt2 <- data.table::data.table(Gnumber = c(rep("DrugA", 9), "DrugB"),
                                DrugName = c(rep("DrugA", 9), "DrugB"),
                                drug_moa = "unknown",
                                drug_moa_2 = "unknown",
                                Gnumber_2 = c(rep("DrugB", 9), "DrugA"),
                                DrugName_2 = c(rep("DrugB", 9), "DrugA"),
                                Concentration = runif(10) + 0.01,
                                Concentration_2 = runif(10) + 0.01,
                                clid = "CL1")
  dt2 <- identify_data_type(dt2)
  split_dt2 <- split_raw_data(dt2)
  expect_true(all(split_dt2[["Gnumber"]] == "DrugA"))
})


test_that("collapse drugs works as expected", {
  idts <- gDRutils::get_env_identifiers()
  cols <- idts[c("drug", "drug2", "drug3",
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

