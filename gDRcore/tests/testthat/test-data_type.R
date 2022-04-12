test_that("identify_data_type and split_raw_data works as expected", {
  df_layout <- merge(cell_lines[7:8, ], drugs[c(4:6), ], by = NULL)
  df_layout <- gDRtestData::add_data_replicates(df_layout)
  df_layout <- gDRtestData::add_concentration(df_layout, concentrations = 10 ^ (seq(-3, .5, .5)))
  
  df_2 <- merge(cell_lines[cell_lines$clid %in% df_layout$clid, ], drugs[c(21, 26), ], by = NULL)
  df_2 <- gDRtestData::add_data_replicates(df_2)
  df_2 <- gDRtestData::add_concentration(df_2, concentrations = 10 ^ (seq(-3, .5, .5)))
  colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")] <- 
    paste0(colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")], "_2")
  
  df_layout_2 <- merge(df_layout, df_2)
  
  df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)
  df <- identify_data_type(df_merged_data)
  expect_equal(ncol(df_merged_data) + 1, ncol(df))
  expect_true("type" %in% names(df))
  
  df_list <- split_raw_data(df)
  expect_true(inherits(df_list, "list"))
  expect_true(all(names(df_list) %in% c("matrix", "single-agent")))
})

test_that("duplicate_single_agent works as expected", {
  df_ <- data.frame(Gnumber=c("DrugA", "vehicle"), Concentration=c(0.3, 0), Gnumber_2=c("vehicle", "DrugA"), Concentration_2=c(0, 0.3))
  expect_equal(duplicate_single_agent2(df_), df_)
  expect_equal(duplicate_single_agent2(df_[1, ]), df_)
  expect_equal(duplicate_single_agent(df_[2, ]), df_)
})
