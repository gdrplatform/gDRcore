data <- "finalMAE_combo_2dose_nonoise2.RDS"
original <- get_synthetic_data(data)

mae <- gDRtestData::generateComboNoNoiseData2(cell_lines, drugs, e_inf, ec50, hill_coef)

test_synthetic_data(original, mae, data)
