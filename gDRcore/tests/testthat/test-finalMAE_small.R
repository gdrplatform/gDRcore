data <- "finalMAE_small.RDS"
original <- get_synthetic_data(data)

mae <- gDRtestData::generateNoiseRawData(cell_lines, drugs, e_inf, ec50, hill_coef)

test_synthetic_data(original, mae, data)
