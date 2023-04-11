test_that("combo_matrix_small: test_synthetic_data", {
  data <- "finalMAE_combo_matrix_small.RDS"
  original <- gDRutils::get_synthetic_data(data)
  
  set.seed(2)
  mae <- purrr::quietly(gDRtestData::generateComboMatrixSmall)(
    cell_lines, drugs, FALSE
  )
  expect_length(mae$warnings, 4)

  test_synthetic_data(original, mae$result, data)
})


input_df <- convert_mae_to_raw_data(mae$result)
input_df$masked <- NULL
mae2 < runDrugResponseProcessingPipeline(input_df)
mae3 <- runDrugResponseProcessingPipeline(df_merged)


avg <- convert_se_assay_to_dt(original[[1]], "Normalized")
avg1 <- convert_se_assay_to_dt(mae$result[[1]], "Normalized")
avg2 <- convert_se_assay_to_dt(mae2[[1]], "Normalized")
avg3 <- convert_se_assay_to_dt(mae3[[1]], "Normalized")

avg <- convert_se_assay_to_dt(original[[1]], "Metrics")
avg1 <- convert_se_assay_to_dt(mae[[1]], "Metrics")
avg2 <- convert_se_assay_to_dt(mae2[[1]], "Metrics")
avg3 <- convert_se_assay_to_dt(mae3[[1]], "Metrics")


summary(avg$x_mean)
summary(avg1$x_mean)
summary(avg2$x_mean)
summary(avg3$x_mean)


library(dplyr)
library(tidyr)


table(input_df$Gnumber)
table(df_merged$Gnumber)

table(input_df$Gnumber_2)
table(df_merged$Gnumber_2)


df1_long <- input_df %>% pivot_longer(cols = everything())
df2_long <- df_merged %>% pivot_longer(cols = everything())

# Sort both data frames by variable and value columns
df1_long_sorted <- df1_long %>% arrange(name, value)
df2_long_sorted <- df2_long %>% arrange(name, value)

# Compare the two data frames
identical(df1_long_sorted, df2_long_sorted)





mae2 <- runDrugResponseProcessingPipeline(as.data.frame(input_df))
