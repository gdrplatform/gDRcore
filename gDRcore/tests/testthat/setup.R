n <- 64
md_df <- data.frame(
  Gnumber = rep(c("vehicle", "untreated", paste0("G", seq(2))), each = 16), 
  DrugName = rep(c("vehicle", "untreated", paste0("GN", seq(2))), each = 16), 
  clid = paste0("C", rep_len(seq(4), n)),
  CellLineName = paste0("N", rep_len(seq(4), n)),
  replicates = rep_len(paste0("R", rep(seq(4), each = 4)), 64),
  drug_moa = "inhibitor",
  ReferenceDivisionTime = rep_len(c(120, 60), n),
  Tissue = "Lung",
  parental_identifier = "CL12345",
  Duration = 160
)

data_df <- data.frame(
  Concentration = rep(c(0, 0, 1, 3), each = 16), 
  ReadoutValue = runif(n, 1000, 5000),
  BackgroundValue = runif(n, 0, 1),
  WellRow = rep_len(LETTERS[1:8], n),
  WellColumn = rep_len(seq(3), n),
  experimenter = "Bob Ross"
)

test_df <- cbind(md_df, data_df)

cell_lines <- gDRtestData::create_synthetic_cell_lines()
drugs <-  gDRtestData::create_synthetic_drugs()
