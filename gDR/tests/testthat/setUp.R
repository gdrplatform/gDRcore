test_df <- data.frame(
  Gnumber = paste0("G", rep(seq(5), 4)), 
  DrugName = paste0("GN", rep(seq(5), 4)),
  drug_moa = "inhibitor",
  clid = paste0("C", rep(seq(2), each = 10)),
  CellLineName = paste0("N", rep(seq(2), each = 10)),
  Tissue = "Lung",
  ReferenceDivisionTime = rep(c(120, 60), each = 10),
  parental_identifier = "CL12345",
  duration = 160, 
  replicates = paste0("R", rep(rep(seq(2), 2), each = 5)),
  Concentration = sample(seq(100), 20), 
  WellRow = rep_len(LETTERS[1:8], 20),
  WellColumn = rep_len(seq(3), 20),
  experimenter = "Bob Ross"
)
