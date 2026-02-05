test_that("map_df works as expected", {
  n <- 64
  md_df <- data.table::data.table(
    Gnumber = rep(c("vehicle", "untreated", paste0("G", seq(2))), each = 16), 
    DrugName = rep(c("vehicle", "untreated", paste0("GN", seq(2))), each = 16), 
    clid = paste0("C", rep_len(seq(4), n)),
    CellLineName = paste0("N", rep_len(seq(4), n)),
    replicates = rep_len(paste0("R", rep(seq(4), each = 4)), 64),
    Duration = 160,
    rn = as.character(seq_len(64))
  )
  
  ref <- md_df$Gnumber %in% c("vehicle", "untreated")
  ref_df <- md_df[ref, ]
  trt_df <- md_df[!ref, ]
  
  ref_cols <- c("clid", "CellLineName", "Duration")
  mapping <- map_df(trt_df, ref_df, ref_cols = ref_cols, ref_type = "untrt_Endpoint")
  
  expect_equal(names(mapping), trt_df$rn)
  expect_equal(length(mapping), nrow(trt_df))
  expect_equal(sort(unique(unname(unlist(mapping)))), sort(ref_df$rn))
})

test_that("Best match is detected for missing controls", {
  trt_dt <- data.table::data.table(
    clid = "C1",
    Duration = 72,
    Tissue = "Lung",
    rn = "T1"
  )
  
  ref_dt <- data.table::data.table(
    clid = c("C1", "C1"),
    Duration = c(72, 72),
    Tissue = c("Lung", "Liver"),
    rn = c("R1", "R2")
  )
  
  ref_cols <- c("clid", "Tissue")
  obs <- map_df(trt_dt, ref_dt, ref_cols = ref_cols, ref_type = "untrt_Endpoint")
  
  expect_equal(obs[["T1"]], "R1")
})

test_that("NAs are returned for missing controls", {
  trt_dt <- data.table::data.table(
    clid = "Cell_A",
    Duration = 72,
    rn = "T1"
  )
  
  ref_dt <- data.table::data.table(
    clid = "Cell_B",
    Duration = 72,
    rn = "R1"
  )
  
  ref_cols <- c("clid")
  obs <- map_df(trt_dt, ref_dt, ref_cols = ref_cols, ref_type = "untrt_Endpoint")
  
  expect_list(obs)
  expect_equal(length(obs[["T1"]]), 0)
})

test_that(".map_references works as expected", {
  mat_elem <- data.table::data.table(
    DrugName = rep(c("untreated", "drugA", "drugB", "untreated"), 2),
    DrugName_2 = rep(c("untreated", "vehicle", "drugA", "drugB"), 2),
    clid = rep(c("C1", "C2"), each = 4)
  )
  
  # Combination data mapping
  obs <- .map_references(mat_elem, rowData_colnames = c("DrugName", "DrugName_2"))
  exp <- list("3" = c("2", "4"), "7" = c("6", "8"))
  expect_equal(obs, exp)
  
  # Single-agent subsetting
  colname <- c("DrugName", "clid")
  mat_elem2 <- mat_elem[, colname, with = FALSE]
  obs2 <- .map_references(mat_elem2, rowData_colnames = "DrugName")
  
  expect_list(obs2)
  expect_true(all(vapply(obs2, is.null, FUN.VALUE = logical(1))))
  
  expect_error(.map_references(mat_elem[, .(DrugName)]), "element of")
})

test_that("map_untreated works as expected", {
  mat_elem <- data.table::data.table(
    DrugName = c("untreated", "drugA", "untreated"),
    DrugName_2 = c("untreated", "untreated", "drugB"),
    DrugName_3 = c("untreated", "untreated", "untreated")
  )
  
  obs <- map_untreated(mat_elem)
  
  expect_logical(obs, any.missing = FALSE, len = 3)
  expect_equal(sum(obs), 1) # Only first row is untreated in all 3 cols
  expect_true(obs[1])
})

test_that(".get_untreated_tag_count helper and assertions work", {
  dt_valid <- data.table::data.table(
    DrugName = c("untreated", "drugA"),
    DrugName_2 = c("untreated", "untreated"),
    clid = "C1"
  )
  
  res <- .get_untreated_tag_count(dt_valid, c("drug_name", "drug_name2"))
  expect_equal(res$ntag, c(2, 1))
  
  expect_error(.get_untreated_tag_count(dt_valid, drug_identifier_keys = "empty"), 
               "Must be element of set")
  expect_error(.get_untreated_tag_count(dt_valid, drug_identifier_keys = "drug"),
               "None of the drug")
})

