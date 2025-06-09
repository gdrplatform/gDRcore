test_that("map_df works as expected", {
  md_df <- unique(md_df)
  md_df$rn <- as.character(seq_len(nrow(md_df)))
  ref <- md_df$Gnumber %in% c("vehicle", "untreated")
  
  ref_df <- md_df[ref, ]
  trt_df <- md_df[!ref, ]
  Keys <- identify_keys(test_df)

  ref_type <- "untrt_Endpoint"
  mapping <- map_df(trt_df, 
                    ref_df, 
                    ref_cols = Keys[[ref_type]],
                    ref_type = ref_type)

  expect_equal(names(mapping), trt_df$rn)
  expect_equal(length(mapping), nrow(trt_df))

  out <- gDRutils::loop(seq(33, 64, 1), function(x) {
    as.character(c((x - 1) %% 16 + 1, 17 + ((x - 1) %% 16)))
    })
  names(out) <- seq(33, 64, 1)

  expect_equal(mapping, out)
  expect_equal(sort(unique(unname(unlist(mapping)))), sort(rownames(ref_df)))
  
  # Test Day0 data with E2
  
  trt_dt <- data.table::data.table(
    clid = c("Cell123456", "Cell123456", "Cell654321", "Cell654321"),
    CellLineName = c("CellLineA", "CellLineA", "CellLineB", "CellLineB"),
    Tissue = c("Liver", "Liver", "Liver", "Liver"),
    parental_identifier = c("ParentA", "ParentA", "ParentB", "ParentB"),
    subtype = c("type1", "type1", "type2", "type2"),
    ReferenceDivisionTime = c(45, 45, 30, 30),
    Gnumber = c("Drug12345678", "Drug12345678", "Drug87654321", "Drug87654321"),
    DrugName = c("DrugX", "DrugX", "DrugY", "DrugY"),
    drug_moa = c("MOA1", "MOA1", "MOA2", "MOA2"),
    Duration = c(120, 120, 120, 120),
    E2 = c("0", "0.0023", "0", "0.0023"),
    rn = c("1", "2", "6", "7")
  )
  
  ref_dt <- data.table::data.table(
    clid = c("Cell123456", "Cell123456", "Cell123456", "Cell654321", "Cell654321", "Cell654321"),
    CellLineName = c("CellLineA", "CellLineA", "CellLineA", "CellLineB", "CellLineB", "CellLineB"),
    Tissue = c("Liver", "Liver", "Liver", "Liver", "Liver", "Liver"),
    parental_identifier = c("ParentA", "ParentA", "ParentA", "ParentB", "ParentB", "ParentB"),
    subtype = c("type1", "type1", "type1", "type2", "type2", "type2"),
    ReferenceDivisionTime = c(NA, NA, NA, 60, 60, 60),
    Gnumber = c("Drug12345678", "Drug12345678", "Drug12345678", "Drug87654321", "Drug87654321", "Drug87654321"),
    DrugName = c("DrugX", "DrugX", "DrugX", "DrugY", "DrugY", "DrugY"),
    drug_moa = c("MOA1", "MOA1", "MOA1", "MOA2", "MOA2", "MOA2"),
    Duration = c(0, 120, 120, 0, 120, 120),
    E2 = c("0", "0", "0.0023", "0", "0", "0.0023"),
    rn = c("3", "4", "5", "8", "9", "10")
  )
  
  ref_cols <- c("CellLineName", "Tissue", "parental_identifier", "subtype", 
                "Barcode", "clid")
  
  ref_type <- "Day0"
  map_override_untrt_controls <- map_df(trt_dt,
                                        ref_dt,
                                        ref_cols = ref_cols,
                                        ref_type = ref_type,
                                        override_untrt_controls = c(E2 = 0.0023))
  expect_list(map_override_untrt_controls)
  expect_length(map_override_untrt_controls, 4)
  expect_equal(as.numeric(unlist(map_override_untrt_controls)), c(3, 3, 8, 8))
  
  
  map_override_untrt_controls2 <- map_df(trt_dt,
                                         ref_dt,
                                         ref_cols = ref_cols,
                                         ref_type = ref_type,
                                         override_untrt_controls = NULL)
  
  expect_list(map_override_untrt_controls2)
  expect_length(map_override_untrt_controls2, 4)
  expect_equal(as.numeric(unlist(map_override_untrt_controls2)), c(3, 3, 8, 8))
})


# TODO: test_that("Best match is detected for missing controls", {}) # nolint

# TODO: test_that("NAs are returned for missing controls", {}) # nolint

mat_elem <- data.table::data.table(DrugName = rep(c("untreated", "drugA", "drugB", "untreated"), 2),
                                   DrugName_2 = rep(c("untreated", "vehicle", "drugA", "drugB"), 2),
                                   clid = rep(c("C1", "C2"), each = 4))

test_that(".map_references works as expected", {
  # Combination data.
  obs <- .map_references(mat_elem, rowData_colnames = c("DrugName", "DrugName_2"))
  exp <- list("3" = c("2", "4"), "7" = c("6", "8"))
  expect_equal(obs, exp)

  # Single-agent data.
  colname <- c("DrugName", "clid")
  mat_elem2 <- mat_elem[, colname, with = FALSE]
  obs2 <- .map_references(mat_elem2, rowData_colnames = c("DrugName", "DrugName_2"))
  exp2 <- list("2" = NULL, "3" = NULL, "6" = NULL, "7" = NULL)
  expect_equal(obs2, exp2)
})

test_that("map_untreated works as expected", {
  # Combination data.
  obs <- map_untreated(mat_elem)
  expect_equal(sum(obs), 2)
})
