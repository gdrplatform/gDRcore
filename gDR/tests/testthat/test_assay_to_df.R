library(testthat)
testthat::context("Test assay_to_df function")

testthat::test_that("Test assay_to_df works", {

seL_file <- system.file(package = "gDR", "testdata/seL.rds")
seL <- readRDS(seL_file)
se <- seL[[3]]
test1 <-  assay_to_df(se, "Avg_Controls")
test2 <-  assay_to_df(se, 4)
expect_equal(test1, test2)

test1 <- test1 %>%
  dplyr::select(-CellLineName, -E2)
test2 <- test2 %>%
  dplyr::select(-CellLineName, -E2)
expect_equal(test1, test2)

experiment_id <- 414
BASE_URL <- "http://reswebappdev301.gene.com:28294"
cell_line_info_tbl <- readRDS(system.file(package = "gDR", "testdata/cell_line_info_tbl.RDS"))
drugs_tbl <- readRDS(system.file(package = "gDR", "testdata/drugs_tbl.RDS"))
## FROM gDRwrapper:
ridTbl <-
  tibble::as_tibble(SummarizedExperiment::rowData(se)) %>% tibble::add_column(rId = rownames(SummarizedExperiment::rowData(se)))
cidTbl <-
  tibble::as_tibble(SummarizedExperiment::colData(se)) %>% tibble::add_column(cId = rownames(SummarizedExperiment::colData(se))) %>%
  dplyr::select_if(!names(.) %in% c("CellLineName", "Tissue", "UserCellLineDesignation"))
combs <-
  expand.grid(ridTbl$rId, cidTbl$cId, stringsAsFactors = FALSE)
colnames(combs) <- c("rId", "cId")
fTbl <- dplyr::left_join(combs, ridTbl, by = "rId") %>%
  dplyr::left_join(cidTbl, by = "cId")
#keys that are defining the unique condtions (usually only 'Duration' and 'clid')
condKeyCols <-
  colnames(fTbl)[!colnames(fTbl) %in% c("rId", "cId", "Gnumber", "DrugName")]
# add condition_ident and treatment_ident (starting from 1)
# TODO: fetch latest 'condition_ident' and 'treatment_ident from db to defined appropriate id
fTbl <-
  seplyr::add_group_indices(fTbl, condKeyCols, "condition_ident") %>%
  dplyr::left_join(drugs_tbl, by = c("Gnumber")) %>%
  dplyr::left_join(cell_line_info_tbl, by = c("clid")) %>%
  tibble::add_column(experiment_id = experiment_id) %>%
  dplyr::mutate(treatment_ident = paste0(experiment_id, "_", condition_ident, "_", seq_len(dplyr::n())))

#get dictTbl
dictTbl <-
  fTbl %>% dplyr::select(rId, cId, condition_ident, treatment_ident) %>%
  dplyr::arrange(condition_ident)

asDf <- gDR::assay_to_df(se, 4)
reqCondCols <- c(
  "experiment_id" = NA,
  "cell_id" = NA,
  "Duration" = NA,
  "UntrtReadout" = NA,
  "Day0Readout" = NA,
  "DivisionTime" = NA,
  "GRbased_on_data" = NA,
  "Untrt_frac_dead_cells" = NA,
  "Day0_frac_dead_cells" = NA
)
finalCondCols <- c(
  "cell_id" = "cell_info_id",
  "Day0_frac_dead_cells" = "day_z_frac_dead_cells",
  "Day0Readout" = "day_z_readout",
  "DivisionTime" = "division_time",
  "GRbased_on_data" = "gr_based_on_data",
  "UntrtReadout" = "untrt_frac_dead_cells",
  "Untrt_frac_dead_cells" = "untrt_readout"
)

condTbl <-
  dplyr::left_join(fTbl, asDf, by = c("rId", "cId")) %>%
  dplyr::select_if(colnames(.) %in% c("condition_ident", names(reqCondCols))) %>%
  dplyr::distinct() %>%
  tibble::add_column(!!!reqCondCols[!names(reqCondCols) %in% names(.)]) %>%
  dplyr::arrange(condition_ident)

treatTbl <-
  fTbl %>% dplyr::select(condition_ident, treatment_ident, drug_id) %>%
  dplyr::group_by(condition_ident) %>%
  dplyr::group_nest(.key = "treatment_metadata")

# complete tibble with data for all the tables from the 'Condition' section
# currently only records for 'condition_metadata' and 'treatment_metadata' added
#TODO: add data for 'condition_codrug', 'condition_annotation' and 'condition_additional_treatments'

### condition_codrug section ###
codrugsTbl <- tibble::tibble()
codrugCols <-
  grep("^DrugName_|^Gnumber_|^Concentration_",
       colnames(fTbl),
       value = TRUE)
if (length(codrugCols)) {
  codrugsTbl <-
    fTbl[, c("rId", "cId", codrugCols)]  %>% dplyr::left_join(dictTbl) %>%
    dplyr::select(-treatment_ident,-rId,-cId) %>% dplyr::distinct()
  # codrugsTbl <-
  #   update_codrugs_tbl(codrugsTbl, base_url = base_url) %>% dplyr::arrange(condition_ident)
}

### condition_additional_treatments section ###
addTrtTbl <- tibble::tibble()
addTrtCols <-
  condKeyCols[!condKeyCols %in% c(codrugCols, "clid", "Duration") &
                condKeyCols %in% colnames(SummarizedExperiment::rowData(se))]
if (length(addTrtCols)) {
  addTrtTbl <-
    fTbl[, c("rId", "cId", addTrtCols)]  %>% dplyr::left_join(dictTbl) %>%
    dplyr::select(-treatment_ident,-rId,-cId) %>% dplyr::distinct() %>% dplyr::arrange(condition_ident)
  
  addTrtTbl <- tidyr::gather(addTrtTbl,
                             key = "metadata_field",
                             value = "metadata_value",
                             -dplyr::one_of("condition_ident")) %>% dplyr::arrange(condition_ident) %>%
    dplyr::group_by(condition_ident) %>%
    dplyr::group_nest(.key = "add_trt")
}

### condition_cell_annots section ###
cellAnnotTbl <- tibble::tibble()
cellAnnotCols <-
  condKeyCols[!condKeyCols %in% c(codrugCols, "clid", "Duration") &
                condKeyCols %in% colnames(SummarizedExperiment::colData(se))]
if (length(cellAnnotCols)) {
  cellAnnotTbl <-
    fTbl[, c("rId", "cId", cellAnnotCols)]  %>% dplyr::left_join(dictTbl) %>%
    dplyr::select(-treatment_ident, -rId, -cId) %>% dplyr::distinct() %>% dplyr::arrange(condition_ident)
  
  cellAnnotTbl <- tidyr::gather(
    cellAnnotTbl,
    key = "annotation_field",
    value = "annotation_value",-dplyr::one_of("condition_ident")
  ) %>% dplyr::arrange(condition_ident) %>%
    dplyr::group_by(condition_ident) %>%
    dplyr::group_nest(.key = "cell_annot")
}

cTbl <-
  dplyr::left_join(condTbl, treatTbl, by = c("condition_ident")) %>% dplyr::select(-condition_ident)  %>%
  dplyr::rename_at(dplyr::vars(names(finalCondCols)), ~ as.character(finalCondCols)) %>%
  dplyr::rename_all(tolower)

# TODO: remove once we agree on data format in condition metadata
# dirty fixes for REST API compatiblity
cTbl$cell_info <-
  lapply(cTbl$cell_info_id, function(x) {
    list(cell_info_id = x)
  })
cTbl$experiment <-
  lapply(cTbl$experiment_id, function(x) {
    list(experiment_id = x)
  })

for (x in seq(length(cTbl$treatment_metadata))) {
  tDrugs <- cTbl$treatment_metadata[[x]]$drug_id
  cTbl$treatment_metadata[[x]]$drug <-
    lapply(tDrugs, function(x) {
      list(drug_id = x)
    })
  cTbl$treatment_metadata[[x]]$drug_id <- NULL
}

cTbl %<>% dplyr::select(-cell_info_id, -experiment_id)
#%>% tibble::as.tibble()


metadataTable2 <- list(dictTbl = dictTbl, condTbl = cTbl, codrugsTbl = codrugsTbl, addTrtTbl = addTrtTbl, cellAnnotTbl = cellAnnotTbl)
checkmate::assert_list(metadataTable2)
})