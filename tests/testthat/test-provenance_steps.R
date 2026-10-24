####
# Provenance of the processing steps: create_SE(), normalize_SE(), average_SE()
####

raw_small <- convert_mae_to_raw_data(gDRutils::get_synthetic_data("finalMAE_small.qs2"))

.run_pipeline <- function(raw) {
  suppressWarnings(suppressMessages(runDrugResponseProcessingPipeline(raw)))[["single-agent"]]
}


test_that("the pipeline records create, normalize and average next to the fit", {
  se <- .run_pipeline(raw_small)
  prov <- gDRutils::get_SE_provenance(se)
  expect_identical(names(prov$steps), c("create", "normalize", "average"))
  expect_identical(names(prov$fits), "Metrics")

  create <- prov$steps$create$settings
  expect_identical(create$readout, "ReadoutValue")
  expect_identical(create$masked$n_wells, 0L)

  normalize <- prov$steps$normalize
  expect_identical(normalize$settings$ndigit_rounding, 4)
  expect_match(normalize$functions$control_mean_fxn$source, "trim = 0.25", fixed = TRUE)

  average <- prov$steps$average
  expect_identical(average$functions$average_FUN$name, "average_FUN")
  expect_identical(average$functions$average_FUN$package, "gDRcore")

  # the steps stay small, and running the pipeline again gives the same object
  expect_lt(length(qs2::qs_serialize(prov)), 4000)
  se_again <- .run_pipeline(raw_small)
  expect_identical(se_again, se)
})


test_that("re-running the pipeline on its own output records the same steps", {
  # reprocessing a processed MAE must give the same object as processing the raw data,
  # so the steps may only record settings, not the form the input came in
  raw <- data.table::copy(raw_small)
  raw$BackgroundValue <- 0
  mae <- suppressWarnings(suppressMessages(runDrugResponseProcessingPipeline(raw)))
  mae_again <- suppressWarnings(suppressMessages(runDrugResponseProcessingPipeline(mae)))
  expect_identical(
    gDRutils::get_SE_provenance(mae_again[["single-agent"]]),
    gDRutils::get_SE_provenance(mae[["single-agent"]])
  )
})


test_that("create_SE counts the masked wells it drops", {
  raw <- data.table::copy(raw_small)
  masked_col <- gDRutils::get_env_identifiers("masked_tag")
  cl_col <- gDRutils::get_env_identifiers("cellline")
  untrt <- gDRutils::get_env_identifiers("untreated_tag")
  drug_col <- gDRutils::get_env_identifiers("drug")
  treated <- which(!raw[[drug_col]] %in% untrt)
  to_mask <- treated[1:7]
  data.table::set(raw, j = masked_col, value = FALSE)
  data.table::set(raw, i = to_mask, j = masked_col, value = TRUE)

  masked <- gDRutils::get_SE_provenance(.run_pipeline(raw))$steps$create$settings$masked
  expect_identical(masked$n_wells, 7L)
  expect_identical(sum(masked$by_cell_line), 7L)
  expect_setequal(names(masked$by_cell_line), unique(raw[[cl_col]][to_mask]))
})


test_that("normalize_SE records which cell lines got no GR value", {
  raw <- data.table::copy(raw_small)
  div_col <- gDRutils::get_env_identifiers("cellline_ref_div_time")
  cl_name_col <- gDRutils::get_env_identifiers("cellline_name")
  first_line <- sort(unique(raw[[cl_name_col]]))[1L]
  data.table::set(raw, i = which(raw[[cl_name_col]] == first_line), j = div_col, value = NA_real_)

  gr <- gDRutils::get_SE_provenance(.run_pipeline(raw))$steps$normalize$settings$gr
  expect_identical(gr$no_gr_cell_lines, first_line)
  expect_false(first_line %in% gr$ref_div_time_cell_lines)
  expect_identical(sum(gr$n_pairs), sum(gr$n_pairs[c("ref_div_time", "none")]))
  expect_gt(gr$n_pairs[["none"]], 0L)
})


test_that("the message about missing GR values stays one short line", {
  expect_identical(.no_gr_message(c("A", "B")), "could not calculate GR values for 2 cell line(s): 'A', 'B'")
  msg <- .no_gr_message(sprintf("CL%03d", 1:300))
  expect_match(msg, "for 300 cell line(s): 'CL001', 'CL002', 'CL003', 'CL004', 'CL005' and 295 more", fixed = TRUE)
})


test_that(".gr_method follows the branches of calculate_GR_value()", {
  expect_identical(.gr_method(c(1, 2), duration = 72, ref_div_time = 30), "day0")
  expect_identical(.gr_method(c(1, NA), duration = 72, ref_div_time = 30), "ref_div_time")
  expect_identical(.gr_method(NA, duration = 72, ref_div_time = NA), "none")
  expect_identical(.gr_method(NA, duration = 72, ref_div_time = 120), "none")
  expect_identical(.gr_method(NA, duration = NA, ref_div_time = 30), "none")
})


test_that("the import record of the input data becomes the first step", {
  raw <- data.table::copy(raw_small)
  files <- data.frame(role = c("manifest", "results"), name = c("manifest.xlsx", "RawData.xlsx"),
                      size = c(100, 200), md5 = c(strrep("a", 32), strrep("b", 32)))
  record <- list(gDRimport = "1.11.9", instrument = "EnVision", files = files)
  data.table::setattr(raw, "gDR_import", record)

  mae <- suppressWarnings(suppressMessages(runDrugResponseProcessingPipeline(raw)))
  prov <- gDRutils::get_SE_provenance(mae[["single-agent"]])
  expect_identical(names(prov$steps), c("import", "create", "normalize", "average"))
  expect_identical(prov$steps$import$settings, record)

  # re-running on the processed MAE keeps the import step it recorded
  mae_again <- suppressWarnings(suppressMessages(runDrugResponseProcessingPipeline(mae)))
  expect_identical(gDRutils::get_SE_provenance(mae_again[["single-agent"]]), prov)
})
