prism_se <- function() {
  mae <- gDRutils::get_synthetic_data("finalMAE_prism")
  MultiAssayExperiment::experiments(mae)[["single-agent"]]
}
# four PRISM curves for three cell lines of the test data, as read by gDRimport::read_PRISM_drc()
prism_curves <- function() {
  data.table::data.table(
    pool_id = c("P1", "P1", "P2", "P1"),
    depmap_id = c("ACH-000004", "ACH-000005", "ACH-000005", "ACH-000006"),
    lua = c("LUA-0001", "LUA-0002", "LUA-0003", "LUA-0004"),
    cell_set = "CS1", pert_name = "G00106", pert_id = "PRC-000000001-001-01", day = 5,
    x_project_id = "TEST_PROJECT", pert_plate = c("PLATE1", "PLATE1", "PLATE2", "PLATE1"),
    fit_name = c("drc_drm_constrained", "dr4pl_initL_constrained", "drc_drm_unconstrained", ""),
    lower_limit = c(0.21, 0.35, 0.29, NA), upper_limit = c(0.98, 1.01, 0.96, NA),
    slope = c(-1.4, -0.9, -1.1, NA), inflection = c(0.85, 2.4, 1.7, NA),
    mse = c(0.011, 0.014, 0.012, NA), mad = c(0.08, 0.09, 0.07, NA),
    frac_var_explained = c(0.91, 0.77, 0.83, NA), successful_fit = c(TRUE, TRUE, TRUE, FALSE),
    auc_riemann = c(0.58, 0.71, 0.66, 0.97), minimum_dose = 0.00046, maximum_dose = 10,
    auc = c(0.6012, 0.7208, 0.6705, NA), log2_auc = c(-0.7341, -0.4723, -0.5767, NA),
    log2_ic50 = c(0.1043, NA, 0.9543, NA), ic50 = c(1.0750, NA, 1.9378, NA), screen = "TEST_SCREEN",
    clid = c("HEL_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
             "HEL9217_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "MONOMAC6_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"),
    Duration = 120
  )
}
prism_rows <- function(se) {
  m <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  m[m$fit_source == "PRISM", ]
}

test_that("add_PRISM_metrics adds every PRISM curve with the gDR curve parameters", {
  se <- prism_se()
  gdr_before <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  curves <- prism_curves()
  out <- add_PRISM_metrics(se, curves)

  prism <- prism_rows(out)
  expect_identical(NROW(prism), NROW(curves))
  data.table::setorder(prism, lua)
  expect_identical(prism$lua, curves$lua)
  expect_identical(unique(prism$normalization_type), "RV")
  expect_equal(prism$x_0, curves$upper_limit)
  expect_equal(prism$x_inf, curves$lower_limit)
  expect_equal(prism$ec50, curves$inflection)
  expect_equal(prism$h, -curves$slope)
  expect_equal(prism$r2, curves$frac_var_explained)
  expect_equal(prism$maxlog10Concentration, log10(curves$maximum_dose))
  expect_equal(prism$prism_auc, curves$auc)
  expect_equal(prism$prism_ic50, curves$ic50)
  expect_identical(prism$fit_type, c(curves$fit_name[1:3], NA))
  expect_true(all(is.na(prism$xc50)) && all(is.na(prism$x_AOC)))
  # two curves of one cell line, told apart by pool, barcode and compound plate
  expect_identical(sum(prism$depmap_id == "ACH-000005"), 2L)
  expect_true(all(prism$N_conc > 0))

  # the gDR fits stay as they were
  gdr_after <- gDRutils::convert_se_assay_to_dt(out, "Metrics")
  gdr_after <- gdr_after[gdr_after$fit_source == "gDR", names(gdr_before), with = FALSE]
  data.table::setorderv(gdr_after, c("rId", "cId", "normalization_type"))
  data.table::setorderv(gdr_before, c("rId", "cId", "normalization_type"))
  expect_equal(gdr_after, gdr_before)
})

test_that("the added parameters give PRISM's curve", {
  curves <- prism_curves()[1]
  p <- prism_rows(add_PRISM_metrics(prism_se(), curves))
  dose <- 10 ^ seq(-3, 1, 0.5)
  prism_curve <- curves$upper_limit +
    (curves$lower_limit - curves$upper_limit) / (1 + (dose / curves$inflection) ^ curves$slope)
  expect_equal(gDRutils::predict_efficacy_from_conc(dose, p$x_inf, p$x_0, p$ec50, p$h), prism_curve)
})

test_that("add_PRISM_metrics replaces the PRISM rows when run again and records the fit", {
  curves <- prism_curves()
  data.table::setattr(curves, "gDR_import", list(
    files = data.frame(role = "prism_drc", name = "drc.csv", md5 = "abc", stringsAsFactors = FALSE)
  ))
  out <- add_PRISM_metrics(prism_se(), curves)
  n <- NROW(gDRutils::convert_se_assay_to_dt(out, "Metrics"))
  again <- add_PRISM_metrics(out, curves)
  expect_identical(NROW(gDRutils::convert_se_assay_to_dt(again, "Metrics")), n)

  prov <- S4Vectors::metadata(again)[["provenance"]]
  expect_true("PRISM" %in% names(prov$fits[["Metrics"]]))
  step <- prov$steps[["prism_drc_import"]]$settings
  expect_identical(step$n_curves, 4L)
  expect_identical(step$n_imported, 4L)
  expect_identical(step$files$md5, "abc")
  expect_identical(step$screen, "TEST_SCREEN")
})

test_that("add_PRISM_metrics leaves out curves without a cell and maps drug identifiers", {
  curves <- prism_curves()
  curves$pert_name <- "BRD-0001"
  extra <- curves[1]
  extra$clid <- "NOT_IN_THE_SE"
  curves <- rbind(curves, extra)

  expect_warning(out <- add_PRISM_metrics(prism_se(), curves, drug_map = c(`BRD-0001` = "G00106")),
                 "1 of 5 PRISM curves have no matching cell")
  expect_identical(NROW(prism_rows(out)), 4L)
  expect_identical(S4Vectors::metadata(out)$provenance$steps$prism_drc_import$settings$n_left_out, 1L)
  # without the map, no drug matches
  expect_error(purrr::quietly(add_PRISM_metrics)(prism_se(), curves), "no PRISM curve matches")
  expect_error(add_PRISM_metrics(prism_se(), curves[, -"inflection"]), "inflection")
})

test_that("add_PRISM_metrics works on the curves read by gDRimport", {
  testthat::skip_if_not(requireNamespace("gDRimport", quietly = TRUE) &&
                          exists("read_PRISM_drc", asNamespace("gDRimport")))
  curves <- gDRimport::read_PRISM_drc(
    system.file("testdata/prism_drc.csv", package = "gDRimport"),
    system.file("testdata/prism_drc_model.csv", package = "gDRimport")
  )
  out <- add_PRISM_metrics(prism_se(), curves)
  expect_identical(NROW(prism_rows(out)), NROW(curves))
  expect_identical(S4Vectors::metadata(out)$provenance$steps$prism_drc_import$settings$files$role,
                   c("prism_drc", "meta_data"))
})
