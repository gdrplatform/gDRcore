#' Add the curve fits delivered by PRISM to the Metrics assay
#'
#' Adds the curves PRISM fitted (read with \code{gDRimport::read_PRISM_drc()}) to the
#' \code{Metrics} assay of processed PRISM data, with \code{fit_source = "PRISM"}, next to the
#' gDR fits. The curves are not refitted, so the values match the PRISM delivery exactly.
#'
#' PRISM fits \code{upper_limit + (lower_limit - upper_limit) / (1 + (dose / inflection)^slope)},
#' which is the gDR curve with \code{x_0 = upper_limit}, \code{x_inf = lower_limit},
#' \code{ec50 = inflection} and \code{h = -slope}. \code{r2} is \code{frac_var_explained} and
#' \code{fit_type} is PRISM's \code{fit_name}. PRISM defines AUC and IC50 differently from gDR
#' (the curve is clamped between 0 and 1, an IC50 outside the tested range is \code{NA}), so
#' those values are kept under their own names (\code{prism_auc}, \code{prism_ic50} and so on)
#' and the gDR columns \code{x_AOC}, \code{xc50}, \code{x_mean} and \code{x_max} stay \code{NA}.
#'
#' A PRISM screen often has several curves for one cell line and compound (different pools,
#' barcodes or compound plates), which gDR averages into one cell. Each of them becomes its
#' own row, with \code{depmap_id}, \code{pert_id}, \code{pool_id}, \code{lua}, \code{cell_set}
#' and \code{pert_plate} to tell them apart. Curves with no matching cell in \code{se} are left
#' out with a warning. Running the function again replaces the PRISM rows of the earlier run.
#'
#' Run it after \code{\link{fit_SE}}, which replaces the whole \code{Metrics} assay.
#'
#' @param se \code{SummarizedExperiment} of single-agent PRISM data, processed by gDR.
#' @param prism_curves \code{data.table} of PRISM curves, as returned by
#'   \code{gDRimport::read_PRISM_drc()}: the columns of the DRC table, the cell line identifier
#'   and Duration.
#' @param drug_map named character vector mapping PRISM's \code{pert_name} to the drug
#'   identifier in \code{se}, for data imported with other drug identifiers. By default the
#'   \code{pert_name} is the drug identifier, as in the PRISM LEVEL5 import.
#' @param normalization_type string, the normalization type of the PRISM curves.
#' @param metrics_assay string, name of the assay with the fits.
#' @param averaged_assay string, name of the assay used to count the concentrations
#'   (\code{N_conc}).
#'
#' @return \code{se} with the PRISM curves in \code{metrics_assay}. The fit is recorded in the
#'   provenance as fit source \code{"PRISM"}, and the step \code{"prism_drc_import"} holds the
#'   input files of \code{prism_curves} (name and md5) and the number of curves imported and
#'   left out. Re-running the recipe needs the same \code{prism_curves}.
#'
#' @examples
#' if (requireNamespace("gDRimport", quietly = TRUE) &&
#'     exists("read_PRISM_drc", asNamespace("gDRimport"))) {
#'   mae <- gDRutils::get_synthetic_data("finalMAE_prism")
#'   se <- MultiAssayExperiment::experiments(mae)[["single-agent"]]
#'   curves <- gDRimport::read_PRISM_drc(
#'     drc_path = system.file("testdata/prism_drc.csv", package = "gDRimport"),
#'     meta_data_path = system.file("testdata/prism_drc_model.csv", package = "gDRimport")
#'   )
#'   se <- add_PRISM_metrics(se, curves)
#'   metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#'   metrics[metrics$fit_source == "PRISM", ]
#' }
#'
#' @keywords runDrugResponseProcessingPipeline
#' @export
add_PRISM_metrics <- function(se,
                              prism_curves,
                              drug_map = NULL,
                              normalization_type = "RV",
                              metrics_assay = "Metrics",
                              averaged_assay = "Averaged") {
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_multi_class(prism_curves, c("data.table", "data.frame"))
  checkmate::assert_character(drug_map, names = "unique", any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_string(normalization_type)
  checkmate::assert_string(metrics_assay)
  checkmate::assert_string(averaged_assay)

  idfs <- gDRutils::get_SE_identifiers(se, simplify = TRUE)
  curves <- data.table::as.data.table(prism_curves)
  checkmate::assert_names(names(curves), must.include = c(
    idfs[["cellline"]], idfs[["duration"]], "pert_name", "fit_name", "lower_limit", "upper_limit",
    "slope", "inflection", "frac_var_explained", "maximum_dose"
  ))
  n_curves <- NROW(curves)
  drug <- if (is.null(drug_map)) curves$pert_name else unname(drug_map[curves$pert_name])

  rdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  rdata$row <- rownames(se)
  cell_lines <- as.character(SummarizedExperiment::colData(se)[[idfs[["cellline"]]]])
  # PRISM curves are single-agent: leave out rows with a second drug
  drug2 <- idfs[["drug2"]]
  if (drug2 %in% names(rdata)) {
    untreated <- gDRutils::get_env_identifiers("untreated_tag")
    rdata <- rdata[is.na(rdata[[drug2]]) | rdata[[drug2]] %in% untreated]
  }
  row_key <- paste(rdata[[idfs[["drug"]]]], as.numeric(rdata[[idfs[["duration"]]]]))
  curves$row <- rdata$row[match(paste(drug, as.numeric(curves[[idfs[["duration"]]]])), row_key)]
  curves$column <- colnames(se)[match(curves[[idfs[["cellline"]]]], cell_lines)]

  matched <- !is.na(curves$row) & !is.na(curves$column)
  if (!any(matched)) {
    stop("no PRISM curve matches a cell in the SE: check the drug identifiers ('drug_map') ",
         "and the cell line identifiers")
  }
  if (!all(matched)) {
    warning(sprintf(
      paste("%d of %d PRISM curves have no matching cell in the SE and are left out",
            "(%d without a cell line match, %d without a drug or Duration match)"),
      sum(!matched), n_curves, sum(is.na(curves$column)), sum(is.na(curves$row))
    ))
  }
  curves <- curves[matched]

  new_dt <- .prism_curves_to_metrics(curves, normalization_type)
  new_dt$N_conc <- .count_concentrations(se, averaged_assay, curves$row, curves$column,
                                         normalization_type, idfs[["concentration"]])
  new_dt$row <- curves$row
  new_dt$column <- curves$column

  se <- gDRutils::persist_fit_assay(se, new_dt, merge = "merge", assay_name = metrics_assay,
                                    row = "row", col = "column",
                                    upsert_key_cols = c("fit_source", "normalization_type"))

  # PRISM fits the curves itself; the recipe adds the same curves again
  se <- gDRutils::record_fit_provenance(
    se, assay = metrics_assay, fit_source = "PRISM", fn = add_PRISM_metrics,
    fit_args = list(drug_map = drug_map, normalization_type = normalization_type,
                    metrics_assay = metrics_assay, averaged_assay = averaged_assay)
  )
  files <- attr(prism_curves, "gDR_import", exact = TRUE)$files
  gDRutils::record_processing_step(
    se, "prism_drc_import",
    settings = list(
      source = "PRISM DRC table",
      files = if (!is.null(files)) as.list(files[, intersect(c("role", "name", "md5"), names(files))]),
      screen = if ("screen" %in% names(curves)) unique(curves$screen),
      n_curves = n_curves,
      n_imported = sum(matched),
      n_left_out = sum(!matched)
    )
  )
}


# One Metrics row per PRISM curve. The curve parameters map onto the gDR model; PRISM's own
# AUC and IC50 keep their names, because their definitions differ from x_AOC and xc50.
#' @keywords internal
.prism_curves_to_metrics <- function(curves, normalization_type) {
  na <- rep(NA_real_, NROW(curves))
  fit_type <- as.character(curves$fit_name)
  fit_type[is.na(fit_type) | !nzchar(fit_type)] <- NA_character_
  out <- data.table::data.table(
    x_mean = na, x_AOC = na, x_AOC_range = na, xc50 = na, x_max = na,
    ec50 = as.numeric(curves$inflection),
    x_inf = as.numeric(curves$lower_limit),
    x_0 = as.numeric(curves$upper_limit),
    h = -as.numeric(curves$slope),
    r2 = as.numeric(curves$frac_var_explained),
    p_value = na, rss = na,
    maxlog10Concentration = log10(as.numeric(curves$maximum_dose)),
    x_sd_avg = na,
    fit_type = fit_type,
    normalization_type = normalization_type,
    fit_source = "PRISM"
  )
  prism_cols <- c(auc = "prism_auc", log2_auc = "prism_log2_auc", auc_riemann = "prism_auc_riemann",
                  ic50 = "prism_ic50", log2_ic50 = "prism_log2_ic50", mse = "prism_mse",
                  mad = "prism_mad", successful_fit = "prism_successful_fit",
                  minimum_dose = "prism_minimum_dose", maximum_dose = "prism_maximum_dose")
  for (col in intersect(names(prism_cols), names(curves))) out[[prism_cols[[col]]]] <- curves[[col]]
  # the curve identity within a gDR cell
  id_cols <- c("depmap_id", "pert_id", "pool_id", "lua", "cell_set", "pert_plate", "x_project_id", "screen")
  for (col in intersect(id_cols, names(curves))) {
    out[[col]] <- as.character(curves[[col]])
  }
  out
}


# Number of concentrations with a value in the averaged data of each (row, column) cell.
#' @keywords internal
.count_concentrations <- function(se, averaged_assay, rows, columns, normalization_type, conc) {
  if (!averaged_assay %in% SummarizedExperiment::assayNames(se)) return(rep(NA_integer_, length(rows)))
  avg <- data.table::as.data.table(BumpyMatrix::unsplitAsDataFrame(
    SummarizedExperiment::assay(se, averaged_assay), row.field = "row", column.field = "column"
  ))
  avg <- avg[avg$normalization_type == normalization_type & !is.na(avg$x) & avg[[conc]] > 0]
  counts <- avg[, list(N_conc = data.table::uniqueN(get(conc))), by = c("row", "column")]
  n <- counts$N_conc[match(paste(rows, columns), paste(counts$row, counts$column))]
  as.integer(n)
}
