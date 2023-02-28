#' If a column called \code{"BackgroundValue"} exists in \code{df_},
#' it will be removed from the \code{readout} column.
#' @rdname runDrugResponseProcessingPipelineFxns
#' @export
#'
create_SE <- function(df_,
                      data_type,
                      readout = "ReadoutValue",
                      control_mean_fxn = function(x) {
                        mean(x, trim = 0.25)
                      },
                      nested_identifiers = NULL,
                      nested_confounders = intersect(names(df_),
                                                     gDRutils::get_env_identifiers("barcode")),
                      override_untrt_controls = NULL) {
  # Assertions:
  stopifnot(any(inherits(df_, "data.frame"), inherits(df_, "DataFrame")))
  checkmate::assert_string(data_type)
  checkmate::assert_string(readout)
  checkmate::assert_function(control_mean_fxn)
  checkmate::assert_character(nested_identifiers, null.ok = TRUE)
  checkmate::assert_character(nested_confounders, null.ok = TRUE)
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  

  if (is.null(nested_identifiers)) {
    nested_identifiers <-
      get_default_nested_identifiers(df_)[[data_model(data_type)]]
  }
  
  if (methods::is(df_, "data.table")) {
    df_ <- S4Vectors::DataFrame(df_)
  }

  nested_keys <- c(nested_identifiers, nested_confounders)
  identifiers <- gDRutils::get_env_identifiers()
  Keys <- identify_keys(df_, nested_keys, override_untrt_controls, identifiers)

  if (!(identifiers$masked_tag %in% colnames(df_))) {
    df_[, identifiers$masked_tag] <- FALSE
  }

  # Remove background value from readout (at least 1e-10 to avoid artefactual normalized values).
  if ("BackgroundValue" %in% colnames(df_)) {
    df_$CorrectedReadout <- pmax(df_[[readout]] - df_$BackgroundValue, 1e-10)
  } else {
    df_$CorrectedReadout <- df_[[readout]]
  }
  
  # overwrite "drug", "drug_name", "drug_moa" with "untreated" if "concentration2" == 0
  if (gDRutils::get_env_identifiers("concentration2") %in% colnames(df_)) {
    single_agent_idx <- df_[[gDRutils::get_env_identifiers("concentration2")]] == 0

    drug_cols <- c("drug", "drug_name", "drug_moa")
    drug2_var <- intersect(unlist(gDRutils::get_env_identifiers(paste0(drug_cols,
"2"), simplify = FALSE)), colnames(df_))

    df_[single_agent_idx, drug2_var] <- gDRutils::get_env_identifiers("untreated_tag")[1]
  }

  # Identify treatments, conditions, and experiment metadata.
  md <- gDRutils::split_SE_components(df_, nested_keys = Keys$nested_keys)
  coldata <- md$condition_md
  rowdata <- md$treatment_md
  exp_md <- md$experiment_md

  mapping_entries <- .create_mapping_factors(rowdata, coldata)
  mapping_entries$groupings <- rownames(mapping_entries)

  refs <- .map_references(mapping_entries)
  trt_conditions <- names(refs)
  sa_conditions <- unique(unname(unlist(refs)))

  treated <- mapping_entries[trt_conditions, ]
  untreated <- mapping_entries[!rownames(mapping_entries) %in% c(trt_conditions, sa_conditions), ]

  ## Map controls.
  controls <- list(untrt_Endpoint = "untrt_Endpoint", Day0 = "Day0")
  ctl_maps <- gDRutils::loop(controls, function(ctl_type) {
    map_df(treated, untreated, override_untrt_controls = override_untrt_controls,
           ref_cols = Keys[[ctl_type]], ref_type = ctl_type)
  })

  ## Combine all controls with respective treatments.
  # Merge raw data back with groupings.
  dfs <- merge(df_, mapping_entries, by = c(colnames(rowdata), colnames(coldata)))

  # Remove all rowdata and coldata. 
  groupings <- dfs$groupings

  ## The mapping_entries contain all exhaustive combinations of treatments and cells.
  ## Not all conditions will actually exist in the data, so filter out those that 
  ## do not exist. 
  treated <- treated[rownames(treated) %in% unique(groupings), ]
  out <- gDRutils::loop(seq_len(nrow(treated)), function(i) {
    trt <- rownames(treated)[i]
    
    trt_df <- dfs[groupings %in% trt, , drop = FALSE]
    refs_df <- dfs[groupings %in% refs[[trt]], , drop = FALSE]
    trt_df <- validate_mapping(trt_df, refs_df, nested_confounders)[, c(md$data_fields, "row_id", "col_id")]
    trt_df$row_id <- unique(dfs[groupings == trt, "row_id"]) # Override the row_id of the references.

    ctl_type <- "untrt_Endpoint"
    untrt_ref <- ctl_maps[[ctl_type]][[trt]]  
    untrt_df <- dfs[groupings %in% untrt_ref, c(md$data_fields, "row_id", "col_id"), drop = FALSE]
    untrt_df <- create_control_df(
      untrt_df, 
      control_cols = Keys[[ctl_type]], 
      control_mean_fxn, 
      out_col_name = "UntrtReadout"
    )

    ctl_type <- "Day0"
    day0_ref <- ctl_maps[[ctl_type]][[trt]]
    day0_df <- dfs[groupings %in% day0_ref, c(md$data_fields, "row_id", "col_id"), drop = FALSE]
    day0_df <- create_control_df(
      day0_df, 
      control_cols = Keys[[ctl_type]],
      control_mean_fxn, 
      out_col_name = "Day0Readout"
    )

    ## Merge all data.frames together.
    # Try to merge by plate, but otherwise just use mean. 
    ctl_df <- untrt_df 
    merge_cols <- intersect(colnames(day0_df), Keys$nested_keys)
    if (nrow(day0_df) > 0L) {
      ctl_df <- merge(untrt_df, day0_df, by = merge_cols, all = TRUE)
    } else {
      ctl_df$Day0Readout <- NA
    } 
    
    row_id <- unique(trt_df$row_id)
    col_id <- unique(trt_df$col_id)
    if (length(row_id) != 1L || length(col_id) != 1L) {
      stop(sprintf("non-unique row_ids: '%s' and col_ids: '%s'", 
        paste0(row_id, collapse = ", "), paste0(col_id, collapse = ", ")))
    }
    ctl_df$row_id <- row_id
    ctl_df$col_id <- col_id
  
    list(ctl_df = ctl_df,
         trt_df = trt_df)
  })

  trt_out <- do.call(rbind, lapply(out, "[[", "trt_df"))
  ctl_out <- do.call(rbind, lapply(out, "[[", "ctl_df"))
  trt_keep <- !colnames(trt_out) %in% c("row_id", "col_id")
  ctl_keep <- !colnames(ctl_out) %in% c("row_id", "col_id")

  trt_mat <- BumpyMatrix::splitAsBumpyMatrix(trt_out[, trt_keep, drop = FALSE],
                                                 row = trt_out$row_id, col = trt_out$col_id)
  ctl_mat <- BumpyMatrix::splitAsBumpyMatrix(ctl_out[, ctl_keep, drop = FALSE],
                                                   row = ctl_out$row_id, col = ctl_out$col_id)
  matsL <- list(RawTreated = trt_mat, Controls = ctl_mat)

  # Filter out to 'treated' conditions only.
  trt_rowdata <- rowdata[rownames(trt_mat), , drop = FALSE] 
  stopifnot(nrow(trt_rowdata) > 0)
  stopifnot(nrow(trt_rowdata) == length(unique(trt_out$row_id)))
 
  se <- SummarizedExperiment::SummarizedExperiment(assays = matsL,
    colData = coldata[match(colnames(trt_mat), rownames(coldata)), ],
    rowData = trt_rowdata,
    metadata = list())

  # Capture important values in experiment metadata.
  se <- gDRutils::set_SE_identifiers(se, identifiers)
  se <- gDRutils::set_SE_experiment_metadata(se, exp_md)
  se <- gDRutils::set_SE_keys(se, Keys)

  se
}

#' @keywords internal
validate_mapping <- function(trt_df, refs_df, nested_confounders) {
  if (!is.null(nested_confounders)) {
    refs_df <- refs_df[unique(trt_df[[nested_confounders]]) == refs_df[[nested_confounders]], ]
  }
  drug_id <- gDRutils::get_env_identifiers("drug")
  drug2_id <- gDRutils::get_env_identifiers("drug2")
  untrt_tag <- gDRutils::get_env_identifiers("untreated_tag")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc2 <- gDRutils::get_env_identifiers("concentration2")
  trt_df <- rbind(trt_df, refs_df[refs_df[[drug_id]] %in% c(unique(c(trt_df[[drug_id]],
                                                                     trt_df[[drug2_id]])), untrt_tag), ])
  # Swap concentrations for single-agent with drug2
  if (conc2 %in% colnames(trt_df)) {
    trt_df[trt_df[[drug_id]] %in% trt_df[[drug2_id]], c(conc, conc2)] <-
      trt_df[trt_df[[drug_id]] %in% trt_df[[drug2_id]], c(conc2, conc)]
  }
  trt_df
}
