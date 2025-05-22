#' get_cell_line_annotation
#'
#' Get cell line annotation data table
#'
#' @param data data.table with cell line identifiers to be matched
#' @param fname string with file name containing the annotation
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @param annotation_package string indicating name of the package containing cell line annotation
#' @return data.table with cell line annotations
#' @keywords annotation
#' @examples
#' data <- data.table::data.table(clid = c("CL1", "CL2", "CL3"))
#' cell_line_annotation <- get_cell_line_annotation(data)
#' @export
#'
get_cell_line_annotation <- function(
    data,
    fname = "cell_lines.csv",
    fill = "unknown",
    annotation_package = if ("gDRinternal" %in% .packages(all.available = TRUE)) {
      "gDRinternal"
    } else {
      "gDRtestData"
    }
) {
  checkmate::assert_data_table(data)
  
  clid <- gDRutils::get_env_identifiers("cellline")
  
  # Read the cell line annotation data
  cell_line_annotation <- if (annotation_package == "gDRinternal") {
    eval(parse(text = paste0(annotation_package, "::get_cell_line_annotations")))(data[[clid]])
  } else {
    data.table::fread(
      system.file("annotation_data", fname, package = annotation_package), header = TRUE
    )
  }
  
  if (NROW(cell_line_annotation) > 0) {
    cmn <- cell_line_annotation[[clid]] %in% data[[clid]]
    cell_line_annotation <- cell_line_annotation[cmn, ]
  }
  
  # Handle missing cell lines
  missing_cell_lines <- setdiff(unique(data[[clid]]), cell_line_annotation[[clid]])
  if (length(missing_cell_lines) > 0) {
    missing_tbl_cell_lines <- data.table::data.table(
      clid = missing_cell_lines,
      CellLineName = missing_cell_lines,
      Tissue = fill,
      ReferenceDivisionTime = NA,
      parental_identifier = fill,
      subtype = fill
    )
    data.table::setnames(missing_tbl_cell_lines, names(cell_line_annotation))
    cell_line_annotation <- rbind(cell_line_annotation, missing_tbl_cell_lines)
  }
  
  # Fill missing values
  ref_div_time <- gDRutils::get_env_identifiers("cellline_ref_div_time")
  
  for (col in setdiff(names(cell_line_annotation), ref_div_time)) {
    cell_line_annotation[[col]][is.na(cell_line_annotation[[col]])] <- fill
  }
  cell_line_annotation[, (ref_div_time) := as.numeric(get(ref_div_time))]
  
  (cell_line_annotation)
}

#' annotate_dt_with_cell_line
#'
#' Annotate cell line data with the provided annotation table
#'
#' @param data data.table with dose-response data
#' @param cell_line_annotation data.table with cell line annotations
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @return data.table with annotated cell lines
#' @keywords annotation
#' @examples
#' data <- data.table::data.table(
#'   clid = c("CL1", "CL2", "CL3"),
#'   Gnumber = c("D1", "D2", "D3")
#' )
#' cell_line_annotation <- get_cell_line_annotation(data)
#' annotated_metadata <- annotate_dt_with_cell_line(data, cell_line_annotation)
#' @export
#'
annotate_dt_with_cell_line <- function(
    data,
    cell_line_annotation,
    fill = "unknown"
) {
  checkmate::assert_data_table(data)
  checkmate::assert_data_table(cell_line_annotation)
  
  cellline <- gDRutils::get_env_identifiers("cellline")
  add_clid <- gDRutils::get_header("add_clid")
  
  checkmate::assert_names(names(cell_line_annotation), must.include = c(cellline, unlist(add_clid)))
  
  # Remove existing annotations if any
  existing_cols <- intersect(unlist(add_clid), names(data))
  if (length(existing_cols) > 0) {
    data[, (existing_cols) := NULL]
  }
  data <- cell_line_annotation[data, on = cellline]
  (data)
}

#' get_drug_annotation
#'
#' Get drug annotation data table
#'
#' @param data data.table with drug identifiers to be matched
#' @param fname string with file name containing the annotation
#' @param fill string indicating how unknown drugs should be filled in the DB
#' @param annotation_package string indicating name of the package containing drug annotation
#' @return data.table with drug annotations
#' @keywords annotation
#' @examples
#' data <- data.table::data.table(Gnumber = c("drug1", "drug2", "drug3"))
#' drug_annotation <- get_drug_annotation(data)
#' @export
#'
get_drug_annotation <- function(
    data,
    fname = "drugs.csv",
    fill = "unknown",
    annotation_package = if ("gDRinternal" %in% .packages(all.available = TRUE)) {
      "gDRinternal"
    } else {
      "gDRtestData"
    }
) {
  checkmate::assert_data_table(data)
  
  drug <- unlist(gDRutils::get_env_identifiers(c("drug", "drug2", "drug3"), simplify = FALSE))
  drug_ann_cols <- unlist(gDRutils::get_env_identifiers(c("drug", "drug_name", "drug_moa"), simplify = FALSE))
  
  
  # Read the drug annotation data
  drug_annotation <- if (annotation_package == "gDRinternal") {
    eval(parse(text = paste0(annotation_package, "::get_drug_annotations")))()
  } else {
    data.table::fread(
      system.file("annotation_data", fname, package = annotation_package), header = TRUE
    )
  }
  
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  all_data_drugs <- setdiff(gDRutils::remove_drug_batch(unlist(data[, intersect(names(data), drug), with = FALSE])),
                            untreated_tag)
  drug_annotation <- drug_annotation[gDRutils::remove_drug_batch(drug_annotation[[drug_ann_cols[["drug"]]]]) %in%
                                       all_data_drugs]
  
  # Handle missing drugs
  missing_drugs <- setdiff(unique(all_data_drugs),
                           gDRutils::remove_drug_batch(drug_annotation[[drug_ann_cols[["drug"]]]]))
  if (length(missing_drugs) > 0) {
    missing_tbl_drugs <- data.table::data.table(
      drug = missing_drugs,
      drug_name = missing_drugs,
      drug_moa = fill
    )
    data.table::setnames(missing_tbl_drugs, drug_ann_cols)
    drug_annotation <- rbind(drug_annotation, missing_tbl_drugs)
  }
  
  (drug_annotation)
}

#' annotate_dt_with_drug
#'
#' Annotate drug data with the provided annotation table
#'
#' @param data data.table with dose-response data
#' @param drug_annotation data.table with drug annotations
#' @param fill string indicating how unknown drugs should be filled in the DB
#' @return data.table with annotated drugs
#' @keywords annotation
#' @examples
#' data <- data.table::data.table(
#'   clid = c("CL1", "CL2", "CL3"),
#'   Gnumber = c("D1", "D2", "D3")
#' )
#' drug_annotation <- get_drug_annotation(data)
#' annotated_metadata <- annotate_dt_with_drug(data, drug_annotation)
#' @export
#'
annotate_dt_with_drug <- function(
    data,
    drug_annotation,
    fill = "unknown"
) {
  checkmate::assert_data_table(data)
  checkmate::assert_data_table(drug_annotation)
  
  drug <- unlist(gDRutils::get_env_identifiers(c("drug", "drug2", "drug3"), simplify = FALSE))
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  drug_name <- unlist(gDRutils::get_env_identifiers(c("drug_name", "drug_name2", "drug_name3"), simplify = FALSE))
  drug_moa <- unlist(gDRutils::get_env_identifiers(c("drug_moa", "drug_moa2", "drug_moa3"), simplify = FALSE))
  drug_idx <- which(drug %in% names(data))
  drug_data <- unlist(data[, intersect(names(data), drug), with = FALSE])
  drug_ann_cols <- unlist(gDRutils::get_env_identifiers(c("drug", "drug_name", "drug_moa"), simplify = FALSE))
  
  # Remove existing annotations if any
  existing_cols <- intersect(c(drug_name, drug_moa), names(data))
  if (length(existing_cols) > 0) {
    data[, (existing_cols) := NULL]
  }
  
  # Extract numeric suffixes from column names
  drug_ids <- sub(".*?(\\d+)$", "\\1", names(drug[drug_idx]))
  drug_ids[grepl("^\\d+$", drug_ids)] <- as.numeric(drug_ids[grepl("^\\d+$", drug_ids)])
  drug_ids[grepl("^\\D+$", drug_ids)] <- 1
  
  drug_identifiers_list <- split(c(drug[drug_idx], drug_name[drug_idx], drug_moa[drug_idx]), drug_ids)
  names(drug_identifiers_list) <- drug[drug_idx]
  drug_identifiers_list <- rev(drug_identifiers_list)
  
  drugs_treated <- setdiff(drug_data, untreated_tag)
  drugs_treated_wo_batch <- gDRutils::remove_drug_batch(drugs_treated) 
  drug_annots_wo_batch <- gDRutils::remove_drug_batch(drug_annotation[[drug[["drug"]]]])
  validated_drugs <- drugs_treated_wo_batch %in% drug_annots_wo_batch
  if (any(!validated_drugs)) {
    missing_tbl_drugs <- data.table::data.table(
      drug = gDRutils::remove_drug_batch(drugs_treated[!validated_drugs]),
      drug_name = drugs_treated[!validated_drugs],
      drug_moa = fill
    )
    data.table::setnames(missing_tbl_drugs, drug_ann_cols)
    drug_annotation <- rbind(drug_annotation, missing_tbl_drugs)
  }
  
  drug_annotation[[drug[["drug"]]]] <- gDRutils::remove_drug_batch(drug_annotation[[drug[["drug"]]]])
  untrt_drug_annotation <- data.table::data.table(
    drug = untreated_tag, drug_name = untreated_tag, drug_moa = untreated_tag)
  data.table::setnames(untrt_drug_annotation, drug_ann_cols)
  drug_annotation <- unique(rbind(
    untrt_drug_annotation,
    drug_annotation
  ))
  
  for (drug_idf in names(drug_identifiers_list)) {
    colnames(drug_annotation) <- drug_identifiers_list[[drug_idf]]
    data$batch <- data[[drug_idf]]
    data[[drug_idf]] <- gDRutils::remove_drug_batch(data[[drug_idf]])
    req_col <- c(drug_idf, setdiff(names(data), names(drug_annotation)))
    data.table::setkeyv(drug_annotation, drug_idf)
    data <- drug_annotation[data[, req_col, with = FALSE], on = drug_idf]
    data[[drug_idf]] <- data$batch
    data$batch <- NULL
  }
  
  (data)
}


#' Retrieve the drug annotation from the annotated dt input
#'
#' @param dt annotated data.table
#'
#' @return data.table with drug annotation
#' @export
#'
#' @keywords annotation
#' @examples
#' dt <- data.table::data.table(Gnumber = "A",
#' DrugName = "drugA",
#' drug_moa = "drug_moa_A")
#' get_drug_annotation_from_dt(dt)
get_drug_annotation_from_dt <- function(dt) {
  checkmate::assert_data_table(dt)
  
  identifiers <- gDRutils::get_env_identifiers()
  drug_ids <- identifiers[grep("drug", names(identifiers))]
  drug_cols <- Filter(function(el) el %in% names(dt), drug_ids)
  
  dt_drug <- dt[, unlist(drug_cols), with = FALSE]
  dt_long <- data.table::melt(dt_drug,
                              measure.vars = patterns(gsub("[0-9]",
                                                           ".",
                                                           paste0("^",
                                                                  unlist(drug_cols[c("drug",
                                                                                     "drug_name",
                                                                                     "drug_moa")])))),
                              value.name = unlist(drug_cols[c("drug",
                                                              "drug_name",
                                                              "drug_moa")]))
  dt_long[, "variable" := NULL]
  unique_dt <- unique(dt_long)
  
  missing_cols <- setdiff(unlist(drug_ids[c("drug", "drug_name", "drug_moa")]), names(unique_dt))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      unique_dt[, (col) := "unknown"]
    }
  }
  
  unique_dt[!unique_dt[[drug_cols[["drug"]]]] %in% gDRutils::get_env_identifiers("untreated_tag"), ]
}


#' Retrieve the cell line annotation from the annotated dt input
#'
#' @param dt annotated data.table
#'
#' @return data.table with cell line annotation
#' @export
#'
#' @keywords annotation
#' @examples
#' dt <- data.table::data.table(Gnumber = "A",
#' clid = "CL123",
#' CellLineName = "cl name",
#' Tissue = "Bone",
#' parental_identifier = "some cl",
#' subtype = "cortical",
#' ReferenceDivisionTime = 5)
#' get_cellline_annotation_from_dt(dt)
get_cellline_annotation_from_dt <- function(dt) {
  checkmate::assert_data_table(dt)
  
  cell_cols <- c(gDRutils::get_env_identifiers("cellline"),
                 gDRutils::get_header("add_clid"))
  
  cell_dt <- dt[, intersect(unlist(cell_cols), names(dt)), with = FALSE]
  
  missing_cols <- setdiff(unlist(cell_cols), names(cell_dt))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      cell_dt[, (col) := "unknown"]
    }
  }
  
  unique(cell_dt)
}


#' annotate_se_with_drug
#'
#' Annotate SummarizedExperiment object with drug annotations
#'
#' @param se SummarizedExperiment object containing dose-response data
#' @param drug_annotation data.table with drug annotations
#' @param fill string indicating how unknown drugs should be filled in the DB
#' @return SummarizedExperiment object with annotated drugs
#' @keywords annotation
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   rowData = data.table::data.table(Gnumber = c("D1", "D2", "D3"))
#' )
#' drug_annotation <- get_drug_annotation(data.table::as.data.table(SummarizedExperiment::rowData(se)))
#' annotated_se <- annotate_se_with_drug(se, drug_annotation)
#' @export
#'
annotate_se_with_drug <- function(
    se,
    drug_annotation,
    fill = "unknown"
) {
  checkmate::assert_class(se, "SummarizedExperiment")
  data <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  annotated_data <- annotate_dt_with_drug(data, drug_annotation, fill)
  SummarizedExperiment::rowData(se) <- annotated_data
  (se)
}

#' annotate_mae_with_drug
#'
#' Annotate MultiAssayExperiment object with drug annotations
#'
#' @param mae MultiAssayExperiment object containing dose-response data
#' @param drug_annotation data.table with drug annotations
#' @param fill string indicating how unknown drugs should be filled in the DB
#' @return MultiAssayExperiment object with annotated drugs
#' @keywords annotation
#' @examples
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'   experiments = list(SummarizedExperiment::SummarizedExperiment(
#'     rowData = data.table::data.table(Gnumber = c("D1", "D2", "D3"))
#'   ))
#' )
#' drug_annotation <- get_drug_annotation(data.table::as.data.table(
#'   SummarizedExperiment::rowData(
#'     MultiAssayExperiment::experiments(mae)[[1]])))
#' annotated_mae <- annotate_mae_with_drug(mae, drug_annotation)
#' @export
#'
annotate_mae_with_drug <- function(
    mae,
    drug_annotation,
    fill = "unknown"
) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  for (i in seq_along(MultiAssayExperiment::experiments(mae))) {
    se <- MultiAssayExperiment::experiments(mae)[[i]]
    MultiAssayExperiment::experiments(mae)[[i]] <- annotate_se_with_drug(se, drug_annotation, fill)
  }
  (mae)
}

#' annotate_se_with_cell_line
#'
#' Annotate SummarizedExperiment object with cell line annotations
#'
#' @param se SummarizedExperiment object containing dose-response data
#' @param cell_line_annotation data.table with cell line annotations
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @return SummarizedExperiment object with annotated cell lines
#' @keywords annotation
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   rowData = data.table::data.table(clid = c("CL1", "CL2", "CL3"))
#' )
#' cell_line_annotation <- get_cell_line_annotation(data.table::as.data.table(SummarizedExperiment::rowData(se)))
#' annotated_se <- annotate_se_with_cell_line(se, cell_line_annotation)
#' @export
#'
annotate_se_with_cell_line <- function(
    se,
    cell_line_annotation,
    fill = "unknown"
) {
  checkmate::assert_class(se, "SummarizedExperiment")
  data <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  annotated_data <- annotate_dt_with_cell_line(data, cell_line_annotation, fill)
  SummarizedExperiment::rowData(se) <- annotated_data
  (se)
}

#' annotate_mae_with_cell_line
#'
#' Annotate MultiAssayExperiment object with cell line annotations
#'
#' @param mae MultiAssayExperiment object containing dose-response data
#' @param cell_line_annotation data.table with cell line annotations
#' @param fill string indicating how unknown cell lines should be filled in the DB
#' @return MultiAssayExperiment object with annotated cell lines
#' @keywords annotation
#' @examples
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'   experiments = list(SummarizedExperiment::SummarizedExperiment(
#'     rowData = data.table::data.table(clid = c("CL1", "CL2", "CL3"))
#'   ))
#' )
#' cell_line_annotation <- get_cell_line_annotation(data.table::as.data.table(
#'   SummarizedExperiment::rowData(
#'     MultiAssayExperiment::experiments(mae)[[1]])))
#' annotated_mae <- annotate_mae_with_cell_line(mae, cell_line_annotation)
#' @export
#'
annotate_mae_with_cell_line <- function(
    mae,
    cell_line_annotation,
    fill = "unknown"
) {
  checkmate::assert_class(mae, "MultiAssayExperiment")
  for (i in seq_along(MultiAssayExperiment::experiments(mae))) {
    se <- MultiAssayExperiment::experiments(mae)[[i]]
    MultiAssayExperiment::experiments(mae)[[i]] <- annotate_se_with_cell_line(se,
                                                                              cell_line_annotation,
                                                                              fill)
  }
  (mae)
}
