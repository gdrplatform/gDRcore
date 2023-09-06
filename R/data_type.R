#' Identify type of data
#'
#' @param df data.table of raw drug response data 
#'           containing both treated and untreated values
#' @param codilution_conc integer of maximum number of concentration ratio 
#'                        of co-treatment to classify as codilution data type;
#'                        defaults to \code{2}
#' @param matrix_conc integer of minimum number of concentration pairs 
#'                    of co-treatment to classify 
#'                    as co-treatment or matrix data type;
#'                    defaults to \code{1}
#' 
#' @examples
#' conc <- rep(seq(0, 0.3, 0.1), 2)
#' ctrl_df <- S4Vectors::DataFrame(
#'   ReadoutValue = c(2, 2, 1, 1, 2, 1),
#'   Concentration = rep(0, 6),
#'   masked = FALSE,
#'   DrugName = rep(c("DRUG_10", "vehicle", "DRUG_8"), 2),
#'   CellLineName = "CELL1"
#' )
#' 
#' trt_df <- S4Vectors::DataFrame(
#'   ReadoutValue = rep(seq(1, 4, 1), 2),
#'   Concentration = conc,
#'   masked = rep(FALSE, 8),
#'   DrugName = c("DRUG_10", "DRUG_8"),
#'   CellLineName = "CELL1"
#' )
#' input_df <- data.table::as.data.table(rbind(ctrl_df, trt_df))
#' input_df$Duration <- 72
#' input_df$CorrectedReadout2 <- input_df$ReadoutValue
#' identify_data_type(input_df)

#' @return 
#' data.table of raw drug response data with additional column \code{type} 
#' with the info of data type for a given row of data.table
#' @export
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
identify_data_type <- function(df,
                               codilution_conc = 2,
                               matrix_conc = 1) {
  
  # find the pairs of drugs with relevant metadata
  drug_ids <- unlist(gDRutils::get_env_identifiers(
    c("drug_name", "drug_name2"), 
    simplify = FALSE
  )) 
  drugs_ids <- drug_ids[which(drug_ids %in% names(df))]
  
  conc_ids <- unlist(gDRutils::get_env_identifiers(
    c("concentration", "concentration2"), 
    simplify = FALSE
  ))
  conc_ids <- conc_ids[which(conc_ids %in% names(df))]
  
  drugs_cotrt_ids <- drugs_ids[!names(drugs_ids) %in% "drug_name"]
  conc_cotrt_ids <- conc_ids[!names(conc_ids) %in% "concentration"]
  
  untreated_tag <- c(gDRutils::get_env_identifiers("untreated_tag"), NA)
  cell <- gDRutils::get_env_identifiers("cellline_name")
  cols_pairs <- intersect(names(df),  c(drug_ids, cell, conc_ids))
  drug_pairs <- unique(df[, cols_pairs, with = FALSE])

  df[, record_id := .I]
  df[, type := NA_character_]
  # loop through the pairs to assess the number of individual 
  # concentration pairs
  # Vectorized operations for calculating 'type'
  zero_concentration <- rowSums(df[, conc_ids, with = FALSE] == 0) == length(conc_ids)
  one_nonzero_concentration <- rowSums(df[, conc_ids, with = FALSE] != 0) == 1
  df$type <- ifelse(zero_concentration, "control",
                    ifelse(one_nonzero_concentration, "single-agent", NA))
  

  # Vectorized operations for calculating 'type' for missing rows
  if (length(conc_ids) > 1) {
    # Identify rows with missing 'type'
    missing_type_rows <- is.na(df$type)
    flat_data <- df[missing_type_rows, cols_pairs, with = FALSE]
    non_zero_concentration2 <- flat_data[[conc_ids[["concentration2"]]]] > 0
    conc_1 <- table(flat_data[[conc_ids[["concentration"]]]])
    conc_2 <- table(flat_data[[conc_ids[["concentration2"]]]])
    n_conc_pairs <- nrow(unique(flat_data[, conc_ids, with = FALSE]))
    conc_ratio <- table(
      round(log10(flat_data[[conc_ids[["concentration"]]]] /
                    flat_data[[conc_ids[["concentration2"]]]]), 2)
    )
    conc_ratio <- conc_ratio[!names(conc_ratio) %in% c("Inf", "-Inf")]
    
    type <- ifelse(length(conc_ratio) <= codilution_conc, "co-dilution", "matrix")
    df$type[missing_type_rows] <- type
  }
  df
}

#' Split raw data into list based on the data types
#'
#' @param df data.table of raw drug response data containing both treated and 
#' untreated values with column specified in `type_col` argument.
#' @param type_col string with column names in `df` with info about data type.
#' Defaults to \code{"type"}.
#' 
#' @examples 
#' cell_lines <- gDRtestData::create_synthetic_cell_lines()
#' drugs <- gDRtestData::create_synthetic_drugs()
#' df_layout <- drugs[4:6, as.list(cell_lines[7:8, ]), names(drugs)]
#' df_layout <- gDRtestData::add_data_replicates(df_layout)
#' df_layout <- gDRtestData::add_concentration(
#'   df_layout,
#'   concentrations = 10 ^ (seq(-3, .5, .5))
#' )
#' 
#' df_2 <- 
#'   drugs[c(21, 26), as.list(cell_lines[which(cell_lines$clid %in% df_layout$clid)]), names(drugs)]
#' df_2 <- gDRtestData::add_data_replicates(df_2)
#' df_2 <- gDRtestData::add_concentration(
#'   df_2,
#'   concentrations = 10 ^ (seq(-3, .5, .5))
#' )
#' colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")] <-
#'   paste0(
#'     colnames(df_2)[colnames(df_2) %in% c(colnames(drugs), "Concentration")],
#'     "_2"
#'   )
#' df_layout_2 <- df_layout[df_2, on = intersect(names(df_layout), names(df_2)), 
#'                         allow.cartesian = TRUE]
#' df_merged_data <- gDRtestData::generate_response_data(df_layout_2, 0)
#' df <- identify_data_type(df_merged_data)
#' split_raw_data(df)
#' 
#' @examples
#' conc <- rep(seq(0, 0.3, 0.1), 2)
#' ctrl_df <- S4Vectors::DataFrame(
#'   ReadoutValue = c(2, 2, 1, 1, 2, 1),
#'   Concentration = rep(0, 6),
#'   masked = FALSE,
#'   DrugName = rep(c("DRUG_10", "vehicle", "DRUG_8"), 2),
#'   CellLineName = "CELL1"
#' )
#' 
#' trt_df <- S4Vectors::DataFrame(
#'   ReadoutValue = rep(seq(1, 4, 1), 2),
#'   Concentration = conc,
#'   masked = rep(FALSE, 8),
#'   DrugName = c("DRUG_10", "DRUG_8"),
#'   CellLineName = "CELL1"
#' )
#' input_df <- data.table::as.data.table(rbind(ctrl_df, trt_df))
#' input_df$Duration <- 72
#' input_df$CorrectedReadout2 <- input_df$ReadoutValue
#' split_df <- identify_data_type(input_df)
#' split_raw_data(split_df)
#'
#' @return list with split data based on its data type
#' @export
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
split_raw_data <- function(df,
                           type_col = "type") {
  
  drug_ids <- unlist(gDRutils::get_env_identifiers(
    c("drug_name", "drug_name2", "drug", "drug2",
      "drug_moa", "drug_moa2", "concentration", "concentration2"), 
    simplify = FALSE)
  ) 
  drug_ids <- drug_ids[which(drug_ids %in% names(df))]
  codrug_ids <- drug_ids[grep("[0-9]", names(drug_ids))]
  conc_idx <- drug_ids[grep("concentration", names(drug_ids))]
  codrug_drug_id <- setdiff(codrug_ids, conc_idx)
  
  cl <- gDRutils::get_env_identifiers("cellline")
  df_list <- split(df, df[[type_col]])
  types <- setdiff(names(df_list), "control")
  control <- df_list[["control"]]
  df_list[["control"]] <- NULL
  cotrt_types <- setdiff(names(df_list), "single-agent")
  
  control_sa_idx <- which(
    rowSums(df[, conc_idx, with = FALSE] == 0) == length(conc_idx)
  )
  control_sa <- df[control_sa_idx, ]
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  
  if (length(cotrt_types) > 0) {
    df_list[cotrt_types] <- gDRutils::loop(cotrt_types, function(x) {
      colnames <- c(cl, drug_ids[["drug_name"]])
      unique_cotrt <- unique(df_list[[x]][, colnames, with = FALSE])
      unique_cotrt_ctrl <- unique(
        control[
          control[[cl]] %in% unique_cotrt[[cl]] &
            control[[drug_ids[["drug_name"]]]] %in% untreated_tag, 
        ][, colnames, with = FALSE])
      cotrt_matching <- rbind(unique_cotrt, unique_cotrt_ctrl)
    df_merged <- rbind(
        df_list[[x]], 
        cotrt_matching[control, on = intersect(names(cotrt_matching), names(control))])
      if (x == "matrix") {
        matrix_data <- rbind(df_merged, df_list[["single-agent"]])
        for (j in conc_idx)
          data.table::set(matrix_data, which(is.na(matrix_data[[j]])), j, 0)
        for (j in codrug_drug_id) {
          data.table::set(matrix_data, which(is.na(matrix_data[[j]])), j, untreated_tag[1])
        }
        matrix_data
      } else {
        df_merged
      }
    })
  }
  
  if ("single-agent" %in% names(df_list)) {
    sa_idx <- gDRutils::loop(
      grep(drug_ids[["concentration"]], drug_ids, value = TRUE), 
      function(x) which(!df_list[["single-agent"]][, x, with = FALSE] == 0)
    )
    sa_idx[["concentration"]] <- NULL
    
    for (codrug in names(sa_idx)) {
      codrug_cols <- grep(
        as.numeric(gsub("\\D", "", codrug)), 
        drug_ids, 
        value = TRUE
      )
      selected_columns <- unname(drug_ids[c("drug_name", "drug", "drug_moa", "concentration")])
      df_list[["single-agent"]][sa_idx[[codrug]], selected_columns] <- 
        df_list[["single-agent"]][sa_idx[[codrug]], codrug_cols, with = FALSE]
    }
    df_list[["single-agent"]][, codrug_ids] <- NULL
    
    selected_columns <- names(df_list[["single-agent"]])
    
    df_list[["single-agent"]] <- rbind(
      df_list[["single-agent"]],
      control_sa[, selected_columns, with = FALSE]
    )
  }
  
  Map(function(x) {
    selected_columns <- which(names(x) != type_col)
    x[, selected_columns, with = FALSE]
  }, df_list)
}
