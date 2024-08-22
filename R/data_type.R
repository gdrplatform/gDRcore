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
#' @keywords data_type
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
#' @keywords data_type
#' @export
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
identify_data_type <- function(df,
                               codilution_conc = 2,
                               matrix_conc = 1) {
  
  # Get drug and concentration identifiers
  drug_ids <- get_relevant_ids(c("drug_name", "drug_name2", "drug_name3"), df)
  conc_ids <- get_relevant_ids(c("concentration", "concentration2", "concentration3"), df)
  
  # Filter out primary drug and concentration identifiers
  drugs_cotrt_ids <- setdiff(drug_ids, "drug_name")
  conc_cotrt_ids <- setdiff(conc_ids, "concentration")
  
  # Get untreated tag and cell line identifiers
  untreated_tag <- c(gDRutils::get_env_identifiers("untreated_tag"), NA)
  cell <- gDRutils::get_env_identifiers("cellline_name")
  
  # Process perturbations
  dt <- process_perturbations(df, drugs_cotrt_ids, conc_cotrt_ids, untreated_tag)
  
  # Double-check the columns after processing perturbations
  drug_ids <- get_relevant_ids(c("drug_name", "drug_name2", "drug_name3"), df)
  conc_ids <- get_relevant_ids(c("concentration", "concentration2", "concentration3"), df)
  
  drugs_cotrt_ids <- setdiff(drug_ids, "drug_name")
  conc_cotrt_ids <- setdiff(conc_ids, "concentration")
  
  
  cols_pairs <- intersect(names(df), c(drug_ids, cell, conc_ids))
  drug_pairs <- unique(df[, cols_pairs, with = FALSE])

  df[, record_id := .I]
  df[, type := NA_character_]
  
  sa_name <- gDRutils::get_supported_experiments("sa")

  controls <- rowSums(df[, conc_ids, with = FALSE] == 0) == length(conc_ids)
  single_agent <- rowSums(df[, conc_ids, with = FALSE] != 0) == 1
  df$type <- ifelse(controls, "control",
                    ifelse(single_agent, sa_name, NA))
  

  if (length(conc_ids) > 1) {
    missing_type_rows <- is.na(df$type)
    flat_data <- df[missing_type_rows, cols_pairs, with = FALSE]
    conc_ratio <- table(
      round(log10(flat_data[[conc_ids[["concentration"]]]] /
                    flat_data[[conc_ids[["concentration2"]]]]), 2)
    )
    conc_ratio <- conc_ratio[!names(conc_ratio) %in% c("Inf", "-Inf")]
    
    type <- ifelse(length(conc_ratio) <= codilution_conc,
                   gDRutils::get_supported_experiments("cd"),
                   gDRutils::get_supported_experiments("combo"))
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
#' @keywords data_type
#' @export
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
split_raw_data <- function(df,
                           type_col = "type") {
  
  df <- collapse_drugs(df)

  drug_ids <- unlist(gDRutils::get_env_identifiers(
    c("drug_name", "drug_name2", "drug_name3", "drug", "drug2", "drug3",
      "drug_moa", "drug_moa2", "drug_moa3", "concentration", "concentration2",
      "concentration3"), 
    simplify = FALSE)
  )
  
  sa_name <- gDRutils::get_supported_experiments("sa")
  
  drug_ids <- drug_ids[which(drug_ids %in% names(df))]
  codrug_ids <- drug_ids[grep("[0-9]", names(drug_ids))]
  conc_idx <- drug_ids[grep("concentration", names(drug_ids))]
  codrug_drug_id <- setdiff(codrug_ids, conc_idx)
  cl <- gDRutils::get_env_identifiers("cellline")

  df_list <- split(df, df[[type_col]])
  types <- setdiff(names(df_list), "control")
  control <- df_list[["control"]]
  df_list[["control"]] <- NULL
  cotrt_types <- setdiff(names(df_list), sa_name)
  
  control_sa_idx <- which(
    rowSums(df[, conc_idx, with = FALSE] == 0) == length(conc_idx)
  )
  control_sa <- df[control_sa_idx, ]
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  
  if (length(cotrt_types) > 0) {
    df_list[cotrt_types] <- gDRutils::loop(cotrt_types, function(x) {
      colnames <- c(cl, drug_ids[["drug_name"]])
      
      df_list[[x]] <- unify_combination_data(df_list[[x]], cl, drug_ids)
      unique_cotrt <- unique(df_list[[x]][, colnames, with = FALSE])
      unique_cotrt_ctrl <- unique(
        control[
          control[[cl]] %in% unique_cotrt[[cl]] &
            control[[drug_ids[["drug_name"]]]] %in% untreated_tag, 
        ][, colnames, with = FALSE])
      
      
      
      cotrt_matching <- unique(rbind(unique_cotrt, unique_cotrt_ctrl))
      df_merged <- rbind(
        df_list[[x]], 
        cotrt_matching[control, on = intersect(names(cotrt_matching), names(control))])
      if (x == gDRutils::get_supported_experiments("combo")) {
        matrix_data <- rbind(df_merged, df_list[[sa_name]])
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
  
  if (any(sa_name == names(df_list))) {
    sa_idx <- gDRutils::loop(
      grep(drug_ids[["concentration"]], drug_ids, value = TRUE), 
      function(x) which(!df_list[[sa_name]][, x, with = FALSE] == 0)
    )
    sa_idx[["concentration"]] <- NULL
    
    for (codrug in names(sa_idx)) {
      codrug_cols <- grep(
        as.numeric(gsub("\\D", "", codrug)), 
        drug_ids, 
        value = TRUE
      )
      selected_columns <- unname(drug_ids[c("drug_name", "drug", "drug_moa", "concentration")])
      df_list[[sa_name]][sa_idx[[codrug]], selected_columns] <- 
        df_list[[sa_name]][sa_idx[[codrug]], codrug_cols, with = FALSE]
    }
    df_list[[sa_name]][, codrug_ids] <- NULL
    
    selected_columns <- names(df_list[[sa_name]])
    
    df_list[[sa_name]] <- rbind(
      df_list[[sa_name]],
      control_sa[, selected_columns, with = FALSE]
    )
  }
  
  Map(function(x) {
    selected_columns <- which(names(x) != type_col)
    x[, selected_columns, with = FALSE]
  }, df_list)
}

#' @keywords internal
#' Function to get relevant identifiers from the environment
#' @param identifiers A character vector of identifier names to fetch from the environment
#' @param df A data frame containing the columns to be checked against the identifiers
#' @return A character vector of relevant identifiers that are present in the data frame
#' @keywords internal
get_relevant_ids <- function(identifiers, df) {
  ids <- unlist(gDRutils::get_env_identifiers(identifiers, simplify = FALSE))
  ids[ids %in% names(df)]
}


#' @keywords internal
#'
#' ensure that the pair of drugs are matching such that there are not data with
#' `drug = drugA` and `drug_2 = drugB`, and `drug = drugB` and `drug_2 = drugA` at the same time
#'
unify_combination_data <- function(df, cl, drug_ids) {
  cotrt_data <- split(df, df[[cl]])
  data.table::rbindlist(lapply(cotrt_data, function(x) {
    duplicated_full_idx <- which(x[[drug_ids[["drug_name"]]]] %chin% x[[drug_ids[["drug_name2"]]]] &
                                   x[[drug_ids[["drug_name2"]]]] %chin% x[[drug_ids[["drug_name"]]]])
    if (length(duplicated_full_idx)) {
      duplicated_data <- x[duplicated_full_idx, ]
      drug_data <- duplicated_data[, drug_ids[c("drug_name", "drug_name2")], with = FALSE]
      unique_rows <- drug_data[!duplicated(t(apply(drug_data, 1, sort))), ]
      idx_duplicated <- duplicated_full_idx[which(is.na(grr_matches(do.call(paste, unique_rows),
                                                                do.call(paste, drug_data))$x))]
      
      primary_drug_idfs <- drug_ids[c("drug", "drug_name", "drug_moa", "concentration")]
      secondary_drug_idfs <- drug_ids[c("drug2", "drug_name2", "drug_moa2", "concentration2")]
      
      x[idx_duplicated,
        c(primary_drug_idfs, secondary_drug_idfs) := .(get(drug_ids[["drug2"]]),
                                                       get(drug_ids[["drug_name2"]]),
                                                       get(drug_ids[["drug_moa2"]]),
                                                       get(drug_ids[["concentration2"]]),
                                                       get(drug_ids[["drug"]]),
                                                       get(drug_ids[["drug_name"]]),
                                                       get(drug_ids[["drug_moa"]]),
                                                       get(drug_ids[["concentration"]])
        )] 
    }
    x
  }))
}


#' @keywords internal
#' avoid `vehicle` in drug columns with low ordinality. 
#' For example, if `drug = vehicle` and `drug_2 = drugA`, change to `drug = drugA` and `drug_2 = vehicle`.
#' propagate to the other columns accordingly (concentration, drug_moa, drug_name)
#'
collapse_drugs <- function(df) {
  
  drug_ids <- unlist(gDRutils::get_env_identifiers(
    c("drug_name", "drug_name2", "drug_name3", "drug", "drug2", "drug3",
      "drug_moa", "drug_moa2", "drug_moa3", "concentration", "concentration2",
      "concentration3"), 
    simplify = FALSE)
  )
  # add '1' to the names for more efficient loop lower
  names(drug_ids)[!grepl("[0-9]", names(drug_ids))] <- paste0(names(drug_ids)[!grepl("[0-9]", names(drug_ids))], "1")
  drug_ids <- drug_ids[which(drug_ids %in% names(df))]
  max_drug <- max(as.numeric(unlist(lapply(names(drug_ids), function(x) substr(x, nchar(x), nchar(x)))),
                             na.rm = TRUE))
  
  if (max_drug > 1) { # collapse treatment iteratively
    for (i in 2:max_drug) {
      for (j in seq_len(max_drug - 1)) {
        idx <- df[[drug_ids[paste0("drug_name", j)]]] %in% gDRutils::get_env_identifiers("untreated_tag") & 
          !(df[[drug_ids[paste0("drug_name", j + 1)]]] %in% gDRutils::get_env_identifiers("untreated_tag"))

        col_idx1 <- drug_ids[grepl(j, names(drug_ids))]
        col_idx2 <- drug_ids[gsub(as.character(j), as.character(j + 1), names(col_idx1))]
        
        df_replace <- df[idx, c(col_idx2, col_idx1), with = FALSE]
        df[idx,
           c(col_idx1, col_idx2) := df_replace]
      }
    }
  }
  df  
}




#' Cleanup additional perturbations in the data.table
#' 
#' This function processes drug and concentration columns in a data.table.
#' It checks if there is only one unique drug (excluding a specified untreated tag)
#' and if there are exactly two doses (one of which is 0). If these conditions are met,
#' it creates a new column named after the drug and fills it with the doses,
#' then removes the original drug and concentration columns.
#'
#' @param dt A data.table containing the data.
#' @param drugs_cotrt_ids A vector of column names related to drugs.
#' @param conc_cotrt_ids A vector of column names related to concentrations.
#' @param untreated_tag A string representing the untreated tag (default is "vehicle").
#' @return A modified data.table with new columns for the drugs and removed original drug and concentration columns.
#' @examples
#' library(data.table)
#' dt <- data.table(
#'   drug1 = c("vehicle", "drugA", "drugA"),
#'   conc1 = c(0, 10, 0),
#'   drug2 = c("vehicle", "drugB", "drugB"),
#'   conc2 = c(0, 20, 0)
#' )
#' drugs_cotrt_ids <- c("drug1", "drug2")
#' conc_cotrt_ids <- c("conc1", "conc2")
#' dt <- process_drug_doses(dt, drugs_cotrt_ids, conc_cotrt_ids)
#' print(dt)
#' @export
process_perturbations <- function(dt,
                                  drugs_cotrt_ids,
                                  conc_cotrt_ids,
                                  untreated_tag = "vehicle") {
  
  # Assertions
  checkmate::assert_data_table(dt)
  checkmate::assert_character(drugs_cotrt_ids, any.missing = FALSE)
  checkmate::assert_character(conc_cotrt_ids, any.missing = FALSE)
  checkmate::assert_true(length(drugs_cotrt_ids) == length(conc_cotrt_ids))
  
  
  # If lengths of drugs_cotrt_ids and conc_cotrt_ids are 0, return dt unchanged
  if (length(drugs_cotrt_ids) == 0 && length(conc_cotrt_ids) == 0) {
    return(dt)
  } else {
    
    # Iterate through each pair of columns in drugs_cotrt_ids and conc_cotrt_ids
    for (i in seq_along(drugs_cotrt_ids)) {
      drug_col <- drugs_cotrt_ids[i]
      conc_col <- conc_cotrt_ids[i]
      
      # Check if there is only one drug in the current drug column (excluding untreated_tag)
      unique_drugs <- unique(dt[[drug_col]])
      unique_drugs <- unique_drugs[!unique_drugs %in% untreated_tag]
      
      if (length(unique_drugs) == 1) {
        # Check if there are only two doses in the current concentration column (0 and another value)
        unique_doses <- unique(dt[[conc_col]])
        
        if (length(unique_doses) == 2 && 0 %in% unique_doses) {
          # Create a new column named after the drug (excluding untreated_tag) and fill it with the doses
          new_column_name <- unique_drugs[1]
          dt[, (new_column_name) := dt[[conc_col]]]
          
          drug_order <- gsub(".*_(\\d+)$|.*", "\\1", drug_col)
          drug_cols <- 
            unique(c(
              get_relevant_ids(paste0(
              c("drug", "drug_name", "drug_moa", "concentration"), drug_order),
              dt),
              drug_col, conc_col))
          
          # Remove the current drug and concentration columns
          dt[, (drug_cols) := NULL]
        }
      }
    }
    return(dt) 
  }
}
