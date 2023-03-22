#' Identify type of data
#'
#' @param df data.frame of raw drug response data containing both treated and untreated values
#' @param codilution_conc integer of maximum number of concentration ratio of co-treatment
#' to classify as codilution data type.
#' Defaults to \code{2}.
#' @param matrix_conc integer of minimum number of concentration pairs of co-treatment
#' to classify as co-treatment or matrix data type.
#' Defaults to \code{1}.
#'
#' @return data.frame of raw drug response data with additional column `type` with the info of
#' data type for a given row of data.frame
#' @details
#' An integer is specified to indicate the maximum number of concentrations to classify data as
#' a cotreatment. Any higher concentration series will be considered as combination matrix data.
#' @export
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
identify_data_type <- function(df,
                               codilution_conc = 2,
                               matrix_conc = 1      # forces any co-treatment as a matrix data model
                               ) {
  
  p_cols <- names(df)
  drug_ids <- unlist(gDRutils::get_env_identifiers(c("drug_name", "drug_name2"), simplify = FALSE)) 
  p_drug_ids <- drug_ids[drug_ids %in% p_cols]
  conc_ids <- unlist(gDRutils::get_env_identifiers(c("concentration",
                                                     "concentration2"), simplify = FALSE))
  p_conc_ids <- conc_ids[conc_ids %in% p_cols]
  
  drugs_cotrt_ids <- setdiff(p_drug_ids, unlist(gDRutils::get_env_identifiers("drug_name")))
  conc_cotrt_ids <- setdiff(p_conc_ids, unlist(gDRutils::get_env_identifiers("concentration")))
  
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  cell <- gDRutils::get_env_identifiers("cellline_name") # TODO: Why do we use cellline_name instead of clid?
  cols_pairs <- intersect(p_cols, c(drug_ids, cell)) # TODO: Can we use p_drug_ids here?
  combinations <- unique(df[, cols_pairs])

  cnt <- seq_len(nrow(df))
  df$type <- NA
  # loop through the pairs to assess the number of individual concentration pairs
  for (idp in seq_len(nrow(combinations))) {
    df_matching <- merge(cbind(df, cnt), combinations[idp, ])
    matching_idx <- df_matching$cnt
    treated <- vapply(gDRutils::loop(df_matching[, p_drug_ids, drop = FALSE],
                             function(x) !x %in% untreated_tag), all, logical(1))
    detect_sa <- sum(treated)
    type <- if (ncol(df[matching_idx, drugs_cotrt_ids, drop = FALSE]) == 0) {
      if (all(df[matching_idx, drug_ids[["drug_name"]]] %in% untreated_tag)) {
        "control"
      } else {
        "single-agent"
      }
    } else if (detect_sa == 1L) {
      "single-agent"
    } else if (detect_sa == 0L) {
      "control" 
    } else {
      NA
    }
    df[matching_idx, "type"]  <- type
    
    if (length(conc_ids) > 1) {
      df[matching_idx, "type"] <- ifelse(rowSums(df[matching_idx, conc_ids, drop = FALSE] == 0) == 1,
                                         "single-agent", df[matching_idx, "type"])
    }
    
    if (all(!is.na(df[matching_idx, "type"]))) {
      next
    }
    flat_data <- df_matching[df_matching[[conc_ids[["concentration2"]]]] > 0, ]
    conc_1 <- table(flat_data[[conc_ids[["concentration"]]]])
    conc_2 <- table(flat_data[[conc_ids[["concentration2"]]]])
    n_conc_pairs <- nrow(unique(flat_data[, conc_ids]))
    conc_ratio <- table(round(log10(flat_data[[conc_ids[["concentration"]]]] /
                                      flat_data[[conc_ids[["concentration2"]]]]), 2))
    conc_ratio <- conc_ratio[!names(conc_ratio) %in% c("Inf", "-Inf")]
    
    type <- 
      if (length(conc_ratio) <= codilution_conc) {
        "co-dilution"
      } else {
        "matrix"
      }
    df[matching_idx, "type"]  <- ifelse(is.na(df[matching_idx, "type"]), type, df[matching_idx, "type"])
  }
  df
}

#' Split raw data into list based on the data types
#'
#' @param df data.frame of raw drug response data containing both treated and untreated values with
#' column specified in `type_col` argument.
#' @param type_col string with column names in `df` with info about data type.
#' Defaults to \code{"type"}.
#'
#' @return list with split data based on its data type
#' @export
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
split_raw_data <- function(df,
                           type_col = "type") {
  
  drug_ids <- unlist(gDRutils::get_env_identifiers(c("drug_name", "drug_name2",
                                                     "drug", "drug2",
                                                     "drug_moa", "drug_moa2",
                                                     "concentration", "concentration2"), simplify = FALSE)) 
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
  control_sa_idx <- which(rowSums(df[, conc_idx, drop = FALSE] == 0) == length(conc_idx))
  control_sa <- df[control_sa_idx, ]
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  
  if (length(cotrt_types) > 0) {
    df_list[cotrt_types] <- gDRutils::loop(cotrt_types, function(x) {
      unique_cotrt <- unique(df_list[[x]][, c(cl, drug_ids[["drug_name"]])])
      unique_cotrt_ctrl <- unique(control[control[[cl]] %in% unique_cotrt[[cl]] &
                                            control[[drug_ids[["drug_name"]]]] %in%
                                            untreated_tag, ][, c(cl, drug_ids[["drug_name"]])])
      cotrt_matching <- rbind(unique_cotrt, unique_cotrt_ctrl)
      df_merged <- rbind(df_list[[x]], merge(cotrt_matching, control))
      if (x == "matrix") {
        matrix_data <- rbind(df_merged, df_list[["single-agent"]])
        matrix_data[, conc_idx][is.na(matrix_data[, conc_idx])] <- 0
        matrix_data[, codrug_drug_id][is.na(matrix_data[, codrug_drug_id])] <- untreated_tag[1]
        matrix_data
      } else {
        df_merged
      }
    })
  }
  
  if ("single-agent" %in% names(df_list)) {
    sa_idx <- gDRutils::loop(grep(drug_ids[["concentration"]], drug_ids, value = TRUE), function(x)
           which(!df_list[["single-agent"]][, x] == 0))
    sa_idx[["concentration"]] <- NULL
    for (codrug in names(sa_idx)) {
      codrug_cols <- grep(as.numeric(gsub("\\D", "", codrug)), drug_ids, value = TRUE)
      df_list[["single-agent"]][sa_idx[[codrug]], drug_ids[c("drug_name", "drug", "drug_moa", "concentration")]] <- 
        df_list[["single-agent"]][sa_idx[[codrug]], codrug_cols]
    }
    df_list[["single-agent"]][, codrug_ids] <- NULL
    df_list[["single-agent"]] <- rbind(df_list[["single-agent"]],
                                       control_sa[, names(df_list[["single-agent"]])])
  }
  
  Map(function(x) {
    x[, names(x) != type_col]
    }, df_list)
}
