#' Identify type of data
#'
#' @param df data.frame of raw drug response data containing both treated and untreated values
#' @param cotreatment_conc integer of maximum number of concentration of co-treatment
#' to classify as cotreatment data type.
#' Defaults to \code{4}.
#' @param codilution_conc integer of maximum number of concentration ratio of co-treatment
#' to classify as codilution data type.
#' Defaults to \code{2}.
#' @param matrix_conc integer of minimum number of concentration pairs of co-treatment
#' to classify as matrix data type.
#' Defaults to \code{4}.
#'
#' @return data.frame of raw drug response data with additional column `type` with the info of
#' data type for a given row of data.frame
#' @export
#'
#' @author Bartosz Czech <bartosz.czech@@contractors.roche.com>
identify_data_type <- function(df,
                               cotreatment_conc = 4,
                               codilution_conc = 2,
                               matrix_conc = 4
                               ) {
  
  # find the pairs of drugs with relevant metadat
  drug_ids <- unlist(gDRutils::get_env_identifiers(c("drugname", "drugname2"), simplify = FALSE)) 
  conc_ids <- unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2"), simplify = FALSE))
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  cell <- gDRutils::get_env_identifiers("cellline")
  cols_pairs <- intersect(names(df),  c(drug_ids, cell))
  drug_pairs <- unique(df[, cols_pairs])
  # drug_pairs_list <- split(drug_pairs[, drug_ids], drug_pairs[[cell]]) #nolint
  
  cnt <- seq_len(nrow(df))
  df$type <- NA
  # loop through the pairs to assess the number of individual concentration pairs
  for (idp in 1:nrow(drug_pairs)) {
    # reverse engineer the type of combination experiment
    df_matching <- merge(cbind(df, cnt), drug_pairs[idp, ])
    df_primary_matching <- merge(cbind(df, cnt), drug_pairs[idp, c(drug_ids[["drugname"]], cell)])
    matching_idx <- df_matching$cnt
    type <- if (is.null(df[matching_idx, conc_ids[["concentration2"]]])) {
      if (all(df[matching_idx, drug_ids[["drugname"]]] %in% untreated_tag)) {
      "control"
    }
    else {
      "single-agent"
    }
      } else if (all(df_primary_matching[conc_ids[["concentration2"]]] == 0) &
                               any(df_primary_matching[conc_ids[["concentration"]]] != 0)) {
        "single-agent"
      } else if (all(df[matching_idx, drug_ids[["drugname"]]] %in% untreated_tag &
                     df[matching_idx, drug_ids[["drugname2"]]] %in% untreated_tag | 
                     df[matching_idx, drug_ids[["drugname"]]] %in% untreated_tag) |
                     all((df[matching_idx, drug_ids[["drugname2"]]] %in% untreated_tag) &
                         any(df_primary_matching[conc_ids[["concentration2"]]] != 0))) {
        "control"
      } else {
        NA
      }
    type
    df[matching_idx, "type"]  <- type
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
      } else if (n_conc_pairs == length(conc_1) * length(conc_2) & length(conc_2) >= matrix_conc) {
        "matrix"
      } else if (length(conc_2) < cotreatment_conc) {
        "cotreatment"
      } else {
        "other"
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
  conc_idx <- intersect(colnames(df), 
                        unlist(gDRutils::get_env_identifiers(c("concentration", "concentration2"), simplify = FALSE)))
  drug <- gDRutils::get_env_identifiers("drug")
  cl <- gDRutils::get_env_identifiers("cellline")
  df_list <- split(df, df[[type_col]])
  types <- setdiff(names(df_list), "control")
  control <- df_list[["control"]]
  df_list[["control"]] <- NULL
  cotrt_types <- setdiff(names(df_list), "single-agent")
  control_sa_idx <- which(rowSums(control[, conc_idx, drop = FALSE] == 0) == length(conc_idx))
  control_sa <- control[control_sa_idx, ]
  untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
  
  if ("single-agent" %in% names(df_list)) {
    df_list[["single-agent"]][, grep("_2", names(df_list[["single-agent"]]))] <- NULL
    df_list[["single-agent"]] <- rbind(df_list[["single-agent"]],
                                       control_sa[, names(df_list[["single-agent"]])])
  }
  if (length(cotrt_types) > 0) {
    df_list[cotrt_types] <- lapply(df_list[cotrt_types], function(x) {
      unique_cotrt <- unique(x[, c(cl, gDRutils::get_env_identifiers("drug"))])
      unique_cotrt_ctrl <- unique(control[control[[cl]] %in% unique_cotrt[[cl]] &
                                            control[[drug]] %in% untreated_tag, ][, c(cl, drug)])
      cotrt_matching <- rbind(unique_cotrt, unique_cotrt_ctrl)
      rbind(x, merge(cotrt_matching, control))
    })
  }
  Map(function(x) {
    x[, names(x) != type_col]
    }, df_list)
}
