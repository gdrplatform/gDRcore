#' Map treated conditions to their respective references.
#'
#' Map treated conditions to their respective Day0, untreated, or single-agent 
#' references using condition metadata.
#'
#' @param trt_md data.table of treated metadata. 
#' @param ref_md data.table of untreated metadata.
#' @param override_untrt_controls named list indicating what treatment metadata 
#' fields should be used as a control. Defaults to \code{NULL}.
#' @param ref_cols character vector of the names of reference columns to 
#' include. Likely obtained from \code{identify_keys()}.
#' @param ref_type string of the reference type to map to.
#' Should be one of \code{c("Day0", "untrt_Endpoint")}.
#'
#' @return named list mapping treated metadata to untreated metadata.
#'
#' @details If \code{override_untrt_controls} is specified, the values in the 
#' named list will supersede the values in \code{trt_md} during the matching 
#' process. This is useful for mapping treatments to specific "standard" 
#' untreated controls.
#' 
#' @seealso identify_keys
#' @keywords map_df
#' @export
#'
map_df <- function(trt_md, 
                   ref_md, 
                   override_untrt_controls = NULL, 
                   ref_cols, 
                   ref_type = c("Day0", "untrt_Endpoint")) {
  
  # Assertions:
  checkmate::assert_data_table(trt_md)
  checkmate::assert_data_table(ref_md)
  checkmate::assert_list(override_untrt_controls, null.ok = TRUE)
  checkmate::assert_character(ref_cols)
  
  ref_type <- match.arg(ref_type)
  
  duration_col <- gDRutils::get_env_identifiers("duration")
  conc_cols <- unlist(gDRutils::get_env_identifiers(
    c("concentration", "concentration2"), 
    simplify = FALSE
  ))
  
  if (ref_type == "Day0") {
    ref_md <- ref_md[get(duration_col) == 0, ]
  }
  
  conc <- cbind(array(0, nrow(ref_md)), 
                ref_md[, intersect(names(ref_md), conc_cols), with = FALSE])
  is_ref_conc <- rowSums(conc == 0) == ncol(conc)
  
  if (ref_type == "Day0") {
    matching_list <- list(T0 = ref_md[[duration_col]] == 0, conc = is_ref_conc)
    matchFactor <- "T0"
  } else if (ref_type == "untrt_Endpoint") {
    matching_list <- list(conc = is_ref_conc)
    matchFactor <- duration_col 
  }
  
  trt_rnames <- trt_md$rn
  ref_rnames <- ref_md$rn
  
  present_ref_cols <- intersect(ref_cols, names(ref_md))
  
  # 1. Exact matches vectorized
  exact_out <- if (length(present_ref_cols) == 0) {
    stats::setNames(replicate(length(trt_rnames), character(0), simplify = FALSE), trt_rnames)
  } else {
    match_l <- grr_matches(
      do.call("paste", trt_md[, ..present_ref_cols]),
      do.call("paste", ref_md[, ..present_ref_cols]),
      all.y = FALSE,
      list = TRUE
    )
    names(match_l) <- trt_rnames
    lapply(match_l, function(x) ref_rnames[sort(x)])
  }
  
  # 2. Search for non-exact matches or overrides
  # We return a list containing both the mapped references and any messages
  res_list <- lapply(seq_along(trt_rnames), function(i) {
    treatment <- trt_rnames[i]
    msg <- NULL
    
    if (length(exact_out[[treatment]]) == 0 ||
        any(is.na(exact_out[[treatment]])) ||
        !is.null(override_untrt_controls)) {
      
      refs <- lapply(present_ref_cols, function(y) {
        ref_md[[y]] == trt_md[rn == treatment, ..y][[1]]
      })
      names(refs) <- present_ref_cols
      
      if (!is.null(override_untrt_controls)) {
        for (overridden in names(override_untrt_controls)) {
          if (overridden %in% names(ref_md)) {
            refs[[overridden]] <- ref_md[[overridden]] == override_untrt_controls[[overridden]]
          }
        }
      }
      
      all_checks <- c(refs, matching_list)
      match_mx <- do.call("rbind", all_checks)
      
      # Calculate scores
      idx <- colSums(match_mx)
      
      # score for metadata columns only (to avoid matching on Duration alone)
      meta_score <- if (length(present_ref_cols) > 0) {
        colSums(match_mx[present_ref_cols, , drop = FALSE])
      } else {
        rep(0, ncol(match_mx))
      }
      
      if (matchFactor %in% rownames(match_mx)) {
        idx <- idx * match_mx[matchFactor, ]
      }
      
      # Identify best matches where at least one metadata field matches
      if (any(idx > 0 & meta_score > 0, na.rm = TRUE)) {
        valid_idx <- which(idx > 0 & meta_score > 0)
        match_idx <- valid_idx[which(idx[valid_idx] == max(idx[valid_idx]))]
        
        msg <- sprintf("Found partial match: ('%s') for treatment: ('%s')", 
                       paste(ref_rnames[match_idx], collapse = ", "), treatment)
        
        return(list(ref = ref_rnames[match_idx], msg = msg))
      } else {
        msg <- sprintf("No partial match found for treatment: ('%s')", treatment)
        return(list(ref = character(0), msg = msg))
      }
    } else {
      return(list(ref = exact_out[[treatment]], msg = NULL))
    }
  })
  
  # Flatten result list into mapping and messages
  out <- lapply(res_list, `[[`, "ref")
  names(out) <- trt_rnames
  
  msgs <- unlist(lapply(res_list, `[[`, "msg"))
  if (length(msgs) > 0) {
    futile.logger::flog.info(paste0(msgs, collapse = "\n"))
  }
  
  out
}

#' Map references
#' 
#' @param mat_elem data.table input
#' @param rowData_colnames character vector of variables for mapping
#' @keywords map_df
#' @return list
#' @export
.map_references <- function(mat_elem, 
                            rowData_colnames = c(gDRutils::get_env_identifiers("duration"), 
                                                 paste0(c("drug", "drug_name", "drug_moa"), "3"))) {
  
  checkmate::assert_data_table(mat_elem)
  checkmate::assert_character(rowData_colnames, null.ok = TRUE)
  
  clid <- gDRutils::get_env_identifiers("cellline")
  checkmate::assert_choice(clid, colnames(mat_elem))
  
  # Avoid recycling code using helper
  tag_info <- .get_untreated_tag_count(mat_elem, c("drug_name", "drug_name2"))
  valid <- tag_info$valid_cols
  
  is_untrt <- tag_info$ntag == tag_info$num_cols
  is_ref   <- tag_info$ntag != 0L & !is_untrt
  
  cotrt_var <- setdiff(rowData_colnames, 
                       gDRutils::get_env_identifiers(
                         c("drug", "drug_name", "drug_moa", paste0(c("drug", "drug_name", "drug_moa"), "2")), 
                         simplify = FALSE
                       )
  )
  cotrt_var <- intersect(cotrt_var, colnames(mat_elem))
  
  mat_elem$rownames <- as.character(seq_len(nrow(mat_elem)))
  trt_elem <- mat_elem[!is_ref & !is_untrt]
  
  out <- vector("list", nrow(trt_elem))
  names(out) <- trt_elem$rownames
  
  if (any(is_ref)) {
    ref_elem <- mat_elem[is_ref]
    
    trtNames <- rep(trt_elem$rownames, length(valid))
    refNames <- rep(ref_elem$rownames, length(valid))
    
    trt <- do.call(paste, do.call(rbind, lapply(valid, function(x) {
      stats::setNames(trt_elem[, c(clid, x), with = FALSE], c(clid, "drug"))
    })))
    
    ref <- do.call(paste, do.call(rbind, lapply(valid, function(x) {
      stats::setNames(ref_elem[, c(clid, x), with = FALSE], c(clid, "drug"))
    })))
    
    matchTrtRef <- grr_matches(trt, ref, list = FALSE, all.y = FALSE)
    matchTrtRef[["x"]] <- trtNames[matchTrtRef[["x"]]]
    matchTrtRef[["y"]] <- refNames[matchTrtRef[["y"]]]
    out <- split(matchTrtRef[["y"]], matchTrtRef[["x"]])
    
    if (length(cotrt_var) > 0) {
      for (i in names(out)) {
        ref_idx <- lapply(na.omit(out[[i]]), function(x) {
          all(ref_elem[rownames == x, ..cotrt_var] == trt_elem[rownames == i, ..cotrt_var])
        })
        out[[i]] <- out[[i]][unlist(ref_idx)]
      }
    }
  }
  out
}

#' Identify untreated rows based on Drug treatment alone
#' 
#' @param mat_elem data.table input
#' @keywords map_df
#' @return logical vector
#' @export
map_untreated <- function(mat_elem) {
  checkmate::assert_data_table(mat_elem)
  tag_info <- .get_untreated_tag_count(mat_elem)
  tag_info$ntag == tag_info$num_cols
}

#' Get the count of untreated tags per row
#' 
#' @param mat_elem data.table input data frame
#' @param drug_identifier_keys character vector of keys to look up identifiers
#'
#' @return list containing ntag, num_cols, and valid_cols
#' @keywords internal
.get_untreated_tag_count <- function(mat_elem, 
                                     drug_identifier_keys = c("drug_name", "drug_name2", "drug_name3")) {
  
  checkmate::assert_data_table(mat_elem)
  valid_cols <- unlist(
    intersect(
      gDRutils::get_env_identifiers(drug_identifier_keys, simplify = FALSE),
      colnames(mat_elem)
    )
  )
  
  if (length(valid_cols) == 0) {
    stop(sprintf("None of the drug identifiers [%s] found", paste(drug_identifier_keys, collapse = ", ")))
  }
  
  untrt_tag <- gDRutils::get_env_identifiers("untreated_tag")
  ntag <- rowSums(mat_elem[, lapply(.SD, `%in%`, untrt_tag), .SDcols = valid_cols])
  
  list(ntag = ntag, num_cols = length(valid_cols), valid_cols = valid_cols)
}