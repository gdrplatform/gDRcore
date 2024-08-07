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
#' Should be one of \code{c("Day0", "untrt_Endpoint", "ref_Endpoint")}.
#' 
#' @examples 
#' n <- 64
#' md_df <- data.table::data.table(
#'   Gnumber = rep(c("vehicle", "untreated", paste0("G", seq(2))), each = 16), 
#'   DrugName = rep(c("vehicle", "untreated", paste0("GN", seq(2))), each = 16), 
#'   clid = paste0("C", rep_len(seq(4), n)),
#'   CellLineName = paste0("N", rep_len(seq(4), n)),
#'   replicates = rep_len(paste0("R", rep(seq(4), each = 4)), 64),
#'   drug_moa = "inhibitor",
#'   ReferenceDivisionTime = rep_len(c(120, 60), n),
#'   Tissue = "Lung",
#'   parental_identifier = "CL12345",
#'   Duration = 160
#' )
#' md_df <- unique(md_df)
#' ref <- md_df$Gnumber %in% c("vehicle", "untreated")
#' ref_df <- md_df[ref, ]
#' trt_df <- md_df[!ref, ]
#' Keys <- identify_keys(trt_df)
#' ref_type <- "untrt_Endpoint"
#' map_df(
#'   trt_df, 
#'   ref_df, 
#'   ref_cols = Keys[[ref_type]],
#'   ref_type = ref_type
#' )
#'
#' @return named list mapping treated metadata to untreated metadata.
#'
#' @details If \code{override_untrt_controls} is specified, 
#' TODO: FILL ME!
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
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  checkmate::assert_character(ref_cols)
  checkmate::assert_character(ref_type)
  
  ref_type <- match.arg(ref_type)
  
  duration_col <- gDRutils::get_env_identifiers("duration")
  conc_cols <- unlist(gDRutils::get_env_identifiers(
    c("concentration", "concentration2"), 
    simplify = FALSE
  ))

  if (ref_type == "Day0") {
    ref_md <- ref_md[get(duration_col) == 0, ]
  }
  
  conc <- cbind(array(0, nrow(ref_md)), # padding to avoid empty df;
                ref_md[, intersect(names(ref_md), conc_cols), drop = FALSE])
  is_ref_conc <- rowSums(conc == 0) == ncol(conc)
  
  
  if (ref_type == "Day0") {
    # Identifying which of the durations have a value of 0.
    matching_list <- list(T0 = ref_md[[duration_col]] == 0, conc = is_ref_conc)
    matchFactor <- "T0"
  } else if (ref_type == "untrt_Endpoint") {
    matching_list <- list(conc = is_ref_conc)
    matchFactor <- duration_col 
  }
  
  trt_rnames <- trt_md$rn
  ref_rnames <- ref_md$rn
  
  # define matrix with matching metadata
  present_ref_cols <- intersect(ref_cols, names(ref_md))
  names(present_ref_cols) <- present_ref_cols
  msgs <- NULL
  
  # 1. there are no matches (present_ref_cols is empty)
  exact_out <- if (length(present_ref_cols) == 0) {
    out <- stats::setNames(replicate(length(trt_rnames), character(0), simplify = FALSE), trt_rnames)
    names(out) <- trt_rnames
    out
  } else {
    # 2. search for exact matches found in the vectorised way
    #    cases with non-exact matches will be returned as NAs
    match_l <-
      grr_matches(
        do.call("paste", trt_md[, present_ref_cols, with = FALSE]),
        do.call("paste", ref_md[, present_ref_cols, with = FALSE]),
        all.y = FALSE,
        list = TRUE
      )
    names(match_l) <- trt_rnames
    lapply(match_l, function(x) {
      ref_rnames[sort(x)]
    })
  }
 
  # 3. only exact matches found 
  out <- if (!anyNA(exact_out) && is.null(override_untrt_controls)) {
    exact_out
    # 4. not all entres are exact matches
    # 4.1 search for potential non-exact matches
    # 4.2 support cases with overriden untreated controls
    # this logic is pretty slow currently 
  } else {

  out <- lapply(seq_along(trt_rnames), function(i) {
    treatment <- trt_rnames[i]
    if (all(is.na(exact_out[[treatment]])) || !is.null(override_untrt_controls)) {
      
      refs <- lapply(present_ref_cols, function(y) {
        unname(unlist(ref_md[, y, with = FALSE]) == unlist(trt_md[which(trt_md$rn == treatment),
                                                                  y, with = FALSE]))
      })
      
      if (!is.null(override_untrt_controls)) {
        for (overridden in names(override_untrt_controls)) {
          refs[[overridden]] <-
            unname(unlist(ref_md[, overridden, with = FALSE]) == override_untrt_controls[[overridden]])
        }
      }
      
      all_checks <- c(refs, matching_list)
      match_mx <- do.call("rbind", all_checks)
      rownames(match_mx) <- names(all_checks)
      match_idx <- which(colSums(match_mx) == ncol(match_mx))
      # No exact match, try to find best match (as many metadata fields as 
      # possible).
      # TODO: rowSums?
      idx <- colSums(match_mx)
      # TODO: Sort this out so that it also takes the average in case multiple 
      # are found.
      idx <- idx * match_mx[matchFactor, ]
      
      if (any(idx > 0, na.rm = TRUE)) {
        match_idx <- which.max(idx)
        msgs <- c(
          msgs, 
          sprintf(
            "Found partial match: ('%s') for treatment: ('%s')",
            rownames(ref_md)[match_idx], treatment
          )
        )
      } else { # failed to find any potential match
        msgs <- c(
          msgs, 
          sprintf("No partial match found for treatment: ('%s')", treatment)
        )
      }
      ref_rnames[match_idx] # TODO: Check that this properly handles 
                            # NAs or NULLs.
    } else {
      exact_out[[treatment]]
    }
  })
  names(out) <- trt_rnames
  out
  }
  
  futile.logger::flog.info(paste0(msgs, collapse = "\n"))
  out
}

#' Map references
#' 
#' @param mat_elem input data frame
#' @param rowData_colnames character vector of variables for the mapping of reference treatments
#'
#' @details
#' Using the given rownames, map the treated and reference conditions.
#' 
#' @keywords map_df
#' @return list
#' 
.map_references <- function(mat_elem, 
                            rowData_colnames = c(gDRutils::get_env_identifiers("duration"), 
                                                paste0(c("drug", "drug_name", "drug_moa"), "3"))) {
  
  clid <- gDRutils::get_env_identifiers("cellline")
  # variables for the mapping of treatments
  valid <- unlist(
    intersect(
      c(
        gDRutils::get_env_identifiers(
          c("drug_name", "drug_name2"), 
          simplify = FALSE
        )
      ),
      colnames(mat_elem)
    )
  )
  # variables for the mapping of reference treatments
  cotrt_var <- setdiff(rowData_colnames, 
     gDRutils::get_env_identifiers(
       c("drug", "drug_name", "drug_moa", paste0(c("drug", "drug_name", "drug_moa"), "2")), 
       simplify = FALSE
     )
  )
  
  drug_cols <- mat_elem[, valid, with = FALSE]

  untrt_tag <- gDRutils::get_env_identifiers("untreated_tag")
  has_tag <- lapply(drug_cols, function(x) x %in% untrt_tag)
  data.table::setDT(has_tag)
  ntag <- rowSums(has_tag)

  is_untrt <- ntag == length(valid)
  is_ref <- ntag != 0L & !is_untrt
  
  mat_elem$rownames <- as.character(seq_len(nrow(mat_elem)))

  # columns with the primary data for the treatment and reference
  trt_elem <- mat_elem[which(!is_ref & !is_untrt)]
  
  out <- vector("list", nrow(trt_elem))
  names(out) <- trt_elem$rownames
  
  if (any(is_ref)) {
    ref_elem <- mat_elem[which(is_ref), ]
    # columns with the matching data for the treatment and reference
    ref_cotrt <- mat_elem[which(is_ref), cotrt_var, with = FALSE]

    # store rownames of trt_elem and ref_elem and replicate them based on the length of 
    # drug columns
    trtNames <- rep(trt_elem$rownames, length(valid))
    refNames <- rep(ref_elem$rownames, length(valid))
    
    # split data.tables to simple model with clid column and drug column
    
    trt <- lapply(valid, function(x) {
      colnames <- c(clid, x) 
      trt_elem[, colnames, with = FALSE] 
    })
    trt <- do.call(
      paste, 
      do.call(
        rbind, 
        lapply(trt, function(x) stats::setNames(x, names(trt[[1]])))
      )
    )
    
    ref <- lapply(valid, function(x) {
      colnames <- c(clid, x) 
      ref_elem[, colnames, with = FALSE] 
    })
    ref <- do.call(
      paste, 
      do.call(
        rbind, 
        lapply(ref, function(x) stats::setNames(x, names(ref[[1]])))
      )
    )
    # match trt and ref
    matchTrtRef <- grr_matches(trt, ref, list = FALSE, all.y = FALSE)
    matchTrtRef[["x"]] <- trtNames[matchTrtRef[["x"]]]
    matchTrtRef[["y"]] <- refNames[matchTrtRef[["y"]]]
    out <- split(matchTrtRef[["y"]], matchTrtRef[["x"]])
    # match the additional variables in the treatment
    if (length(ref_cotrt) && length(cotrt_var)) {
      for (i in names(out)) {
        # matching the ref_elem to the trt_elem for the cotrt_var
        ref_idx <- lapply(na.omit(out[[i]]), function(x) {
          ref_elem[rn == x, cotrt_var, with = FALSE] ==
              trt_elem[rn == i, cotrt_var, with = FALSE]
          })
        out[[i]] <- out[[i]][unlist(lapply(ref_idx, all))]
      }
    }
    out
  } else {
    out
  }
}


#' Identify untreated rows based on Drug treatment alone
#' 
#' @param mat_elem input data frame
#'
#' @details
#' Using the given rownames, map the untreated conditions
#' 
#' @keywords map_df
#' @return list
#' 
map_untreated <- function(mat_elem) {
  # TODO: avoid recycling code of map_references above
  #
  clid <- gDRutils::get_env_identifiers("cellline")
  valid <- unlist(
    intersect(
      c(
        gDRutils::get_env_identifiers(
          c("drug_name", "drug_name2", "drug_name3"), 
          simplify = FALSE
        )
      ),
      colnames(mat_elem)
    )
  )
  
  drug_cols <- mat_elem[, valid, with = FALSE]

  untrt_tag <- gDRutils::get_env_identifiers("untreated_tag")
  ntag <- rowSums(drug_cols[,
                            lapply(.SD, `%in%`, untrt_tag),
                            .SDcols = names(drug_cols)])
  ntag == length(valid)
}
