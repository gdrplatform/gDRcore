#' Map treated conditions to their respective references.
#'
#' Map treated conditions to their respective Day0, untreated, or single-agent references using condition metadata.
#'
#' @param trt_md data.frame of treated metadata. 
#' @param ref_md data.frame of untreated metadata.
#' @param override_untrt_controls named list indicating what treatment metadata fields should be used as a control.
#' Defaults to \code{NULL}.
#' @param ref_cols character vector of the names of reference columns to include.
#' Likely obtained from \code{identify_keys()}.
#' @param ref_type string of the reference type to map to.
#' Should be one of \code{c("Day0", "untrt_Endpoint", "ref_Endpoint")}.
#'
#' @return named list mapping treated metadata to untreated metadata.
#'
#' @details If \code{override_untrt_controls} is specified, 
#' TODO: FILL ME!
#' 
#' @seealso identify_keys
#' @export
#'
map_df <- function(trt_md, 
                   ref_md, 
                   override_untrt_controls = NULL, 
                   ref_cols, 
                   ref_type = c("Day0", "untrt_Endpoint")) {
  # Assertions:
  checkmate::assert_class(trt_md, "data.frame")
  checkmate::assert_class(ref_md, "data.frame")
  ref_type <- match.arg(ref_type)

  duration_col <- gDRutils::get_env_identifiers("duration")

  conc <- cbind(array(0, nrow(ref_md)), # padding to avoid empty df;
    ref_md[, agrep("Concentration", colnames(ref_md)), drop = FALSE])
  is_ref_conc <- apply(conc, 1, function(z) {
    all(z == 0)
    })

  if (ref_type == "Day0") {
    # Identifying which of the durations have a value of 0.
    matching_list <- list(T0 = ref_md[, duration_col] == 0, conc = is_ref_conc)
    matchFactor <- "T0"
  } else if (ref_type == "untrt_Endpoint") {
    matching_list <- list(conc = is_ref_conc)
    matchFactor <- duration_col 
  }

  trt_rnames <- rownames(trt_md)

  # define matrix with matching metadata
  present_ref_cols <- intersect(ref_cols, names(ref_md))
  names(present_ref_cols) <- present_ref_cols

  out <- list("vector", length(trt_rnames))
  msgs <- NULL
  for (i in seq_len(length(trt_rnames))) {
    treatment <- trt_rnames[i]
    refs <- lapply(present_ref_cols, function(y) {
      ref_md[, y] == trt_md[treatment, y]
      })

    if (!is.null(override_untrt_controls)) {
        for (overridden in names(override_untrt_controls)) {
            refs[[overridden]] <- ref_md[, overridden] == override_untrt_controls[[overridden]]
    }}

    all_checks <- c(refs, matching_list)
    match_mx <- do.call("rbind", all_checks)
    rownames(match_mx) <- names(all_checks)
    match_idx <- which(apply(match_mx, 2, all)) # test matching conditions
    if (length(match_idx) == 0) {
      # No exact match, try to find best match (as many metadata fields as possible).
      # TODO: rowSums?
      idx <- apply(match_mx, 2, function(y) sum(y, na.rm = TRUE)) 
      # TODO: Sort this out so that it also takes the average in case multiple are found.
      idx <- idx * match_mx[matchFactor, ]
      
      if (any(idx > 0)) {
        match_idx <- which.max(idx)
        msgs <- c(msgs, sprintf("Found partial match: ('%s') for treatment: ('%s')",
          rownames(ref_md)[match_idx], treatment))
      } else { # failed to find any potential match
        msgs <- c(msgs, sprintf("No partial match found for treatment: ('%s')", treatment))
      }
    }
    out[[i]] <- rownames(ref_md)[match_idx] # TODO: Check that this properly handles NAs. 
  }
  futile.logger::flog.info(paste0(msgs, collapse = "\n"))
  names(out) <- trt_rnames
  out
}


.map_cotreatments <- function() {
  # TODO: Temporarily remove the cotreatments endpoint from create_SE,
  # since I think we will want to be able to do that in fit_SE later anyway.
  # creates another list for the co-treatment end points that are missing

      ref_type <- "ref_Endpoint"
      cotrt_ref <- ref_maps[[ref_type]][[trt]]  
      if (length(cotrt_ref) > 0L) {
        cotrt_df <- dfs[groupings %in% cotrt_ref, , drop = FALSE]
        cotrt_df <- create_control_df(
          cotrt_df, 
          control_cols = Keys[[ref_type]], 
          control_mean_fxn, 
          out_col_name = "RefReadout"
        )
      } else if (length(ref_maps[[paste0("cotrt_", ref_type)]][[trt]]) > 0L) {
        cotrt_ref <- ref_maps[[paste0("cotrt_", ref_type)]][[trt]]
        cotrt_df <- dfs[groupings %in% cotrt_ref, , drop = FALSE]

        if (any(cotrt_df$Concentration == treated$Concentration_2[treated$groupings %in% trt])) {
            cotrt_df <- create_control_df(
                cotrt_df[cotrt_df$Concentration == treated$Concentration_2[treated$groupings %in% trt], ], 
                control_cols = Keys[[ref_type]], 
                control_mean_fxn, 
                out_col_name = "RefReadout"
            )
        } else {
          cotrt_df <- infer_control_df(
            cotrt_df,
            treated$Concentration_2[treated$groupings %in% trt],
            control_cols = Keys[[ref_type]],
            control_mean_fxn,
            out_col_name = "RefReadout"
            )
        }
      } else {
        # Set the cotrt reference to NA if not found 
        cotrt_df <- untrt_df 
        cotrt_df$UntrtReadout <- NA
        colnames(cotrt_df)[grepl("UntrtReadout", colnames(cotrt_df))] <- "RefReadout"
      }
   
  ref_maps[["cotrt_ref_Endpoint"]] <- NULL
  # focus on cases where the reference may be as primary drug (common in co-treatment experiments)
  if (paste0(identifiers$drug, "_2") %in% colnames(treated)) {
    
    # NOTE: may have to deal with override_untrt_controls 

    ref_type <- "ref_Endpoint"
    missing_cotrt <- vapply(ref_maps[[ref_type]], function(x) {
      length(x) == 0L
      }, TRUE)
    
    # Then look amongst the treated to fill any missing cotrt references.
    if (any(missing_cotrt)) {
        # try to find the co-treated reference among treated data (with Drug/Drug_2 swap)    
        pseudo_untreated <- treated[treated$Concentration_2 == 0, ]
        # remove Concentration as is will have to be matched with the Concentration
        pseudo_untreated$Concentration_2 <- NULL 
        
        # swap columns related to drug and drug_2
        idx_1 <- which(colnames(pseudo_untreated) %in% 
            c(identifiers$drug, 
              identifiers$drugname,
              identifiers$drug_moa))
        idx_2 <- which(colnames(pseudo_untreated) %in% 
            paste0(c(identifiers$drug, 
                identifiers$drugname,
                identifiers$drug_moa), "_2"))
        colnames(pseudo_untreated)[idx_1] <- paste0(colnames(pseudo_untreated)[idx_1], "_2")
        colnames(pseudo_untreated)[idx_2] <- gsub("_2", "", colnames(pseudo_untreated)[idx_2])

        ref_maps[["cotrt_ref_Endpoint"]] <- map_df(treated[missing_cotrt, ], pseudo_untreated, 
            override_untrt_controls = override_untrt_controls, ref_cols = Keys[[ref_type]], ref_type = ref_type)

    } # we may be able to extend to other cases if applicable
  }

  ## TODO: Check for failed cotreatment mappings. 
}


#' @details
#' Using the given rownames, map the treated and reference conditions.
.map_references <- function(mat_elem) {
  clid <- get_env_identifiers("cellline")
  valid <- intersect(c(get_env_identifiers(c("drugname", "drugname2"))), colnames(mat_elem))
  drug_cols <- mat_elem[valid]

  untrt_tag <- get_env_identifiers("untreated_tag")
  pattern <- paste0(sprintf("^%s$", untrt_tag), collapse = "|")
  pattern <- "vehicle|untreated"
  has_tag <- as.data.frame(lapply(drug_cols, function(x) grepl(pattern, x)))
  ntag <- rowSums(has_tag)

  is_untrt <- ntag == length(valid)
  is_ref <- ntag != 0L & !is_untrt

  trt <- mat_elem[!is_ref & !is_untrt, ]
  ref <- mat_elem[is_ref, ]

  out <- vector("list", nrow(trt))
  names(out) <- rownames(trt)

  if (any(is_ref)) {
    compare_cols <- c(valid, clid)
    
    for (t in rownames(trt)) {
      refs <- NULL
      for (r in rownames(ref)) {
        if (all(setdiff(as.character(ref[r, compare_cols]), as.character(trt[t, compare_cols])) %in% untrt_tag)) {
          refs <- c(refs, r)
        }
      }
      out[[t]] <- refs
    }
  }
  out
}
