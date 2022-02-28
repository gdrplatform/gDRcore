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
  checkmate::assert_vector(override_untrt_controls, null.ok = TRUE)
  checkmate::assert_character(ref_cols)
  checkmate::assert_character(ref_type)
  
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
  ref_rnames <- rownames(ref_md)
  
  # define matrix with matching metadata
  present_ref_cols <- intersect(ref_cols, names(ref_md))
  names(present_ref_cols) <- present_ref_cols
  
  out <- list("vector", length(trt_rnames))
  msgs <- NULL
  
  # Parallel computing
  clusters <- parallel::makeCluster(detect_cores(), type = "FORK")
  doParallel::registerDoParallel(clusters)
  
  out <- foreach::foreach(i = seq_along(trt_rnames)) %dopar% {
    treatment <- trt_rnames[i]
    refs <- lapply(present_ref_cols, function(y) {
      ref_md[, y] == trt_md[treatment, y]
    })
    
    if (!is.null(override_untrt_controls)) {
      for (overridden in names(override_untrt_controls)) {
        refs[[overridden]] <- ref_md[, overridden] == override_untrt_controls[[overridden]]
      }
    }
    
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
    ref_rnames[match_idx] # TODO: Check that this properly handles NAs or NULLs. 
  }
  parallel::stopCluster(clusters)
  
  futile.logger::flog.info(paste0(msgs, collapse = "\n"))
  names(out) <- trt_rnames
  out
}


#' @details
#' Using the given rownames, map the treated and reference conditions.
.map_references <- function(mat_elem) {
  clid <- gDRutils::get_env_identifiers("cellline")
  valid <- unlist(intersect(c(gDRutils::get_env_identifiers(c("drug_name", "drug_name2"), simplify = FALSE)),
                     colnames(mat_elem)))
  drug_cols <- mat_elem[valid]

  untrt_tag <- gDRutils::get_env_identifiers("untreated_tag")
  mat_elem[mat_elem == untrt_tag[[2]]] <- untrt_tag[[1]]
  pattern <- paste0(sprintf("^%s$", untrt_tag), collapse = "|")
  has_tag <- as.data.frame(lapply(drug_cols, function(x) grepl(pattern, x)))
  ntag <- rowSums(has_tag)

  is_untrt <- ntag == length(valid)
  is_ref <- ntag != 0L & !is_untrt

  trt <- mat_elem[!is_ref & !is_untrt, ]
  ref <- mat_elem[is_ref, ]

  out <- vector("list", nrow(trt))
  names(out) <- rownames(trt)
  
  if (any(is_ref)) {
    # store rownames of trt and ref and replicate them based on the length of drug columns
    trtNames <- rep(rownames(trt), length(valid))
    refNames <- rep(rownames(ref), length(valid))
    
    # split data.frames to simple model with clid column and drug column
    trt <- lapply(valid, function(x) trt[, c(clid, x)])
    trt <- do.call(paste, do.call(rbind, lapply(trt, function(x) setNames(x, names(trt[[1]])))))
    
    ref <- lapply(valid, function(x) ref[, c(clid, x)])
    ref <- do.call(paste, do.call(rbind, lapply(ref, function(x) setNames(x, names(ref[[1]])))))
    
    # match trt and ref
    matchTrtRef <- grr::matches(trt, ref, list = FALSE, all.y = FALSE)
    matchTrtRef[["x"]] <- trtNames[matchTrtRef[["x"]]]
    matchTrtRef[["y"]] <- refNames[matchTrtRef[["y"]]]
    
    out <- split(matchTrtRef[["y"]], matchTrtRef[["x"]])
    out
  } else {
    out
  }
}
