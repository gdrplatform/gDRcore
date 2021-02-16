#' @export
#'
calculate_reference_values_SE <- function(se) {
  fitting_refs <- SummarizedExperiment::assays(se)[["Control"]] 
  for (i in rownames(se)) {
    for (j in colnames(se)) {
      ref_df <- refs[i, j]
      ## Calculate the reference metrics to be used in downstreaming curve fitting to calculate metrics. 
      ## These values will later serve as the asymptotic values for the fit_curves function.
      ref_df$RefRelativeViability <- round(ref_df$RefReadout/ref_df$UntrtReadout, nDigits_rounding)
      ref_df$RefGRvalue <- calculate_GR_value(ref_df$RefReadout, ref_df$UntrtReadout, day0readout = ref_df$Day0Readout, nDigits_rounding)
      ref_df$DivisionTime <- round(rdata[i, duration_col] / log2(ref_df$UntrtReadout/ref_df$Day0Readout), nDigits_rounding)
      fitting_refs[i, j] <- ref_df
    }
  }
  assays(se)[["Control"]] <- refs
  se
}


####################
# Helper functions
####################

#' Create a control dataframe for a treatment-cell line combination.
#'
#' @param df_
#' @param key
#' @param control_mean_fxn
#' @param out_col_name string of the output readout that will replace \code{CorrectedReadout}.
#'
#'
#' @export
#'
create_control_df <- function(df_, Keys, key, control_mean_fxn, out_col_name) {
  if (length(df_) != 0L) {
    # Rename CorrectedReadout.
    df_ <- df_[, c("CorrectedReadout", intersect(Keys[[key]], colnames(df_)))]
    colnames(df_)[grepl("CorrectedReadout", colnames(df_))] <- out_col_name

    # Aggregate by all non-readout data (the metadata).
    df_ <- aggregate(df_[, out_col_name, drop = FALSE], 
                     by = as.list(df_[, colnames(df_) != out_col_name, drop = FALSE]),
	             function(x) control_mean_fxn(x))
  }
  df_
}


#' @export
#'
extrapolate_references <- function(ref_df, ref_conc) {
  if (length(unique(ref_df$Concentration[!ref_df$masked])) > 3) {
    drc_fit <- drc::drm(
      CorrectedReadout ~ Concentration,
      data = ref_df[!ref_df$masked,],
      fct = drc::LL.4(), # para = c(Hill, x_inf, x0, c50)
      start = c(2, min(ref_df$CorrectedReadout), max(ref_df$CorrectedReadout),
		  median(ref_df$Concentration)),
      lowerl = c(1e-5, min(ref_df$CorrectedReadout) * 0.8, # wide range
		  min(ref_df$CorrectedReadout) * 0.9,
		  min(ref_df$Concentration)/1e3),
      upperl =  c(12, max(ref_df$CorrectedReadout) * 1.1,
		  max(ref_df$CorrectedReadout) * 1.2,
		  max(ref_df$Concentration) * 1e3)
    )
    corrected_readout <- predict(drc_fit, data.frame(Concentration = ref_conc))
  } else {
    corrected_readout <- NA
  }
  corrected_readout
}


#' @export
#'
identify_treatment_references <- function(row_maps_cotrt, trt_rdata, ref_rdata) {
  # Reassess the cases without a match to find equivalent drug and concentration (only 2 drugs).
  # Test if one can use one of the treatment as a reference.
  if ('Gnumber_2' %in% colnames(trt_rdata)) {
    no_matches <- names(row_maps_cotrt)[sapply(row_maps_cotrt, length) == 0L]
    for (rnames in no_matches) {
      # loop through the rows without co-treatment matched
      ref_metadata_idx <- setdiff(intersect(Keys$ref_Endpoint, names(ref_rdata)),
				  c('Gnumber_2', "DrugName_2", 'Concentration_2'))
      names(ref_metadata_idx) <- ref_metadata_idx

      ref_match <- apply(as.matrix(c((
	lapply(ref_metadata_idx, function(y) # matching the metadata
	      unlist(trt_rdata[, y, drop = FALSE] ==
		  (trt_rdata[rnames, y, drop = FALSE]))
	    )), # matching the drugs with mapping from Gnumber to Gnumber_2
	    list(Gnumber = trt_rdata$Gnumber ==
	      trt_rdata[rnames,'Gnumber_2']),
	    list(Gnumber_2 = trt_rdata$Gnumber_2 %in% 
	      gDRutils::get_identifier('untreated_tag')))),
	    2, all)
      if (any(ref_match)) {
	row_maps_cotrt[rnames] <- rownames(normSE)[ref_match]
      }
    }
  }
}


######################################
# Below has not been evaluated yet
######################################

calculate_reference_values_SE_part2 <- function() {
  # Match the reference endpoint with the same co-treatment.
  row_maps_cotrt <- map_df(split_list$treated, split_list$treated, ref_type = "ref_Endpoint")
  row_maps_cotrt <- identify_treatment_references(row_maps_cotrt, split_list$treated, split_list$untreated)

  # reference co-treatment is not always present
  if (i %in% names(row_maps_cotrt) && length(row_maps_cotrt[[i]]) > 0) {
    if (all(row_maps_cotrt[[i]] %in% rownames(ctrlSE))) {
      # get all the co-treatment reference endpoint data
      df_ref <- do.call(rbind, lapply(row_maps_cotrt[[i]], function(x) ctrl_original[x, col_maps[j]][[1]]))
      df_ref <- df_ref[, c("CorrectedReadout", intersect(Keys$ref_Endpoint, colnames(df_ref))), drop = FALSE]
      colnames(df_ref)[grepl("CorrectedReadout", colnames(df_ref))] <- "RefReadout"
      # Aggregate by all non-readout data (the metadata).
      df_ref <- aggregate(df_ref[, "RefReadout", drop = FALSE], 
        by = as.list(df_ref[, which(grepl("RefReadout", colnames(df_ref))), drop = FALSE]), 
        function(x) control_mean_fct(x))

      # check if all control have matching co-treated wells are on the same plate
      df_end <- merge(df_end, df_ref, by = intersect(colnames(df_end), c('Barcode', Keys$discard_keys)), all = TRUE)
      if (any(is.na(df_ref$Barcode))) {
	futile.logger::flog.warn("Control data for the drug are propagated to other plates with co-drug controls. Treatment Id: '%s' Cell_line Id: '%s'", i, j)

	# propagate average values to the other plates
	df_end$UntrtReadout[is.na(df_end$UntrtReadout)] <- mean(df_end$UntrtReadout, na.rm = TRUE) # TODO: Do we want the control_mean_fxn here? 
	df_end$RefReadout[is.na(df_end$RefReadout)] <- mean(df_end$RefReadout, na.rm = TRUE)
      }

      #gladkia: assert for control data
      if (nrow(df_end) == 0L) {
	stop(sprintf("Control dataframe failed. Treatment Id: '%s' Cell_line Id: '%s'", i, j))
      }
    } else if (all(row_maps_cotrt[[i]] %in% rownames(se))) {
      # case of the reference being with Gnumber == Gnumber_2
      ref_conc <- rdata[i, 'Concentration_2']
      df_ref <- do.call(rbind,
	lapply(row_maps_cotrt[[i]], function(x) {
	  if (any(se_original[x, col_maps[j]][[1]]$Concentration == ref_conc)) {
	    # the reference value with same concentration is found
	    se_original[x, col_maps[j]][
	      se_original[x, col_maps[j]][[1]]$Concentration == ref_conc, ]
	  } else {
	    # Infer the reference with proper concentration through a fit.
	    corrected_readout <- extrapolate_references(ref_df, ref_conc)
	    df_ref <- data.frame(Concentration = ref_conc, CorrectedReadout = corrected_readout)
	  }
	})
      )
      
      df_ref <- df_ref[, c("CorrectedReadout", intersect(Keys$ref_Endpoint, colnames(df_ref))), drop = FALSE]
      colnames(df_ref)[grepl("CorrectedReadout", colnames(df_ref))] <- "RefReadout"
      df_ref <- aggregate(df_ref[, "RefReadout", drop = FALSE], 
	by = as.list(df_ref[, which(grepl("RefReadout", colnames(df_ref))), drop = FALSE]), 
	function(x) control_mean_fct(x))

      # check if all control have matching co-treated wells are on the same plate
      if (all(df_end$Barcode %in% df_ref$Barcode) && all(df_ref$Barcode %in% df_end$Barcode)) {
	df_end <- merge(df_end, df_ref, by = intersect(colnames(df_end), c('Barcode', Keys$discard_keys)))
      } else {
	futile.logger::flog.warn("Control data for the drug are propagated to other plates with co-drug controls. Treatment Id: '%s' Cell_line Id: '%s'", i, j)
	# propagate average values to the other plates
	df_end <- merge(df_end, df_ref, by = "Barcode", all = TRUE)
	
	df_end$UntrtReadout[is.na(df_end$UntrtReadout)] <- mean(df_end$UntrtReadout, na.rm = TRUE)
	df_end$RefReadout[is.na(df_end$RefReadout)] <- mean(df_end$RefReadout, na.rm = TRUE)
      }
    } else {
      stop(sprintf("Reference failed. Treatment Id: '%s' Cell_line Id: '%s'", i, j))
    }
  } else if (i %in% names(row_maps_cotrt) && length(row_maps_cotrt[[i]]) == 0L) {
    futile.logger::flog.warn("No reference condition found for Treatment Id: '%s' Cell_line Id: '%s'", i, j)
    df_end$RefReadout <- NA
  } else {
    df_end$RefReadout <- df_end$UntrtReadout
  }

  df_end <- df_end[, c("CorrectedReadout", intersect(Keys$untrt_Endpoint, colnames(df_end))), drop = FALSE]
  colnames(df_end)[1] <- "UntrtReadout"
  if (ncol(df_end) > 1L) {
    # I think one example of this would be the 'Replicate' field, but I'm not sure what the other cases would be.
    df_end <- aggregate(df_end[, 1, drop = FALSE],
      by = as.list(df_end[, -1, drop = FALSE]),
      function(x) control_mean_fct(x))
  } else { # There are no untrt_Endpoint keys in the data.frame. 
    df_end <- DataFrame(UntrtReadout = control_mean_fct(df_end$UntrtReadout))
  }

  # TODO: I think this is basically just trying to figure out if any overriding key_values exist. 
  # Come back to this.
  if (!is.null(key_values) & length(key_values) > 0) {
    for (i in which(names(key_values) %in% names(ctrl_rdata))) {
      if (is.numeric(key_values[i])) {
	row_endpoint_value_filter <- row_endpoint_value_filter &
	    (ctrl_rdata[, names(key_values)[i]] == key_values[i] &
		    !is.na(ctrl_rdata[, names(key_values)[i]]))
      } else {
	row_endpoint_value_filter <- row_endpoint_value_filter &
	    (ctrl_rdata[, names(key_values)[i] ] %in% key_values[i])
      }
    }
  }

  # Match the columns. 
  # TODO: Wait a minute. Do we need to do the splitting for hte columns as well? 
  # mapping for columns; 1 to 1 unless overridden by key_values
#  col_maps <- array(colnames(ctrlSE), dimnames = list(colnames(normSE)))
#  if (any(names(key_values) %in% names(norm_cdata))) {
#    col_maps[] <- colnames(ctrlSE)[
#      which(key_values[names(key_values) %in% names(norm_cdata)] ==
#	  SummarizedExperiment::colData(ctrlSE)[, names(SummarizedExperiment::colData(ctrlSE)) %in% names(key_values)])]
#  }


  #gladkia: assert for merged study/control data
  ctrl_bcodes <- sort(unique(df_ctrl$Barcode))
  trt_bcodes <- sort(unique(se_original[i, j][[1]]$Barcode))
  # check if all treated values have matching controls on the same plate
  if (!all(trt_bcodes %in% ctrl_bcodes)) {
    # if not, propagate to all plates
    futile.logger::flog.warn("Control data are averaged and propagated to treatment plates. Treatment Id: '%s' (plates '%s') Control plates: '%s'", i, paste(trt_bcodes, collapse = ", "), paste(ctrl_bcodes, collapse = ", ")
      )
    data.table::setDF(data.table::rbindlist(list(df_ctrl, cbind(data.frame(Barcode = setdiff(trt_bcodes, ctrl_bcodes)),
      t(colMeans(df_ctrl[, setdiff(colnames(df_ctrl), "Barcode")])))), fill = TRUE))
  }

  # works with by = character(0) but changes the order of rows
  df_merged <- merge(data.frame(se_original[i, j][[1]]), data.frame(df_ctrl),
    by = intersect(colnames(df_ctrl), c('Barcode', Keys$discard_keys)),
    all.x = TRUE)

  ctrl_original <- SummarizedExperiment::assay(gDR::aapply(ctrlSE, function(x) x[!x$masked,]))
  # need to keep original data for the case in which reference is such that Gnumber == Gnumber_2
}
