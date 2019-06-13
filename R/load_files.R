
# openxlsx skip the first emprty rows and cannot be overridden
#' @import readxl 
#' @import readr
#' @import stringr

#' @export
load_data = function(manifest_file, df_template_files, results_file, log_str) {

    log_str = c(log_str, '', 'load_merge_data')

    if (is.data.frame(df_template_files)) {# for the shiny app
        template_file = df_template_files$datapath
        template_filename = df_template_files$name
    } else {
        template_filename = template_file
    }

    manifest = load_manifest(manifest_file, log_str)
    treatments = load_templates(df_template_files, log_str)
    data = load_results(results_file, log_str)

    # check the all template files are available
    if (!all(unique(manifest$Template[manifest$Barcode %in% data$Barcode])
                    %in% basename(template_filename))) {
        ErrorMsg = paste('Some template files are missing:',
                    paste(setdiff(unique(manifest$Template[manifest$Barcode %in% data$Barcode]),
                                     basename(template_filename)), collapse = ' ; '))
        log_str = c(log_str, 'Error in load_merge_data:')
        log_str = c(log_str, ErrorMsg)
        writeLines(log_str)
        stop(ErrorMsg)
    }
    return( list(manifest = manifest, treatments = treatments, data = data) )
}



#' @export
load_manifest = function (manifest_file, log_str) {
    # manifest_file is a string or a vector of strings

    log_str = c(log_str, '', 'load_manifest')

    # read files
    manifest_data = lapply(manifest_file, function(x) {
        df = tryCatch( { read_excel(x, col_names = T) },
            error = function(e) { read_tsv(x, col_names=T, skip_empty_rows=T) }
            )
    } )

    # check default headers are in each df
    dump = sapply(lapply(1:length(manifest_file),
            function(x) c(manifest_file[x], manifest_data[x])),
                function(x) check_metadata_names(colnames(x[[2]]), log_str,
                    df_name=x[[1]], df_type='manifest'))

    cat_manifest_data = bind_rows(manifest_data)
    colnames(cat_manifest_data) = check_metadata_names(colnames(cat_manifest_data),
                    log_str, 'manifest')
    # check that barcodes are unique
    stopifnot(dim(cat_manifest_data)[1] == length(unique(cat_manifest_data$Barcode)))
    # add error message

    cat_manifest_data$Template = basename(cat_manifest_data$Template)

    print('Manifest loaded:')
    print(dim(cat_manifest_data))
    return(cat_manifest_data)
}



#' @export
load_templates = function (df_template_files, log_str) {
    # template_file is a string or a vector of strings
    log_str = c(log_str, '', 'load_templates')

    if (is.data.frame(df_template_files)) {# for the shiny app
        template_file = df_template_files$datapath
        template_filename = df_template_files$name
    } else {
        template_file = df_template_files
        template_filename = basename(template_file)
    }

    # read sheets in files
    template_sheets = lapply(template_file, excel_sheets)
    # check Gnumber is present in each df
    sapply(lapply(1:length(template_file), function(x) c(template_file[x], template_sheets[x])),
        function(x) check_metadata_names(x[[2]], log_str, df_name=x[[1]], df_type='template'))

    all_templates = data.frame()
    for (iF in 1:length(template_file)) {
        print(paste('Loading', template_sheets[[iF]]))
        # first check that the sheet names are ok
        # identify Gnumber sheet (case insensitive)
        Gnumber_idx = grep('Gnumber$', template_sheets[[iF]], ignore.case = T)
        # case of untreated plate
        if (length(template_sheets[[iF]])==1) {
            if(length(Gnumber_idx)==0 || Gnumber_idx!=1) {
                ErrorMsg = paste('In untreated template file', template_file[[iF]],
                    ', sheet name must be Gnumber')
                log_str = c(log_str, 'Error in load_templates:')
                log_str = c(log_str, ErrorMsg)
                writeLines(log_str)
                stop(ErrorMsg)
            }
            df = read_excel(template_file[[iF]], sheet = template_sheets[[iF]][1],
                    col_names = paste0('x', 1:48), range = 'A1:AV32')
            if ( !(all(toupper(unlist(df)[!is.na(unlist(df))]) %in% c('UNTREATED', 'VEHICLE')))) {
                    ErrorMsg = paste('In untreated template file', template_file[[iF]],
                        ', entries mush be Vehicle or Untreated')
                    log_str = c(log_str, 'Error in load_templates:')
                    log_str = c(log_str, ErrorMsg)
                    writeLines(log_str)
                    stop(ErrorMsg)
            }
        } else {
        # normal case
            check_metadata_names(template_sheets[[iF]], log_str,
                    df_name=template_filename[iF],
                    df_type='template_treatment')
        }
        # read the different sheets and check for plate size
        # enforce range to avoid skipping empty rows at the beginning
        df = read_excel(template_file[[iF]], sheet = template_sheets[[iF]][Gnumber_idx],
                col_names = paste0('x', 1:48), range = 'A1:AV32', col_types="text")
        # get the plate size
        n_row = 2**ceiling(log2(max(which(apply(!is.na(df), 1, any)))))
        n_col = 1.5*2**ceiling(log2(max(which(apply(!is.na(df), 2, any)))/1.5))
        n_row = max(n_row, n_col/1.5)
        n_col = max(1.5*n_row, n_col)
        plate_range = ifelse(n_col<26, paste0('A1:', LETTERS[n_col], n_row), 'A1:AV32')

        # need to adapt for 1536 well plates
        df_template = expand.grid(WellRow = LETTERS[1:n_row], WellColumn = 1:n_col)

        for (iS in template_sheets[[iF]]) {
            df = as.data.frame(read_excel(template_file[[iF]], sheet = iS,
                            col_names = paste0('x', 1:n_col), range = plate_range))
            df$WellRow = LETTERS[1:n_row]
            df_melted = reshape2::melt(df, id.vars='WellRow')
            colnames(df_melted)[3] = iS
            colnames(df_melted)[colnames(df_melted) == 'variable'] = 'WellColumn'
            df_melted$WellColumn = gsub('x', '', df_melted$WellColumn)
            df_template = merge(df_template, df_melted, by=c('WellRow', 'WellColumn'))
        }
        df_template$Template = template_filename[iF]
        colnames(df_template) = check_metadata_names(colnames(df_template), log_str,
                            df_name=template_filename[iF])
        all_templates = bind_rows(all_templates, df_template)

    }
    print('Templates loaded:')
    print(dim(all_templates))
    return(all_templates)
}



#' @export
load_results = function(df_results_files, log_str) {
    # results_file is a string or a vector of strings
    log_str = c(log_str, '', 'load_results')

    if (is.data.frame(df_results_files)) {# for the shiny app
        results_file = df_results_files$datapath
        results_filename = df_results_files$name
    } else {
        results_file = df_results_files
        results_filename = basename(results_file)
    }

    stopifnot(sapply(results_file, file.exists))

    # test if the result files are .tsv or .xls(x) files
    isExcel = sapply(results_file, function(x) tryCatch(
        { excel_sheets(x)
            return(TRUE) },
            error = function(x) { return(FALSE) }
    ))

    # read sheets in files; warning if more than one sheet (unexpected but can be handled)
    results_sheets = vector('list', length(results_file))
    results_sheets[!isExcel] = 0
    results_sheets[isExcel] = lapply(results_file[isExcel], excel_sheets)
    if (any(lapply(results_sheets, length)>1)) {
        WarnMsg = paste('multiple sheets in result file:',
                results_file[lapply(results_sheets, length)>1])
        log_str = c(log_str, paste('Warning in ', match.call()[[1]]))
        log_str = c(log_str, WarnMsg)
        warning(WarnMsg)
    }

    # read all files and sheets
    all_results = data.frame()
    for (iF in 1:length(results_file)) {
        for (iS in results_sheets[[iF]]) {
            print(paste('Reading file', results_file[[iF]], '; sheet', iS))
            if (iS == 0) {
                df = read_tsv(results_file[[iF]], col_names=F, skip_empty_rows=T)
                # skip_empty_rows flag needs to be TRUE even if it ends up not skipping empty rows
                if (dim(df)[2] == 1) { # likely a csv file
                    df = read_csv(results_file[[iF]], col_names=F, skip_empty_rows=T)
                }
            } else { # expect an Excel spreadsheet
                if (length(results_sheets[[iF]])>1) {# if multiple sheets, assume 1 plate per sheet
                    df = read_excel(results_file[[iF]], sheet = iS, #col_names = F,
                            col_names = paste0('x', 1:48), range = 'A1:AV32')
                } else {
                    df = read_excel(results_file[[iF]], sheet = iS, col_names = F)
                    colnames(df) = col_names = paste0('x', 1:dim(df)[2])
                }
                df = df[, !apply(df[1:48,],2,function(x) all(is.na(x)))] # remove extra columns
                                        # limit to first 48 rows in case Protocol information is
                                        # exported which generate craps at the end of the file
            }
            full_rows = !apply(df,1,function(x) all(is.na(x))) # not empty rows
            # if big gap, delete what is at the bottom (Protocol information)
            gaps = min(which(full_rows)[ (diff(which(full_rows))>20) ]+1, dim(df)[1])
            df = df[ which(full_rows)[which(full_rows) <= gaps], ] # remove extra rows
            df = df[ , !apply(df,2,function(x) all(is.na(x)))] # remove empty columns
            # get the plate size
            n_col = 1.5*2**ceiling(log2(dim(df)[2]/1.5))
            n_row = n_col/1.5

            # get the barcode(s) in the sheet; expected in column C (third one)
            Barcode_idx = which(as.data.frame(df)[,3] %in% 'Barcode')
            # run through all plates
            for (iB in Barcode_idx) {
                # two type of format depending on where Background information is placed
                if (any(as.data.frame(df)[iB+(1:4),1] %in% 'Background information')) {
                    ref_bckgrd = which(as.data.frame(df)[iB+(1:4),1] %in% 'Background information')
                    readout_offset = 2 + ref_bckgrd
                    stopifnot(as.character(df[iB+ref_bckgrd+1,4]) %in% 'Signal')
                    BackgroundValue = as.numeric(df[iB+ref_bckgrd+2,4])
                } else {
                    # export without background information
                    # case of " Exported with EnVision Workstation version 1.13.3009.1409 "
                    readout_offset = 1
                    BackgroundValue = 0
                }

                # check the structure of file is ok
                check_values = as.matrix(df[iB+readout_offset+c(0,1, n_row, n_row+1), n_col])
                if (any(c(is.na(check_values[2:3]), !is.na(check_values[c(1,4)])))) {
                    ErrorMsg = paste('In result file', results_filename[[iF]], '(sheet', iS,
                        ') readout values are misplaced for plate', as.character(df[iB+1,3]))
                    log_str = c(log_str, 'Error in load_results:')
                    log_str = c(log_str, ErrorMsg)
                    writeLines(log_str)
                    stop(ErrorMsg)
                }

                readout = as.matrix(df[iB+readout_offset+(1:n_row),1:n_col])

                # check that the plate size is consistent and contains values
                if (any(is.na(readout))) {
                    ErrorMsg = paste('In result file', results_filename[[iF]], '(sheet', iS,
                        ') readout values are missing for plate', as.character(df[iB+1,3]))
                    log_str = c(log_str, 'Error in load_results:')
                    log_str = c(log_str, ErrorMsg)
                    writeLines(log_str)
                    stop(ErrorMsg)
                }

                df_results = data.frame(
                    Barcode = as.character(df[iB+1,3]),
                    WellRow = LETTERS[1:n_row],
                    WellColumn = as.vector(t(matrix(1:n_col,n_col,n_row))),
                    ReadoutValue = as.numeric(as.vector(readout)),
                    BackgroundValue = BackgroundValue
                )
                print(paste('Plate', as.character(df[iB+1,3]),
                        'read;', dim(df_results)[1], 'wells'))
                all_results = rbind(all_results, df_results)
            }
            print('File done')
        }
    }
    return(all_results)
}



#' @export
check_metadata_names = function(col_df, log_str, df_name = '', df_type = NULL) {

    log_str = c(log_str, '   check_metadata_names')
    # first check for required column names
    if (!is.null(df_type)) {
        if (df_type == 'manifest') {
                expected_headers = c('Barcode', 'Time', 'Template')
        }
        if (df_type == 'template') {
                expected_headers = c('Gnumber')
        }
        if (df_type == 'template_treatment') {
                expected_headers = c('Gnumber', 'Concentration')
        }

        headersOK = toupper(expected_headers) %in% toupper(col_df)
        if (any(!headersOK)) {
            ErrorMsg = paste(df_name,
                'does not contains all expected headers for a', df_type, '; ',
                paste(expected_headers[ !(expected_headers %in% col_df) ], collpase = ' ; '),
                ' required')
            log_str = c(log_str, 'Error in check_metadata_names:')
            log_str = c(log_str, ErrorMsg)
            writeLines(log_str)
            stop(ErrorMsg)
        }
        if (df_type == 'template_treatment') {
            # assess if multiple drugs and proper pairing
            n_drug = agrep('Gnumber', col_df, ignore.case = T)
            n_conc = agrep('Concentration', col_df, ignore.case = T)
            if (length(n_drug) != length(n_conc)) {
                ErrorMsg = paste('Treatment template', df_name,
                    'does not contains the same number of Gnumber_* and Concentration_* sheets')
                log_str = c(log_str, 'Error in check_metadata_names:')
                log_str = c(log_str, ErrorMsg)
                writeLines(log_str)
                stop(ErrorMsg)
            }
            if (length(n_drug)>1) {
                trt_sheets = c(paste0('Gnumber_', 2:length(n_drug)),
                    paste0('Concentration_', 2:length(n_conc)))
                if (!(all(toupper(trt_sheets) %in% toupper(col_df)))) {
                    ErrorMsg = paste('Treatment template', df_name,
                        'does not contains: ',
                        paste(trt_sheets[!(toupper(trt_sheets) %in% toupper(col_df))],
                            collapse = ' ; '))
                    log_str = c(log_str, 'Error in check_metadata_names:')
                    log_str = c(log_str, ErrorMsg)
                    writeLines(log_str)
                    stop(ErrorMsg)
                }
            }
        }
    }

    corrected_names = col_df

    # remove spaces and convert to WordUppercase
    names_spaces = regexpr('\\s', corrected_names)>0
    if (any(names_spaces)) {
        for (i in which(names_spaces)) {
            s <- strsplit(corrected_names[i], " ")[[1]]
            corrected_names[i] = paste(toupper(substring(s, 1, 1)), substring(s, 2),
              sep = "", collapse = "")
        }

        WarnMsg = paste('Metadata field names for', df_name,
            'cannot contain spaces --> corrected to: ',
                paste(corrected_names[names_spaces], collapse = ' ; '))
        log_str = c(log_str, 'Warning in check_metadata_names:')
        log_str = c(log_str, WarnMsg)

        warning(WarnMsg)
    }

    # check for wrong metadata field names (including dash, starting with number, ... )
    bad_names = regexpr('\\W', corrected_names)>0 | regexpr('\\d', corrected_names)==1
    if (any(bad_names)) {
        ErrorMsg = paste('Metadata field names for', df_name,
            'cannot contain special characters or start with a number: ',
                paste(corrected_names[bad_names], collapse = ' ; '))
        log_str = c(log_str, 'Error in check_metadata_names:')
        log_str = c(log_str, ErrorMsg)
        writeLines(log_str)
        stop(ErrorMsg)
    }

    # common headers that are written in a specific way
    # throw warning if close match and correct upper/lower case for consistency
    controlled_headers=c('CLID', 'Time', 'Barcode', 'Template', 'Gnumber', 'Concentration',
                paste0('Gnumber_', 2:10), paste0('Concentration_', 2:10))
    for (i in 1:length(controlled_headers)) {
        case_match = setdiff(
            grep(paste0(controlled_headers[i],'$'), corrected_names, ignore.case = T),
                        grep(paste0(controlled_headers[i],'$'), corrected_names))
        if (length(case_match)>0){
            WarnMsg = paste('Header', corrected_names[case_match], 'in', df_name,
                                    'corrected to', controlled_headers[i])
            corrected_names[case_match] = controlled_headers[i]
            log_str = c(log_str, 'Warning in check_metadata_names:')
            log_str = c(log_str, WarnMsg)
            warning(WarnMsg)
        }

        # more noise than actually useful
        # fuzzy_match = setdiff(agrep(controlled_headers[i], corrected_names),
        #                         grep(controlled_headers[i], corrected_names))
        # if (length(fuzzy_match)>0){
        #     WarnMsg = paste('Header', corrected_names[fuzzy_match], 'in', df_name,
        #                     'looks similar to', controlled_headers[i], '; Please check for typo')
        #     log_str = c(log_str, 'Warning in check_metadata_names:')
        #     log_str = c(log_str, WarnMsg)
        #     warning(WarnMsg)
        # }
    }

    # check for headers that are reserved for downstream analyses
    ReservedHeaders = c('CellLineName', 'Tissue', 'ReferenceDivisionTime', 'DrugName',
                    paste0('DrugName_', 2:10), "ReadoutValue", "BackgroundValue",
                    "CorrectedReadout", "Day0Readout", 'GRvalue', 'RelativeViability',
                    'DivisionTime')
    if (any(corrected_names %in% ReservedHeaders)) {
        ErrorMsg = paste('Metadata field name: ',
            paste(intersect(ReservedHeaders, corrected_names), collapse = ' ; '),
            ' in', df_name, 'is not valid (reserved for output)')
        log_str = c(log_str, 'Error in check_metadata_names:')
        log_str = c(log_str, ErrorMsg)
        writeLines(log_str)
        stop(ErrorMsg)
    }

    return(corrected_names)
}
