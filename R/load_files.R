library(readxl) # openxlsx skip the first emprty rows and cannot be overridden
library(readr)
library(stringr)

load_manifest = function (manifest_file, log_file) {
    # manifest_file is a string or a vector of strings
    # log_file is an open file to sink errors and warnings

    # read files
    manifest_data = lapply(manifest_file, function(x) {
        df = tryCatch( { read_excel(x, col_names = T) },
            error = function(e) { read_tsv(x, col_names=T, skip_empty_rows=T) }
            )
    } )

    # check default headers are in each df
    dump = sapply(lapply(1:length(manifest_file),
            function(x) c(manifest_file[x], manifest_data[x])),
                function(x) check_metadata_names(colnames(x[[2]]), log_file,
                    df_name=x[[1]], df_type='manifest'))

    cat_manifest_data = bind_rows(manifest_data)
    colnames(cat_manifest_data) = check_metadata_names(colnames(cat_manifest_data),
                    log_file, 'manifest')
    # check that barcodes are unique
    stopifnot(dim(cat_manifest_data)[1] == length(unique(cat_manifest_data$Barcode)))
    cat_manifest_data$Template = basename(cat_manifest_data$Template)

    return(cat_manifest_data)
}


load_templates = function (template_file, log_file) {
    # template_file is a string or a vector of strings
    # log_file is an open file to sink errors and warnings

    # read sheets in files
    template_sheets = lapply(template_file, excel_sheets)
    # check Gnumber is present in each df
    sapply(lapply(1:2, function(x) c(template_file[x], template_sheets[x])),
        function(x) check_metadata_names(x[[2]], log_file, df_name=x[[1]], df_type='template'))

    all_templates = data.frame()
    for (iF in 1:length(template_file)) {
        # first check that the sheet names are ok
        # case of untreated plate
        if (length(template_sheets[[iF]])==1) {
            stopifnot(template_sheets[[iF]] == 'Gnumber')
            df = read_excel(template_file[[iF]], sheet = 'Gnumber',
                    col_names = paste0('x', 1:48), range = 'A1:AV32')
            stopifnot(all(toupper(unlist(df)[!is.na(unlist(df))]) %in% c('UNTREATED', 'VEHICLE')))
        } else {
        # normal case
            check_metadata_names(template_sheets[[iF]], log_file, df_name=template_file[iF],
                    df_type='template_treatment')
        }

        # read the different sheets and check for plate size
        # enforce range to avoid skipping empty rows at the beginning
        df = read_excel(template_file[[iF]], sheet = 'Gnumber',
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
            df = read_excel(template_file[[iF]], sheet = iS,
                            col_names = paste0('x', 1:n_col), range = plate_range)
            df$WellRow = LETTERS[1:n_row]
            df_melted = melt(df, id='WellRow', value.name = iS)
            colnames(df_melted)[colnames(df_melted) == 'variable'] = 'WellColumn'
            df_melted$WellColumn = gsub('x', '', df_melted$WellColumn)
            df_template = merge(df_template, df_melted, by=c('WellRow', 'WellColumn'))
        }
        df_template$Template =  basename(template_file[[iF]])
        all_templates = bind_rows(all_templates, df_template)
    }
    return(all_templates)
}

load_results = function (results_file, log_file) {
    # results_file is a string or a vector of strings
    # log_file is an open file to sink errors and warnings

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
        writeLines(paste('Warning in ', match.call()[[1]]), log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }

    # read all files and sheets
    all_results = data.frame()
    for (iF in 1:length(results_file)) {
        for (iS in results_sheets[[iF]]) {
            if (iS == 0) {
                df = read_tsv(results_file[[iF]], col_names=F, skip_empty_rows=T)
                # skip_empty_rows flag needs to be TRUE even if it ends up not skipping empty rows
            } else { # expect an Excel spreadsheet
                df = read_excel(results_file[[iF]], sheet = iS, #col_names = F,
                        col_names = paste0('x', 1:48), range = 'A1:AV32')
                df = df[, !apply(df,2,function(x) all(is.na(x)))]
            }
            # get the plate size
            n_col = 1.5*2**ceiling(log2(dim(df)[2]/1.5))
            n_row = n_col/1.5

            # get the barcode(s) in the sheet; expected in column C (third one)
            Barcode_idx = which(as.data.frame(df)[,3] %in% 'Barcode')
            # run through all plates
            for (iB in Barcode_idx) {
                # two type of format depending on where Background information is placed
                if (df[iB+3,1] %in% 'Background information') {
                    readout_offset = 6
                    stopifnot(as.character(df[iB+4,4]) %in% 'Signal')
                    BackgroundValue = as.numeric(df[iB+5,4])
                } else {
                    # export without background information
                    # case of " Exported with EnVision Workstation version 1.13.3009.1409 "
                    readout_offset = 2
                    BackgroundValue = 0
                }

                readout = as.matrix(df[iB+readout_offset+(1:n_row),1:n_col])

                # check that the plate size is consistent and contains values
                stopifnot(all(!is.na(readout)))

                df_results = data.frame(
                    Barcode = as.character(df[iB+1,3]),
                    WellRow = LETTERS[1:n_row],
                    WellColumn = as.vector(t(matrix(1:n_col,n_col,n_row))),
                    ReadoutValue = as.numeric(as.vector(readout)),
                    BackgroundValue = BackgroundValue
                )
                all_results = rbind(all_results, df_results)
            }
        }
    }
    return(all_results)
}


check_metadata_names = function(col_df, log_file, df_name = '', df_type = NULL) {

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

        headersOK = expected_headers %in% col_df
        if (any(!headersOK)) {
            ErrorMsg = paste(df_name,
                'does not contains all expected headers for a', df_type, '; ',
                paste(expected_headers[ !(expected_headers %in% col_df) ], collpase = ' ; '),
                ' required')
            writeLines('Error in check_metadata_names:', log_file)
            writeLines(ErrorMsg, log_file)
            close(log_file)
            stop(ErrorMsg)
        }
        if (df_type == 'template_treatment') {
            # assess if multiple drugs and proper pairing
            n_drug = agrep('Gnumber', col_df)
            n_conc = agrep('Concentration', col_df)
            if (length(n_drug) != length(n_conc)) {
                ErrorMsg = paste('Treatment template', df_name,
                    'does not contains the same number of Gnumber_* and Concentration_* sheets')
                writeLines('Error in check_metadata_names:', log_file)
                writeLines(ErrorMsg, log_file)
                close(log_file)
                stop(ErrorMsg)
            }
            if (length(n_drug)>1) {
                trt_sheets = c(paste0('Gnumber_', 2:length(n_drug)),
                    paste0('Concentration_', 2:length(n_conc)))
                if (!(all(trt_sheets %in% col_df))) {
                    ErrorMsg = paste('Treatment template', df_name,
                        'does not contains: ',
                        paste(trt_sheets[!(trt_sheets %in% col_df)], collapse = ' ; '))
                    writeLines('Error in check_metadata_names:', log_file)
                    writeLines(ErrorMsg, log_file)
                    close(log_file)
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
        writeLines('Warning in check_metadata_names:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }

    # check for wrong metadata field names (including dash, starting with number, ... )
    bad_names = regexpr('\\W', corrected_names)>0 | regexpr('\\d', corrected_names)==1
    if (any(bad_names)) {
        ErrorMsg = paste('Metadata field names for', df_name,
            'cannot contain special characters or start with a number: ',
                paste(corrected_names[bad_names], collapse = ' ; '))
        writeLines('Error in check_metadata_names:', log_file)
        writeLines(ErrorMsg, log_file)
        close(log_file)
        stop(ErrorMsg)
    }

    # common headers that are written in a specific way
    # throw warning if close match and correct upper/lower case for consistency
    controlled_headers = c('CLid', 'Media', 'Ligand')
    for (i in 1:length(controlled_headers)) {
        case_match = setdiff(grep(controlled_headers[i], corrected_names, ignore.case = T),
                                grep(controlled_headers[i], corrected_names))
        if (length(case_match)>0){
            corrected_names[case_match] = controlled_headers[i]
            WarnMsg = paste('Header', corrected_names[case_match], 'in', df_name,
                                    'corrected to', controlled_headers[i])
            writeLines('Warning in check_metadata_names:', log_file)
            writeLines(WarnMsg, log_file)
            warning(WarnMsg)
        }

        fuzzy_match = setdiff(agrep(controlled_headers[i], corrected_names),
                                grep(controlled_headers[i], corrected_names))
        if (length(fuzzy_match)>0){
            WarnMsg = paste('Header', corrected_names[fuzzy_match], 'in', df_name,
                            'looks similar to', controlled_headers[i], '; Please check for typo')
            writeLines('Warning in check_metadata_names:', log_file)
            writeLines(WarnMsg, log_file)
            warning(WarnMsg)
        }
    }

    # check for headers that are reserved for downstream analyses
    ReservedHeaders = c('CellLineName', 'Tissue', 'ReferenceDivisionTime', 'DrugName',
                    paste0('DrugName_', 2:10), "ReadoutValue", "BackgroundValue",
                    "CorrectedReadout", "Day0Readout", 'GRvalue', 'RelViability', 'DivisionTime')
    if (any(corrected_names %in% ReservedHeaders)) {
        ErrorMsg = paste('Metadata field name: ',
            paste(intersect(ReservedHeaders, corrected_names), collapse = ' ; '),
            ' in', df_name, 'is not valid (reserved for output)')
        writeLines('Error in check_metadata_names:', log_file)
        writeLines(ErrorMsg, log_file)
        close(log_file)
        stop(ErrorMsg)
    }

    return(corrected_names)
}
