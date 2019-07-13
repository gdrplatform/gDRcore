source("analyze_data.R") # to get the function identify_keys

# DB structure is described in https://drive.google.com/open?id=1gX5ja_dSdygr2KYTqYUiENKWxu9HkOEz

# mimic the mySQL db
gdr_projects = data.frame()
condition_data = data.frame()
condition_codrug = data.frame()
condition_additional_treatment = data.frame()

project = 1

if (project == 1) {

    project_raw_files = data.frame(project_number = project,
                file_type = c('manifest', 'template', 'template', 'results', 'results'),
                file_uri = c('test_Kyle1/Manifest file_0077&Palbo.tsv',
                            'test_Kyle1/Template_BT474_untreated.xlsx',
                            'test_Kyle1/Template_BT474_7daytreated_KE.xlsx',
                            'test_Kyle1/cyquant_BT474_day0.xlsx',
                            'test_Kyle1/cyquant_BT474_day7.xlsx'))

    project_processed_files = data.frame(project_number = project,
                file_type = c('normalized', 'averaged', 'metrics'),
                file_uri = c('test_Kyle1/ref_df_normalized.tsv',
                            'test_Kyle1/ref_df_averaged.tsv',
                            'test_Kyle1/ref_df_metrics.tsv'))

    df_normalized = read.table(paste0('../inst/testdata/',project_processed_files[1,3]), header = T)
    df_averaged = read.table(paste0('../inst/testdata/', project_processed_files[2,3]), header = T)
    df_metrics = read.table(paste0('../inst/testdata/', project_processed_files[3,3]), header = T)

    project_data = c(project_number = project,
                        username = 'kyle',
                        description = 'Kyle1',
                        approved = NULL,
                        date_processed = date())
}

# need to update:
#  - CLID --> clid
#  - add DivisionTime in the df_metrics table

keys = identify_keys(df_metrics)$DoseResp
condition_keys = setdiff(keys, c("DrugName", 'Concentration', 'Gnumber', "CellLineName", "Tissue",
                                keys[grep('DrugName', keys)]))
sub_condition_data = df_metrics[, condition_keys]   # , 'DivisionTime')]
stopifnot(dim(sub_condition_data)[1] == dim(unique(sub_condition_data))[1])
sub_condition_data$project_number = project
sub_condition_data$condition_number = dim(condition_data)[1] + (1:dim(sub_condition_data)[1])

# get the Gnumbers in their own table if any
N_add_Drugs = condition_keys[grep('Gnumber', condition_keys)]
sub_condition_codrug = data.frame()
for (d in N_add_Drugs) {
    codrug = sub_condition_data[, c('condition_number', d,
                gsub('Gnumber','Concentration', d))]
    colnames(codrug)[-1] = c('Gnumber','Concentration')
    sub_condition_codrug = rbind(sub_condition_codrug, codrug)
}
condition_codrug = rbind(condition_codrug, sub_condition_codrug)

# get the other condition treatments
add_treatments = setdiff(condition_keys[ c(-agrep('Gnumber', condition_keys),
                                    -agrep('Concentration', condition_keys))],
                                c('Time', 'CLID'))
sub_condition_additional_treatment = data.frame()
for (d in add_treatments) {
    cotrt = cbind(sub_condition_data[, c('condition_number', d)], metadata_field = d)
    colnames(cotrt)[2] = 'metadata_value'
    sub_condition_additional_treatment = rbind(sub_condition_additional_treatment, cotrt)
}
condition_additional_treatment = rbind(condition_additional_treatment,
                                    sub_condition_additional_treatment)


df_metrics[,keys]
