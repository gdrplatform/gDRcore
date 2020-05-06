devtools::load_all('../')

# drug combination: 2 x 3 drugs (9 doses by 1 doses + 0) in triplicates + single agent for all
df_normalized = data.frame( CellLineName = as.vector(t(matrix(c('COV318', 'HCC2218'), 2, 270))),
                    Tissue = as.vector(t(matrix(c('Ovary', 'Breast'), 6, 90))),
                    Duration = 72,
                    clid = as.vector(t(matrix(c('CL586056', 'CL131828'), 6, 90))),
                    Barcode = as.vector(t(matrix(paste0('Plate', 1:6), 6, 90))),
                    Gnumber = as.vector(t(matrix(
                      c('G02001876', 'G02442104', 'G00505032', 'G02967907', 'G00033829',
                        'G02001876', 'G02442104', 'G02001876', 'G02442104'), 54, 10))),
                    Gnumber_2 = as.vector(t(matrix(
                      c('Vehicle', 'Vehicle', 'Vehicle', 'Vehicle', 'Vehicle',
                      'G00505032', 'G00033829','G00033829', 'G00505032'), 54, 10))),
                    Concentration = as.vector(matrix(c(0, 10**(seq(-2.5,1,.5))), 540, 1)),
                    Concentration_2 = 0,
                    RelativeViability = 1,
                    DivisionTime = as.vector(t(matrix(c(58,90), 2, 270))),
                    ReferenceDivisionTime = as.vector(t(matrix(c(60,85), 2, 270))),
                    BackgroundValue = 1)

df_normalized$Concentration_2[ df_normalized$Gnumber_2 != 'Vehicle'] = 1

# apply the drug effect
ec50 = matrix(c(.08, 2, .5, .2, .03, 100,  .2, 2, .02, 3, .5, 100),6,2)
colnames(ec50) = c('COV318', 'HCC2218')
rownames(ec50) = c('G02001876', 'G02442104', 'G00505032', 'G02967907', 'G00033829', 'Vehicle')

e_max = matrix(c(.62, .9, .5, .6, .75, 1,  .85, 1, .7, .3, .8, 1),6,2)
colnames(e_max) = c('COV318', 'HCC2218')
rownames(e_max) = c('G02001876', 'G02442104','G00505032', 'G02967907', 'G00033829', 'Vehicle')

df_normalized$ReadoutValue = 100 *
    apply( df_normalized, 1, function(x)
      (e_max[x['Gnumber'],x['CellLineName']] + (1-e_max[x['Gnumber'],x['CellLineName']])*
        (ec50[x['Gnumber'],x['CellLineName']]**2 / (as.numeric(x['Concentration'])**2 +
          ec50[x['Gnumber'],x['CellLineName']]**2))) *
      (e_max[x['Gnumber_2'],x['CellLineName']] + (1-e_max[x['Gnumber_2'],x['CellLineName']])*
        (ec50[x['Gnumber_2'],x['CellLineName']]**2 / (as.numeric(x['Concentration_2'])**2 +
          ec50[x['Gnumber_2'],x['CellLineName']]**2)))
        )

# add some synergy/antagonism
idx = df_normalized$CellLineName == 'HCC2218' & df_normalized$Gnumber == 'G02442104' &
                    df_normalized$Gnumber_2 == 'G00505032'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] +
                            .2 * (100 - df_normalized$ReadoutValue[idx])

idx = df_normalized$CellLineName == 'HCC2218' & df_normalized$Gnumber == 'Trametinib' &
                                df_normalized$Gnumber_2 == 'G00033829'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] -
                            .1 * (100 - df_normalized$ReadoutValue[idx])

idx = df_normalized$CellLineName == 'COV318' & df_normalized$Gnumber == 'Palbociclib'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] +
                            .1 * (100 - df_normalized$ReadoutValue[idx])

idx = df_normalized$CellLineName == 'COV318' & df_normalized$Gnumber_2 != 'G00033829'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] -
                            .15 * (100 - df_normalized$ReadoutValue[idx])

df_normalized$ReadoutValue = pmin(df_normalized$ReadoutValue, 102)

# add some noise
set.seed(1)
df_normalized$ReadoutValue = df_normalized$ReadoutValue + runif(dim(df_normalized)[1], -4, 3)

levels(df_normalized$Gnumber) = c(levels(df_normalized$Gnumber), 'Vehicle')
df_normalized$Gnumber[df_normalized$Concentration == 0] = 'Vehicle'
levels(df_normalized$Gnumber_2) = c(levels(df_normalized$Gnumber_2), 'Vehicle')
df_normalized$Gnumber_2[df_normalized$Concentration_2 == 0] = 'Vehicle'

swap_idx = df_normalized$Concentration == 0 & df_normalized$Gnumber_2 != 'Vehicle'
df_normalized[swap_idx, c('Gnumber', 'Concentration')] =
  df_normalized[swap_idx, c('Gnumber_2', 'Concentration_2')]
df_normalized$Gnumber_2[swap_idx] = 'Vehicle'
df_normalized$Concentration_2[swap_idx] = 0

# decompose de table into original files
df_raw_data = df_normalized[order(df_normalized$Barcode),]
df_raw_data$Template = 'Template_trt.tsv'
df_raw_data$WellRow = sort(rep(LETTERS[1:8],12))[-6:-1]
df_raw_data$WellColumn = rep(1:12,8)[-6:-1]
levels(df_raw_data$Barcode) = paste0('Plate',1:6)

df_manifest = unique(df_raw_data[, c('Barcode', 'Template', 'Duration', 'clid')])
df_data = df_raw_data[,c('Barcode', 'WellRow', 'WellColumn', 'ReadoutValue', 'BackgroundValue')]
df_treatment = unique(df_raw_data[, c('Template', 'WellRow', 'WellColumn', 'Gnumber',
                'Gnumber_2', 'Concentration', 'Concentration_2')])

dir.create('../inst/testdata/data12')
write.table(df_manifest, '../inst/testdata/data12/Manifest_data12.tsv',
        sep='\t', quote=F, row.names=F)
write.table(df_data, '../inst/testdata/data12/RawData_data12.tsv',
        sep='\t', quote=F, row.names=F)
for (trt_f in unique(df_treatment$Template)) {
        write.table(df_treatment[df_treatment$Template==trt_f,
                setdiff(colnames(df_treatment), 'Template')],
            paste0('../inst/testdata/data12/', trt_f),
                sep='\t', quote=F, row.names=F)
}

write.table(df_normalized, '../inst/testdata/data12/calculated_normalized_data.tsv',
            sep='\t', quote=F, row.names=F)
