
devtools::load_all('~/Rpackages/GeneDataScreenR/') # will be accessible on servers but currently working on branch   gDR
devtools::load_all('../gDR/')


qcs <- GeneDataScreenR::QCSession('QCS-29233') # projet 39
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG1"), "/QCSobject.RDS"))

qcs <- GeneDataScreenR::QCSession('QCS-22668') # project 17
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG2"), "/QCSobject.RDS"))

qcs <- GeneDataScreenR::QCSession('QCS-31246')
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG3"), "/QCSobject.RDS"))

qcs <- GeneDataScreenR::QCSession('QCS-30210')
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG4"), "/QCSobject.RDS"))

qcs <- GeneDataScreenR::QCSession('QCS-26438') # --> co-dilution series
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG5"), "/QCSobject.RDS"))

qcs <- GeneDataScreenR::QCSession('QCS-32933') # --> cotreatment matrix
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG6"), "/QCSobject.RDS"))

qcs <- GeneDataScreenR::QCSession('QCS-21554') # --> codilution and fixed concentration
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG7"), "/QCSobject.RDS"))
