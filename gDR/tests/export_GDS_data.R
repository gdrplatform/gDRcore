
devtools::load_all('~/Rpackages/GeneDataScreenR/') # will be accessible on servers but currently working on branch   gDR
devtools::load_all('../gDR/')


qcs <- QCSession('QCS-29233') # projet 39
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG1"), "/QCSobject.RDS"))

qcs <- QCSession('QCS-22668') # project 17
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG2"), "/QCSobject.RDS"))

qcs <- QCSession('QCS-31246')
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG3"), "/QCSobject.RDS"))

qcs <- QCSession('QCS-30210')
saveRDS(qcs, paste0(system.file(package = "gDR", "testdata", "dataG4"), "/QCSobject.RDS"))
