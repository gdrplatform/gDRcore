
# additional function to check if all paths in vector are readable
is.readables <- function(paths){
  missing_path_string <- paste(paths[as.logical(-file.access(paths, 4))], collapse = ', ', sep = '   ')
  message <- paste0("Following path(s) with no read permission found: '", missing_path_string,"'")
  assertthat::assert_that(sum(file.access(paths,4)) == 0, msg = message)
}

