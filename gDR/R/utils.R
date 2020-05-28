
DICTIONARY <- list(
  Untreated = "untreated",
  Vehicle = "vehicle",
  vehcle = "vehicle"
)

standardize_record_values <- function(x, dictionary = DICTIONARY){
  for (i in 1:length(dictionary)) {
    x[x == names(dictionary[i])] <- dictionary[[i]]
  }
  x
}
