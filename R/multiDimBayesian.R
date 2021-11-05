

reliabilityMultiDimBayesian <- function(jaspResults, dataset, options) {

  dataset <- .readData(dataset, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }

  return()

}
