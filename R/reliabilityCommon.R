

.reliabilityReadData <- function(dataset, options) {

  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.numeric = variables, columns.as.factor = NULL, exclude.na.listwise = NULL)
  }
  return(dataset)
}

.reliabilityCheckErrors <- function(dataset, options) {

  # check for existing inverse
  .checkInverse <- function() {
    if (length(options[["variables"]]) > 2) {
      use.cases <- "everything"
      if (anyNA(dataset)) {
        if (options[["missingValues"]] == "excludeCasesPairwise")
          use.cases <- "pairwise.complete.obs"
        else if (options[["missingValues"]] == "excludeCasesListwise")
          use.cases <- "complete.obs"
      }
      if (isTryError(try(solve(cov(dataset, use = use.cases)),silent=TRUE))) {
        return(gettext("The covariance matrix of the data is not invertible"))
      }
    }
    return(NULL)
  }

  .hasErrors(dataset = dataset, options = options, perform = "run",
             type = c("infinity", "variance", "observations", "varCovData"),
             observations.amount = " < 3",
             varCovData.corFun = function(x) cor(x, use = "pairwise.complete.obs"),
             custom = .checkInverse,
             exitAnalysisIfErrors = TRUE)

}


.reliabilityCheckLoadings <- function(dataset, variables) {
  if (ncol(dataset > 2)) {
    prin <- psych::principal(dataset)
    idx <- prin[["loadings"]] < 0
    sidx <- sum(idx)
    if (sidx == 0) {
      footnote <- ""
    } else {
      footnote <- sprintf(ngettext(length(variables[idx]),
                                   "The following item correlated negatively with the scale: %s. ",
                                   "The following items correlated negatively with the scale: %s. "),
                          paste(variables[idx], collapse = ", "))
    }
  } else {
    return("Please enter at least 3 Variables to do an analysis")
  }
}

.reverseScoreItems <- function(dataset, options) {
  dataset_rev <- as.matrix(dataset) # fails for string factors!
  cols <- match(unlist(options[["reverseScaledItems"]]), .unv(colnames(dataset)))
  total <- apply(as.matrix(dataset[, cols]), 2, min, na.rm = T) + apply(as.matrix(dataset[, cols]), 2, max, na.rm = T)
  dataset_rev[ ,cols] <- matrix(rep(total, nrow(dataset)), nrow(dataset), length(cols), byrow=T) - dataset[ ,cols]
  return(dataset_rev)
}

.getStateContainer <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies=c("variables", "reverseScaledItems", "noSamples", "missingValues", "bootType",
                  "setSeed", "seedValue", "intervalOn")
    )

    return(jaspResults[["stateContainer"]])
}