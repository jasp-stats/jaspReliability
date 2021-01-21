

.getStateContainerF <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainerF"]]))
    return(jaspResults[["stateContainerF"]])

  jaspResults[["stateContainerF"]] <- createJaspContainer(dependencies=c("variables", "reverseScaledItems", "noSamples",
                                                                         "missingValues", "bootType","setSeed",
                                                                         "seedValue", "intervalOn")
  )

  return(jaspResults[["stateContainerF"]])
}

.getStateContainerB <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainerB"]]))
    return(jaspResults[["stateContainerB"]])

  jaspResults[["stateContainerB"]] <- createJaspContainer(dependencies=c("variables", "reverseScaledItems", "noSamples",
                                                                         "noBurnin", "noThin", "noChains",
                                                                         "missingValues","setSeed", "seedValue")
  )

  return(jaspResults[["stateContainerB"]])
}
