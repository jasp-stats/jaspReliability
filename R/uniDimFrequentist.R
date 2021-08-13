
reliabilityUniDimFrequentist <- function(jaspResults, dataset, options) {

  dataset <- .readData(dataset, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }

  .checkErrors(dataset, options)


  model <- .frequentistPreCalc(jaspResults, dataset, options)
  options <- .scaleItemBoxAlign(options)

  model[["derivedOptions"]] <- .frequentistDerivedOptions(options)
  model[["omegaScale"]] <- .frequentistOmegaScale(jaspResults, dataset, options, model)
  model[["omegaItem"]] <- .frequentistOmegaItem(jaspResults, dataset, options, model)
  model[["alphaScale"]] <- .frequentistAlphaScale(jaspResults, dataset, options, model)
  model[["alphaItem"]] <- .frequentistAlphaItem(jaspResults, dataset, options, model)
  model[["lambda2Scale"]] <- .frequentistLambda2Scale(jaspResults, dataset, options, model)
  model[["lambda2Item"]] <- .frequentistLambda2Item(jaspResults, dataset, options, model)
  model[["lambda6Scale"]] <- .frequentistLambda6Scale(jaspResults, dataset, options, model)
  model[["lambda6Item"]] <- .frequentistLambda6Item(jaspResults, dataset, options, model)
  model[["glbScale"]] <- .frequentistGlbScale(jaspResults, dataset, options, model)
  model[["glbItem"]] <- .frequentistGlbItem(jaspResults, dataset, options, model)

  model[["averageInterItemCor"]] <- .frequentistAverageCor(jaspResults, dataset, options, model)
  model[["meanScale"]] <- .frequentistMean(jaspResults, dataset, options, model)
  model[["sdScale"]] <- .frequentistStdDev(jaspResults, dataset, options, model)
  model[["itemRestCor"]] <- .frequentistItemRestCor(jaspResults, dataset, options, model)
  model[["meanItem"]] <- .frequentistMeanItem(jaspResults, dataset, options, model)
  model[["sdItem"]] <- .frequentistSdItem(jaspResults, dataset, options, model)

  model[["scaleResults"]] <- .frequentistComputeScaleResults(jaspResults, dataset, options, model)

  .frequentistScaleTable(jaspResults, model, options)
  .frequentistItemTable(jaspResults, model, options)
  .frequentistSingleFactorFitTable(jaspResults, model, options)


  return()

}

.frequentistDerivedOptions <- function(options) {

  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale",
                                           "glbScale", "averageInterItemCor", "meanScale", "sdScale")]),
    itemDroppedSelected = unlist(options[c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item",
                                           "glbItem", "itemRestCor", "meanItem", "sdItem")]),
    namesEstimators     = list(
      tables = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                 "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"),
      tables_item = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                      gettext("Greatest Lower Bound"), gettext("Item-rest correlation"), gettext("mean"), gettext("sd")),
      coefficients = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                       gettext("Greatest Lower Bound")))
  )
  return(derivedOptions)
}


.getStateContainerF <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = c("variables", "reverseScaledItems", "noSamples",
                                                                          "missingValues", "bootType", "setSeed",
                                                                          "seed", "intervalOn", "disableSampleSave")
  )

  return(jaspResults[["stateContainer"]])
}

.freqItemDroppedStats <- function(Cov, f = function(){}) {

  out <- numeric(ncol(Cov))
  for (i in seq_len(ncol(Cov))) {
    out[i] <- f(Cov[-i, -i])
  }
  return(out)
}
