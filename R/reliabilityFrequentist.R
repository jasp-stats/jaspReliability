
reliabilityFrequentist <- function(jaspResults, dataset, options) {

  dataset <- .readData(dataset, options)
  .checkErrors(dataset, options)

  model <- .frequentistPreCalc(jaspResults, dataset, options)
  model[["itemDroppedCovs"]] <- .frequentistItemDroppedMats(jaspResults, dataset, options, model)
  model[["derivedOptions"]] <- .frequentistDerivedOptions(options)
  model[["omega"]] <- .frequentistOmega(jaspResults, dataset, options, model)
  model[["alpha"]] <- .frequentistAlpha(jaspResults, dataset, options, model)
  model[["lambda2"]] <- .frequentistLambda2(jaspResults, dataset, options, model)
  model[["lambda6"]] <- .frequentistLambda6(jaspResults, dataset, options, model)
  model[["glb"]] <- .frequentistGlb(jaspResults, dataset, options, model)
  model[["average"]] <- .frequentistAverageCor(jaspResults, dataset, options, model)
  model[["mean"]] <- .frequentistMean(jaspResults, dataset, options, model)
  model[["sd"]] <- .frequentistStdDev(jaspResults, dataset, options, model)
  model[["itemRestCor"]] <- .frequentistItemRestCor(jaspResults, dataset, options, model)
  model[["itemMean"]] <- .frequentistItemMean(jaspResults, dataset, options, model)
  model[["itemSd"]] <- .frequentistItemSd(jaspResults, dataset, options, model)

  .frequentistScaleTable(         jaspResults, model, options)
  .frequentistItemTable(          jaspResults, model, options)
  .freqentistSingleFactorFitTable(jaspResults, model, options)

  return()

}

.frequentistDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("omegaScale","alphaScale", "lambda2Scale", "lambda6Scale",
                                            "glbScale", "averageInterItemCor", "meanScale", "sdScale")]),
    itemDroppedSelected = unlist(options[c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item",
                                            "glbItem", "itemRestCor", "itemMean", "itemSd")]),
    namesEstimators     = list(
      tables = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                 "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"),
      tables_item = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                      gettext("Greatest Lower Bound"), gettext("Item-rest correlation"), gettext("mean"), gettext("sd")),
      coefficients = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                       gettext("Greatest Lower Bound"))
    )
  )
  return(derivedOptions)
}

