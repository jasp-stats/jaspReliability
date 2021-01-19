

reliabilityBayesian <- function(jaspResults, dataset, options) {

  # sink(file="~/Downloads/logbay.txt")
  # on.exit(sink(NULL))

  dataset <- .readData(dataset, options)
  .checkErrors(dataset, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }



  model <- .BayesianPreCalc(jaspResults, dataset, options)
  options <- .scaleItemBoxAlign(options)
  model[["itemDroppedCovs"]] <- .BayesianItemDroppedMats(jaspResults, dataset, options, model)
  model[["derivedOptions"]] <- .BayesianDerivedOptions(options)
  model[["omega"]] <- .BayesianOmega(jaspResults, dataset, options, model)
  model[["alpha"]] <- .BayesianAlpha(jaspResults, dataset, options, model)
  model[["lambda2"]] <- .BayesianLambda2(jaspResults, dataset, options, model)
  model[["lambda6"]] <- .BayesianLambda6(jaspResults, dataset, options, model)
  model[["glb"]] <- .BayesianGlb(jaspResults, dataset, options, model)
  model[["average"]] <- .BayesianAverageCor(jaspResults, dataset, options, model)
  model[["mean"]] <- .BayesianMean(jaspResults, dataset, options, model)
  model[["sd"]] <- .BayesianStdDev(jaspResults, dataset, options, model)
  model[["itemRestCor"]] <- .BayesianItemRestCor(jaspResults, dataset, options, model)
  model[["itemMean"]] <- .BayesianItemMean(jaspResults, dataset, options, model)
  model[["itemSd"]] <- .BayesianItemSd(jaspResults, dataset, options, model)

  .BayesianScaleTable(         jaspResults, model, options)
  .BayesianItemTable(          jaspResults, model, options)
  .BayesianProbTable(          jaspResults, model, options)
  .BayesianPosteriorPlot(      jaspResults, model, options)
  .BayesianIfItemPlot(         jaspResults, model, options)
  .BayesianPosteriorPredictive(jaspResults, model, options)
  .BayesianTracePlot(          jaspResults, model, options)

  return()

}

.BayesianDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                           "averageInterItemCor", "meanScale", "sdScale")]),
    selectedEstimatorsPlots  = unlist(options[c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale",
                                                "glbScale")]),
    itemDroppedSelected = unlist(options[c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item","glbItem",
                                           "itemRestCor", "itemMean", "itemSd")]),
    itemDroppedSelectedItem = unlist(options[c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item",
                                               "glbItem")]),

    namesEstimators     = list(
      tables = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                 "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"),
      tables_item = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                      gettext("Greatest Lower Bound"), gettext("Item-rest correlation"), gettext("mean"), gettext("sd")),
      coefficients = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                       gettext("Greatest Lower Bound"), gettext("Item-rest correlation")),
      plots = list(expression("McDonald's"~omega), expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]),
                   expression("Guttman's"~lambda[6]), gettext("Greatest Lower Bound")),
      plotsNoGreek = c("omega", "alpha", "lambda2", "lambda6", "glb")
    )

  )
  return(derivedOptions)
}


