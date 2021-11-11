

reliabilityUniDimBayesian <- function(jaspResults, dataset, options) {

  options <- .parseAndStoreFormulaOptions(jaspResults, options, "iwScale")

  dataset <- .readData(dataset, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }


  .checkErrors(dataset, options, Bayes = TRUE)

  model <- .BayesianPreCalc(jaspResults, dataset, options)

  options <- .scaleItemBoxAlign(options)


  model[["derivedOptions"]] <- .BayesianDerivedOptions(options)
  model[["gibbsCor"]] <- .BayesianStdCov(jaspResults, dataset, options, model)
  model[["omegaScale"]] <- .BayesianOmegaScale(jaspResults, dataset, options, model)
  model[["omegaItem"]] <- .BayesianOmegaItem(jaspResults, dataset, options, model)
  model[["omegaScaleStd"]] <- .BayesianOmegaScaleStd(jaspResults, dataset, options, model)
  model[["omegaItemStd"]] <- .BayesianOmegaItemStd(jaspResults, dataset, options, model)
  model[["alphaScale"]] <- .BayesianAlphaScale(jaspResults, dataset, options, model)
  model[["alphaItem"]] <- .BayesianAlphaItem(jaspResults, dataset, options, model)
  model[["lambda2Scale"]] <- .BayesianLambda2Scale(jaspResults, dataset, options, model)
  model[["lambda2Item"]] <- .BayesianLambda2Item(jaspResults, dataset, options, model)
  model[["lambda6Scale"]] <- .BayesianLambda6Scale(jaspResults, dataset, options, model)
  model[["lambda6Item"]] <- .BayesianLambda6Item(jaspResults, dataset, options, model)
  model[["glbScale"]] <- .BayesianGlbScale(jaspResults, dataset, options, model)
  model[["glbItem"]] <- .BayesianGlbItem(jaspResults, dataset, options, model)

  model[["averageInterItemCor"]] <- .BayesianAverageCor(jaspResults, dataset, options, model)
  model[["meanScale"]] <- .BayesianMean(jaspResults, dataset, options, model)
  model[["sdScale"]] <- .BayesianStdDev(jaspResults, dataset, options, model)
  model[["itemRestCor"]] <- .BayesianItemRestCor(jaspResults, dataset, options, model)
  model[["meanItem"]] <- .BayesianMeanItem(jaspResults, dataset, options, model)
  model[["sdItem"]] <- .BayesianSdItem(jaspResults, dataset, options, model)

  model[["scaleResults"]] <- .BayesianComputeScaleResults(jaspResults, options, model)
  model[["itemResults"]] <- .BayesianComputeItemResults(jaspResults, options, model)

  .BayesianScaleTable(jaspResults, model, options)
  .BayesianItemTable(jaspResults, model, options)
  .BayesianProbTable(jaspResults, model, options)
  .BayesianLoadingsTable(jaspResults, model, options)
  .BayesianPosteriorPlot(jaspResults, model, options)
  .BayesianIfItemPlot(jaspResults, model, options)
  .omegaPosteriorPredictive(jaspResults, model, options)
  .BayesianTracePlot(jaspResults, model, options)
  return()

}

.BayesianDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                           "averageInterItemCor", "meanScale", "sdScale")]),
    selectedEstimatorsPlots  = unlist(options[c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale",
                                                "glbScale")]),
    itemDroppedSelected = unlist(options[c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item", "glbItem",
                                           "itemRestCor", "meanItem", "sdItem")]),
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

