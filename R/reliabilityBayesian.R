

reliabilityBayesian<- function(jaspResults, dataset, options) {

  sink(file="~/Downloads/log.txt")
  on.exit(sink(NULL))

  dataset <- .reliabilityReadData(dataset, options)
  .reliabilityCheckErrors(dataset, options)

  model <- .reliabilityBayesianPreCalc(jaspResults, dataset, options)
  model[["alpha"]] <- .reliabilityBayesianAlpha(jaspResults, dataset, options, model)
  # lambda2 <- .reliabilityFrequentistLambda2(jaspResults, dataset, options)
  # lambda6 <- .reliabilityFrequentistLambda6(jaspResults, dataset, options)
  # glb <- .reliabilityFrequentistGlb(jaspResults, dataset, options)
  # omega <- .reliabilityFrequentistOmega(jaspResults, dataset, options)
  # scaleStats <- .reliabilityFrequentistScaleStats(jaspResults, dataset, options)


  # .frequentistReliabilityScaleTable(         jaspResults, model, options)
  # .frequentistReliabilityItemTable(          jaspResults, model, options)
  # .freqentistReliabilitySingleFactorFitTable(jaspResults, model, options)

  return()

}

.BayesianReliabilityDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("mcDonaldScale", "alphaScale", "guttman2Scale", "guttman6Scale", "glbScale",
                                           "averageInterItemCor", "meanScale", "sdScale")]),
    selectedEstimatorsPlots  = unlist(options[c("mcDonaldScale", "alphaScale", "guttman2Scale", "guttman6Scale",
                                                "glbScale")]),
    itemDroppedSelected = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "guttman6Item","glbItem",
                                           "itemRestCor", "meanItem", "sdItem")]),
    itemDroppedSelectedItem = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "guttman6Item",
                                               "glbItem")]),
    namesEstimators     = list(
      tables = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                 "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"),
      tables_item = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                      gettext("Greatest Lower Bound"), gettext("Item-rest correlation"), gettext("mean"), gettext("sd")),
      coefficients = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                       gettext("Greatest Lower Bound"), gettext("Item-rest correlation")),
      plots = list(expression("McDonald's"~omega), expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]),
                   expression("Guttman's"~lambda[6]), gettext("Greatest Lower Bound"))
    ),
    order_end = c(5, 1, 2, 3, 4) # order for plots and such, put omega to the front
  )
  return(derivedOptions)
}


.reliabilityBayesianPreCalc <- function(jaspResults, dataset, options) {


}

