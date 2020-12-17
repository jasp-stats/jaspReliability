
reliabilityFrequentist <- function(jaspResults, dataset, options) {

  # sink(file="~/Downloads/log.txt")
  # on.exit(sink(NULL))

  dataset <- .reliabilityReadData(dataset, options)
  .reliabilityCheckErrors(dataset, options)

  model <- .reliabilityFrequentistPreCalc(jaspResults, dataset, options)
  # model[["alpha"]] <- .reliabilityFrequentistAlpha(jaspResults, dataset, options, model)
  model[["guttman2"]] <- .reliabilityFrequentistGuttman2(jaspResults, dataset, options, model)

  # lambda6 <- .reliabilityFrequentistLambda6(jaspResults, dataset, options)
  # glb <- .reliabilityFrequentistGlb(jaspResults, dataset, options)
  # omega <- .reliabilityFrequentistOmega(jaspResults, dataset, options)
  # scaleStats <- .reliabilityFrequentistScaleStats(jaspResults, dataset, options)


  # .frequentistReliabilityScaleTable(         jaspResults, model, options)
  # .frequentistReliabilityItemTable(          jaspResults, model, options)
  # .freqentistReliabilitySingleFactorFitTable(jaspResults, model, options)

  return()

}

.reliabilityFrequentistDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("mcDonaldScale","alphaScale", "guttman2Scale", "guttman6Scale",
                                            "glbScale", "averageInterItemCor", "meanScale", "sdScale")]),
    itemDroppedSelected = unlist(options[c("mcDonaldItem", "alphaItem", "guttman2Item", "guttman6Item",
                                            "glbItem", "itemRestCor", "meanItem", "sdItem")]),
    namesEstimators     = list(
      tables = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                 "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"),
      tables_item = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                      gettext("Greatest Lower Bound"), gettext("Item-rest correlation"), gettext("mean"), gettext("sd")),
      coefficients = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                       gettext("Greatest Lower Bound")),
      plots = list(expression("McDonald's"~omega), expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]),
                   expression("Guttman's"~lambda[6]), gettext("Greatest Lower Bound"))
    )
  )

  return(derivedOptions)
}


.reliabilityFrequentistPreCalc <- function(jaspResults, dataset, options) {

  if (!is.null(jaspResults[["stateContainer"]][["modelObj"]]$object))
    return(jaspResults[["stateContainer"]][["modelObj"]]$object)

  derivedOptions <- .reliabilityFrequentistDerivedOptions(options)

  # what if no coefficient boxes are checked?
  if(!any(derivedOptions[["selectedEstimators"]])) {
    variables <- options[["variables"]]
    if (length(options[["reverseScaledItems"]]) > 0L) {
      dataset <- .reverseScoreItems(dataset, options)
    }
    model <- list()
    model[["footnote"]] <- .reliabilityCheckLoadings(dataset, variables)
    return(model)
  }

  model <- jaspResults[["modelObj"]]$object

  if (is.null(model)) {
    model <- list()

    # check if any items correlate negatively with the scale
    model[["footnote"]] <- .reliabilityCheckLoadings(dataset, options[["variables"]])

    # check for missings and determine the missing handling
    if (anyNA(dataset)) {
      if (options[["missingValues"]] == "excludeCasesPairwise") {
        missing <- "pairwise"
        use.cases <- "pairwise.complete.obs"
        model[["footnote"]] <- gettextf("%s Of the observations, pairwise complete cases were used. ",
                                        model[["footnote"]])
      } else {
        pos <- which(is.na(dataset), arr.ind = TRUE)[, 1]
        dataset <- dataset[-pos, ]
        use.cases <- "complete.obs"
        model[["footnote"]] <- gettextf("%s Of the observations, %1.f complete cases were used. ",
                                        model[["footnote"]], nrow(dataset))
      }
    } else {
      use.cases <- "everything"
    }

    # check for inverse scored items
    if (length(options[["reverseScaledItems"]]) > 0L) {
      dataset <- .reverseScoreItems(dataset, options)
    }
    p <- ncol(dataset); n <- nrow(dataset)
    cc <- cov(dataset, use = use.cases)
    model[["data_cov"]] <- cc
  }

  # check if interval is checked and bootstrapped covariance sample has to be generated
  if (options[["intervalOn"]] && is.null(model[["bootsamp"]]) &&
      (
        (options[["mcDonaldScale"]] && (options[["omegaEst"]] == "pfa")) ||
        (options[["alphaScale"]] && !(options[["alphaInterval"]] == "alphaAnalytic")) ||
        options[["guttman2Scale"]] || options[["guttman6Scale"]] || options[["glbScale"]]
        )
      ) {

    if (options[["setSeed"]])
      set.seed(options[["seedValue"]])

    boot_cov <- array(0, c(options[["noSamples"]], p, p))

    if (options[["bootType"]] == "bootNonPara") {
      for (i in 1:options[["noSamples"]]){
        boot_data <- as.matrix(dataset[sample.int(n, size = n, replace = TRUE), ])
        boot_cov[i, , ] <- cov(boot_data, use = use.cases)
      }

    } else if (options[["bootType"]] == "bootPara") {
      for (i in 1:options[["noSamples"]]){
        boot_data <- MASS::mvrnorm(n, colMeans(dataset, na.rm = T), cc)
        boot_cov[i, , ] <- cov(boot_data)
      }
    }
    model[["bootsamp"]] <- boot_cov
  }

  # check for the covariance matrices of the item deleted cases
  if (is.null(model[["itemDroppedCovs"]]) &&
      options[["mcDonaldItem"]] || options[["alphaItem"]] || options[["guttman2Item"]] ||
      options[["guttman6Item"]] || options[["glbItem"]]) {

    Ctmp <- array(0, c(p, p - 1, p - 1))
    for (i in 1:p){
      Ctmp[i, , ] <- cc[-i, -i]
    }
    model[["itemDroppedCovs"]] <- Ctmp
  }

  if (is.null(jaspResults[["stateContainer"]])) {
    jaspResults[["stateContainer"]] <- createJaspContainer()
    jaspResults[["stateContainer"]]$dependOn(
      options = c("variables", "reverseScaledItems", "noSamples", "missingValues", "bootType",
                  "setSeed", "seedValue", "intervalOn")
    )
  }

  modelState <- createJaspState(model)
  jaspResults[["stateContainer"]] <- modelState

  return(model)
}

