
#' @importFrom jaspBase createJaspContainer createJaspHtml createJaspPlot createJaspQmlSource createJaspState createJaspTable
#' progressbarTick .quitAnalysis .readDataSetToEnd startProgressbar

#' @importFrom stats approxfun cor cov cov2cor density integrate median na.omit pnorm qchisq qnorm quantile rgamma rnorm sd var

#' @export
unidimensionalReliabilityBayesian <- function(jaspResults, dataset, options) {

  options <- jaspBase::.parseAndStoreFormulaOptions(jaspResults, options, "inverseWishartPriorScale")

  # check for listwise deletion
  datasetOld <- dataset
  dataset <- .handleData(datasetOld, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }

  .checkErrors(dataset, options, Bayes = TRUE)

  model <- .BayesianPreCalc(jaspResults, dataset, options, datasetOld)

  options <- .scaleItemBoxAlignB(options)


  model[["derivedOptions"]] <- .BayesianDerivedOptions(options)
  model[["gibbsCor"]] <- .BayesianStdCov(jaspResults, dataset, options, model)
  model[["singleFactor"]] <- .BayesianSingleFactorModel(jaspResults, dataset, options, model)
  model[["singleFactorItem"]] <- .BayesianSingleFactorModelItem(jaspResults, dataset, options, model)
  model[["scaleOmega"]] <- .BayesianOmegaScale(jaspResults, dataset, options, model)
  model[["itemDeletedOmega"]] <- .BayesianOmegaItem(jaspResults, dataset, options, model)
  model[["scaleAlpha"]] <- .BayesianAlphaScale(jaspResults, dataset, options, model)
  model[["itemDeletedAlpha"]] <- .BayesianAlphaItem(jaspResults, dataset, options, model)
  model[["scaleLambda2"]] <- .BayesianLambda2Scale(jaspResults, dataset, options, model)
  model[["itemDeletedLambda2"]] <- .BayesianLambda2Item(jaspResults, dataset, options, model)
  model[["scaleSplithalf"]] <- .BayesianSplithalfScale(jaspResults, dataset, options, model)
  model[["itemDeletedSplithalf"]] <- .BayesianSplithalfItem(jaspResults, dataset, options, model)
  model[["averageInterItemCorrelation"]] <- .BayesianAverageCor(jaspResults, dataset, options, model)
  model[["scaleMean"]] <- .BayesianMean(jaspResults, dataset, options, model)
  model[["scaleVar"]] <- .BayesianVar(jaspResults, dataset, options, model)
  model[["scaleSd"]] <- .BayesianStdDev(jaspResults, dataset, options, model)
  model[["itemRestCorrelation"]] <- .BayesianItemRestCorrelation(jaspResults, dataset, options, model)
  model[["itemMean"]] <- .BayesianMeanItem(jaspResults, dataset, options, model)
  model[["itemVar"]] <- .BayesianVarItem(jaspResults, dataset, options, model)
  model[["itemSd"]] <- .BayesianSdItem(jaspResults, dataset, options, model)

  model[["scaleResults"]] <- .BayesianComputeScaleResults(jaspResults, options, model)
  model[["itemResults"]] <- .BayesianComputeItemResults(jaspResults, options, model)

  model[["omegaFitMeasures"]] <- .BayesianOmegaFitMeasures(jaspResults, dataset, options, model)

  .BayesianScaleTable(jaspResults, model, options)
  .BayesianItemTable(jaspResults, model, options)
  .BayesianProbTable(jaspResults, model, options)
  .BayesianLoadingsTable(jaspResults, model, options)
  .BayesianOmegaFitMeasuresTable(jaspResults, model, options)
  .BayesianPosteriorPlot(jaspResults, model, options)
  .BayesianIfItemPlot(jaspResults, model, options)
  .omegaPosteriorPredictive(jaspResults, model, options)
  .BayesianTracePlot(jaspResults, model, options)
  return()

}

.BayesianDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("scaleOmega", "scaleAlpha", "scaleLambda2", "scaleSplithalf",
                                           "averageInterItemCorrelation", "scaleMean", "scaleVar", "scaleSd")]),
    selectedEstimatorsPlots  = unlist(options[c("scaleOmega", "scaleAlpha", "scaleLambda2", "scaleSplithalf")]),
    itemDroppedSelected = unlist(options[c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2", "itemDeletedSplithalf",
                                           "itemRestCorrelation", "itemMean", "itemVar", "itemSd")]),
    itemDroppedSelectedItem = unlist(options[c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2", "itemDeletedSplithalf")]),

    namesEstimators     = list(
      tables = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2", "Split-half coefficient",
                 "Average interitem correlation", "Mean", "Variance", "SD"),
      tables_item = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2", gettext("Split-half coefficient"),
                      gettext("Item-rest correlation"), gettext("Mean"), gettext("Variance"), gettext("SD")),
      coefficients = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2", gettext("Split-half coefficient"),
                       gettext("Item-rest correlation")),
      plots = list(expression("Coefficient"~omega), expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]),
                   gettext("Split-half coefficient")),
      plotsNoGreek = c("omega", "alpha", "lambda2", "splithalf")
    )

  )
  return(derivedOptions)
}


#### Precalulate results ####

.BayesianPreCalc <- function(jaspResults, dataset, options, datasetOld) {

  if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object)) {
    if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object$gibbsSamp))
      return(.getStateContainerB(jaspResults)[["modelObj"]]$object)
  }


  derivedOptions <- .BayesianDerivedOptions(options)

  # what if no coefficient boxes are checked?
  if (!any(derivedOptions[["selectedEstimators"]]) && !any(derivedOptions[["itemDroppedSelected"]])) {
    variables <- options[["variables"]]
    empty <-  TRUE
    model <- list(empty = empty)
    model[["footnote"]] <- .checkLoadings(dataset, variables)
    return(model)
  }


  # what if too few variables are entered:
  if (length(options[["variables"]]) < 2) {
    empty <-  TRUE
    model <- list(empty = empty)
    model[["footnote"]] <- .atLeast2Variables()
    model[["itemsDropped"]] <- colnames(dataset)
    return(model)
  }

  model <- jaspResults[["modelObj"]]$object

  if (is.null(model)) {
    model <- list()

    model[["footnote"]] <- ""

    # check for missings and determine the missing handling

    # check for missings and determine the missing handling
    if (options[["naAction"]] == "listwise" && nrow(datasetOld) > nrow(dataset)) { # this indicates listwise deletion
      model[["use.cases"]] <- "complete.obs"
      model[["pairwise"]] <- FALSE
      model[["footnote"]] <- gettextf("%s Of the observations, %1.f complete cases were used. ",
                                      model[["footnote"]], nrow(dataset))

    } else if (anyNA(dataset)) { # when pairwise deletion
      model[["use.cases"]] <- "pairwise.complete.obs"
      model[["pairwise"]] <- TRUE
      model[["footnote"]] <- gettextf("%s Of the observations, pairwise complete cases were used. ",
                                      model[["footnote"]])

    } else {
      model[["use.cases"]] <- "everything"
      model[["pairwise"]] <- FALSE
    }

    cc <- cov(dataset, use = model[["use.cases"]])
    model[["data_cov"]] <- cc
    model[["k"]] <- ncol(dataset)
    model[["n"]] <- nrow(dataset)

    # check if any items correlate negatively with the scale
    model[["footnote"]] <- sprintf("%s%s", model[["footnote"]], .checkLoadings(dataset, options[["variables"]]))
  }

  jaspBase::.setSeedJASP(options)

  # check if posterior cov sample already exists and any of the relevant coefficients are checked
  if (is.null(model[["gibbsSamp"]]) &&
      (options[["scaleAlpha"]] || options[["scaleLambda2"]] || options[["scaleSplithalf"]] ||
       options[["averageInterItemCorrelation"]])
  ) {

    startProgressbar(options[["samples"]] * options[["chains"]])
    dataset <- scale(dataset, scale = FALSE)
    c_out <- try(Bayesrel:::covSamp(dataset, options[["samples"]], options[["burnin"]],
                                    options[["thinning"]], options[["chains"]],
                                    model[["pairwise"]], progressbarTick, k0 = options[["inverseWishartPriorScale"]],
                                    df0 = options[["inverseWishartPriorDf"]]), silent = TRUE)
    if (inherits(c_out, "try-error")) {
      if (model[["pairwise"]]) {
        .quitAnalysis(gettext("Sampling the posterior covariance matrix for either one of [alpha, lambda2, lambda6, glb] failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
      } else {
        .quitAnalysis(gettext("Sampling the posterior covariance matrix for either one of [alpha, lambda2, lambda6, glb] failed."))
      }
    }

    model[["gibbsSamp"]] <- c_out$cov_mat

  }

  model[["progressbarLength"]] <- options[["chains"]] *
    length(seq(1, options[["samples"]] - options[["burnin"]], options[["thinning"]]))

  model[["itemsDropped"]] <- colnames(dataset)

  if (options[["samplesSavingDisabled"]])
    return(model)

  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model, dependencies = c("inverseWishartPriorScale", "inverseWishartPriorDf"))

  return(model)
}


# standardize covariance matrix sample
.BayesianStdCov <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["gibbsCor"]]$object))
    return(.getStateContainerB(jaspResults)[["gibbsCor"]]$object)

  if (is.null(model[["gibbsSamp"]])) {
    return()
  } else {
    if (options[["coefficientType"]] == "standardized" || options[["averageInterItemCorrelation"]] || options[["scaleSplithalf"]]) {
      out <- model[["gibbsSamp"]]
      startProgressbar(model[["progressbarLength"]])
      for (i in seq_len(nrow(model[["gibbsSamp"]]))) {
        for (j in seq_len(ncol(model[["gibbsSamp"]]))) {
          out[i, j, , ] <- .cov2cor.callback(model[["gibbsSamp"]][i, j, , ], progressbarTick)
        }
      }
      stateContainer <- .getStateContainerB(jaspResults)
      stateContainer[["gibbsCor"]] <- createJaspState(out, dependencies = c("coefficientType", "averageInterItemCorrelation", "scaleSplithalf"))
    } else {
      return()
    }
  }
  return(out)
}


.BayesianSingleFactorModel <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["singleFactorObj"]]$object))
    return(.getStateContainerB(jaspResults)[["singleFactorObj"]]$object)

  out <- model[["singleFactor"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleOmega"]] && is.null(model[["empty"]])) {

    if (is.null(out[["loadings"]])) {
      startProgressbar(options[["samples"]] * options[["chains"]])

      dataset <- scale(dataset, scale = FALSE)

      jaspBase::.setSeedJASP(options)

      tmp_out <- try(Bayesrel:::omegaSampler(dataset, options[["samples"]], options[["burnin"]],
                                             options[["thinning"]], options[["chains"]],
                                             model[["pairwise"]], progressbarTick,
                                             a0 = options[["inverseGammaPriorShape"]], b0 = options[["inverseGammaPriorScale"]],
                                             m0 = options[["normalPriorMean"]]), silent = TRUE)
      if (model[["pairwise"]] && inherits(tmp_out, "try-error")) {
        .quitAnalysis(gettext("Sampling the posterior factor model for omega failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
      }

      out[["loadings"]] <- tmp_out$lambda
      out[["residuals"]] <- tmp_out$psi
      out[["factor_var"]] <- tmp_out$phi
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["singleFactorObj"]] <- createJaspState(out, dependencies = c("scaleOmega", "inverseGammaPriorShape", "inverseGammaPriorScale",
                                                                                 "normalPriorMean"))
  }

  return(out)
}


.BayesianSingleFactorModelItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["singleFactorItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["singleFactorItemObj"]]$object)

  out <- model[["singleFactorItem"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedOmega"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      return()
    }

    if (is.null(out[["itemLoadings"]])) {

      startProgressbar(options[["samples"]] * options[["chains"]] * ncol(dataset))

      dataset <- scale(dataset, scale = FALSE)

      jaspBase::.setSeedJASP(options)

      out[["itemLoadings"]] <- array(0, c(options[["chains"]],
                                          length(seq(1, options[["samples"]] - options[["burnin"]], options[["thinning"]])),
                                          ncol(dataset), ncol(dataset) - 1))
      out[["itemResiduals"]] <- array(0, c(options[["chains"]],
                                           length(seq(1, options[["samples"]] - options[["burnin"]], options[["thinning"]])),
                                           ncol(dataset), ncol(dataset) - 1))

      for (i in seq_len(ncol(dataset))) {
        tmp_out <- Bayesrel:::omegaSampler(dataset[, -i],
                                           options[["samples"]], options[["burnin"]], options[["thinning"]],
                                           options[["chains"]], model[["pairwise"]], progressbarTick,
                                           a0 = options[["inverseGammaPriorShape"]], b0 = options[["inverseGammaPriorScale"]],
                                           m0 = options[["normalPriorMean"]])
        out[["itemLoadings"]][, , i, ] <- tmp_out$lambda
        out[["itemResiduals"]][, , i, ] <- tmp_out$psi
      }

      dd <- dim(out[["itemLoadings"]])
      out[["itemLoadings"]] <- array(out[["itemLoadings"]], c(dd[1] * dd[2], ncol(dataset), ncol(dataset) - 1))
      out[["itemResiduals"]] <- array(out[["itemResiduals"]], c(dd[1] * dd[2], ncol(dataset), ncol(dataset) - 1))

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["singleFactorItemObj"]] <- createJaspState(out, dependencies = c("itemDeletedOmega", "inverseGammaPriorShape", "inverseGammaPriorScale",
                                                                                     "normalPriorMean"))
  }

  return(out)
}




#### Calculate coefficients ####

.BayesianOmegaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["scaleOmegaObj"]]$object))
    return(.getStateContainerB(jaspResults)[["scaleOmegaObj"]]$object)

  out <- model[["scaleOmega"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleOmega"]] && is.null(model[["empty"]]) && !is.null(model[["singleFactor"]])) {

    if (is.null(out[["samp"]])) {

      startProgressbar(model[["progressbarLength"]])

      if (options[["coefficientType"]] == "unstandardized") {

        out[["samp"]] <- coda::mcmc(.omegaOnArray(model[["singleFactor"]][["loadings"]],
                                                  model[["singleFactor"]][["residuals"]],
                                                  progressbarTick))
      }

      if (options[["coefficientType"]] == "standardized" || options[["standardizedLoadings"]]) {

        lstd <- .stdFactorLoads(model[["singleFactor"]][["loadings"]],
                                model[["singleFactor"]][["residuals"]])
        estd <- 1 - lstd^2

        if (options[["coefficientType"]] == "standardized") {
          out[["samp"]] <- coda::mcmc(.omegaOnArray(lstd, estd, progressbarTick))
        }

        if (options[["standardizedLoadings"]]) {
          out[["loadingsStdSamp"]] <- lstd
        }
      }

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleOmegaObj"]] <- createJaspState(out, dependencies = c("scaleOmega", "coefficientType",
                                                                               "standardizedLoadings", "inverseGammaPriorShape", "inverseGammaPriorScale",
                                                                               "normalPriorMean"))
  }

  return(out)
}


.BayesianOmegaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemDeletedOmegaObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemDeletedOmegaObj"]]$object)

  out <- model[["itemDeletedOmega"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedOmega"]] && is.null(model[["empty"]]) && !is.null(model[["singleFactorItem"]])) {

    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")
      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      if (options[["coefficientType"]] == "unstandardized") {

        out[["itemSamp"]] <- coda::mcmc(.omegaOnArray(model[["singleFactorItem"]][["itemLoadings"]],
                                                      model[["singleFactorItem"]][["itemResiduals"]],
                                                      progressbarTick))
      } else { # standardized
        lstd <- .stdFactorLoads(model[["singleFactorItem"]][["itemLoadings"]],
                                model[["singleFactorItem"]][["itemResiduals"]])
        estd <- 1 - lstd^2

        out[["itemSamp"]] <- coda::mcmc(.omegaOnArray(lstd, estd, progressbarTick))
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemDeletedOmegaObj"]] <- createJaspState(out, dependencies = c("itemDeletedOmega", "coefficientType",
                                                                              "inverseGammaPriorShape", "inverseGammaPriorScale", "normalPriorMean"))
  }

  return(out)
}


.BayesianAlphaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["scaleAlphaObj"]]$object))
    return(.getStateContainerB(jaspResults)[["scaleAlphaObj"]]$object)

  out <- model[["scaleAlpha"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleAlpha"]] && is.null(model[["empty"]])) {

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])

      out[["samp"]] <- coda::mcmc(apply(model[[if (options[["coefficientType"]] == "unstandardized") "gibbsSamp" else "gibbsCor"]],
                                        MARGIN = c(1, 2), Bayesrel:::applyalpha, progressbarTick))

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleAlphaObj"]] <- createJaspState(out, dependencies = c("scaleAlpha", "inverseWishartPriorScale", "inverseWishartPriorDf",
                                                                               "coefficientType"))
  }

  return(out)
}

.BayesianAlphaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemDeletedAlphaObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemDeletedAlphaObj"]]$object)

  out <- model[["itemDeletedAlpha"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedAlpha"]] && is.null(model[["empty"]])) {
    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")
      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesianItemDroppedStats(model[[if (options[["coefficientType"]] == "unstandardized") "gibbsSamp" else "gibbsCor"]],
                                                  Bayesrel:::applyalpha, progressbarTick)


    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemDeletedAlphaObj"]] <- createJaspState(out, dependencies = c("itemDeletedAlpha", "inverseWishartPriorScale", "inverseWishartPriorDf",
                                                                              "coefficientType"))
  }

  return(out)
}



.BayesianLambda2Scale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["scaleLambda2Obj"]]$object))
    return(.getStateContainerB(jaspResults)[["scaleLambda2Obj"]]$object)

  out <- model[["scaleLambda2"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleLambda2"]] && is.null(model[["empty"]])) {

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])

      out[["samp"]] <- coda::mcmc(apply(model[[if (options[["coefficientType"]] == "unstandardized") "gibbsSamp" else "gibbsCor"]],
                                        MARGIN = c(1, 2), Bayesrel:::applylambda2, progressbarTick))

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleLambda2Obj"]] <- createJaspState(out, dependencies = c("scaleLambda2", "inverseWishartPriorScale", "inverseWishartPriorDf",
                                                                                 "coefficientType"))
  }

  return(out)
}

.BayesianLambda2Item <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemDeletedLambda2Obj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemDeletedLambda2Obj"]]$object)

  out <- model[["itemDeletedLambda2"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedLambda2"]] && is.null(model[["empty"]])) {
    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")

      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesianItemDroppedStats(model[[if (options[["coefficientType"]] == "unstandardized") "gibbsSamp" else "gibbsCor"]],
                                                  Bayesrel:::applylambda2, progressbarTick)
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemDeletedLambda2Obj"]] <- createJaspState(out, dependencies = c("itemDeletedLambda2", "inverseWishartPriorScale", "inverseWishartPriorDf",
                                                                                "coefficientType"))
  }

  return(out)
}


.BayesianSplithalfScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["scaleSplithalfObj"]]$object))
    return(.getStateContainerB(jaspResults)[["scaleSplithalfObj"]]$object)

  out <- model[["scaleSplithalf"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleSplithalf"]] && is.null(model[["empty"]]) && !is.null(model[["gibbsCor"]])) {

    nit <- ncol(dataset)
    splits <- split(seq_len(nit), 1:2)
    startProgressbar(model[["progressbarLength"]])
    out[["samp"]] <- matrix(NA, nrow(model[["gibbsCor"]]), ncol(model[["gibbsCor"]]))
    for (i in seq_len(nrow(model[["gibbsCor"]]))) {
      for (j in seq_len(ncol(model[["gibbsCor"]]))) {
        out[["samp"]][i, j] <- .splithalfCor(model[["gibbsCor"]][i, j, , ], splits, progressbarTick)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleSplithalfObj"]] <- createJaspState(out, dependencies = c("scaleSplithalf", "inverseWishartPriorScale", "inverseWishartPriorDf"))
  }

  return(out)
}

.BayesianSplithalfItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemDeletedSplithalfObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemDeletedSplithalfObj"]]$object)

  out <- model[["itemDeletedSplithalf"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedSplithalf"]] && is.null(model[["empty"]])) {
    if (ncol(dataset) < 3) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")

      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesianItemDroppedStats(model[["gibbsCor"]], .splithalfCor, progressbarTick, splithalf = TRUE)
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemDeletedSplithalfObj"]] <- createJaspState(out, dependencies = c("itemDeletedSplithalf", "inverseWishartPriorScale", "inverseWishartPriorDf",
                                                                                       "coefficientType"))
  }

  return(out)
}

.BayesianAverageCor <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["avgCorObj"]]$object))
    return(.getStateContainerB(jaspResults)[["avgCorObj"]]$object)

  out <- model[["average"]]
  if (is.null(out))
    out <- list()

  if (options[["averageInterItemCorrelation"]] && !is.null(model[["gibbsSamp"]])) {

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])
      out[["samp"]] <- matrix(0, nrow(model[["gibbsSamp"]]), ncol(model[["gibbsSamp"]]))
      lowerTriangleIndex <- which(lower.tri(model[["gibbsSamp"]][1, 1, , ]))
      for (i in seq_len(nrow(model[["gibbsSamp"]]))) {
        for (j in seq_len(ncol(model[["gibbsSamp"]]))) {
          corm <- .cov2cor.callback(model[["gibbsSamp"]][i, j, , ], progressbarTick)
          out[["samp"]][i, j] <- mean(corm[lowerTriangleIndex])
        }
      }
      out[["samp"]] <- coda::mcmc(out[["samp"]])

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["avgCorObj"]] <- createJaspState(out, dependencies = c("averageInterItemCorrelation", "inverseWishartPriorScale", "inverseWishartPriorDf"))
  }

  return(out)
}


.BayesianMean <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["meanObj"]]$object))
    return(.getStateContainerB(jaspResults)[["meanObj"]]$object)

  out <- model[["mean"]]
  if (is.null(out))
    out <- list()
  if (options[["scaleMean"]] && is.null(model[["empty"]])) {
    out[["est"]] <- if (options[["meanSdScoresMethod"]] == "sumScores")
      mean(rowSums(dataset, na.rm = TRUE))
    else
      mean(rowMeans(dataset, na.rm = TRUE))

    out[["cred"]] <- c(NA_real_, NA_real_)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("scaleMean", "meanSdScoresMethod"))
  }
  return(out)
}

.BayesianVar <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["varObj"]]$object))
    return(.getStateContainerB(jaspResults)[["varObj"]]$object)

  out <- model[["var"]]
  if (is.null(out))
    out <- list()
  if (options[["scaleVar"]] && is.null(model[["empty"]])) {
    out[["est"]] <- if (options[["meanSdScoresMethod"]] == "sumScores")
      var(rowSums(dataset, na.rm = TRUE))
    else
      var(rowMeans(dataset, na.rm = TRUE))

    out[["cred"]] <- c(NA_real_, NA_real_)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["varObj"]] <- createJaspState(out, dependencies = c("scaleVar", "meanSdScoresMethod"))
  }
  return(out)
}


.BayesianStdDev <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["sdObj"]]$object))
    return(.getStateContainerB(jaspResults)[["sdObj"]]$object)

  out <- model[["sd"]]
  if (is.null(out))
    out <- list()
  if (options[["scaleSd"]] && is.null(model[["empty"]])) {
    out[["est"]] <- if (options[["meanSdScoresMethod"]] == "sumScores")
      sd(rowSums(dataset, na.rm = TRUE))
    else
      sd(rowMeans(dataset, na.rm = TRUE))

    out[["cred"]] <- c(NA_real_, NA_real_)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["sdObj"]] <- createJaspState(out, dependencies = c("scaleSd", "meanSdScoresMethod"))
  }
  return(out)
}



.BayesianItemRestCorrelation <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemRestObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemRestObj"]]$object)

  out <- model[["itemRestCorrelation"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemRestCorrelation"]] && is.null(model[["empty"]])) {

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(options[["samples"]] * options[["chains"]] * ncol(dataset))

      jaspBase::.setSeedJASP(options)
      out[["itemSamp"]] <- .itemRestCorrelation(dataset, options[["samples"]], options[["burnin"]],
                                        options[["thinning"]], options[["chains"]], model[["pairwise"]],
                                        callback = progressbarTick, options[["inverseWishartPriorScale"]],
                                        options[["inverseWishartPriorDf"]])
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = c("itemRestCorrelation", "inverseWishartPriorScale", "inverseWishartPriorDf"))
  }

  return(out)
}

.BayesianMeanItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemMeanObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemMeanObj"]]$object)

  out <- model[["itemMean"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemMean"]] && is.null(model[["empty"]])) {

    out[["itemEst"]] <- colMeans(dataset, na.rm = TRUE)

    out[["itemCred"]] <- matrix(NA_real_, ncol(dataset), 2)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemMeanObj"]] <- createJaspState(out, dependencies = "itemMean")
  }
  return(out)
}

.BayesianVarItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemVarObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemVarObj"]]$object)

  out <- model[["itemVar"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemVar"]] && is.null(model[["empty"]])) {

    out[["itemEst"]] <- apply(dataset, 2, var, na.rm = TRUE)

    out[["itemCred"]] <- matrix(NA_real_, ncol(dataset), 2)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemVarObj"]] <- createJaspState(out, dependencies = "itemVar")
  }
  return(out)
}

.BayesianSdItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemSdObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemSdObj"]]$object)

  out <- model[["itemSd"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemSd"]] && is.null(model[["empty"]])) {

    out[["itemEst"]] <- apply(dataset, 2, sd, na.rm = TRUE)
    out[["itemCred"]] <- matrix(NA_real_, ncol(dataset), 2)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemSdObj"]] <- createJaspState(out, dependencies = "itemSd")
  }
  return(out)
}



.BayesianComputeScaleResults <- function(jaspResults, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["scaleResultsObj"]]$object))
    return(.getStateContainerB(jaspResults)[["scaleResultsObj"]]$object)

  out <- model[["scaleResults"]]
  if (is.null(out))
    out <- list()

  if (is.null(model[["empty"]])) {

    ciValue <- options[["scaleCiLevel"]]

    selected <- names(which(model[["derivedOptions"]][["selectedEstimators"]]))
    # first the coefficients with samples
    sampellist <- model[selected]
    samps <- .sampleListHelper(sampellist, "samp")

    out[["est"]] <- lapply(samps, .getPointEstFun(options[["pointEstimate"]]))
    out[["cred"]] <- lapply(samps, function(x) coda::HPDinterval(coda::mcmc(c(x)), prob = ciValue))

    if (options[["rHat"]]) {
      out[["rHat"]] <- lapply(samps, function(x) {
        # for each samp, (1) convert the rows to a coda object, (2) convert to mcmc.list, so that (3) gelman.diag is happy.
        coda::gelman.diag(coda::as.mcmc.list(lapply(seq_len(nrow(x)), function(i) coda::mcmc(x[i, ]))))[["psrf"]][, 1]
      })
    }
    if (options[["effectiveSampleSize"]]) {
      out[["effectiveSampleSize"]] <- lapply(samps, function(x) {
        sum(coda::effectiveSize(t(x)))
      })
    }

    if (options[["standardizedLoadings"]]) {
      out[["loadingsStd"]] <- apply(model[["scaleOmega"]][["loadingsStdSamp"]], 3, .getPointEstFun(options[["pointEstimate"]]))
    }

    # check for mean var and sd
    if ("scaleMean" %in% selected) {
      out[["est"]][["scaleMean"]] <- model[["scaleMean"]][["est"]]
      out[["cred"]][["scaleMean"]] <- model[["scaleMean"]][["cred"]]
      if (options[["rHat"]]) {
        out[["rHat"]][["scaleMean"]] <- NA_real_
      }
      if (options[["effectiveSampleSize"]]) {
        out[["effectiveSampleSize"]][["scaleMean"]] <- NA_real_
      }
    }
    if ("scaleVar" %in% selected) {
      out[["est"]][["scaleVar"]] <- model[["scaleVar"]][["est"]]
      out[["cred"]][["scaleVar"]] <- model[["scaleVar"]][["cred"]]
      if (options[["rHat"]]) {
        out[["rHat"]][["scaleVar"]] <- NA_real_
      }
      if (options[["effectiveSampleSize"]]) {
        out[["effectiveSampleSize"]][["scaleVar"]] <- NA_real_
      }
    }
    if ("scaleSd" %in% selected) {
      out[["est"]][["scaleSd"]] <- model[["scaleSd"]][["est"]]
      out[["cred"]][["scaleSd"]] <- model[["scaleSd"]][["cred"]]
      if (options[["rHat"]]) {
        out[["rHat"]][["scaleSd"]] <- NA_real_
      }
      if (options[["effectiveSampleSize"]]) {
        out[["effectiveSampleSize"]][["scaleSd"]] <- NA_real_
      }
    }

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleResultsObj"]] <- createJaspState(out, dependencies = c("scaleCiLevel",
                                                                                 "scaleMean", "scaleSd", "scaleVar", "rHat",
                                                                                 "scaleAlpha", "scaleOmega",
                                                                                 "scaleLambda2", "averageInterItemCorrelation",
                                                                                 "meanSdScoresMethod", "inverseWishartPriorScale", "inverseWishartPriorDf",
                                                                                 "inverseGammaPriorShape", "inverseGammaPriorScale",
                                                                                 "coefficientType", "pointEstimate",
                                                                                 "normalPriorMean", "standardizedLoadings",
                                                                                 "effectiveSampleSize"))

  }

  return(out)

}


.BayesianComputeItemResults <- function(jaspResults, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemResultsObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemResultsObj"]]$object)

  out <- model[["itemResults"]]
  if (is.null(out))
    out <- list()

  chosen <- any(model[["derivedOptions"]][["itemDroppedSelected"]])

  if (is.null(model[["empty"]]) && chosen) {

    ciValue <- options[["itemCiLevel"]]

    selected <- names(which(model[["derivedOptions"]][["itemDroppedSelected"]]))
    # first the coefficients with samples
    sampellist <- model[selected]
    samps <- .sampleListHelper(sampellist, "itemSamp")

    if (options[["pointEstimate"]] == "mean") {
      out[["est"]] <- lapply(samps, colMeans)
    } else { # median
      out[["est"]] <- lapply(samps, .colMedians)
    }
    out[["cred"]] <- lapply(samps, function(x) coda::HPDinterval(coda::mcmc(x), prob = ciValue))

    # check for mean and sd
    if ("itemMean" %in% selected) {
      out[["est"]][["itemMean"]] <- model[["itemMean"]][["itemEst"]]
      out[["cred"]][["itemMean"]] <- model[["itemMean"]][["itemCred"]]
    }
    if ("itemVar" %in% selected) {
      out[["est"]][["itemVar"]] <- model[["itemVar"]][["itemEst"]]
      out[["cred"]][["itemVar"]] <- model[["itemVar"]][["itemCred"]]
    }
    if ("itemSd" %in% selected) {
      out[["est"]][["itemSd"]] <- model[["itemSd"]][["itemEst"]]
      out[["cred"]][["itemSd"]] <- model[["itemSd"]][["itemCred"]]
    }

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemResultsObj"]] <- createJaspState(out, dependencies = c("itemDeletedOmega",  "itemDeletedAlpha",
                                                                                "itemDeletedLambda2", "itemCiLevel",
                                                                                "itemRestCorrelation", "itemMean", "itemSd", "itemVar",
                                                                                "inverseWishartPriorScale", "inverseWishartPriorDf", "inverseGammaPriorShape", "inverseGammaPriorScale",
                                                                                "coefficientType", "pointEstimate", "normalPriorMean"))

  }

  return(out)
}


# see https://www.rdocumentation.org/packages/blavaan/versions/0.3-17/topics/blavFitIndices

.BayesianOmegaFitMeasures <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["omegaFitMeasuresObj"]]$object))
    return(.getStateContainerB(jaspResults)[["omegaFitMeasuresObj"]]$object)

  out <- model[["omegaFitMeasures"]]
  if (is.null(out)) {
    out <- list()
  }

  if (options[["scaleOmega"]] && options[["omegaFitMeasures"]] && is.null(model[["empty"]])) {

    n <- nrow(dataset)
    k <- ncol(dataset)
    pstar <- k * (k + 1) / 2 # unique elements in the covariance matrix, variances + covariances
    # cdat <- .rescale(cov(dataset), n) # rescale to n, not n-1
    dataset <- scale(dataset, scale = FALSE)
    cdat <- cov(dataset, use = model[["use.cases"]])
    implieds <- .implCovs(model[["singleFactor"]][["loadings"]],
                          model[["singleFactor"]][["residuals"]],
                          model[["singleFactor"]][["factor_var"]])
    dd <- dim(implieds)
    startProgressbar(dd[1] * dd[2] + # for model implied covariance matrices
                       options[["samples"]] * options[["chains"]] + # for null model sampling
                       dd[1] * dd[2]) # for null model implied covariance matrices

    ### Chisqs ###
    LL1 <- sum(.dmultinorm(dataset, cdat)) # loglikelihood saturated model
    LR_obs <- apply(implieds, c(1, 2), .LRblav, data = dataset, basell = LL1, callback = progressbarTick) # loglikelihoods tested model

    lsm <- apply(model[["singleFactor"]][["loadings"]], 3, mean)
    psm <- apply(model[["singleFactor"]][["residuals"]], 3, mean)
    implM <- tcrossprod(lsm) + diag(psm) # mean implied covariance matrix
    Dtm <- .LRblav(dataset, implM, LL1) # deviance of the mean model implied cov matrix
    out[["lr"]] <- Dtm

    out[["srmr"]] <- .SRMR(cdat, implM)
    ### BRMSEA ###
    Dm <- mean(LR_obs) # mean deviance of the model implied cov matrices
    pD <- Dm - Dtm # effective number of parameters (free parameters)
    out[["rmsea"]] <- .BRMSEA(LR_obs, pstar, pD, n)
    # rmsea_m <- sqrt((Dtm - (pstar - pD)) / ((pstar - pD) * n))

    ### CFI and TLI need a nullmodel:
    res_null <- try(Bayesrel:::omegaSamplerNull(dataset, options[["samples"]], options[["burnin"]],
                                                options[["thinning"]], options[["chains"]],
                                                model[["pairwise"]], progressbarTick,
                                                a0 = options[["inverseGammaPriorShape"]], b0 = options[["inverseGammaPriorScale"]],
                                                m0 = options[["normalPriorMean"]]), silent = TRUE)

    if (model[["pairwise"]] && inherits(res_null, "try-error")) {
      .quitAnalysis(gettext("Sampling the posterior null model failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
    }

    ps_null <- res_null$psi
    ps_null <- apply(ps_null, 3, as.vector)
    impl_null <- array(0, c(nrow(ps_null), k, k))
    for (i in 1:nrow(ps_null)) {
      impl_null[i, , ] <- diag(ps_null[i, ])
    }

    LR_null <- apply(impl_null, 1, .LRblav, data = dataset, basell = LL1, callback = progressbarTick)
    Dm_null <- mean(LR_null)

    psm_null <- colMeans(ps_null)
    implM_null <- diag(psm_null)
    Dtm_null <- .LRblav(dataset, implM_null, LL1)
    pD_null <- Dm_null - Dtm_null

    cfis <- 1 - ((LR_obs - pstar) / (LR_null - pstar))
    tlis <- (((LR_null - pD_null) / (pstar - pD_null)) - ((LR_obs - pD) / (pstar - pD))) /
      (((LR_null - pD_null) / (pstar - pD_null)) - 1)

    cfis[cfis < 0] <- 0
    tlis[tlis < 0] <- 0
    cfis[cfis > 1] <- 1
    tlis[tlis > 1] <- 1

    out[["cfi"]] <- cfis
    out[["tli"]] <- tlis


    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["omegaFitMeasuresObj"]] <- createJaspState(out, dependencies = c("scaleOmega", "inverseGammaPriorShape", "inverseGammaPriorScale",
                                                                                "normalPriorMean", "omegaFitMeasures"))
  }

  return(out)
}





#### Output tables ####

.BayesianScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Bayesian Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("scaleCiLevel", "scaleMean", "scaleSd", "scaleVar", "rHat",
                                  "scaleAlpha", "scaleOmega", "scaleLambda2",
                                  "averageInterItemCorrelation", "meanSdScoresMethod", "effectiveSampleSize"))

  pointEstimate <- gettextf("Posterior %s", options[["pointEstimate"]])

  scaleTable$addColumnInfo(name = "coefficient", title = gettext("Coefficient"), type = "string")
  scaleTable$addColumnInfo(name = "estimate", title = pointEstimate, type = "number")
  ci <- format(100 * options[["scaleCiLevel"]], digits = 3, drop0trailing = TRUE)
  scaleTable$addColumnInfo(name = "lower", title = "Lower", type = "number", overtitle = gettextf("%s%% CI", ci))
  scaleTable$addColumnInfo(name = "upper", title = "Upper", type = "number", overtitle = gettextf("%s%% CI", ci))

  if (options[["rHat"]]) {
    scaleTable$addColumnInfo(name = "rHat", title = gettext("R-hat"), type = "number")
  }
  if (options[["effectiveSampleSize"]]) {
    scaleTable$addColumnInfo(name = "ess", title = gettext("ESS"), type = "number")
  }

  scaleTable$position <- 1
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idxSelected <- which(selected)

  dt <- data.frame(matrix(".", nrow = length(idxSelected), ncol = 0))
  if (any(selected)) {
    for (i in 1:length(idxSelected)) {
      dt$coefficient[i] <- opts[idxSelected[i]]
    }
  }

  # if no coefficients selected or not enough variables:
  if (!.is.empty(model)) {
    dt$estimate <- unlist(model[["scaleResults"]][["est"]], use.names = FALSE)
    dt$lower <- sapply(model[["scaleResults"]][["cred"]], function(x) x[1])
    dt$upper <- sapply(model[["scaleResults"]][["cred"]], function(x) x[2])

    if (options[["rHat"]]) {
      dt$rHat <- unlist(model[["scaleResults"]][["rHat"]], use.names = FALSE)
    }
    if (options[["effectiveSampleSize"]]) {
      dt$ess <- unlist(model[["scaleResults"]][["effectiveSampleSize"]], use.names = FALSE)
    }
  }

  scaleTable$setData(dt)

  if (model[["footnote"]] != "") {
    scaleTable$addFootnote(model[["footnote"]])
  }

  return()
}


.BayesianItemTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["itemTable"]]$object) ||
      !any(model[["derivedOptions"]][["itemDroppedSelected"]]))
    return()

  derivedOptions <- model[["derivedOptions"]]

  itemTable <- createJaspTable(gettext("Bayesian Individual Item Reliability Statistics"))

  itemTable$dependOn(options = c("itemDeletedOmega",  "itemDeletedAlpha",  "itemDeletedLambda2",
                                 "itemCiLevel", "itemRestCorrelation", "itemMean", "itemSd",
                                 "scaleAlpha", "scaleOmega", "scaleLambda2"))

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemTable$position <- 2
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  overTitles <- format(derivedOptions[["namesEstimators"]][["tables_item"]], digits = 3, drop0trailing = TRUE)
  overTitles <- gettextf("%s (if item dropped)", overTitles)
  cred <- format(100 * options[["itemCiLevel"]], digits = 3, drop0trailing = TRUE)

  selected <- derivedOptions[["itemDroppedSelected"]]
  idxSelected <- which(selected)
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]]
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]

  itemTable[["variable"]] <- model[["itemsDropped"]]

  footnote <- ""
  if (!is.null(unlist(options[["reverseScaledItems"]])))
    footnote <- .addFootnoteReverseScaledItems(options)

  pointEstimate <- gettextf("Posterior %s", options[["pointEstimate"]])

  fewItemProblem <- FALSE
  for (i in idxSelected) {
    if (estimators[i] %in% coefficients) {
      if (estimators[i] == "Item-rest correlation") { # no item deleted for item rest cor
        itemTable$addColumnInfo(name = paste0("postMean", i), title = pointEstimate, type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
      } else {
        itemTable$addColumnInfo(name = paste0("postMean", i), title = pointEstimate, type = "number",
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", cred), type = "number",
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", cred), type = "number",
                                overtitle = overTitles[i])
        fewItemProblem <- TRUE
      }
    } else {
      itemTable$addColumnInfo(name = paste0("postMean", i), title = estimators[i], type = "number")
    }
  }


  if (is.null(model[["empty"]])) {
    tb <- data.frame(variable = model[["itemsDropped"]])
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (estimators[i] %in% coefficients) { # reliability coefficient
        rows <- length(options[["variables"]])
        if (rows < 3 && nm != "itemRestCorrelation") {
          newtb <- cbind(postMean = rep(NaN, rows), matrix(NaN, rows, 2, dimnames = list(NULL, c("lower", "upper"))))
        } else {
          newtb <- cbind(postMean = model[["itemResults"]][["est"]][[nm]], model[["itemResults"]][["cred"]][[nm]])
        }
      } else {
        newtb <- cbind(postMean = model[["itemResults"]][["est"]][[nm]])
      }
      colnames(newtb) <- paste0(colnames(newtb), i)
      tb <- cbind(tb, newtb)
    }
    itemTable$setData(tb)

    if (length(options[["variables"]]) < 3 && fewItemProblem) {
      footnote <- gettextf("%s Please enter at least 3 variables for the if item dropped statistics.", footnote)
    }
  }

  if (footnote != "") {
    itemTable$addFootnote(footnote)
  }

  itemTable$position <- 2
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  return()
}



.BayesianProbTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["probabilityTable"]]$object))
    return()

  # check if the values are in the proper order, meaning the lower field value is smaller than the upper field
  if (options[["probabilityTableLowerBound"]] > options[["probabilityTableUpperBound"]]) {
    low <- options[["probabilityTableUpperBound"]]
    high <- options[["probabilityTableLowerBound"]]
    footnote <- gettext("The bounds you entered have been rearranged in increasing order to provide meaningful results.")
  } else {
    low <- options[["probabilityTableLowerBound"]]
    high <- options[["probabilityTableUpperBound"]]
    footnote <- ""
  }
  probabilityTable <- createJaspTable(
    gettextf("Probability that Reliability Coefficient is Larger than %.2f and Smaller than %.2f", low, high))
  probabilityTable$dependOn(options = c("probabilityTableLowerBound", "probabilityTable", "probabilityTableUpperBound"))

  overTitle <- gettext("Probability")
  probabilityTable$addColumnInfo(name = "statistic", title = gettext("Coefficient"), type = "string")
  probabilityTable$addColumnInfo(name = "prior", title = gettext("Prior"), type = "number", overtitle = overTitle)
  probabilityTable$addColumnInfo(name = "posterior", title = gettext("Posterior"), type = "number", overtitle = overTitle)

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  idxSelected  <- which(selected)

  if (options[["probabilityTable"]] && is.null(model[["empty"]])) {

    n.item <- model[["k"]]

    end <- 512
    xx <- seq(0, 1, length.out = 512)
    poslow <- end - sum(xx > low)
    poshigh <- end - sum(xx > high)

    probsPost <- numeric(sum(selected))
    probsPrior <- numeric(sum(selected))

    for (i in seq_len(length(idxSelected))) {
      nm <- names(idxSelected[i])

      samp_tmp <- as.vector(model[[nm]][["samp"]])
      probsPost[i] <- mean(samp_tmp > low) - mean(samp_tmp > high)

      if (nm == "scaleOmega") {
        startProgressbar(2e3)
      } else {
        startProgressbar(4e3)
      }
      prior <- .samplePrior(n.item, nm, progressbarTick, options[["inverseWishartPriorScale"]], options[["inverseWishartPriorDf"]],
                            options[["inverseGammaPriorShape"]], options[["inverseGammaPriorScale"]], options[["normalPriorMean"]])

      probsPrior[i] <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]]) -
        sum(prior[["y"]][poshigh:end]) / sum(prior[["y"]])

    }

    df <- data.frame(statistic = opts[idxSelected], prior = probsPrior, posterior = probsPost)
    probabilityTable$setData(df)

    if (footnote != "") {
      probabilityTable$addFootnote(footnote)
    }

    probabilityTable$position <- 3
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["probabilityTable"]] <- probabilityTable
  }

  return()
}


# use posterior mean and median in the table, also flip the table

.BayesianOmegaFitMeasuresTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["fitTable"]]$object))
    return()

  fitTable <- createJaspTable(gettextf("Fit Measures for the Single-Factor Model"))

  fitTable$dependOn(options = c("scaleOmega", "omegaFitMeasures", "omegaFitMeasuresCutoffRmsea", "omegaFitMeasuresCutoffCfiTli", "pointEstimate",
                                "omegaFitMeasuresCiLevel"))

  fitTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  fitTable$addColumnInfo(name = "lr", title = "B-LR", type = "number")
  fitTable$addColumnInfo(name = "srmr", title = "B-SRMR", type = "number")
  fitTable$addColumnInfo(name = "rmsea", title = "B-RMSEA", type = "number")
  fitTable$addColumnInfo(name = "cfi", title = "B-CFI", type = "number")
  fitTable$addColumnInfo(name = "tli", title = "B-TLI", type = "number")

  cred <- gettextf("%s%% CI",
                   format(100 * options[["omegaFitMeasuresCiLevel"]], digits = 3, drop0trailing = TRUE))
  intervalLow <- gettextf("%s lower bound", cred)
  intervalUp <- gettextf("%s upper bound", cred)

  pointEstimate <- gettextf("Posterior %s", options[["pointEstimate"]])

  if (options[["scaleOmega"]] && options[["omegaFitMeasures"]] && is.null(model[["empty"]])) {

    pointEstimates <- vapply(model[["omegaFitMeasures"]], .getPointEstFun(options[["pointEstimate"]]), numeric(1))

    creds <- vapply(model[["omegaFitMeasures"]][3:5], function(x) coda::HPDinterval(coda::mcmc(c(x))), numeric(2))
    creds <- cbind(matrix(NA_real_, 2, 2), creds) # no entry for LR and SRMR

    cutoffs_saturated <- vapply(model[["omegaFitMeasures"]][3], function(x) mean(x < options[["omegaFitMeasuresCutoffRmsea"]]),
                                numeric(1))
    cutoffs_null <- vapply(model[["omegaFitMeasures"]][-(1:3)], function(x) mean(x > options[["omegaFitMeasuresCutoffCfiTli"]]),
                           numeric(1))
    cutoffs <- c(NA_real_, NA_real_, cutoffs_saturated, cutoffs_null) # no entry for LR and SRMR

    df <- data.frame(estimate = c(pointEstimate, intervalLow, intervalUp, "Relative to cutoff"))
    dfnums <- rbind(pointEstimates, creds, cutoffs)
    df <- cbind(df, dfnums)
    fitTable$setData(df)

    fitTable$addFootnote(gettext("'Relative to cutoff'-row denotes the probability that the B-RMSEA is smaller than
                         the corresponding cutoff and the probabilities that the B-CFI/TLI are larger than the
                         corresponding cutoff."))

    fitTable$position <- 4
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["fitTable"]] <- fitTable
  }

  return()
}



.BayesianLoadingsTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["loadTable"]]$object))
    return()

  loadTable <- createJaspTable(gettext("Standardized Loadings of the Single-Factor Model"))

  loadTable$dependOn(options = c("scaleOmega", "standardizedLoadings"))

  loadTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")
  loadTable$addColumnInfo(name = "loadings", title = gettext("Standardized loading"), type = "number")

  derivedOptions <- model[["derivedOptions"]]

  if (options[["scaleOmega"]] && options[["standardizedLoadings"]] && is.null(model[["empty"]])) {



    df <- data.frame(
      variable = options[["variables"]],
      loadings = model[["scaleResults"]][["loadingsStd"]])
    loadTable$setData(df)

    loadTable$position <- 5
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["loadTable"]] <- loadTable
  }

  return()
}



#### Output plots ####

.BayesianPosteriorPlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainer"]]$object))
    return()

  plotContainer <- createJaspContainer(gettext("Posterior Plots"))
  plotContainer$dependOn(options = c("posteriorPlot", "posteriorPlotShaded", "probabilityTable", "probabilityTableLowerBound",
                                     "probabilityTableUpperBound", "posteriorPlotFixedRange", "posteriorPlotPriorDisplayed", "scaleCiLevel",
                                     "scaleAlpha", "scaleOmega", "scaleLambda2"))

  derivedOptions <- model[["derivedOptions"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  idxSelected   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

  if (options[["posteriorPlotShaded"]] && options[["probabilityTable"]]) {
    if (options[["probabilityTableLowerBound"]] > options[["probabilityTableUpperBound"]]) {
      low <- options[["probabilityTableUpperBound"]]
      high <- options[["probabilityTableLowerBound"]]
    } else {
      low <- options[["probabilityTableLowerBound"]]
      high <- options[["probabilityTableUpperBound"]]
    }
    posteriorPlotShaded <- c(low, high)
  } else {
    posteriorPlotShaded <- NULL
  }

  if (options[["posteriorPlot"]] && is.null(model[["empty"]])) {
    n.item <- model[["k"]]

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (is.null(plotContainer[[nmsObjsNoGreek[i]]])) {
        if (options[["posteriorPlotPriorDisplayed"]]) {

          if (nm == "scaleOmega") {
            startProgressbar(2e3)
          } else {
            startProgressbar(4e3)
          }
          prior <- .samplePrior(n.item, nm, progressbarTick, options[["inverseWishartPriorScale"]], options[["inverseWishartPriorDf"]],
                                options[["inverseGammaPriorShape"]], options[["inverseGammaPriorScale"]], options[["normalPriorMean"]])
        } else {
          prior <- NULL
        }

        p <- .makeSinglePosteriorPlot(model[[nm]], model[["scaleResults"]][["cred"]][[nm]], nmsLabs[[i]],
                                      options[["posteriorPlotFixedRange"]], posteriorPlotShaded,
                                      options[["posteriorPlotPriorDisplayed"]], prior)
        plotObj <- createJaspPlot(plot = p, title = nmsObjs[i])
        plotObj$dependOn(options = names(idxSelected[i]))
        plotObj$position <- i
        plotContainer[[nmsObjsNoGreek[i]]] <- plotObj

      }
    }
    plotContainer$position <- 6
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["plotContainer"]] <- plotContainer
  }

  return()
}

.makeSinglePosteriorPlot <- function(coefList, coefResults, nms, posteriorPlotFixedRange, shade = NULL, priorTrue, priorSample) {

  # TODO: consider precomputing all densities (maybe with kernsmooth?) and reducing memory that way

  samp_tmp <- as.vector(coefList[["samp"]])
  if (posteriorPlotFixedRange) {
    d <- stats::density(samp_tmp, from = 0, to = 1, n = 2^10)
  } else {
    d <- stats::density(samp_tmp, n = 2^10)
  }
  datDens <- data.frame(x = d$x, y = d$y)
  datPrior <- data.frame(x = priorSample$x, y = priorSample$y)


  xBreaks <- jaspGraphs::getPrettyAxisBreaks(datDens$x)

  # max height posterior is at 90% of plot area; remainder is for credible interval
  ymax <- max(d$y) / .9
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, ymax))
  ymax <- max(yBreaks)
  scaleCriRound <- round(coefResults, 3)
  datCri <- data.frame(xmin = scaleCriRound[1L], xmax = scaleCriRound[2L], y = .925 * ymax)
  height <- (ymax - .925 * ymax) / 2
  if (posteriorPlotFixedRange) {
    if (datCri$xmin == datCri$xmax) { # if two zeros, the interval is merged together
      datTxt <- data.frame(x = c(datCri$xmin, datCri$xmax),
                           y = 0.985 * ymax,
                           label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
                           stringsAsFactors = FALSE)
    } else {
      datTxt <- data.frame(x = c(datCri$xmin - .08, datCri$xmax + .08),
                           y = 0.985 * ymax,
                           label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
                           stringsAsFactors = FALSE)
    }

  } else {
    datTxt <- data.frame(x = c(datCri$xmin, datCri$xmax),
                         y = 0.985 * ymax,
                         label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
                         stringsAsFactors = FALSE)
  }


  if (datCri$xmin[1L] < 0) {
    datCri$xmin[1L] <- 0
    datTxt$x[1L] <- 0
    datTxt$label[1L] <- "< 0"
  } else if (datCri$xmin[1L] == 0) {
    datCri$xmin[1L] <- 0
    datTxt$x[1L] <- 0
    datTxt$label[1L] <- "0"
  }

  # if the bounds are less than 0.05 away from 0 or 1, expand the axis by 0.1 so the credible interval text does not
  # get chopped off.
  xExpand <- .1 * ((c(0, 1) - datTxt$x) <= 0.05)
  if (posteriorPlotFixedRange && max(datDens$y, na.rm = T) >= 1000) xExpand <- c(xExpand[1], .05)
  # with large numebrs on the y-axis, the x-axis labels to the right get cut off sometimes, when the range is fixed

  g <- ggplot2::ggplot(data = datDens, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(size = .85) +
    ggplot2::geom_errorbarh(data = datCri, mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y = y),
                            height = height, inherit.aes = FALSE) +
    ggplot2::geom_text(data = datTxt, mapping = ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
                       size = 5) +
    ggplot2::scale_y_continuous(name = gettext("Density"), breaks = yBreaks, limits = range(yBreaks)) +
    ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand, limits = range(xBreaks))

  if (!is.null(shade)) {
    datFilter <- datDens[datDens[["x"]] >= shade[1] & datDens[["x"]] <= shade[2], ]
    if (length(datFilter$x) == 0)
      datFilter <- data.frame(x = 0, y = 0)
    g <- g + ggplot2::geom_ribbon(data = datFilter, mapping = ggplot2::aes(ymin = 0, ymax = y),
                                  fill = "grey", alpha = 0.95) +
      ggplot2::geom_line(size = .85)
  }

  if (priorTrue) {
    g <- g + ggplot2::geom_line(data = datPrior, mapping = ggplot2::aes(x = x, y = y),
                                linetype = "dashed", size = .85) +
      ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, limits = range(xBreaks),
                                  expand = xExpand)

  }

  return(jaspGraphs::themeJasp(g))

}



.BayesianIfItemPlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainerItem"]]$object))
    return()

  plotContainerItem <- createJaspContainer(gettext("If Item Dropped Posterior Plots"))
  plotContainerItem$dependOn(options = c("variables", "itemDeletedPlot",
                                         "itemCiLevel", "itemDeletedPlotOrderedType", "itemDeletedPlotOrdered",
                                         "itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2"))

  derivedOptions <- model[["derivedOptions"]]

  selected <- derivedOptions[["itemDroppedSelectedItem"]]
  idxSelected   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables_item"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

  if (options[["itemDeletedPlotOrdered"]]) {
    ordering <- options[["itemDeletedPlotOrderedType"]]
  } else {
    ordering <- NULL
  }

  if (is.null(model[["empty"]]) && options[["itemDeletedPlot"]]) {
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (is.null(plotContainerItem[[nmsObjsNoGreek[i]]])) {

        if (length(options[["variables"]]) < 3) {
          plotObjItem <- createJaspPlot(plot = NULL, title = nmsObjs[i], width = 400)
          plotObjItem$setError(gettext("Please enter at least 3 variables for the if item dropped plot"))
        } else {
          # use the coefficient item object in the model list, and the object directly before it,
          # which is always the corresponding scale object
          prevNumber <- which(names(model) == nm) - 1
          name <- unlist(strsplit(nm, "Deleted"))[2]
          coefPos <- grep(name, names(model[["scaleResults"]][["est"]]))

          p <- .makeIfItemPlot(model[[nm]], model[[prevNumber]],
                               model[["itemResults"]][["est"]][[nm]],
                               model[["scaleResults"]][["est"]][[coefPos]],
                               nmsLabs[[i]],
                               options[["itemCiLevel"]],
                               ordering = ordering, model[["itemsDropped"]])
          plotObjItem <- createJaspPlot(plot = p, title = nmsObjs[i], width = 400)
          plotObjItem$dependOn(options = names(idxSelected[i]))
          plotObjItem$position <- i
          if (is.null(p)) {
            plotObjItem$setError(gettext("KLD ordering failed because two variables have perfect correlation"))
          }
        }
        plotContainerItem[[nmsObjsNoGreek[i]]] <- plotObjItem

      }
    }
    plotContainerItem$position <- 7
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["plotContainerItem"]] <- plotContainerItem
  }

  return()
}

.makeIfItemPlot <- function(coefItem, coefScale, coefItemEst, coefScaleEst, nms, int, ordering, variables) {
  n_row <- length(variables)
  lower <- (1 - int) / 2
  upper <- int + (1 - int) / 2

  samp_tmp <- as.vector(coefScale[["samp"]])
  dat <- data.frame(as.matrix(samp_tmp), row.names =  NULL)
  names(dat) <- "value"
  dat$colos <- "1"
  dat$var <- "original"

  dat_del <- t(as.matrix(as.data.frame(coefItem[["itemSamp"]])))
  names <- jaspBase::decodeColNames(variables)

  for (i in n_row:1) {
    tmp <- as.data.frame(dat_del[i, ])
    colnames(tmp) <- "value"
    tmp$var <- names[i]
    tmp$colos <- "2"
    dat <- rbind(dat, tmp)
  }
  dat$var <- factor(dat$var, levels = unique(dat$var))

  if (!is.null(ordering)) {
    est <- as.data.frame(coefItemEst)
    est[n_row + 1, ] <- 1
    colnames(est) <- "value"
    est$name <- c(names, "original")

    if (ordering == "mean") {
      dists <- abs(coefScaleEst - coefItemEst)
      dists[length(dists) + 1] <- 0
      est <- est[order(dists, decreasing = FALSE), ]
      dat$var <- factor(dat$var, levels = c(est$name))

    } else if (ordering == "kullbackLeibler") {
      samps <- coefItem[["itemSamp"]]
      og_samp <- samp_tmp

      dists <- try(apply(samps, 2, .KLD.statistic, y = og_samp)) # kl divergence
      ### when there are only three variables and two of them have almost perfect correlation, KLD fails
      if (any(round(cor(samps)[lower.tri(cor(samps))], 3) == 1) && inherits(dists, "try-error")) {
        return(NULL)
      } else {
        dists[length(dists) + 1] <- 0
        est <- est[order(dists), ]
        dat$var <- factor(dat$var, levels = c(est$name))
      }

    } else if (ordering == "kolmogorovSmirnov") {
      samps <- coefItem[["itemSamp"]]
      og_samp <- samp_tmp
      dists <- apply(samps, 2, .ks.test.statistic, y = og_samp) # ks distance
      dists[length(dists) + 1] <- 0
      est <- est[order(dists), ]
      dat$var <- factor(dat$var, levels = c(est$name))
    }
  }

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = value, y = var, fill = colos)) +
    ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = c(lower, 0.5, upper),
                                  alpha = .85, show.legend = FALSE, scale = 1) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(colour = "black")) +
    ggplot2::ylab(gettext("Item Dropped")) +
    ggplot2::xlab(nms) +
    ggplot2::scale_fill_grey() +
    ggplot2::scale_y_discrete(expand = ggplot2::expand_scale(mult = c(0.1, 0.25)))

  return(jaspGraphs::themeJasp(g))

}



.omegaPosteriorPredictive <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["omegaPPC"]]$object))
    return()

  if (options[["omegaPosteriorPredictiveCheck"]] && options[["scaleOmega"]] && is.null(model[["empty"]])) {

    ll <- apply(model[["singleFactor"]][["loadings"]], 3, as.vector)
    rr <- apply(model[["singleFactor"]][["residuals"]], 3, as.vector)
    phi <- c(model[["singleFactor"]][["factor_var"]])

    cobs <- model[["data_cov"]]

    k <- ncol(cobs)
    nsamp <- nrow(ll)
    ee_impl <- matrix(0, nsamp, k)
    for (i in seq_len(nsamp)) {
      ctmp <- ll[i, ] %*% t(phi[i]) %*% t(ll[i, ]) + diag(rr[i, ])
      dtmp <- MASS::mvrnorm(model[["n"]], rep(0, k), ctmp)
      ee_impl[i, ] <- eigen(cov(dtmp), only.values = TRUE)$values
    }

    eframe <- data.frame(number = seq(1, k), eigen_value = eigen(cobs)$values)
    eframe$eigen_sim_low <- apply(ee_impl, 2, quantile, prob = .025)
    eframe$eigen_sim_up <- apply(ee_impl, 2, quantile, prob = .975)

    yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(eframe$eigen_sim_up)))

    g <- ggplot2::ggplot(eframe, mapping = ggplot2::aes(x = number, y = eigen_value)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = eigen_sim_low, ymax = eigen_sim_up), color = "grey55", width = 0.2,
                             size = 1) +
      ggplot2::geom_point(size = 2.25) +
      ggplot2::xlim(c(1, k)) +
      ggplot2::scale_y_continuous(name = gettext("Eigenvalue"), breaks = yBreaks, limits = range(yBreaks)) +
      ggplot2::scale_x_continuous(name = gettext("Eigenvalue No."),
                                  breaks = seq(1, k),
                                  expand = ggplot2::expand_scale(mult = c(.1, .1)))

    g <- jaspGraphs::themeJasp(g)

    plot <- createJaspPlot(plot = g, title = "Posterior Predictive Check Omega", width = 350)
    plot$dependOn(options = c("omegaPosteriorPredictiveCheck", "scaleOmega", "coefficientType"))

    plot$position <- 8
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["omegaPPC"]] <- plot
  }

  return()
}



.BayesianTracePlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainerTP"]]$object))
    return()

  if (is.null(model[["empty"]]) && options[["tracePlot"]]) {

    plotContainerTP <- createJaspContainer(gettext("Convergence Traceplot"))
    plotContainerTP$dependOn(options = c("tracePlot", "scaleAlpha", "scaleOmega", "scaleLambda2"))

    derivedOptions <- model[["derivedOptions"]]
    selected <- derivedOptions[["selectedEstimatorsPlots"]]
    idxSelected   <- which(selected)
    nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
    nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]]
    nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (is.null(plotContainerTP[[nmsObjsNoGreek[i]]])) {

        p <- .makeTracePlot(model[[nm]], nmsLabs[[i]])
        plotObjTP <- createJaspPlot(plot = p, title = nmsObjs[i], width = 400)
        plotObjTP$dependOn(options = names(idxSelected[i]))
        plotObjTP$position <- i
        plotContainerTP[[nmsObjsNoGreek[i]]] <- plotObjTP

      }
    }

    plotContainerTP$position <- 9
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["plotContainerTP"]] <- plotContainerTP

  }

  return()
}


.makeTracePlot <- function(coefList, nms) {

  dd <- coefList[["samp"]]
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, length(dd[1, ])))

  dv <- cbind(dd[1, ], 1, seq(1, ncol(dd)))
  for (j in 2:nrow(dd)) {
    dv <- rbind(dv, cbind(dd[j, ], j, seq(1, ncol(dd))))
  }
  dat <- data.frame(dv)
  colnames(dat) <- c("Value", "chain", "Iterations")
  dat$chain <- as.factor(dat$chain)

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = Iterations, y = Value, colour = chain)) +
    ggplot2::geom_line(size = .3) +
    ggplot2::ylab(nms) +
    ggplot2::scale_x_continuous(name = gettext("Iterations"), breaks = xBreaks,
                                limits = range(xBreaks),
                                expand = ggplot2::expand_scale(mult = c(0.05, 0.1)))

  return(jaspGraphs::themeJasp(g))

}





#### Help functions ####


.getStateContainerB <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = c("variables", "reverseScaledItems",
                                                                          "samples", "burnin", "thinning",
                                                                          "chains", "naAction", "setSeed",
                                                                          "seed", "samplesSavingDisabled")
  )

  return(jaspResults[["stateContainer"]])
}


.summarizePosteriorItems <- function(samples, ciValue) {
  return(list(
    colMeans(samples),
    coda::HPDinterval(coda::mcmc(samples), prob = ciValue)
  ))
}


.samplePrior <- function(k, estimate, callback = function(){}, k0, df0, a0, b0, m0) {

  n_samp <- 2e3

  if (estimate == "scaleOmega") {
    H0 <- 1 # prior multiplier matrix for lambdas variance
    l0k <- rep(m0, k) # prior lambdas
    a0k <- a0 # prior gamma function for psis
    b0k <- b0 # prior gamma for psi
    prioromega <- numeric(n_samp)
    for (i in seq_len(n_samp)) {
      invpsi <- rgamma(k, a0k, b0k)
      psi <- 1 / invpsi
      lambda <- rnorm(k, l0k, sqrt(psi * H0))
      prioromega[i] <- Bayesrel:::omegaBasic(lambda, psi)
      callback()
    }
    out <- density(prioromega, from = 0, to = 1, n = 512)
    return(out)
  }

  v0 <- df0
  k0 <- k0
  t <- diag(k)
  T0 <- solve(t / k0)
  m <- array(0, c(n_samp, k, k))

  for (i in seq_len(n_samp)) {
    m[i, , ] <- LaplacesDemon::rinvwishart(v0, T0)
    callback()
  }

  if (estimate == "scaleAlpha") {
    prioralpha <- apply(m, MARGIN = 1, Bayesrel:::applyalpha, callback)
    out <- density(prioralpha, from = 0, to = 1, n = 512)
    return(out)

  }
  if (estimate == "scaleLambda2") {
    priorlambda2 <- apply(m, MARGIN = 1, Bayesrel:::applylambda2, callback)
    out <- density(priorlambda2, from = 0, to = 1, n = 512)
    return(out)

  }

  if (estimate == "scaleSplithalf") {
    nit <- k
    splits <- split(seq_len(nit), 1:2)
    priorsplithalf<- apply(m, MARGIN = 1, .splithalfCor, splits = splits, callback)
    out <- density(priorsplithalf, from = 0, to = 1, n = 512)
    return(out)

  }

}

.BayesianItemDroppedStats <- function(cov_samp, f1 = function(){}, callback = function(){}, splithalf = FALSE) {

  dd <- dim(cov_samp)
  nit <- dd[3] - 1
  splits <- split(seq_len(nit), 1:2)
  out <- matrix(0, dd[1] * dd[2], dd[3])
  cov_samp <- array(cov_samp, c(dd[1] * dd[2], dd[3], dd[3]))
  if (splithalf) {
    for (i in seq_len(dd[3])) {
      out[, i] <- apply(cov_samp[, -i, -i], c(1), f1, callback = callback, splits = splits)
    }
  } else {
    for (i in seq_len(dd[3])) {
      out[, i] <- apply(cov_samp[, -i, -i], c(1), f1, callback = callback)
    }
  }

  return(out)
}


.itemRestCorrelation <- function(dataset, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0, df0) {

  ircor_samp <- matrix(0, n.chains * length(seq(1, n.iter - n.burnin, thin)), ncol(dataset))
  for (i in seq(ncol(dataset))) {
    help_dat <- cbind(as.matrix(dataset[, i]), rowSums(as.matrix(dataset[, -i]), na.rm = TRUE))
    ircor_samp[, i] <- .WishartCorTransform(help_dat, n.iter = n.iter, n.burnin = n.burnin, thin = thin,
                                            n.chains = n.chains, pairwise = pairwise, callback = callback, k0, df0)
  }

  return(ircor_samp)
}

.WishartCorTransform <- function(x, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0, df0) {
  tmp_cov <- Bayesrel:::covSamp(x, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0 = k0, df0 = df0)$cov_mat
  dd <- dim(tmp_cov)
  tmp_cov <- array(tmp_cov, c(dd[1] * dd[2], dd[3], dd[4]))
  tmp_cor <- apply(tmp_cov, c(1), cov2cor)
  out <- tmp_cor[2, ]
  callback()
  return(out)
}

.colMedians <- function(x) {
  return(apply(x, 2, median))
}


.stdFactorLoads <- function(ll, ee, callback = function(){}) {
  ds <- dim(ll)
  out <- ll
  for (i in seq_len(ds[1])) {
    for (j in seq_len(ds[2])) {
      implV <- diag(ll[i, j, ] %*% t(ll[i, j, ]) + diag(ee[i, j, ]))
      out[i, j, ] <- ll[i, j, ] / sqrt(implV)
      callback()
    }
  }
  return(out)
}

.omegaOnArray <- function(ll, ee, callback = function(){}) {
  ds <- dim(ll)
  out <- matrix(0, ds[1], ds[2])
  for (i in seq_len(ds[1])) {
    for (j in seq_len(ds[2])) {
      out[i, j] <- Bayesrel:::omegaBasic(ll[i, j, ], ee[i, j, ])
      callback()
    }
  }
  return(out)
}

.implCovs <- function(ll, ee, pp, callback = function(){}) {
  ds <- dim(ll)
  out <- array(0, c(ds[1], ds[2], ds[3], ds[3]))
  for (i in seq_len(ds[1])) {
    for (j in seq_len(ds[2])) {
      out[i, j, , ] <- ll[i, j, ] %*% t(pp[i, j]) %*% t(ll[i, j, ]) + diag(ee[i, j, ])
      callback()
    }
  }
  return(out)
}


.SRMR <- function(cdat, impl) {
  nvar <- ncol(cdat)
  e <- (nvar * (nvar + 1)) / 2
  sqrt.d <- 1 / sqrt(diag(cdat))
  D <- diag(sqrt.d, ncol = length(sqrt.d))
  R <- D %*% (cdat - impl) %*% D
  srmr <- sqrt(sum(R[lower.tri(R, diag = TRUE)]^2) / e)
  return(srmr)
}


.LRblav <- function(data, cmat, basell, callback = function(){}) {
  tmpll <- .dmultinorm(data, cmat)
  out <- 2 * (basell - sum(tmpll))
  callback()
  return(out)
}


.BRMSEA <- function(chisq, p, pD, n) {
  dChisq <- (chisq - p)
  dChisq[dChisq < 0] <- 0
  rmsea <- sqrt(dChisq / ((p - pD) * n))
  return(rmsea)
}


# borrowed that from mnormt package
.dmultinorm <- function(x, varcov, mm = 0) {
  d <- ncol(varcov)
  X <- t(x - mm)
  varcov <- (varcov + t(varcov))/2
  u <- chol(varcov, pivot = FALSE)
  inv <- chol2inv(u)
  logdet <- 2 * sum(log(diag(u)))
  Q <- colSums((inv %*% X) * X)
  Q <- Q[!is.na(Q)]
  logPDF <- as.vector(Q + d * logb(2 * pi) + logdet)/(-2)
  return(logPDF)
}

.getPointEstFun <- function(pointEstimateFunString) {
  return(switch(pointEstimateFunString,
                "mean"   = mean,
                "median" = median,
                stop("getPointEstFun has no case for value: ", pointEstimateFunString))
  )
}

# change options when scale box is unchecked
.scaleItemBoxAlignB <- function(options) {
  opts <- options
  if (!options[["scaleOmega"]])
    opts[["itemDeletedOmega"]] <- FALSE
  if (!options[["scaleAlpha"]])
    opts[["itemDeletedAlpha"]] <- FALSE
  if (!options[["scaleLambda2"]])
    opts[["itemDeletedLambda2"]] <- FALSE

  return(opts)

}
