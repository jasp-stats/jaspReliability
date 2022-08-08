
#' @importFrom jaspBase createJaspContainer createJaspHtml createJaspPlot createJaspQmlSource createJaspState createJaspTable
#' progressbarTick .quitAnalysis .readDataSetToEnd startProgressbar

#' @importFrom stats approxfun cor cov cov2cor density integrate median na.omit pnorm qchisq qnorm quantile rgamma rnorm sd var

#' @export
unidimensionalReliabilityBayesian <- function(jaspResults, dataset, options) {

  options <- jaspBase::.parseAndStoreFormulaOptions(jaspResults, options, "iwScale")

  dataset <- .readData(dataset, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }


  .checkErrors(dataset, options, Bayes = TRUE)

  model <- .BayesianPreCalc(jaspResults, dataset, options)

  options <- .scaleItemBoxAlign(options)


  model[["derivedOptions"]] <- .BayesianDerivedOptions(options)
  model[["gibbsCor"]] <- .BayesianStdCov(jaspResults, dataset, options, model)
  model[["singleFactor"]] <- .BayesianSingleFactorModel(jaspResults, dataset, options, model)
  model[["singleFactorItem"]] <- .BayesianSingleFactorModelItem(jaspResults, dataset, options, model)
  model[["omegaScale"]] <- .BayesianOmegaScale(jaspResults, dataset, options, model)
  model[["omegaItem"]] <- .BayesianOmegaItem(jaspResults, dataset, options, model)
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

  model[["fitMeasures"]] <- .BayesianFitMeasures(jaspResults, dataset, options, model)

  .BayesianScaleTable(jaspResults, model, options)
  .BayesianItemTable(jaspResults, model, options)
  .BayesianProbTable(jaspResults, model, options)
  .BayesianLoadingsTable(jaspResults, model, options)
  .BayesianFitMeasuresTable(jaspResults, model, options)
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


# -------------------------------------------
#       Bayesian precalulate results
# -------------------------------------------


.BayesianPreCalc <- function(jaspResults, dataset, options) {

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

    # when listwise deletion is chosen the values are deleted upon reading in the data,
    # before entering this whole analysis, without this we would never know the initial
    # size of the data
    tmp <- .readDataSetToEnd(columns.as.numeric = unlist(options[["variables"]]))
    old_n <- nrow(tmp)

    if (nrow(dataset) < old_n) { # this indicates listwise deletion
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
      (options[["alphaScale"]] || options[["lambda2Scale"]] || options[["lambda6Scale"]] || options[["glbScale"]] ||
       options[["averageInterItemCor"]])
  ) {

    startProgressbar(options[["noSamples"]] * options[["noChains"]])
    dataset <- scale(dataset, scale = FALSE)
    c_out <- try(Bayesrel:::covSamp(dataset, options[["noSamples"]], options[["noBurnin"]],
                                    options[["noThin"]], options[["noChains"]],
                                    model[["pairwise"]], progressbarTick, k0 = options[["iwScale"]],
                                    df0 = options[["iwDf"]]), silent = TRUE)
    if (inherits(c_out, "try-error")) {
      if (model[["pairwise"]]) {
        .quitAnalysis(gettext("Sampling the posterior covariance matrix for either one of [alpha, lambda2, lambda6, glb] failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
      } else {
        .quitAnalysis(gettext("Sampling the posterior covariance matrix for either one of [alpha, lambda2, lambda6, glb] failed."))
      }
    }

    model[["gibbsSamp"]] <- c_out$cov_mat

  }

  model[["progressbarLength"]] <- options[["noChains"]] *
    length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]]))

  model[["itemsDropped"]] <- colnames(dataset)

  if (options[["disableSampleSave"]])
    return(model)

  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model, dependencies = c("iwScale", "iwDf"))

  return(model)
}


# standardize covariance matrix sample
.BayesianStdCov <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["gibbsCor"]]$object))
    return(.getStateContainerB(jaspResults)[["gibbsCor"]]$object)

  if (is.null(model[["gibbsSamp"]])) {
    return()
  } else {
    if (options[["stdCoeffs"]] == "unstand") {
      return()
    } else { # standardized
      out <- model[["gibbsSamp"]]
      startProgressbar(model[["progressbarLength"]])
      for (i in seq_len(nrow(model[["gibbsSamp"]]))) {
        for (j in seq_len(ncol(model[["gibbsSamp"]]))) {
          out[i, j, , ] <- .cov2cor.callback(model[["gibbsSamp"]][i, j, , ], progressbarTick)
        }
      }
      stateContainer <- .getStateContainerB(jaspResults)
      stateContainer[["gibbsCorObj"]] <- createJaspState(out, dependencies = "stdCoeffs")
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

  if (options[["omegaScale"]] && is.null(model[["empty"]])) {

    if (is.null(out[["loadings"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]])

      dataset <- scale(dataset, scale = FALSE)

      jaspBase::.setSeedJASP(options)

      tmp_out <- try(Bayesrel:::omegaSampler(dataset, options[["noSamples"]], options[["noBurnin"]],
                                             options[["noThin"]], options[["noChains"]],
                                             model[["pairwise"]], progressbarTick,
                                             a0 = options[["igShape"]], b0 = options[["igScale"]],
                                             m0 = options[["loadMean"]]), silent = TRUE)
      if (model[["pairwise"]] && inherits(tmp_out, "try-error")) {
        .quitAnalysis(gettext("Sampling the posterior factor model for omega failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
      }

      out[["loadings"]] <- tmp_out$lambda
      out[["residuals"]] <- tmp_out$psi
      out[["factor_var"]] <- tmp_out$phi
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["singleFactorObj"]] <- createJaspState(out, dependencies = c("omegaScale", "igShape", "igScale",
                                                                                 "loadMean"))
  }

  return(out)
}


.BayesianSingleFactorModelItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["singleFactorItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["singleFactorItemObj"]]$object)

  out <- model[["singleFactorItem"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaItem"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      return()
    }

    if (is.null(out[["itemLoadings"]])) {

      startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

      dataset <- scale(dataset, scale = FALSE)

      jaspBase::.setSeedJASP(options)

      out[["itemLoadings"]] <- array(0, c(options[["noChains"]],
                                          length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]])),
                                          ncol(dataset), ncol(dataset) - 1))
      out[["itemResiduals"]] <- array(0, c(options[["noChains"]],
                                           length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]])),
                                           ncol(dataset), ncol(dataset) - 1))

      for (i in seq_len(ncol(dataset))) {
        tmp_out <- Bayesrel:::omegaSampler(dataset[, -i],
                                           options[["noSamples"]], options[["noBurnin"]], options[["noThin"]],
                                           options[["noChains"]], model[["pairwise"]], progressbarTick,
                                           a0 = options[["igShape"]], b0 = options[["igScale"]],
                                           m0 = options[["loadMean"]])
        out[["itemLoadings"]][, , i, ] <- tmp_out$lambda
        out[["itemResiduals"]][, , i, ] <- tmp_out$psi
      }

      dd <- dim(out[["itemLoadings"]])
      out[["itemLoadings"]] <- array(out[["itemLoadings"]], c(dd[1] * dd[2], ncol(dataset), ncol(dataset) - 1))
      out[["itemResiduals"]] <- array(out[["itemResiduals"]], c(dd[1] * dd[2], ncol(dataset), ncol(dataset) - 1))

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["singleFactorItemObj"]] <- createJaspState(out, dependencies = c("omegaItem", "igShape", "igScale",
                                                                                     "loadMean"))
  }

  return(out)
}




# -------------------------------------------
#       Bayesian calculate coefficients
# -------------------------------------------


.BayesianOmegaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["omegaScaleObj"]]$object))
    return(.getStateContainerB(jaspResults)[["omegaScaleObj"]]$object)

  out <- model[["omegaScale"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaScale"]] && is.null(model[["empty"]]) && !is.null(model[["singleFactor"]])) {

    if (is.null(out[["samp"]])) {

      startProgressbar(model[["progressbarLength"]])

      if (options[["stdCoeffs"]] == "unstand") {

        out[["samp"]] <- coda::mcmc(.omegaOnArray(model[["singleFactor"]][["loadings"]],
                                                  model[["singleFactor"]][["residuals"]],
                                                  progressbarTick))
      }

      if (options[["stdCoeffs"]] == "stand" || options[["dispLoadings"]]) {

        lstd <- .stdFactorLoads(model[["singleFactor"]][["loadings"]],
                                model[["singleFactor"]][["residuals"]])
        estd <- 1 - lstd^2

        if (options[["stdCoeffs"]] == "stand") {
          out[["samp"]] <- coda::mcmc(.omegaOnArray(lstd, estd, progressbarTick))
        }

        if (options[["dispLoadings"]]) {
          out[["loadingsStdSamp"]] <- lstd
        }
      }

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["omegaScaleObj"]] <- createJaspState(out, dependencies = c("omegaScale", "stdCoeffs",
                                                                               "dispLoadings", "igShape", "igScale",
                                                                               "loadMean"))
  }

  return(out)
}


.BayesianOmegaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["omegaItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["omegaItemObj"]]$object)

  out <- model[["omegaItem"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaItem"]] && is.null(model[["empty"]]) && !is.null(model[["singleFactorItem"]])) {

    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")
      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      if (options[["stdCoeffs"]] == "unstand") {

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

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["omegaItemObj"]] <- createJaspState(out, dependencies = c("omegaItem", "stdCoeffs",
                                                                              "igShape", "igScale", "loadMean"))
  }

  return(out)
}


.BayesianAlphaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["alphaScaleObj"]]$object))
    return(.getStateContainerB(jaspResults)[["alphaScaleObj"]]$object)

  out <- model[["alphaScale"]]
  if (is.null(out))
    out <- list()

  if (options[["alphaScale"]] && is.null(model[["empty"]])) {

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])

      out[["samp"]] <- coda::mcmc(apply(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]],
                                        MARGIN = c(1, 2), Bayesrel:::applyalpha, progressbarTick))

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["alphaScaleObj"]] <- createJaspState(out, dependencies = c("alphaScale", "iwScale", "iwDf",
                                                                               "stdCoeffs"))
  }

  return(out)
}

.BayesianAlphaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["alphaItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["alphaItemObj"]]$object)

  out <- model[["alphaItem"]]
  if (is.null(out))
    out <- list()

  if (options[["alphaItem"]] && is.null(model[["empty"]])) {
    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")
      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesItemDroppedStats(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]],
                                                  Bayesrel:::applyalpha, progressbarTick)


    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["alphaItemObj"]] <- createJaspState(out, dependencies = c("alphaItem", "iwScale", "iwDf",
                                                                              "stdCoeffs"))
  }

  return(out)
}



.BayesianLambda2Scale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["lambda2ScaleObj"]]$object))
    return(.getStateContainerB(jaspResults)[["lambda2ScaleObj"]]$object)

  out <- model[["lambda2Scale"]]
  if (is.null(out))
    out <- list()

  if (options[["lambda2Scale"]] && is.null(model[["empty"]])) {

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])

      out[["samp"]] <- coda::mcmc(apply(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]],
                                        MARGIN = c(1, 2), Bayesrel:::applylambda2, progressbarTick))

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda2ScaleObj"]] <- createJaspState(out, dependencies = c("lambda2Scale", "iwScale", "iwDf",
                                                                                 "stdCoeffs"))
  }

  return(out)
}

.BayesianLambda2Item <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["lambda2ItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["lambda2ItemObj"]]$object)

  out <- model[["lambda2Item"]]
  if (is.null(out))
    out <- list()

  if (options[["lambda2Item"]] && is.null(model[["empty"]])) {
    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")

      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesItemDroppedStats(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]],
                                                  Bayesrel:::applylambda2, progressbarTick)
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda2ItemObj"]] <- createJaspState(out, dependencies = c("lambda2Item", "iwScale", "iwDf",
                                                                                "stdCoeffs"))
  }

  return(out)
}



.BayesianLambda6Scale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["lambda6ScaleObj"]]$object))
    return(.getStateContainerB(jaspResults)[["lambda6ScaleObj"]]$object)

  out <- model[["lambda6Scale"]]
  if (is.null(out))
    out <- list()

  if (options[["lambda6Scale"]] && is.null(model[["empty"]])) {

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])

      out[["samp"]] <- coda::mcmc(apply(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]],
                                        MARGIN = c(1, 2), Bayesrel:::applylambda6, progressbarTick))
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda6ScaleObj"]] <- createJaspState(out, dependencies = c("lambda6Scale", "iwScale", "iwDf",
                                                                                 "stdCoeffs"))
  }

  return(out)
}

.BayesianLambda6Item <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["lambda6ItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["lambda6ItemObj"]]$object)

  out <- model[["lambda6Item"]]
  if (is.null(out))
    out <- list()

  if (options[["lambda6Item"]] && is.null(model[["empty"]])) {
    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")

      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesItemDroppedStats(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]],
                                                  Bayesrel:::applylambda6, progressbarTick)
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda6ItemObj"]] <- createJaspState(out, dependencies = c("lambda6Item", "iwScale", "iwDf",
                                                                                "stdCoeffs"))
  }

  return(out)
}


.BayesianGlbScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["glbScaleObj"]]$object))
    return(.getStateContainerB(jaspResults)[["glbScaleObj"]]$object)

  out <- model[["glbScale"]]
  if (is.null(out))
    out <- list()

  if (options[["glbScale"]] && is.null(model[["empty"]])) {

    if (is.null(out[["samp"]])) {

      dd <- dim(model[["gibbsSamp"]])
      out[["samp"]] <- matrix(0, dd[1], dd[2])

      startProgressbar(dd[1] * 3)

      for (i in seq_len(dd[1])) {
        out[["samp"]][i, ] <- Bayesrel:::glbOnArrayCustom(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]][i, , , ],
                                                          callback = progressbarTick)
      }



    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["glbScaleObj"]] <- createJaspState(out, dependencies = c("glbScale", "iwScale", "iwDf",
                                                                             "stdCoeffs"))
  }

  return(out)
}

.BayesianGlbItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["glbItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["glbItemObj"]]$object)

  out <- model[["glbItem"]]
  if (is.null(out))
    out <- list()

  if (options[["glbItem"]] && is.null(model[["empty"]])) {
    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")

      return(out)
    }

    if (is.null(out[["itemSamp"]])) {
      # special case glb, because it works with arrays not only matrices, small speedup...
      dd <- dim(model[["gibbsSamp"]])
      out[["itemSamp"]] <- matrix(0, dd[1] * dd[2], dd[3])

      startProgressbar(3 * ncol(dataset))

      cov_samp <- array(model[[if (options[["stdCoeffs"]] == "unstand") "gibbsSamp" else "gibbsCor"]],
                        c(dd[1] * dd[2], dd[3], dd[3]))
      for (i in seq_len(dd[3])) {
        out[["itemSamp"]][, i] <- Bayesrel:::glbOnArrayCustom(cov_samp[, -i, -i], callback = progressbarTick)
      }

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["glbItemObj"]] <- createJaspState(out, dependencies = c("glbItem", "iwScale", "iwDf", "stdCoeffs"))
  }

  return(out)
}


.BayesianAverageCor <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["avgCorObj"]]$object))
    return(.getStateContainerB(jaspResults)[["avgCorObj"]]$object)

  out <- model[["average"]]
  if (is.null(out))
    out <- list()

  if (options[["averageInterItemCor"]] && !is.null(model[["gibbsSamp"]])) {

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

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["avgCorObj"]] <- createJaspState(out, dependencies = c("averageInterItemCor", "iwScale", "iwDf"))
  }

  return(out)
}


.BayesianMean <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["meanObj"]]$object))
    return(.getStateContainerB(jaspResults)[["meanObj"]]$object)

  out <- model[["mean"]]
  if (is.null(out))
    out <- list()
  if (options[["meanScale"]] && is.null(model[["empty"]])) {
    out[["est"]] <- if (options[["scoresMethod"]] == "sumScores")
      mean(rowSums(dataset, na.rm = TRUE))
    else
      mean(rowMeans(dataset, na.rm = TRUE))

    out[["cred"]] <- c(NA_real_, NA_real_)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("meanScale", "scoresMethod"))
  }
  return(out)
}


.BayesianStdDev <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["sdObj"]]$object))
    return(.getStateContainerB(jaspResults)[["sdObj"]]$object)

  out <- model[["sd"]]
  if (is.null(out))
    out <- list()
  if (options[["sdScale"]] && is.null(model[["empty"]])) {
    out[["est"]] <- if (options[["scoresMethod"]] == "sumScores")
      sd(rowSums(dataset, na.rm = TRUE))
    else
      sd(rowMeans(dataset, na.rm = TRUE))

    out[["cred"]] <- c(NA_real_, NA_real_)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["sdObj"]] <- createJaspState(out, dependencies = c("sdScale", "scoresMethod"))
  }
  return(out)
}



.BayesianItemRestCor <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemRestObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemRestObj"]]$object)

  out <- model[["itemRestCor"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemRestCor"]] && is.null(model[["empty"]])) {

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

      jaspBase::.setSeedJASP(options)
      out[["itemSamp"]] <- .itemRestCor(dataset, options[["noSamples"]], options[["noBurnin"]],
                                        options[["noThin"]], options[["noChains"]], model[["pairwise"]],
                                        callback = progressbarTick, options[["iwScale"]])
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = c("itemRestCor", "iwScale", "iwDf"))
  }

  return(out)
}

.BayesianMeanItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["meanItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["meanItemObj"]]$object)

  out <- model[["meanItem"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["meanItem"]] && is.null(model[["empty"]])) {

    out[["itemEst"]] <- colMeans(dataset, na.rm = TRUE)

    out[["itemCred"]] <- matrix(NA_real_, ncol(dataset), 2)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["meanItemObj"]] <- createJaspState(out, dependencies = "meanItem")
  }
  return(out)
}

.BayesianSdItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["sdItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["sdItemObj"]]$object)

  out <- model[["sdItem"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["sdItem"]] && is.null(model[["empty"]])) {

    out[["itemEst"]] <- apply(dataset, 2, sd, na.rm = TRUE)
    out[["itemCred"]] <- matrix(NA_real_, ncol(dataset), 2)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["sdItemObj"]] <- createJaspState(out, dependencies = "sdItem")
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

    ciValue <- options[["credibleIntervalValueScale"]]

    selected <- names(which(model[["derivedOptions"]][["selectedEstimators"]]))
    # first the coefficients with samples
    sampellist <- model[selected]
    samps <- .sampleListHelper(sampellist, "samp")

    out[["est"]] <- lapply(samps, .getPointEstFun(options[["pointEst"]]))
    out[["cred"]] <- lapply(samps, function(x) coda::HPDinterval(coda::mcmc(c(x)), prob = ciValue))

    if (options[["rHat"]]) {
      out[["rHat"]] <- lapply(samps, function(x) {
        # for each samp, (1) convert the rows to a coda object, (2) convert to mcmc.list, so that (3) gelman.diag is happy.
        coda::gelman.diag(coda::as.mcmc.list(lapply(seq_len(nrow(x)), function(i) coda::mcmc(x[i, ]))))[["psrf"]][, 1]
      })
    }

    if (options[["dispLoadings"]]) {
      out[["loadingsStd"]] <- apply(model[["omegaScale"]][["loadingsStdSamp"]], 3, .getPointEstFun(options[["pointEst"]]))
    }

    # check for mean and sd
    if ("meanScale" %in% selected) {
      out[["est"]][["meanScale"]] <- model[["meanScale"]][["est"]]
      out[["cred"]][["meanScale"]] <- model[["meanScale"]][["cred"]]
      if (options[["rHat"]]) {
        out[["rHat"]][["meanScale"]] <- NA_real_
      }
    }
    if ("sdScale" %in% selected) {
      out[["est"]][["sdScale"]] <- model[["sdScale"]][["est"]]
      out[["cred"]][["sdScale"]] <- model[["sdScale"]][["cred"]]
      if (options[["rHat"]]) {
        out[["rHat"]][["sdScale"]] <- NA_real_
      }
    }

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleResultsObj"]] <- createJaspState(out, dependencies = c("credibleIntervalValueScale",
                                                                                 "meanScale", "sdScale", "rHat",
                                                                                 "alphaScale", "omegaScale",
                                                                                 "lambda2Scale", "lambda6Scale",
                                                                                 "glbScale","averageInterItemCor",
                                                                                 "scoresMethod", "iwScale", "iwDf",
                                                                                 "igShape", "igScale",
                                                                                 "stdCoeffs", "pointEst",
                                                                                 "loadMean", "dispLoadings"))

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

    ciValue <- options[["credibleIntervalValueItem"]]

    selected <- names(which(model[["derivedOptions"]][["itemDroppedSelected"]]))
    # first the coefficients with samples
    sampellist <- model[selected]
    samps <- .sampleListHelper(sampellist, "itemSamp")

    if (options[["pointEst"]] == "mean") {
      out[["est"]] <- lapply(samps, colMeans)
    } else { # median
      out[["est"]] <- lapply(samps, .colMedians)
    }
    out[["cred"]] <- lapply(samps, function(x) coda::HPDinterval(coda::mcmc(x), prob = ciValue))

    # check for mean and sd
    if ("meanItem" %in% selected) {
      out[["est"]][["meanItem"]] <- model[["meanItem"]][["itemEst"]]
      out[["cred"]][["meanItem"]] <- model[["meanItem"]][["itemCred"]]

    }
    if ("sdItem" %in% selected) {
      out[["est"]][["sdItem"]] <- model[["sdItem"]][["itemEst"]]
      out[["cred"]][["sdItem"]] <- model[["sdItem"]][["itemCred"]]

    }

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemResultsObj"]] <- createJaspState(out, dependencies = c("omegaItem",  "alphaItem",
                                                                                "lambda2Item",  "lambda6Item",
                                                                                "glbItem","credibleIntervalValueItem",
                                                                                "itemRestCor", "meanItem", "sdItem",
                                                                                "iwScale", "iwDf", "igShape", "igScale",
                                                                                "stdCoeffs", "pointEst", "loadMean"))

  }

  return(out)
}


# see https://www.rdocumentation.org/packages/blavaan/versions/0.3-17/topics/blavFitIndices

.BayesianFitMeasures <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["fitMeasuresObj"]]$object))
    return(.getStateContainerB(jaspResults)[["fitMeasuresObj"]]$object)

  out <- model[["fitMeasures"]]
  if (is.null(out)) {
    out <- list()
  }

  if (options[["omegaScale"]] && options[["fitMeasures"]] && is.null(model[["empty"]])) {

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
                       options[["noSamples"]] * options[["noChains"]] + # for null model sampling
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
    res_null <- try(Bayesrel:::omegaSamplerNull(dataset, options[["noSamples"]], options[["noBurnin"]],
                                                options[["noThin"]], options[["noChains"]],
                                                model[["pairwise"]], progressbarTick,
                                                a0 = options[["igShape"]], b0 = options[["igScale"]],
                                                m0 = options[["loadMean"]]), silent = TRUE)

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
    stateContainer[["fitMeasuresObj"]] <- createJaspState(out, dependencies = c("omegaScale", "igShape", "igScale",
                                                                                "loadMean", "fitMeasures"))
  }

  return(out)
}




# ------------------------------------------- #
#       Bayesian output tables
# ------------------------------------------- #


.BayesianScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Bayesian Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("credibleIntervalValueScale", "meanScale", "sdScale", "rHat",
                                  "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                  "averageInterItemCor", "scoresMethod"))

  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  interval <- gettextf("%s%% CI",
                       format(100 * options[["credibleIntervalValueScale"]], digits = 3, drop0trailing = TRUE))
  intervalLow <- gettextf("%s lower bound", interval)
  intervalUp <- gettextf("%s upper bound", interval)

  pointEst <- gettextf("Posterior %s", options[["pointEst"]])

  if (options[["rHat"]]) {
    allData <- data.frame(estimate = c(pointEst, intervalLow, intervalUp, "R-hat"))
  } else {
    allData <- data.frame(estimate = c(pointEst, intervalLow, intervalUp))
  }

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idxSelected <- which(selected)

  for (i in idxSelected) {
    scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
  }

  if (.is.empty(model)) {
    scaleTable$setData(allData)
    if (model[["footnote"]] != "") {
      scaleTable$addFootnote(model[["footnote"]])
    }
    scaleTable$position <- 1
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleTable"]] <- scaleTable
    return()
  }

  for (j in seq_along(idxSelected)) {
    i <- idxSelected[j]
    nm <- names(idxSelected[j])

    newData <- data.frame(est = c(unlist(model[["scaleResults"]][["est"]][[nm]], use.names = FALSE),
                                  unlist(model[["scaleResults"]][["cred"]][[nm]], use.names = FALSE)))

    if (options[["rHat"]]) {
      newData <- rbind(newData, model[["scaleResults"]][["rHat"]][[nm]])
    }
    colnames(newData) <- paste0(colnames(newData), i)
    allData <- cbind(allData, newData)
  }

  scaleTable$setData(allData)

  if (model[["footnote"]] != "") {
    scaleTable$addFootnote(model[["footnote"]])
  }

  scaleTable$position <- 1
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  return()
}


.BayesianItemTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["itemTable"]]$object) ||
      !any(model[["derivedOptions"]][["itemDroppedSelected"]]))
    return()

  derivedOptions <- model[["derivedOptions"]]

  itemTable <- createJaspTable(gettext("Bayesian Individual Item Reliability Statistics"))

  itemTable$dependOn(options = c("omegaItem",  "alphaItem",  "lambda2Item",  "lambda6Item", "glbItem",
                                 "credibleIntervalValueItem", "itemRestCor", "meanItem", "sdItem",
                                 "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale"))

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemTable$position <- 2
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  overTitles <- format(derivedOptions[["namesEstimators"]][["tables_item"]], digits = 3, drop0trailing = TRUE)
  overTitles <- gettextf("%s (if item dropped)", overTitles)
  cred <- format(100 * options[["credibleIntervalValueItem"]], digits = 3, drop0trailing = TRUE)

  selected <- derivedOptions[["itemDroppedSelected"]]
  idxSelected <- which(selected)
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]]
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]


  itemTable[["variable"]] <- model[["itemsDropped"]]

  footnote <- ""
  if (!is.null(unlist(options[["reverseScaledItems"]])))
    footnote <- .addFootnoteReverseScaledItems(options)

  pointEst <- gettextf("Posterior %s", options[["pointEst"]])

  fewItemProblem <- FALSE
  for (i in idxSelected) {
    if (estimators[i] %in% coefficients) {
      if (estimators[i] == "Item-rest correlation") { # no item deleted for item rest cor
        itemTable$addColumnInfo(name = paste0("postMean", i), title = pointEst, type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
      } else {
        itemTable$addColumnInfo(name = paste0("postMean", i), title = pointEst, type = "number",
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

      if (i %in% 1:6) {
        rows <- length(options[["variables"]])
        if (rows < 3 && nm != "itemRestCor") {
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

  if (!is.null(.getStateContainerB(jaspResults)[["probTable"]]$object))
    return()

  # check if the values are in the proper order, meaning the lower field value is smaller than the upper field
  if (options[["probTableValueLow"]] > options[["probTableValueHigh"]]) {
    low <- options[["probTableValueHigh"]]
    high <- options[["probTableValueLow"]]
    footnote <- gettext("The bounds you entered have been rearranged in increasing order to provide meaningful results.")
  } else {
    low <- options[["probTableValueLow"]]
    high <- options[["probTableValueHigh"]]
    footnote <- ""
  }
  probTable <- createJaspTable(
    gettextf("Probability that Reliability Statistic is Larger than %.2f and Smaller than %.2f", low, high))
  probTable$dependOn(options = c("probTableValueLow", "probTable", "probTableValueHigh"))

  overTitle <- gettext("Probability")
  probTable$addColumnInfo(name = "statistic", title = gettext("Statistic"), type = "string")
  probTable$addColumnInfo(name = "prior", title = gettext("Prior"), type = "number", overtitle = overTitle)
  probTable$addColumnInfo(name = "posterior", title = gettext("Posterior"), type = "number", overtitle = overTitle)

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  idxSelected  <- which(selected)

  if (options[["probTable"]] && is.null(model[["empty"]])) {

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

      if (nm == "omegaScale") {
        startProgressbar(2e3)
      } else {
        startProgressbar(4e3)
      }

      prior <- .samplePrior(n.item, nm, progressbarTick, options[["iwScale"]], options[["iwDf"]],
                            options[["igShape"]], options[["igScale"]], options[["loadMean"]])

      probsPrior[i] <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]]) -
        sum(prior[["y"]][poshigh:end]) / sum(prior[["y"]])

    }

    df <- data.frame(statistic = opts[idxSelected], prior = probsPrior, posterior = probsPost)
    probTable$setData(df)

    if (footnote != "") {
      probTable$addFootnote(footnote)
    }

    probTable$position <- 3
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["probTable"]] <- probTable
  }

  return()
}


# use posterior mean and median in the table, also flip the table

.BayesianFitMeasuresTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["fitTable"]]$object))
    return()

  fitTable <- createJaspTable(gettextf("Fit Measures for the Single-Factor Model"))

  fitTable$dependOn(options = c("omegaScale", "fitMeasures", "fitCutoffSat", "fitCutoffNull", "pointEst",
                                "credibleIntervalValueFitMeasures"))

  fitTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  fitTable$addColumnInfo(name = "lr", title = "B-LR", type = "number")
  fitTable$addColumnInfo(name = "srmr", title = "B-SRMR", type = "number")
  fitTable$addColumnInfo(name = "rmsea", title = "B-RMSEA", type = "number")
  fitTable$addColumnInfo(name = "cfi", title = "B-CFI", type = "number")
  fitTable$addColumnInfo(name = "tli", title = "B-TLI", type = "number")

  cred <- gettextf("%s%% CI",
                   format(100 * options[["credibleIntervalValueFitMeasures"]], digits = 3, drop0trailing = TRUE))
  intervalLow <- gettextf("%s lower bound", cred)
  intervalUp <- gettextf("%s upper bound", cred)

  pointEst <- gettextf("Posterior %s", options[["pointEst"]])

  if (options[["omegaScale"]] && options[["fitMeasures"]] && is.null(model[["empty"]])) {

    pointEsts <- vapply(model[["fitMeasures"]], .getPointEstFun(options[["pointEst"]]), numeric(1))

    creds <- vapply(model[["fitMeasures"]][3:5], function(x) coda::HPDinterval(coda::mcmc(c(x))), numeric(2))
    creds <- cbind(matrix(NA_real_, 2, 2), creds) # no entry for LR and SRMR

    cutoffs_saturated <- vapply(model[["fitMeasures"]][3], function(x) mean(x < options[["fitCutoffSat"]]),
                                numeric(1))
    cutoffs_null <- vapply(model[["fitMeasures"]][-(1:3)], function(x) mean(x > options[["fitCutoffNull"]]),
                           numeric(1))
    cutoffs <- c(NA_real_, NA_real_, cutoffs_saturated, cutoffs_null) # no entry for LR and SRMR

    df <- data.frame(estimate = c(pointEst, intervalLow, intervalUp, "Relative to cutoff"))
    dfnums <- rbind(pointEsts, creds, cutoffs)
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

  loadTable$dependOn(options = c("omegaScale", "dispLoadings"))

  loadTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")
  loadTable$addColumnInfo(name = "loadings", title = gettext("Standardized loading"), type = "number")

  derivedOptions <- model[["derivedOptions"]]

  if (options[["omegaScale"]] && options[["dispLoadings"]] && is.null(model[["empty"]])) {



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




# -------------------------------------------
#       Bayesian output plots
# -------------------------------------------


.BayesianPosteriorPlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainer"]]$object))
    return()

  plotContainer <- createJaspContainer(gettext("Posterior Plots"))
  plotContainer$dependOn(options = c("plotPosterior", "shadePlots", "probTable", "probTableValueLow",
                                     "probTableValueHigh", "fixXRange", "dispPrior", "credibleIntervalValueScale",
                                     "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale"))

  derivedOptions <- model[["derivedOptions"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  idxSelected   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

  if (options[["shadePlots"]] && options[["probTable"]]) {
    if (options[["probTableValueLow"]] > options[["probTableValueHigh"]]) {
      low <- options[["probTableValueHigh"]]
      high <- options[["probTableValueLow"]]
    } else {
      low <- options[["probTableValueLow"]]
      high <- options[["probTableValueHigh"]]
    }
    shadePlots <- c(low, high)
  } else {
    shadePlots <- NULL
  }

  if (options[["plotPosterior"]] && is.null(model[["empty"]])) {
    n.item <- model[["k"]]

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (is.null(plotContainer[[nmsObjsNoGreek[i]]])) {
        if (options[["dispPrior"]]) {

          if (nm == "omegaScale") {
            startProgressbar(2e3)
          } else {
            startProgressbar(4e3)
          }
          prior <- .samplePrior(n.item, nm, progressbarTick, options[["iwScale"]], options[["iwDf"]],
                                options[["igShape"]], options[["igScale"]], options[["loadMean"]])
        } else {
          prior <- NULL
        }

        p <- .makeSinglePosteriorPlot(model[[nm]], model[["scaleResults"]][["cred"]][[nm]], nmsLabs[[i]],
                                      options[["fixXRange"]], shadePlots,
                                      options[["dispPrior"]], prior)
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

.makeSinglePosteriorPlot <- function(coefList, coefResults, nms, fixXRange, shade = NULL, priorTrue, priorSample) {

  # TODO: consider precomputing all densities (maybe with kernsmooth?) and reducing memory that way

  samp_tmp <- as.vector(coefList[["samp"]])
  if (fixXRange) {
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
  if (fixXRange) {
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
  if (fixXRange && max(datDens$y, na.rm = T) >= 1000) xExpand <- c(xExpand[1], .05)
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
  plotContainerItem$dependOn(options = c("variables", "plotItem",
                                         "credibleIntervalValueItem", "orderType", "orderItem",
                                         "omegaItem", "alphaItem", "lambda2Item", "lambda6Item", "glbItem"))

  derivedOptions <- model[["derivedOptions"]]

  selected <- derivedOptions[["itemDroppedSelectedItem"]]
  idxSelected   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables_item"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

  if (options[["orderItem"]]) {
    ordering <- options[["orderType"]]
  } else {
    ordering <- NULL
  }

  if (is.null(model[["empty"]]) && options[["plotItem"]]) {
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
          name <- unlist(strsplit(nm, "Item"))
          coefPos <- grep(name, names(model[["scaleResults"]][["est"]]))

          p <- .makeIfItemPlot(model[[nm]], model[[prevNumber]],
                               model[["itemResults"]][["est"]][[nm]],
                               model[["scaleResults"]][["est"]][[coefPos]],
                               nmsLabs[[i]],
                               options[["credibleIntervalValueItem"]],
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

    if (ordering == "orderItemMean") {
      dists <- abs(coefScaleEst - coefItemEst)
      dists[length(dists) + 1] <- 0
      est <- est[order(dists, decreasing = FALSE), ]
      dat$var <- factor(dat$var, levels = c(est$name))

    } else if (ordering == "orderItemKL") {
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

    } else if (ordering == "orderItemKS") {
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

  if (options[["dispPPC"]] && options[["omegaScale"]] && is.null(model[["empty"]])) {

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
    plot$dependOn(options = c("dispPPC", "omegaScale", "stdCoeffs"))

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
    plotContainerTP$dependOn(options = c("tracePlot", "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale",
                                         "glbScale"))

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





# -------------------------------------------
#       Bayesian help functions
# -------------------------------------------


.getStateContainerB <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = c("variables", "reverseScaledItems",
                                                                          "noSamples", "noBurnin", "noThin",
                                                                          "noChains", "missingValues", "setSeed",
                                                                          "seed", "disableSampleSave")
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

  if (estimate == "omegaScale") {
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

  if (estimate == "alphaScale") {
    prioralpha <- apply(m, MARGIN = 1, Bayesrel:::applyalpha, callback)
    out <- density(prioralpha, from = 0, to = 1, n = 512)
    return(out)

  }
  if (estimate == "lambda2Scale") {
    priorlambda2 <- apply(m, MARGIN = 1, Bayesrel:::applylambda2, callback)
    out <- density(priorlambda2, from = 0, to = 1, n = 512)
    return(out)

  }
  if (estimate == "lambda6Scale") {
    priorlambda6 <- apply(m, MARGIN = 1, Bayesrel:::applylambda6, callback)
    out <- density(priorlambda6, from = 0, to = 1, n = 512)
    return(out)

  }
  if (estimate == "glbScale") {
    priorglb <- Bayesrel:::glbOnArrayCustom(m, callback)
    out <- density(priorglb, from = 0, to = 1, n = 512)
    return(out)

  }

}

.BayesItemDroppedStats <- function(cov_samp, f1 = function(){}, callback = function(){}) {

  dd <- dim(cov_samp)
  out <- matrix(0, dd[1] * dd[2], dd[3])
  cov_samp <- array(cov_samp, c(dd[1] * dd[2], dd[3], dd[3]))
  for (i in seq_len(dd[3])) {
    out[, i] <- apply(cov_samp[, -i, -i], c(1), f1, callback)
  }

  return(out)
}



.itemRestCor <- function(dataset, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0) {

  ircor_samp <- matrix(0, n.chains * length(seq(1, n.iter - n.burnin, thin)), ncol(dataset))
  for (i in seq(ncol(dataset))) {
    help_dat <- cbind(as.matrix(dataset[, i]), rowMeans(as.matrix(dataset[, -i]), na.rm = TRUE))
    ircor_samp[, i] <- .WishartCorTransform(help_dat, n.iter = n.iter, n.burnin = n.burnin, thin = thin,
                                            n.chains = n.chains, pairwise = pairwise, callback = callback, k0)
  }

  return(ircor_samp)
}

.WishartCorTransform <- function(x, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0) {
  tmp_cov <- Bayesrel:::covSamp(x, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0 = k0, df0 = NULL)$cov_mat
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

.getPointEstFun <- function(pointEstFunString) {
  return(switch(pointEstFunString,
                "mean"   = mean,
                "median" = median,
                stop("getPointEstFun has no case for value: ", pointEstFunString))
  )
}
