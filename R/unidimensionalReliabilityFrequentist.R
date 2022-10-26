


#' @export
unidimensionalReliabilityFrequentist <- function(jaspResults, dataset, options) {

  dataset <- .readData(dataset, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }

  .checkErrors(dataset, options)


  model <- .frequentistPreCalc(jaspResults, dataset, options)
  options <- .scaleItemBoxAlign(options)

  model[["derivedOptions"]] <- .frequentistDerivedOptions(options)
  model[["scaleOmega"]] <- .frequentistOmegaScale(jaspResults, dataset, options, model)
  model[["itemDeletedOmega"]] <- .frequentistOmegaItem(jaspResults, dataset, options, model)
  model[["scaleAlpha"]] <- .frequentistAlphaScale(jaspResults, dataset, options, model)
  model[["itemDeletedAlpha"]] <- .frequentistAlphaItem(jaspResults, dataset, options, model)
  model[["scaleLambda2"]] <- .frequentistLambda2Scale(jaspResults, dataset, options, model)
  model[["itemDeletedLambda2"]] <- .frequentistLambda2Item(jaspResults, dataset, options, model)
  model[["scaleLambda6"]] <- .frequentistLambda6Scale(jaspResults, dataset, options, model)
  model[["itemDeletedLambda6"]] <- .frequentistLambda6Item(jaspResults, dataset, options, model)
  model[["scaleGreatestLowerBound"]] <- .frequentistGlbScale(jaspResults, dataset, options, model)
  model[["itemDeletedGreatestLowerBound"]] <- .frequentistGlbItem(jaspResults, dataset, options, model)

  model[["averageInterItemCorrelation"]] <- .frequentistAverageCor(jaspResults, dataset, options, model)
  model[["scaleMean"]] <- .frequentistMean(jaspResults, dataset, options, model)
  model[["scaleSd"]] <- .frequentistStdDev(jaspResults, dataset, options, model)
  model[["itemRestCorrelation"]] <- .frequentistItemRestCor(jaspResults, dataset, options, model)
  model[["itemMean"]] <- .frequentistMeanItem(jaspResults, dataset, options, model)
  model[["itemSd"]] <- .frequentistSdItem(jaspResults, dataset, options, model)

  model[["scaleResults"]] <- .frequentistComputeScaleResults(jaspResults, dataset, options, model)

  .frequentistScaleTable(jaspResults, model, options)
  .frequentistItemTable(jaspResults, model, options)
  .frequentistSingleFactorFitTable(jaspResults, model, options)


  return()

}

.frequentistDerivedOptions <- function(options) {

  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("scaleOmega", "scaleAlpha", "scaleLambda2", "scaleLambda6",
                                           "scaleGreatestLowerBound", "averageInterItemCorrelation", "scaleMean", "scaleSd")]),
    itemDroppedSelected = unlist(options[c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2", "itemDeletedLambda6",
                                           "itemDeletedGreatestLowerBound", "itemRestCorrelation", "itemMean", "itemSd")]),
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

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = c("variables", "reverseScaledItems", "bootstrapSamples",
                                                                          "naAction", "bootstrapType", "setSeed",
                                                                          "seed", "ci", "samplesSavingDisabled")
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




# -------------------------------------------
#       Frequentist precalculate results
# -------------------------------------------

.frequentistPreCalc <- function(jaspResults, dataset, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["modelObj"]]$object)) {
    if (!is.null(.getStateContainerF(jaspResults)[["modelObj"]]$object$bootSamp))
      return(.getStateContainerF(jaspResults)[["modelObj"]]$object)
  }



  derivedOptions <- .frequentistDerivedOptions(options)

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

    k <- ncol(dataset); n <- nrow(dataset)
    model[["k"]] <- k
    model[["n"]] <- n
    cc <- cov(dataset, use = model[["use.cases"]])
    model[["data_cov"]] <- cc
    model[["data_cor"]] <- cor(dataset, use = model[["use.cases"]])
    model[["itemsDropped"]] <- colnames(dataset)


    # check if any items correlate negatively with the scale
    model[["footnote"]] <- sprintf("%s%s", model[["footnote"]], .checkLoadings(dataset, options[["variables"]]))
  }

  jaspBase::.setSeedJASP(options)

  # check if interval is checked and bootstrapped covariance sample has to be generated
  if (is.null(model[["bootSamp"]]) &&
      options[["ci"]] &&
      (
        (options[["scaleOmega"]] && options[["omegaEstimationMethod"]] == "pfa") ||
        (options[["scaleAlpha"]] && !(options[["alphaIntervalMethod"]] == "analytic")) ||
        options[["scaleLambda2"]] ||
        options[["scaleLambda6"]] ||
        options[["scaleGreatestLowerBound"]] ||
        options[["averageInterItemCorrelation"]]
      )
  ) {


    boot_cov <- array(0, c(options[["bootstrapSamples"]], k, k))

    startProgressbar(options[["bootstrapSamples"]])

    if (options[["bootstrapType"]] == "parametric") {
      model[["parametric"]] <- TRUE
      for (i in seq_len(options[["bootstrapSamples"]])) {
        boot_data <- MASS::mvrnorm(n, colMeans(dataset, na.rm = TRUE), cc)
        boot_cov[i, , ] <- cov(boot_data)
        progressbarTick()
      }
    } else {
      model[["parametric"]] <- FALSE
      for (i in seq_len(options[["bootstrapSamples"]])) {
        boot_data <- as.matrix(dataset[sample.int(n, size = n, replace = TRUE), ])
        boot_cov[i, , ] <- cov(boot_data, use = model[["use.cases"]])
        progressbarTick()
      }
    }
    model[["bootSamp"]] <- boot_cov
  }

  if (options[["samplesSavingDisabled"]])
    return(model)

  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model)

  return(model)
}




# -------------------------------------------
#       Frequentist calculate coefficients
# -------------------------------------------


.frequentistOmegaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["scaleOmegaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["scaleOmegaObj"]]$object)

  out <- model[["scaleOmega"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleOmega"]] && is.null(model[["empty"]]) && options[["ci"]]) {

    ciValue <- options[["ciLevel"]]
    if (options[["omegaEstimationMethod"]] == "cfa") {
      if (options[["omegaIntervalMethod"]] == "bootstrapped") {
        parametric <- options[["bootstrapType"]] == "parametric"
        omegaboot <- out[["samp"]]
        if (is.null(omegaboot)) {
          startProgressbar(options[["bootstrapSamples"]])
          jaspBase::.setSeedJASP(options)
          omegaboot <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = FALSE,
                                                pairwise = model[["pairwise"]], parametric = parametric,
                                                n.boot = options[["bootstrapSamples"]], callback = progressbarTick)
          out[["samp"]] <- omegaboot[["omega_boot"]]
        }
      }
    } else { # omega with pfa
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["bootstrapSamples"]])
        out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applyomegaPFA, callback = progressbarTick)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleOmegaObj"]] <- createJaspState(out, dependencies = c("scaleOmega", "omegaEstimationMethod",
                                                                               "omegaIntervalMethod"))
  }

  return(out)
}

.frequentistOmegaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedOmegaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedOmegaObj"]]$object)

  out <- model[["itemDeletedOmega"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedOmega"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    if (options[["omegaEstimationMethod"]] == "cfa") {

      dataset <- scale(dataset, scale = FALSE)
      # do we have to compute item dropped values
      if (is.null(out[["itemDropped"]])) {
        out[["itemDropped"]] <- numeric(ncol(dataset))
        for (i in seq_len(ncol(dataset))) {
          out[["itemDropped"]][i] <- Bayesrel:::applyomegaCFAData(dataset[, -i], interval = .95,
                                                                  pairwise = model[["pairwise"]])
        }
      }
      if (anyNA(out[["itemDropped"]])) {
        out[["error"]] <- gettext("Omega item dropped statistics with CFA failed.")
        out[["itemDropped"]] <- rep(NaN, ncol(dataset))
      }

    } else { # omega with pfa
      # do we have to compute item dropped values
      if (is.null(out[["itemDropped"]]))
        out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applyomegaPFA)

      if (anyNA(out[["itemDropped"]]))
        out[["error"]] <- gettext("Omega item dropped statistics with PFA failed.")

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedOmegaObj"]] <- createJaspState(out, dependencies = c("itemDeletedOmega", "omegaEstimationMethod"))
  }

  return(out)
}


.frequentistAlphaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["scaleAlphaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["scaleAlphaObj"]]$object)

  out <- model[["scaleAlpha"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleAlpha"]] && is.null(model[["empty"]])) {

    # alpha unstandardized
    if (options[["alphaType"]] == "unstandardized") {
      # do we need an interval estimate?
      if (options[["ci"]]) {
        # should the interval be bootstrapped
        if (options[["alphaIntervalMethod"]] == "bootstrapped") {
          if (is.null(out[["samp"]])) {
            startProgressbar(options[["bootstrapSamples"]])
            out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applyalpha, callback = progressbarTick)
          }
        }
      }
    } else { # alpha standardized
      # do we need an interval estimate?
      if (options[["ci"]]) {
        # should the interval be bootstrapped
        if (options[["alphaIntervalMethod"]] == "bootstrapped") {
          if (is.null(out[["sampCor"]])) {
            out[["sampCor"]] <- numeric(options[["bootstrapSamples"]])
            for (i in seq_len(options[["bootstrapSamples"]])) {
              out[["sampCor"]][i] <- Bayesrel:::applyalpha(cov2cor(model[["bootSamp"]][i, , ]))
            }
          }
        }
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleAlphaObj"]] <- createJaspState(out, dependencies = c("scaleAlpha", "alphaType",
                                                                               "alphaIntervalMethod"))
  }
  return(out)
}

.frequentistAlphaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedAlphaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedAlphaObj"]]$object)

  out <- model[["itemDeletedAlpha"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedAlpha"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    if (options[["alphaType"]] == "unstandardized") { # alpha unstandardized
      # do we have to compute item dropped values
      if (is.null(out[["itemDropped"]]))
        out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applyalpha)

    } else { # alpha standardized
      if (is.null(out[["itemDropped"]])) {
        out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cor"]], Bayesrel:::applyalpha)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedAlphaObj"]] <- createJaspState(out, dependencies = c("itemDeletedAlpha", "alphaType"))
  }
  return(out)
}


.frequentistLambda2Scale <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleLambda2Obj"]]$object))
    return(.getStateContainerF(jaspResults)[["scaleLambda2Obj"]]$object)

  out <- model[["scaleLambda2"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["scaleLambda2"]]  && is.null(model[["empty"]])) {
    # do we need an interval estimate?
    if (options[["ci"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["bootstrapSamples"]])
        out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applylambda2, callback = progressbarTick)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleLambda2Obj"]] <- createJaspState(out, dependencies = "scaleLambda2")
  }
  return(out)
}

.frequentistLambda2Item <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedLambda2Obj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedLambda2Obj"]]$object)

  out <- model[["itemDeletedLambda2"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["itemDeletedLambda2"]]  && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    if (is.null(out[["itemDropped"]]))
      out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applylambda2)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedLambda2Obj"]] <- createJaspState(out, dependencies = "itemDeletedLambda2")
  }
  return(out)
}


.frequentistLambda6Scale <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleLambda6Obj"]]$object))
    return(.getStateContainerF(jaspResults)[["scaleLambda6Obj"]]$object)

  out <- model[["scaleLambda6"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["scaleLambda6"]]  && is.null(model[["empty"]])) {
    if (options[["ci"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["bootstrapSamples"]])
        out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applylambda6, callback = progressbarTick)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleLambda6Obj"]] <- createJaspState(out, dependencies = "scaleLambda6")
  }
  return(out)
}

.frequentistLambda6Item <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedLambda6Obj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedLambda6Obj"]]$object)

  out <- model[["itemDeletedLambda6"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["itemDeletedLambda6"]]  && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    if (is.null(out[["itemDropped"]]))
      out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applylambda6)
    if (anyNA(out[["itemDropped"]])) {
      out[["error"]] <- gettext("Lambda6 item dropped statistics failed.")
      out[["itemDropped"]] <- rep(NaN, ncol(dataset))
    }
    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedLambda6Obj"]] <- createJaspState(out, dependencies = "itemDeletedLambda6")
  }
  return(out)
}


# check the error handling of the glb !!!!!!!!
.frequentistGlbScale <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleGreatestLowerBoundObj"]]$object))
    return(.getStateContainerF(jaspResults)[["scaleGreatestLowerBoundObj"]]$object)

  out <- model[["scaleGreatestLowerBound"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["scaleGreatestLowerBound"]]  && is.null(model[["empty"]])) {
    # do we need an interval estimate?
    if (options[["ci"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(4)
        out[["samp"]] <- Bayesrel:::glbOnArrayCustom(model[["bootSamp"]], callback = progressbarTick)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleGreatestLowerBoundObj"]] <- createJaspState(out, dependencies = "scaleGreatestLowerBound")
  }
  return(out)
}

.frequentistGlbItem <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedGreatestLowerBoundObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedGreatestLowerBoundObj"]]$object)

  out <- model[["itemDeletedGreatestLowerBound"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["itemDeletedGreatestLowerBound"]]  && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    # do we have to compute item dropped values
    if (is.null(out[["itemDropped"]])) {
      # special case glb since it has build in array functionality, but it might be only slightly faster
      itemDroppedCovs <- array(0, c(model[["k"]], model[["k"]]-1, model[["k"]]-1))
      for (i in seq_len(model[["k"]])) {
        itemDroppedCovs[i, , ] <- model[["data_cov"]][-i, -i]
      }
      out[["itemDropped"]] <- c(Bayesrel:::glbOnArrayCustom(itemDroppedCovs))
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedGreatestLowerBoundObj"]] <- createJaspState(out, dependencies = "itemDeletedGreatestLowerBound")
  }
  return(out)
}


.frequentistAverageCor <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["avgObj"]]$object))
    return(.getStateContainerF(jaspResults)[["avgObj"]]$object)

  out <- model[["averageInterItemCorrelation"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["averageInterItemCorrelation"]] && is.null(model[["empty"]])) {
    if (options[["ci"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["bootstrapSamples"]])
        out[["samp"]] <- numeric(options[["bootstrapSamples"]])
        for (i in seq_len(options[["bootstrapSamples"]])) {
          corm <- .cov2cor.callback(model[["bootSamp"]][i, , ], progressbarTick)
          out[["samp"]][i] <- mean(corm[lower.tri(corm)])
        }
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["avgObj"]] <- createJaspState(out,
                                                  dependencies = c("averageInterItemCorrelation"))
  }
  return(out)
}

.frequentistMean <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["meanObj"]]$object))
    return(.getStateContainerF(jaspResults)[["meanObj"]]$object)

  out <- model[["scaleMean"]]
  if (is.null(out))
    out <- list()
  if (options[["scaleMean"]] && is.null(model[["empty"]])) {
    ciValue <- options[["ciLevel"]]

    if (options[["meanSdScoresMethod"]] == "sumScores") {
      out[["est"]] <- mean(rowSums(dataset, na.rm = TRUE))
      sdmean <- sd(rowSums(dataset, na.rm = TRUE))
    } else {
      out[["est"]] <- mean(rowMeans(dataset, na.rm = TRUE))
      sdmean <- sd(rowMeans(dataset, na.rm = TRUE))
    }

    if (options[["ci"]]) {
      zz <- qnorm(1 - (1 - ciValue) / 2)
      out[["conf"]] <- c(out[["est"]] - zz * (sdmean / sqrt(model[["n"]])),
                         out[["est"]] + zz * (sdmean / sqrt(model[["n"]])))
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("scaleMean", "meanSdScoresMethod"))
  }
  return(out)
}

.frequentistStdDev <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["sdObj"]]$object))
    return(.getStateContainerF(jaspResults)[["sdObj"]]$object)

  out <- model[["scaleSd"]]
  if (is.null(out))
    out <- list()
  if (options[["scaleSd"]] && is.null(model[["empty"]])) {
    ciValue <- options[["ciLevel"]]

    out[["est"]] <- if (options[["meanSdScoresMethod"]] == "sumScores")
      sd(rowSums(dataset, na.rm = TRUE))
    else
      sd(rowMeans(dataset, na.rm = TRUE))

    if (options[["ci"]]) {
      chiValueLow <- qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
      out[["conf"]] <- c(sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueLow),
                         sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueHigh))
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["sdObj"]] <- createJaspState(out, dependencies = c("scaleSd", "meanSdScoresMethod"))
  }
  return(out)
}


.frequentistItemRestCor <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemRestObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemRestObj"]]$object)

  out <- model[["itemRestCorrelation"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemRestCorrelation"]]  && is.null(model[["empty"]])) {
    out[["itemDropped"]] <- numeric(ncol(dataset))
    for (i in seq_len(ncol(dataset))) {
      out[["itemDropped"]][i] <- cor(as.matrix(dataset[, i]), rowMeans(as.matrix(dataset[, -i]), na.rm = TRUE),
                                     use = model[["use.cases"]])
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = "itemRestCorrelation")
  }
  return(out)
}

.frequentistMeanItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemMeanObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemMeanObj"]]$object)

  out <- model[["itemMean"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemMean"]]  && is.null(model[["empty"]])) {
    out[["itemDropped"]] <- colMeans(dataset, na.rm = TRUE)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemMeanObj"]] <- createJaspState(out, dependencies = "itemMean")
  }
  return(out)
}

.frequentistSdItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemSdObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemSdObj"]]$object)

  out <- model[["itemSd"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemSd"]]  && is.null(model[["empty"]])) {

    out[["itemDropped"]] <- apply(dataset, 2, sd, na.rm = TRUE)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemSdObj"]] <- createJaspState(out, dependencies = "itemSd")
  }
  return(out)
}




.frequentistComputeScaleResults <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["scaleResultsObj"]]$object))
    return(.getStateContainerF(jaspResults)[["scaleResultsObj"]]$object)

  out <- model[["scaleResults"]]
  if (is.null(out))
    out <- list()

  if (is.null(model[["empty"]])) {

    ciValue <- options[["ciLevel"]]

    # go one coefficient at a time, because there are too many special options for a generic solution
    # #### omega ####
    if (options[["scaleOmega"]]) {
      if (options[["omegaEstimationMethod"]] == "cfa") {
        dataset <- scale(dataset, scale = FALSE)
        omegaO <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = TRUE,
                                           pairwise = model[["pairwise"]])
        if (is.na(omegaO[["omega"]])) {
          out[["error"]][["scaleOmega"]] <- gettext("Omega calculation with CFA failed.
                                                    Try changing to PFA in Advanced Options")
          out[["est"]][["scaleOmega"]] <- NaN
        } else {
          out[["fit"]][["scaleOmega"]] <- omegaO[["indices"]]
          out[["est"]][["scaleOmega"]] <- omegaO[["omega"]]

          if (options[["ci"]]) {
            if (options[["omegaIntervalMethod"]] == "analytic") {
              out[["conf"]][["scaleOmega"]] <- c(omegaO$omega_lower, omegaO$omega_upper)
            } else { # omega interval bootstrapped
              if (!is.null(model[["scaleOmega"]][["samp"]])) {
                if (sum(!is.na(model[["scaleOmega"]][["samp"]])) >= 2) {
                  out[["conf"]][["scaleOmega"]] <- quantile(model[["scaleOmega"]][["samp"]],
                                                            probs = c((1 - ciValue)/2, 1 - (1 - ciValue) / 2),
                                                            na.rm = TRUE)
                } else {
                  out[["error"]][["scaleOmega"]] <- gettext("Omega bootstrapped interval calculation with CFA failed.
                                                            Try changing to PFA in 'Advanced Options'")
                  out[["conf"]][["scaleOmega"]] <- NaN
                }
              }
            }
          }
        }
      } else { # omega method is pfa
        out[["est"]][["scaleOmega"]] <- Bayesrel:::applyomegaPFA(model[["data_cov"]])
        if (is.na(out[["est"]][["scaleOmega"]])) {
          out[["error"]][["scaleOmega"]] <- gettext("Omega calculation with PFA failed.")
          out[["est"]][["scaleOmega"]] <- NaN
        } else {
          if (options[["ci"]]) {
            if (!is.null(model[["scaleOmega"]][["samp"]])) {
              if (sum(!is.na(model[["scaleOmega"]][["samp"]])) >= 2) {
                out[["conf"]][["scaleOmega"]] <- quantile(model[["scaleOmega"]][["samp"]],
                                                          probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                          na.rm = TRUE)
              } else {
                out[["error"]][["scaleOmega"]] <- gettext("Omega interval calculation with PFA failed.")
                out[["conf"]][["scaleOmega"]] <- NaN
              }
            }
          }
        }
      }
    }

    # #### alpha ####
    if (options[["scaleAlpha"]]) {
      # alpha unstandardized
      if (options[["alphaType"]] == "unstandardized") {
        out[["est"]][["scaleAlpha"]] <- Bayesrel:::applyalpha(model[["data_cov"]])
        if (options[["ci"]]) {
          # should the interval be analytic
          if (options[["alphaIntervalMethod"]] == "analytic") {
            out[["conf"]][["scaleAlpha"]] <- Bayesrel:::ciAlpha(1 - ciValue, model[["n"]], model[["data_cov"]])
          } else { # alpha interval bootstrapped
            if (!is.null(model[["scaleAlpha"]][["samp"]])) {
              out[["conf"]][["scaleAlpha"]] <- quantile(model[["scaleAlpha"]][["samp"]],
                                                        probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                        na.rm = TRUE)
            }
          }
        }
      } else { # alpha standardized
        ccor <- model[["data_cor"]]
        out[["est"]][["scaleAlpha"]] <- Bayesrel:::applyalpha(ccor)
        if (options[["ci"]]) {
          # should the interval be analytic
          if (options[["alphaIntervalMethod"]] == "analytic") {
            out[["conf"]][["scaleAlpha"]] <- Bayesrel:::ciAlpha(1 - ciValue, model[["n"]], ccor)
          } else {
            if (!is.null(model[["scaleAlpha"]][["sampCor"]])) {
              out[["conf"]][["scaleAlpha"]] <- quantile(model[["scaleAlpha"]][["sampCor"]],
                                                        probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                        na.rm = TRUE)
            }
          }
        }
      }
    }

    # #### lambda 6 ####
    if (options[["scaleLambda6"]]) {
      out[["est"]][["scaleLambda6"]] <- Bayesrel:::applylambda6(model[["data_cov"]])
      if (is.na(out[["est"]][["scaleLambda6"]])) {
        out[["error"]][["scaleLambda6"]] <- gettext("Lambda6 calculation failed.")
      }
      if (options[["ci"]]) {
        if (sum(!is.na(model[["scaleLambda6"]][["samp"]])) >= 2) {
          out[["conf"]][["scaleLambda6"]] <- quantile(model[["scaleLambda6"]][["samp"]],
                                                      probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
        } else {
          out[["error"]][["scaleLambda6"]] <- gettext("Lambda6 interval calculation failed.")
        }

      }
    }

    if (options[["ci"]]) {
      # the intervals for coefficients with bootstrap samples can be done generic
      bootCoeffs <- c("scaleLambda2", "scaleGreatestLowerBound", "averageInterItemCorrelation")
      selected <- bootCoeffs[bootCoeffs %in% names(which(model[["derivedOptions"]][["selectedEstimators"]]))]

      sampellist <- model[selected]
      samps <- .sampleListHelper(sampellist, "samp")
      out[["conf"]][selected] <- lapply(samps, function(x) {quantile(x, probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                                     na.rm = TRUE)})
    }

    # point estimates
    if (options[["scaleLambda2"]]) {
      out[["est"]][["scaleLambda2"]] <- Bayesrel:::applylambda2(model[["data_cov"]])
    }
    if (options[["scaleGreatestLowerBound"]]) {
      out[["est"]][["scaleGreatestLowerBound"]] <- Bayesrel:::glbOnArrayCustom(model[["data_cov"]])
    }
    if (options[["averageInterItemCorrelation"]]) {
      out[["est"]][["averageInterItemCorrelation"]] <- mean(model[["data_cor"]][lower.tri(model[["data_cor"]])])
    }

    # just copying for mean and sd
    if (options[["scaleMean"]]) {
      out[["est"]][["scaleMean"]] <- model[["scaleMean"]][["est"]]
      out[["conf"]][["scaleMean"]] <- model[["scaleMean"]][["conf"]]
    }
    if (options[["scaleSd"]]) {
      out[["est"]][["scaleSd"]] <- model[["scaleSd"]][["est"]]
      out[["conf"]][["scaleSd"]] <- model[["scaleSd"]][["conf"]]
    }


    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleResultsObj"]] <- createJaspState(out, dependencies = c("ciLevel",
                                                                                 "scaleMean", "scaleSd",
                                                                                 "scaleAlpha", "scaleOmega",
                                                                                 "scaleLambda2", "scaleLambda6",
                                                                                 "scaleGreatestLowerBound","averageInterItemCorrelation",
                                                                                 "meanSdScoresMethod",
                                                                                 "omegaEstimationMethod", "omegaIntervalMethod",
                                                                                 "alphaIntervalMethod", "alphaType"))


  }

  return(out)

}




# -------------------------------------------
#       Frequentist output tables
# -------------------------------------------


.frequentistScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Frequentist Scale Reliability Statistics"))

  scaleTable$dependOn(options = c("scaleOmega", "scaleAlpha", "scaleLambda2", "scaleLambda6", "scaleGreatestLowerBound",
                                  "averageInterItemCorrelation", "scaleMean", "scaleSd", "meanSdScoresMethod",
                                  "omegaEstimationMethod", "omegaIntervalMethod", "alphaType", "alphaIntervalMethod"))
  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  if (options[["ci"]]) {
    interval <- gettextf("%s%% CI",
                         format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE))
    intervalLow <- gettextf("%s lower bound", interval)
    intervalUp <- gettextf("%s upper bound", interval)
    allData <- data.frame(estimate = c(gettext("Point estimate"), intervalLow, intervalUp))
  } else {
    allData <- data.frame(estimate = c(gettext("Point estimate")))
  }

  scaleTable$position <- 1
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  if (!is.null(model[["scaleResults"]][["error"]][["scaleOmega"]])) {
    model[["footnote"]] <- paste(model[["footnote"]], model[["scaleResults"]][["error"]][["scaleOmega"]])
  }
  if (!is.null(model[["scaleResults"]][["error"]][["scaleLambda6"]])) {
    model[["footnote"]] <- paste(model[["footnote"]], model[["scaleResults"]][["error"]][["scaleLambda6"]])
  }

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idxSelected <- which(selected)


  # if no coefficients selected:
  if (.is.empty(model)) {
    scaleTable$setData(allData)


    for (i in idxSelected) {
      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
    }

    if (model[["footnote"]] != "") {
      scaleTable$addFootnote(model[["footnote"]])
    }

    return()
  }

  if (options[["ci"]]) {
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
      newData <- data.frame(est = c(unlist(model[["scaleResults"]][["est"]][[nm]], use.names = FALSE),
                                    unlist(model[["scaleResults"]][["conf"]][[nm]], use.names = FALSE)))
      colnames(newData) <- paste0(colnames(newData), i)
      allData <- cbind(allData, newData)
    }
  } else {
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
      newData <- data.frame(est = c(unlist(model[["scaleResults"]][["est"]][[nm]], use.names = FALSE)))
      colnames(newData) <- paste0(colnames(newData), i)
      allData <- cbind(allData, newData)
    }
  }

  scaleTable$setData(allData)

  if (model[["footnote"]] != "") {
    scaleTable$addFootnote(model[["footnote"]])
  }
  return()
}



.frequentistItemTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["itemTable"]]$object) ||
      !any(model[["derivedOptions"]][["itemDroppedSelected"]])) {
    return()
  }

  derivedOptions <- model[["derivedOptions"]]

  itemTable <- createJaspTable(gettext("Frequentist Individual Item Reliability Statistics"))
  itemTable$dependOn(options = c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2", "itemDeletedLambda6", "itemDeletedGreatestLowerBound",
                                 "itemMean", "itemRestCorrelation", "itemSd",
                                 "scaleOmega", "scaleAlpha", "scaleLambda2", "scaleLambda6", "scaleGreatestLowerBound"))
  # adding the scale options fixes a bug, where the item table would remain displayed
  # after one had checked a scale coefficient box and the item coefficient box and then unchecked the scale coeff box

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemTable$position <- 2
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  footnote <- ""

  if (!is.null(model[["itemDeletedOmega"]][["error"]])) {
    footnote <- paste(footnote, model[["itemDeletedOmega"]][["error"]])
  }
  if (!is.null(model[["itemDeletedLambda6"]][["error"]])) {
    footnote <- paste(footnote, model[["itemDeletedOmega"]][["error"]])
  }

  selected <- derivedOptions[["itemDroppedSelected"]]
  coefficientsTable <- derivedOptions[["namesEstimators"]][["tables_item"]]
  overTitle <- gettext("If item dropped")
  idxSelected <- which(selected)
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]


  if (length(model[["itemsDropped"]]) > 0) {
    itemTable[["variable"]] <- model[["itemsDropped"]]

    if (!is.null(unlist(options[["reverseScaledItems"]])))
      footnote <- .addFootnoteReverseScaledItems(options)
  }

  fewItemProblem <- FALSE
  for (i in idxSelected) {
    if (coefficientsTable[i] %in% coefficients) {
      itemTable$addColumnInfo(name = paste0("pointEstimate", i), title = coefficientsTable[i], type = "number",
                              overtitle = overTitle)
      fewItemProblem <- TRUE
    } else {
      itemTable$addColumnInfo(name = paste0("pointEstimate", i), title = coefficientsTable[i], type = "number")
    }
  }

  if (is.null(model[["empty"]])) {
    tb <- data.frame(variable = model[["itemsDropped"]])

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])
      newtb <- cbind(pointEstimate = model[[nm]][["itemDropped"]])
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

  return()
}



# once the package is updated check this again and apply:
.frequentistSingleFactorFitTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["fitTable"]]$object) ||
      is.null(model[["scaleResults"]][["fit"]][["scaleOmega"]]) || !options[["omegaFitMeasures"]])
    return()

  fitTable <- createJaspTable(gettextf("Fit Measures of Single Factor Model Fit"))

  fitTable$addColumnInfo(name = "measure", title = gettext("Fit measure"),   type = "string")
  fitTable$addColumnInfo(name = "value",  title = gettext("Value"), type = "number")

  opts <- c("Chi-Square", "df", "p.value", "RMSEA", "Lower 90% CI RMSEA", "Upper 90% CI RMSEA", "SRMR")
  allData <- data.frame(
    measure = opts,
    value = as.vector(model[["scaleResults"]][["fit"]][["scaleOmega"]])
  )
  fitTable$setData(allData)


  fitTable$dependOn(options = c("variables", "scaleOmega", "reverseScaledItems", "omegaFitMeasures", "naAction",
                                "omegaEstimationMethod"))
  fitTable$position <- 3
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["fitTable"]] <- fitTable

  return()

}



# ----------------------------------------------
#       Common unidim. reliability functions
#       put here because there is space
# ----------------------------------------------

# This is a temporary fix
# TODO: remove it when R will solve this problem!
gettextf <- function(fmt, ..., domain = NULL)  {
  return(sprintf(gettext(fmt, domain = domain), ...))
}

.readData <- function(dataset, options) {

  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(
      columns.as.numeric  = variables,
      exclude.na.listwise = if (options[["naAction"]] == "listwise") variables else NULL
    )
  }
  return(dataset)
}


.checkErrors <- function(dataset, options, Bayes = FALSE) {

  if (Bayes) {
    # check for sensible MCMC values
    .checkMCMCBounds <- function() {
      # are the remaining samples after burnin and thinning larger than 2?
      if (ceiling((options[["samples"]] - options[["burnin"]]) / options[["thinning"]]) < 2) {
        return(gettext("Too few MCMC samples will remain after the removal of the burnin samples and thinning,
                     please increase No.samples."))
      }
      return(NULL)
    }


  } else {
    .checkMCMCBounds <- NULL
  }

  jaspBase::.hasErrors(dataset = dataset,
                       type = c("infinity", "variance", "observations"),
                       observations.amount = " < 3",
                       infinity.target = options$variables,
                       variance.target = options$variables,
                       observations.target = options$variables,
                       custom = .checkMCMCBounds,
                       exitAnalysisIfErrors = TRUE)
}


.checkLoadings <- function(dataset, variables) {
  if (ncol(dataset >= 2)) {
    # check for negative loadings:
    prin <- psych::principal(dataset)
    idx <- prin[["loadings"]] < 0
    sidx <- sum(idx)
    if (sidx == 0) {
      footnote <- ""
    } else {
      footnote <- sprintf(ngettext(sidx,
                                   "The following item correlated negatively with the scale: %s. ",
                                   "The following items correlated negatively with the scale: %s. "),
                          paste(variables[idx], collapse = ", "))
    }

    # check for perfect correlations:
    cr <- cor(dataset, use = "pairwise.complete.obs")
    cr[lower.tri(cr, diag = TRUE)] <- 0
    pos <- which(round(cr, 3) == 1, arr.ind = TRUE)
    if (length(pos) == 0) {
      footnote <- gettextf("%s", footnote)
    } else {
      for (i in seq_len(nrow(pos))) {
        footnote <- gettextf("%1$s Variables %2$s and %3$s correlated perfectly. ",
                             footnote, variables[pos[i, 1]], variables[pos[i, 2]])
      }
    }

    return(footnote)

  } else {
    return(.atLeast2Variables())
  }
}

.reverseScoreItems <- function(dataset, options) {
  dataset_rev <- as.matrix(dataset) # fails for string factors!
  cols <- match(unlist(options[["reverseScaledItems"]]), colnames(dataset))
  total <- apply(as.matrix(dataset[, cols]), 2, min, na.rm = TRUE) +
    apply(as.matrix(dataset[, cols]), 2, max, na.rm = TRUE)
  dataset_rev[, cols] <- matrix(rep(total, nrow(dataset)), nrow(dataset), length(cols), byrow = TRUE) -
    as.matrix(dataset[, cols])
  return(as.data.frame(dataset_rev))
}


.cov2cor.callback <- function(C, callback) {
  callback()
  return(cov2cor(C))
}

# calculate the kullback leibler distance between two samples
.KLD.statistic <- function(x, y) {
  # transform the samples to PDFs:
  xdf <- .get_approx_density(x)
  ydf <- .get_approx_density(y)

  xx <- seq(0, 1, length.out = 1e4)
  t <- LaplacesDemon::KLD(xdf(xx), ydf(xx))
  t$sum.KLD.py.px
}

# calculate the kolomogorov smirnov distances between some samples and the original sample
.ks.test.statistic <- function(x, y) {
  t <- stats::ks.test(x, y)
  t$statistic
}

# konvert empirical samples to cumulative density functions
.get_approx_density <- function(x) {
  d <- density(x, n = 2^12)
  f <- approxfun(d$x, d$y, yleft = 0, yright = 0)
  c <- integrate(f, 0, 1)$value
  return(
    function(x) {
      return(f(x) / c)
    }
  )
}



# change options when scale box is unchecked
.scaleItemBoxAlign <- function(options) {
  opts <- options
  if (!options[["scaleOmega"]])
    opts[["itemDeletedOmega"]] <- FALSE
  if (!options[["scaleAlpha"]])
    opts[["itemDeletedAlpha"]] <- FALSE
  if (!options[["scaleLambda2"]])
    opts[["itemDeletedLambda2"]] <- FALSE
  if (!options[["scaleLambda6"]])
    opts[["itemDeletedLambda6"]] <- FALSE
  if (!options[["scaleGreatestLowerBound"]])
    opts[["itemDeletedGreatestLowerBound"]] <- FALSE

  return(opts)

}

.addFootnoteReverseScaledItems <- function(options) {
  out <- sprintf(ngettext(length(options[["reverseScaledItems"]]),
                          "The following item was reverse scaled: %s. ",
                          "The following items were reverse scaled: %s. "),
                 paste(options[["reverseScaledItems"]], collapse = ", "))
  return(out)
}

.atLeast2Variables <- function() {
  return(gettext("Please enter at least 2 variables to do an analysis"))
}

.is.empty <- function(model) {
  !is.null(model[["empty"]])
}


.sampleListHelper <- function(ll, llName) {
  samps <- lapply(ll, `[[`, llName)
  samps[lengths(samps) == 0] <- NULL
  return(samps)
}

