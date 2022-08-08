


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
      options[["intervalOn"]] &&
      (
        (options[["omegaScale"]] && options[["omegaMethod"]] == "pfa") ||
        (options[["alphaScale"]] && !(options[["alphaInterval"]] == "alphaAnalytic")) ||
        options[["lambda2Scale"]] ||
        options[["lambda6Scale"]] ||
        options[["glbScale"]] ||
        options[["averageInterItemCor"]]
      )
  ) {


    boot_cov <- array(0, c(options[["noSamples"]], k, k))

    startProgressbar(options[["noSamples"]])

    if (options[["bootType"]] == "parametric") {
      model[["parametric"]] <- TRUE
      for (i in seq_len(options[["noSamples"]])) {
        boot_data <- MASS::mvrnorm(n, colMeans(dataset, na.rm = TRUE), cc)
        boot_cov[i, , ] <- cov(boot_data)
        progressbarTick()
      }
    } else {
      model[["parametric"]] <- FALSE
      for (i in seq_len(options[["noSamples"]])) {
        boot_data <- as.matrix(dataset[sample.int(n, size = n, replace = TRUE), ])
        boot_cov[i, , ] <- cov(boot_data, use = model[["use.cases"]])
        progressbarTick()
      }
    }
    model[["bootSamp"]] <- boot_cov
  }

  if (options[["disableSampleSave"]])
    return(model)

  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model)

  return(model)
}




# -------------------------------------------
#       Frequentist calculate coefficients
# -------------------------------------------


.frequentistOmegaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["omegaScaleObj"]]$object))
    return(.getStateContainerF(jaspResults)[["omegaScaleObj"]]$object)

  out <- model[["omegaScale"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaScale"]] && is.null(model[["empty"]]) && options[["intervalOn"]]) {

    ciValue <- options[["confidenceIntervalValue"]]
    if (options[["omegaMethod"]] == "cfa") {
      if (options[["omegaInterval"]] == "omegaBoot") {
        parametric <- options[["bootType"]] == "parametric"
        omegaboot <- out[["samp"]]
        if (is.null(omegaboot)) {
          startProgressbar(options[["noSamples"]])
          jaspBase::.setSeedJASP(options)
          omegaboot <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = FALSE,
                                                pairwise = model[["pairwise"]], parametric = parametric,
                                                n.boot = options[["noSamples"]], callback = progressbarTick)
          out[["samp"]] <- omegaboot[["omega_boot"]]
        }
      }
    } else { # omega with pfa
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["noSamples"]])
        out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applyomegaPFA, callback = progressbarTick)
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["omegaScaleObj"]] <- createJaspState(out, dependencies = c("omegaScale", "omegaMethod",
                                                                               "omegaInterval"))
  }

  return(out)
}

.frequentistOmegaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["omegaItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["omegaItemObj"]]$object)

  out <- model[["omegaItem"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaItem"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    if (options[["omegaMethod"]] == "cfa") {

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

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["omegaItemObj"]] <- createJaspState(out, dependencies = c("omegaItem", "omegaMethod"))
  }

  return(out)
}


.frequentistAlphaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["alphaScaleObj"]]$object))
    return(.getStateContainerF(jaspResults)[["alphaScaleObj"]]$object)

  out <- model[["alphaScale"]]
  if (is.null(out))
    out <- list()

  if (options[["alphaScale"]] && is.null(model[["empty"]])) {

    # alpha unstandardized
    if (options[["alphaMethod"]] == "alphaUnstand") {
      # do we need an interval estimate?
      if (options[["intervalOn"]]) {
        # should the interval be bootstrapped
        if (options[["alphaInterval"]] == "alphaBoot") {
          if (is.null(out[["samp"]])) {
            startProgressbar(options[["noSamples"]])
            out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applyalpha, callback = progressbarTick)
          }
        }
      }
    } else { # alpha standardized
      # do we need an interval estimate?
      if (options[["intervalOn"]]) {
        # should the interval be bootstrapped
        if (options[["alphaInterval"]] == "alphaBoot") {
          if (is.null(out[["sampCor"]])) {
            out[["sampCor"]] <- numeric(options[["noSamples"]])
            for (i in seq_len(options[["noSamples"]])) {
              out[["sampCor"]][i] <- Bayesrel:::applyalpha(cov2cor(model[["bootSamp"]][i, , ]))
            }
          }
        }
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["alphaScaleObj"]] <- createJaspState(out, dependencies = c("alphaScale", "alphaMethod",
                                                                               "alphaInterval"))
  }
  return(out)
}

.frequentistAlphaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["alphaItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["alphaItemObj"]]$object)

  out <- model[["alphaItem"]]
  if (is.null(out))
    out <- list()

  if (options[["alphaItem"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    if (options[["alphaMethod"]] == "alphaUnstand") { # alpha unstandardized
      # do we have to compute item dropped values
      if (is.null(out[["itemDropped"]]))
        out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applyalpha)

    } else { # alpha standardized
      if (is.null(out[["itemDropped"]])) {
        out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cor"]], Bayesrel:::applyalpha)
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["alphaItemObj"]] <- createJaspState(out, dependencies = c("alphaItem", "alphaMethod"))
  }
  return(out)
}


.frequentistLambda2Scale <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["lambda2ScaleObj"]]$object))
    return(.getStateContainerF(jaspResults)[["lambda2ScaleObj"]]$object)

  out <- model[["lambda2Scale"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["lambda2Scale"]]  && is.null(model[["empty"]])) {
    # do we need an interval estimate?
    if (options[["intervalOn"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["noSamples"]])
        out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applylambda2, callback = progressbarTick)
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda2ScaleObj"]] <- createJaspState(out, dependencies = "lambda2Scale")
  }
  return(out)
}

.frequentistLambda2Item <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["lambda2ItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["lambda2ItemObj"]]$object)

  out <- model[["lambda2Item"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["lambda2Item"]]  && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NaN, NaN)
      return(out)
    }

    if (is.null(out[["itemDropped"]]))
      out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applylambda2)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda2ItemObj"]] <- createJaspState(out, dependencies = "lambda2Item")
  }
  return(out)
}


.frequentistLambda6Scale <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["lambda6ScaleObj"]]$object))
    return(.getStateContainerF(jaspResults)[["lambda6ScaleObj"]]$object)

  out <- model[["lambda6Scale"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["lambda6Scale"]]  && is.null(model[["empty"]])) {
    if (options[["intervalOn"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["noSamples"]])
        out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applylambda6, callback = progressbarTick)
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda6ScaleObj"]] <- createJaspState(out, dependencies = "lambda6Scale")
  }
  return(out)
}

.frequentistLambda6Item <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["lambda6ItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["lambda6ItemObj"]]$object)

  out <- model[["lambda6Item"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["lambda6Item"]]  && is.null(model[["empty"]])) {

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
    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda6ItemObj"]] <- createJaspState(out, dependencies = "lambda6Item")
  }
  return(out)
}


# check the error handling of the glb !!!!!!!!
.frequentistGlbScale <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["glbScaleObj"]]$object))
    return(.getStateContainerF(jaspResults)[["glbScaleObj"]]$object)

  out <- model[["glbScale"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["glbScale"]]  && is.null(model[["empty"]])) {
    # do we need an interval estimate?
    if (options[["intervalOn"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(4)
        out[["samp"]] <- Bayesrel:::glbOnArrayCustom(model[["bootSamp"]], callback = progressbarTick)
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["glbScaleObj"]] <- createJaspState(out, dependencies = "glbScale")
  }
  return(out)
}

.frequentistGlbItem <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["glbItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["glbItemObj"]]$object)

  out <- model[["glbItem"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["glbItem"]]  && is.null(model[["empty"]])) {

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

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["glbItemObj"]] <- createJaspState(out, dependencies = "glbItem")
  }
  return(out)
}


.frequentistAverageCor <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["avgObj"]]$object))
    return(.getStateContainerF(jaspResults)[["avgObj"]]$object)

  out <- model[["averageInterItemCor"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["averageInterItemCor"]] && is.null(model[["empty"]])) {
    if (options[["intervalOn"]]) {
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["noSamples"]])
        out[["samp"]] <- numeric(options[["noSamples"]])
        for (i in seq_len(options[["noSamples"]])) {
          corm <- .cov2cor.callback(model[["bootSamp"]][i, , ], progressbarTick)
          out[["samp"]][i] <- mean(corm[lower.tri(corm)])
        }
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["avgObj"]] <- createJaspState(out,
                                                  dependencies = c("averageInterItemCor"))
  }
  return(out)
}

.frequentistMean <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["meanObj"]]$object))
    return(.getStateContainerF(jaspResults)[["meanObj"]]$object)

  out <- model[["meanScale"]]
  if (is.null(out))
    out <- list()
  if (options[["meanScale"]] && is.null(model[["empty"]])) {
    ciValue <- options[["confidenceIntervalValue"]]

    if (options[["scoresMethod"]] == "sumScores") {
      out[["est"]] <- mean(rowSums(dataset, na.rm = TRUE))
      sdmean <- sd(rowSums(dataset, na.rm = TRUE))
    } else {
      out[["est"]] <- mean(rowMeans(dataset, na.rm = TRUE))
      sdmean <- sd(rowMeans(dataset, na.rm = TRUE))
    }

    if (options[["intervalOn"]]) {
      zz <- qnorm(1 - (1 - ciValue) / 2)
      out[["conf"]] <- c(out[["est"]] - zz * (sdmean / sqrt(model[["n"]])),
                         out[["est"]] + zz * (sdmean / sqrt(model[["n"]])))
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("meanScale", "scoresMethod"))
  }
  return(out)
}

.frequentistStdDev <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["sdObj"]]$object))
    return(.getStateContainerF(jaspResults)[["sdObj"]]$object)

  out <- model[["sdScale"]]
  if (is.null(out))
    out <- list()
  if (options[["sdScale"]] && is.null(model[["empty"]])) {
    ciValue <- options[["confidenceIntervalValue"]]

    out[["est"]] <- if (options[["scoresMethod"]] == "sumScores")
      sd(rowSums(dataset, na.rm = TRUE))
    else
      sd(rowMeans(dataset, na.rm = TRUE))

    if (options[["intervalOn"]]) {
      chiValueLow <- qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
      out[["conf"]] <- c(sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueLow),
                         sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueHigh))
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["sdObj"]] <- createJaspState(out, dependencies = c("sdScale", "scoresMethod"))
  }
  return(out)
}


.frequentistItemRestCor <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemRestObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemRestObj"]]$object)

  out <- model[["itemRestCor"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemRestCor"]]  && is.null(model[["empty"]])) {
    out[["itemDropped"]] <- numeric(ncol(dataset))
    for (i in seq_len(ncol(dataset))) {
      out[["itemDropped"]][i] <- cor(as.matrix(dataset[, i]), rowMeans(as.matrix(dataset[, -i]), na.rm = TRUE),
                                     use = model[["use.cases"]])
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = "itemRestCor")
  }
  return(out)
}

.frequentistMeanItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["meanItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["meanItemObj"]]$object)

  out <- model[["meanItem"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["meanItem"]]  && is.null(model[["empty"]])) {
    out[["itemDropped"]] <- colMeans(dataset, na.rm = TRUE)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["meanItemObj"]] <- createJaspState(out, dependencies = "meanItem")
  }
  return(out)
}

.frequentistSdItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["sdItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["sdItemObj"]]$object)

  out <- model[["sdItem"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["sdItem"]]  && is.null(model[["empty"]])) {

    out[["itemDropped"]] <- apply(dataset, 2, sd, na.rm = TRUE)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["sdItemObj"]] <- createJaspState(out, dependencies = "sdItem")
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

    ciValue <- options[["confidenceIntervalValue"]]

    # go one coefficient at a time, because there are too many special options for a generic solution
    # #### omega ####
    if (options[["omegaScale"]]) {
      if (options[["omegaMethod"]] == "cfa") {
        dataset <- scale(dataset, scale = FALSE)
        omegaO <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = TRUE,
                                           pairwise = model[["pairwise"]])
        if (is.na(omegaO[["omega"]])) {
          out[["error"]][["omegaScale"]] <- gettext("Omega calculation with CFA failed.
                                                    Try changing to PFA in Advanced Options")
          out[["est"]][["omegaScale"]] <- NaN
        } else {
          out[["fit"]][["omegaScale"]] <- omegaO[["indices"]]
          out[["est"]][["omegaScale"]] <- omegaO[["omega"]]

          if (options[["intervalOn"]]) {
            if (options[["omegaInterval"]] == "omegaAnalytic") {
              out[["conf"]][["omegaScale"]] <- c(omegaO$omega_lower, omegaO$omega_upper)
            } else { # omega interval bootstrapped
              if (!is.null(model[["omegaScale"]][["samp"]])) {
                if (sum(!is.na(model[["omegaScale"]][["samp"]])) >= 2) {
                  out[["conf"]][["omegaScale"]] <- quantile(model[["omegaScale"]][["samp"]],
                                                            probs = c((1 - ciValue)/2, 1 - (1 - ciValue) / 2),
                                                            na.rm = TRUE)
                } else {
                  out[["error"]][["omegaScale"]] <- gettext("Omega bootstrapped interval calculation with CFA failed.
                                                            Try changing to PFA in 'Advanced Options'")
                  out[["conf"]][["omegaScale"]] <- NaN
                }
              }
            }
          }
        }
      } else { # omega method is pfa
        out[["est"]][["omegaScale"]] <- Bayesrel:::applyomegaPFA(model[["data_cov"]])
        if (is.na(out[["est"]][["omegaScale"]])) {
          out[["error"]][["omegaScale"]] <- gettext("Omega calculation with PFA failed.")
          out[["est"]][["omegaScale"]] <- NaN
        } else {
          if (options[["intervalOn"]]) {
            if (!is.null(model[["omegaScale"]][["samp"]])) {
              if (sum(!is.na(model[["omegaScale"]][["samp"]])) >= 2) {
                out[["conf"]][["omegaScale"]] <- quantile(model[["omegaScale"]][["samp"]],
                                                          probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                          na.rm = TRUE)
              } else {
                out[["error"]][["omegaScale"]] <- gettext("Omega interval calculation with PFA failed.")
                out[["conf"]][["omegaScale"]] <- NaN
              }
            }
          }
        }
      }
    }

    # #### alpha ####
    if (options[["alphaScale"]]) {
      # alpha unstandardized
      if (options[["alphaMethod"]] == "alphaUnstand") {
        out[["est"]][["alphaScale"]] <- Bayesrel:::applyalpha(model[["data_cov"]])
        if (options[["intervalOn"]]) {
          # should the interval be analytic
          if (options[["alphaInterval"]] == "alphaAnalytic") {
            out[["conf"]][["alphaScale"]] <- Bayesrel:::ciAlpha(1 - ciValue, model[["n"]], model[["data_cov"]])
          } else { # alpha interval bootstrapped
            if (!is.null(model[["alphaScale"]][["samp"]])) {
              out[["conf"]][["alphaScale"]] <- quantile(model[["alphaScale"]][["samp"]],
                                                        probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                        na.rm = TRUE)
            }
          }
        }
      } else { # alpha standardized
        ccor <- model[["data_cor"]]
        out[["est"]][["alphaScale"]] <- Bayesrel:::applyalpha(ccor)
        if (options[["intervalOn"]]) {
          # should the interval be analytic
          if (options[["alphaInterval"]] == "alphaAnalytic") {
            out[["conf"]][["alphaScale"]] <- Bayesrel:::ciAlpha(1 - ciValue, model[["n"]], ccor)
          } else {
            if (!is.null(model[["alphaScale"]][["sampCor"]])) {
              out[["conf"]][["alphaScale"]] <- quantile(model[["alphaScale"]][["sampCor"]],
                                                        probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                        na.rm = TRUE)
            }
          }
        }
      }
    }

    # #### lambda 6 ####
    if (options[["lambda6Scale"]]) {
      out[["est"]][["lambda6Scale"]] <- Bayesrel:::applylambda6(model[["data_cov"]])
      if (is.na(out[["est"]][["lambda6Scale"]])) {
        out[["error"]][["lambda6Scale"]] <- gettext("Lambda6 calculation failed.")
      }
      if (options[["intervalOn"]]) {
        if (sum(!is.na(model[["lambda6Scale"]][["samp"]])) >= 2) {
          out[["conf"]][["lambda6Scale"]] <- quantile(model[["lambda6Scale"]][["samp"]],
                                                      probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
        } else {
          out[["error"]][["lambda6Scale"]] <- gettext("Lambda6 interval calculation failed.")
        }

      }
    }

    if (options[["intervalOn"]]) {
      # the intervals for coefficients with bootstrap samples can be done generic
      bootCoeffs <- c("lambda2Scale", "glbScale", "averageInterItemCor")
      selected <- bootCoeffs[bootCoeffs %in% names(which(model[["derivedOptions"]][["selectedEstimators"]]))]

      sampellist <- model[selected]
      samps <- .sampleListHelper(sampellist, "samp")
      out[["conf"]][selected] <- lapply(samps, function(x) {quantile(x, probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                                     na.rm = TRUE)})
    }

    # point estimates
    if (options[["lambda2Scale"]]) {
      out[["est"]][["lambda2Scale"]] <- Bayesrel:::applylambda2(model[["data_cov"]])
    }
    if (options[["glbScale"]]) {
      out[["est"]][["glbScale"]] <- Bayesrel:::glbOnArrayCustom(model[["data_cov"]])
    }
    if (options[["averageInterItemCor"]]) {
      out[["est"]][["averageInterItemCor"]] <- mean(model[["data_cor"]][lower.tri(model[["data_cor"]])])
    }

    # just copying for mean and sd
    if (options[["meanScale"]]) {
      out[["est"]][["meanScale"]] <- model[["meanScale"]][["est"]]
      out[["conf"]][["meanScale"]] <- model[["meanScale"]][["conf"]]
    }
    if (options[["sdScale"]]) {
      out[["est"]][["sdScale"]] <- model[["sdScale"]][["est"]]
      out[["conf"]][["sdScale"]] <- model[["sdScale"]][["conf"]]
    }


    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleResultsObj"]] <- createJaspState(out, dependencies = c("confidenceIntervalValue",
                                                                                 "meanScale", "sdScale",
                                                                                 "alphaScale", "omegaScale",
                                                                                 "lambda2Scale", "lambda6Scale",
                                                                                 "glbScale","averageInterItemCor",
                                                                                 "scoresMethod",
                                                                                 "omegaMethod", "omegaInterval",
                                                                                 "alphaInterval", "alphaMethod"))


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

  scaleTable$dependOn(options = c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                  "averageInterItemCor", "meanScale", "sdScale", "scoresMethod",
                                  "omegaMethod", "omegaInterval", "alphaMethod", "alphaInterval"))
  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  if (options[["intervalOn"]]) {
    interval <- gettextf("%s%% CI",
                         format(100 * options[["confidenceIntervalValue"]], digits = 3, drop0trailing = TRUE))
    intervalLow <- gettextf("%s lower bound", interval)
    intervalUp <- gettextf("%s upper bound", interval)
    allData <- data.frame(estimate = c(gettext("Point estimate"), intervalLow, intervalUp))
  } else {
    allData <- data.frame(estimate = c(gettext("Point estimate")))
  }

  scaleTable$position <- 1
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  if (!is.null(model[["scaleResults"]][["error"]][["omegaScale"]])) {
    model[["footnote"]] <- paste(model[["footnote"]], model[["scaleResults"]][["error"]][["omegaScale"]])
  }
  if (!is.null(model[["scaleResults"]][["error"]][["lambda6Scale"]])) {
    model[["footnote"]] <- paste(model[["footnote"]], model[["scaleResults"]][["error"]][["lambda6Scale"]])
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

  if (options[["intervalOn"]]) {
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
  itemTable$dependOn(options = c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item", "glbItem",
                                 "meanItem", "itemRestCor", "sdItem",
                                 "omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale"))
  # adding the scale options fixes a bug, where the item table would remain displayed
  # after one had checked a scale coefficient box and the item coefficient box and then unchecked the scale coeff box

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemTable$position <- 2
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  footnote <- ""

  if (!is.null(model[["omegaItem"]][["error"]])) {
    footnote <- paste(footnote, model[["omegaItem"]][["error"]])
  }
  if (!is.null(model[["lambda6Item"]][["error"]])) {
    footnote <- paste(footnote, model[["omegaItem"]][["error"]])
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
      itemTable$addColumnInfo(name = paste0("pointEst", i), title = coefficientsTable[i], type = "number",
                              overtitle = overTitle)
      fewItemProblem <- TRUE
    } else {
      itemTable$addColumnInfo(name = paste0("pointEst", i), title = coefficientsTable[i], type = "number")
    }
  }

  if (is.null(model[["empty"]])) {
    tb <- data.frame(variable = model[["itemsDropped"]])

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])
      newtb <- cbind(pointEst = model[[nm]][["itemDropped"]])
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
      is.null(model[["scaleResults"]][["fit"]][["omegaScale"]]) || !options[["fitMeasures"]])
    return()

  fitTable <- createJaspTable(gettextf("Fit Measures of Single Factor Model Fit"))

  fitTable$addColumnInfo(name = "measure", title = gettext("Fit measure"),   type = "string")
  fitTable$addColumnInfo(name = "value",  title = gettext("Value"), type = "number")

  opts <- c("Chi-Square", "df", "p.value", "RMSEA", "Lower 90% CI RMSEA", "Upper 90% CI RMSEA", "SRMR")
  allData <- data.frame(
    measure = opts,
    value = as.vector(model[["scaleResults"]][["fit"]][["omegaScale"]])
  )
  fitTable$setData(allData)


  fitTable$dependOn(options = c("variables", "omegaScale", "reverseScaledItems", "fitMeasures", "missingValues",
                                "omegaMethod"))
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
      exclude.na.listwise = if (options[["missingValues"]] == "excludeCasesListwise") variables else NULL
    )
  }
  return(dataset)
}


.checkErrors <- function(dataset, options, Bayes = FALSE) {

  if (Bayes) {
    # check for sensible MCMC values
    .checkMCMCBounds <- function() {
      # are the remaining samples after burnin and thinning larger than 2?
      if (ceiling((options[["noSamples"]] - options[["noBurnin"]]) / options[["noThin"]]) < 2) {
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
        footnote <- gettextf("%s Variables %s and %s correlated perfectly. ",
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
  if (!options[["omegaScale"]])
    opts[["omegaItem"]] <- FALSE
  if (!options[["alphaScale"]])
    opts[["alphaItem"]] <- FALSE
  if (!options[["lambda2Scale"]])
    opts[["lambda2Item"]] <- FALSE
  if (!options[["lambda6Scale"]])
    opts[["lambda6Item"]] <- FALSE
  if (!options[["glbScale"]])
    opts[["glbItem"]] <- FALSE

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

