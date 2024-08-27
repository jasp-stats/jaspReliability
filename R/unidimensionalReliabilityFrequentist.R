


#' @export
unidimensionalReliabilityFrequentist <- function(jaspResults, dataset, options) {

  sink(file = "~/Downloads/log.txt")
  on.exit(sink(NULL))

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

  model[["averageInterItemCorrelation"]] <- .frequentistAverageCor(jaspResults, dataset, options, model)
  model[["scaleMean"]] <- .frequentistMean(jaspResults, dataset, options, model)
  model[["scaleVar"]] <- .frequentistVar(jaspResults, dataset, options, model)
  model[["scaleSd"]] <- .frequentistStdDev(jaspResults, dataset, options, model)
  model[["itemRestCorrelation"]] <- .frequentistItemRestCor(jaspResults, dataset, options, model)
  model[["itemMean"]] <- .frequentistMeanItem(jaspResults, dataset, options, model)
  model[["itemVar"]] <- .frequentistVarItem(jaspResults, dataset, options, model)
  model[["itemSd"]] <- .frequentistSdItem(jaspResults, dataset, options, model)

  model[["scaleResults"]] <- .frequentistComputeScaleResults(jaspResults, dataset, options, model)

  .frequentistScaleTable(jaspResults, model, options)
  .frequentistItemTable(jaspResults, model, options)
  .frequentistSingleFactorFitTable(jaspResults, model, options)
  .frequentistLoadingsTable(jaspResults, model, options)


  return()

}

.frequentistDerivedOptions <- function(options) {

  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("scaleOmega", "scaleAlpha", "scaleLambda2",
                                           "averageInterItemCorrelation", "scaleMean", "scaleVar", "scaleSd")]),
    itemDroppedSelected = unlist(options[c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2",
                                           "itemRestCorrelation", "itemMean", "itemVar", "itemSd")]),
    namesEstimators     = list(
      tables = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2",
                 gettext("Average interitem correlation"), gettext("Mean"), gettext("Variance"), gettext("SD")),
      tables_item = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2",
                      gettext("Item-rest correlation"), gettext("Mean"), gettext("Variance"), gettext("SD")),
      coefficients = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2"))
  )
  return(derivedOptions)
}


.getStateContainerF <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = c("variables", "reverseScaledItems", "bootstrapSamples",
                                                                          "naAction", "bootstrapType", "setSeed",
                                                                          "seed", "ci", "samplesSavingDisabled", "intervalMethod")
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




###### Precalculate results #####

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
      (
        (options[["scaleOmega"]] && options[["omegaEstimationMethod"]] == "pfa") ||
        (options[["scaleAlpha"]] && options[["intervalMethod"]] == "bootstrapped") ||
        options[["scaleLambda2"]] ||
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




# ##### Frequentist calculate coefficients #####
#'
#' So the structure here is to first calculcate the bootstrap samples per coefficient
#' And in the computeScaleResults function to summarize the results


.frequentistOmegaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["scaleOmegaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["scaleOmegaObj"]]$object)

  out <- model[["scaleOmega"]]
  if (is.null(out))
    out <- list()

  if (options[["scaleOmega"]] && is.null(model[["empty"]])) {

    ciValue <- options[["ciLevel"]]
    if (options[["omegaEstimationMethod"]] == "cfa") {
      if (options[["intervalMethod"]] == "bootstrapped") {
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
    stateContainer[["scaleOmegaObj"]] <- createJaspState(out, dependencies = c("scaleOmega", "omegaEstimationMethod"))
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
      out[["itemDropped"]] <- c(NA, NA)
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
        out[["itemDropped"]] <- rep(NA, ncol(dataset))
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
      # should the interval be bootstrapped
      if (options[["intervalMethod"]] == "bootstrapped") {
        if (is.null(out[["samp"]])) {
          startProgressbar(options[["bootstrapSamples"]])
          out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applyalpha, callback = progressbarTick)
        }
      }

    } else { # alpha standardized
      # should the interval be bootstrapped
      if (options[["intervalMethod"]] == "bootstrapped") {
        if (is.null(out[["sampCor"]])) {
          out[["sampCor"]] <- numeric(options[["bootstrapSamples"]])
          for (i in seq_len(options[["bootstrapSamples"]])) {
            out[["sampCor"]][i] <- Bayesrel:::applyalpha(cov2cor(model[["bootSamp"]][i, , ]))
          }
        }
      }

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleAlphaObj"]] <- createJaspState(out, dependencies = c("scaleAlpha", "alphaType"))
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
      out[["itemDropped"]] <- c(NA, NA)
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

    if (is.null(out[["samp"]])) {
      startProgressbar(options[["bootstrapSamples"]])
      out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applylambda2, callback = progressbarTick)
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
      out[["itemDropped"]] <- c(NA, NA)
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


.frequentistAverageCor <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["avgObj"]]$object))
    return(.getStateContainerF(jaspResults)[["avgObj"]]$object)

  out <- model[["averageInterItemCorrelation"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["averageInterItemCorrelation"]] && is.null(model[["empty"]])) {
    if (is.null(out[["samp"]])) {
      startProgressbar(options[["bootstrapSamples"]])
      out[["samp"]] <- numeric(options[["bootstrapSamples"]])
      for (i in seq_len(options[["bootstrapSamples"]])) {
        corm <- .cov2cor.callback(model[["bootSamp"]][i, , ], progressbarTick)
        out[["samp"]][i] <- mean(corm[lower.tri(corm)])
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

    zz <- qnorm(1 - (1 - ciValue) / 2)
    out[["se"]] <- sdmean / sqrt(model[["n"]])
    out[["conf"]] <- c(out[["est"]] - zz * out[["se"]],
                       out[["est"]] + zz * out[["se"]])


    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("scaleMean", "meanSdScoresMethod"))
  }
  return(out)
}


.frequentistVar <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["varObj"]]$object))
    return(.getStateContainerF(jaspResults)[["varObj"]]$object)

  out <- model[["scaleVar"]]
  if (is.null(out))
    out <- list()
  if (options[["scaleVar"]] && is.null(model[["empty"]])) {
    ciValue <- options[["ciLevel"]]

    xx <- if (options[["meanSdScoresMethod"]] == "sumScores") {
      rowSums(dataset, na.rm = TRUE)
    } else {
      rowMeans(dataset, na.rm = TRUE)
    }

    out[["est"]] <- var(xx)
    out[["se"]] <- .seVar(x = xx)

    if (options[["intervalMethodVar"]] == "chisq") {
      chiValueLow <- qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
      out[["conf"]] <- c(((model[["n"]] - 1) * out[["est"]]) / chiValueLow,
                         ((model[["n"]] - 1) * out[["est"]]) / chiValueHigh)
    } else { # wald interval
      out[["conf"]] <- out[["est"]] + c(-1, 1) * qnorm(1 - (1 - ciValue) / 2) * out[["se"]]
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["varObj"]] <- createJaspState(out, dependencies = c("scaleVar", "meanSdScoresMethod", "intervalMethodVar"))
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

    xx <- if (options[["meanSdScoresMethod"]] == "sumScores") {
      rowSums(dataset, na.rm = TRUE)
    } else {
      rowMeans(dataset, na.rm = TRUE)
    }
    out[["est"]] <- sd(xx)
    out[["se"]] <- .seSd(x = xx)

    if (options[["intervalMethodVar"]] == "chisq") {
      chiValueLow <- qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
      out[["conf"]] <- c(sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueLow),
                         sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueHigh))
    } else {
      out[["conf"]] <- out[["est"]] + c(-1, 1) * qnorm(1 - (1 - ciValue) / 2) * out[["se"]]
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

.frequentistVarItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemVarObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemVarObj"]]$object)

  out <- model[["itemVar"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemVar"]]  && is.null(model[["empty"]])) {

    out[["itemDropped"]] <- apply(dataset, 2, var, na.rm = TRUE)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemVarObj"]] <- createJaspState(out, dependencies = "itemVar")
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
    # omega
    if (options[["scaleOmega"]]) {
      if (options[["omegaEstimationMethod"]] == "cfa") {
        datasetOmega <- scale(dataset, scale = FALSE)
        omegaO <- Bayesrel:::omegaFreqData(datasetOmega, interval = ciValue, omega.int.analytic = TRUE,
                                           pairwise = model[["pairwise"]])
        if (is.na(omegaO[["omega"]])) {
          out[["error"]][["scaleOmega"]] <- gettext("Omega calculation with CFA failed.
                                                    Try changing to PFA in Advanced Options")
          out[["est"]][["scaleOmega"]] <- NA
        } else {
          out[["fit"]][["scaleOmega"]] <- omegaO[["indices"]]
          out[["est"]][["scaleOmega"]] <- omegaO[["omega"]]

          if (options[["intervalMethod"]] == "analytic") {
            pp <- lavaan::parameterEstimates(omegaO$fit.object)
            out[["se"]][["scaleOmega"]] <- pp[pp$label == "omega", "se"]
            out[["conf"]][["scaleOmega"]] <- c(omegaO$omega_lower, omegaO$omega_upper)
          } else { # omega interval bootstrapped
            if (!is.null(model[["scaleOmega"]][["samp"]])) {
              if (sum(!is.na(model[["scaleOmega"]][["samp"]])) >= 2) {
                out[["conf"]][["scaleOmega"]] <- quantile(model[["scaleOmega"]][["samp"]],
                                                          probs = c((1 - ciValue)/2, 1 - (1 - ciValue) / 2),
                                                          na.rm = TRUE)
                out[["se"]][["scaleOmega"]] <- sd(model[["scaleOmega"]][["samp"]], na.rm = TRUE)
              } else {
                out[["error"]][["scaleOmega"]] <- gettext("Omega bootstrapped interval calculation with CFA failed.
                                                          Try changing to PFA in 'Advanced Options'")
                out[["conf"]][["scaleOmega"]] <- NA
              }
            }
          }

          if (options[["standardizedLoadings"]]) {

            lavObj <- omegaO$fit.object
            loads <- lavaan::inspect(lavObj, what = "std")$lambda
            out[["loadings"]] <- loads
          }
        }
      } else { # omega method is pfa
        omOut <- applyomegaPFA(model[["data_cov"]], loadings = TRUE)
        out[["est"]][["scaleOmega"]] <- omOut[["om"]]
        if (is.na(omOut[["om"]])) {
          out[["error"]][["scaleOmega"]] <- gettext("Omega calculation with PFA failed.")
          out[["est"]][["scaleOmega"]] <- NA
        } else {
          if (!is.null(model[["scaleOmega"]][["samp"]])) {
            if (sum(!is.na(model[["scaleOmega"]][["samp"]])) >= 2) {
              if (options[["intervalMethod"]] == "analytic") {
                out[["conf"]][["scaleOmega"]] <- c(NA, NA)
                out[["error"]][["scaleOmega"]] <- gettext("The analytic confidence interval is not available for coefficient omega obtained with PFA.")
              } else {
                out[["conf"]][["scaleOmega"]] <- quantile(model[["scaleOmega"]][["samp"]],
                                                          probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                          na.rm = TRUE)
                out[["se"]][["scaleOmega"]] <- sd(model[["scaleOmega"]][["samp"]], na.rm = TRUE)
              }

            } else {
              out[["error"]][["scaleOmega"]] <- gettext("Omega interval calculation with PFA failed.")
              out[["conf"]][["scaleOmega"]] <- NA
            }
          }

          if (options[["standardizedLoadings"]]) {
            out[["loadings"]] <- omOut[["loadings"]]
          }
        }
      }
    }

    # alpha
    if (options[["scaleAlpha"]]) {
      # alpha unstandardized
      if (options[["alphaType"]] == "unstandardized") {
        out[["est"]][["scaleAlpha"]] <- Bayesrel:::applyalpha(model[["data_cov"]])

        # should the interval be analytic
        if (options[["intervalMethod"]] == "analytic") {
          out[["se"]][["scaleAlpha"]] <- .seLambda3(dataset)
          out[["conf"]][["scaleAlpha"]] <- out[["est"]][["scaleAlpha"]] + c(-1, 1) * out[["se"]][["scaleAlpha"]] * qnorm(1 - (1 - ciValue) / 2)
        } else { # alpha interval bootstrapped
          if (!is.null(model[["scaleAlpha"]][["samp"]])) {
            out[["conf"]][["scaleAlpha"]] <- quantile(model[["scaleAlpha"]][["samp"]],
                                                      probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                      na.rm = TRUE)
            out[["se"]][["scaleAlpha"]] <- sd(model[["scaleAlpha"]][["samp"]], na.rm = TRUE)
          }
        }

      } else { # alpha standardized
        ccor <- model[["data_cor"]]
        out[["est"]][["scaleAlpha"]] <- Bayesrel:::applyalpha(ccor)

        # should the interval be analytic
        if (options[["intervalMethod"]] == "analytic") {

          # TODO: need to scale the dataset, but for now this leads to an error,
          # however, there is another PR with standardized coefficients anyways, so leave this for now.
          # out[["se"]][["scaleAlpha"]] <- .seLambda3(scale(dataset))
          out[["se"]][["scaleAlpha"]] <- NA
          out[["conf"]][["scaleAlpha"]] <- out[["est"]][["scaleAlpha"]] + c(-1, 1) * out[["se"]][["scaleAlpha"]] * qnorm(1 - (1 - ciValue) / 2)

        } else {
          if (!is.null(model[["scaleAlpha"]][["sampCor"]])) {
            out[["conf"]][["scaleAlpha"]] <- quantile(model[["scaleAlpha"]][["sampCor"]],
                                                      probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                      na.rm = TRUE)
            out[["se"]][["scaleAlpha"]] <- sd(model[["scaleAlpha"]][["sampCor"]], na.rm = TRUE)
          }
        }

      }
    }


    # lambda 2
    if (options[["scaleLambda2"]]) {
      out[["est"]][["scaleLambda2"]] <- Bayesrel:::applylambda2(model[["data_cov"]])
      if (options[["intervalMethod"]] == "bootstrapped") {
        samp <- model[["scaleLambda2"]][["samp"]]
        out[["conf"]][["scaleLambda2"]] <- quantile(samp, probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
        out[["se"]][["scaleLambda2"]] <- sd(samp, na.rm = TRUE)
      } else { # interval analytic
        out[["se"]][["scaleLambda2"]] <- .seLambda2(dataset)
        out[["conf"]][["scaleLambda2"]] <- out[["est"]][["scaleLambda2"]] + c(-1, 1) * out[["se"]][["scaleLambda2"]] * qnorm(1 - (1 - ciValue) / 2)
      }
    }

    # average interitem correlation
    if (options[["averageInterItemCorrelation"]]) {
      out[["est"]][["averageInterItemCorrelation"]] <- mean(model[["data_cor"]][lower.tri(model[["data_cor"]])])
      if (options[["intervalMethod"]] == "bootstrapped") {
        samp <- model[["averageInterItemCorrelation"]][["samp"]]
        out[["conf"]][["averageInterItemCorrelation"]] <- quantile(samp, probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
        out[["se"]][["averageInterItemCorrelation"]] <- sd(samp, na.rm = TRUE)
      } else { # interval analytic
        # TODO: what is the SE of the average interitem correlation?
        out[["se"]][["averageInterItemCorrelation"]] <- NA
        out[["conf"]][["averageInterItemCorrelation"]] <- out[["est"]][["averageInterItemCorrelation"]] + c(-1, 1) * out[["se"]][["averageInterItemCorrelation"]] * qnorm(1 - (1 - ciValue) / 2)
      }
    }

    # just copying for mean and sd
    if (options[["scaleMean"]]) {
      out[["est"]][["scaleMean"]] <- model[["scaleMean"]][["est"]]
      out[["se"]][["scaleMean"]] <- model[["scaleMean"]][["se"]]
      out[["conf"]][["scaleMean"]] <- model[["scaleMean"]][["conf"]]
    }

    # just copying for mean and sd
    if (options[["scaleVar"]]) {
      out[["est"]][["scaleVar"]] <- model[["scaleVar"]][["est"]]
      out[["se"]][["scaleVar"]] <- model[["scaleVar"]][["se"]]
      out[["conf"]][["scaleVar"]] <- model[["scaleVar"]][["conf"]]
    }

    if (options[["scaleSd"]]) {
      out[["est"]][["scaleSd"]] <- model[["scaleSd"]][["est"]]
      out[["se"]][["scaleSd"]] <- model[["scaleSd"]][["se"]]
      out[["conf"]][["scaleSd"]] <- model[["scaleSd"]][["conf"]]
    }


    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleResultsObj"]] <- createJaspState(out, dependencies = c("ciLevel",
                                                                                 "scaleMean", "scaleSd", "scaleVar",
                                                                                 "scaleAlpha", "scaleOmega",
                                                                                 "scaleLambda2",
                                                                                 "averageInterItemCorrelation",
                                                                                 "meanSdScoresMethod",
                                                                                 "omegaEstimationMethod", "intervalMethod",
                                                                                 "intervalMethodVar", "alphaType"))


  }

  return(out)

}




##### Output tables ####

.frequentistScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Frequentist Scale Reliability Statistics"))

  scaleTable$dependOn(options = c("scaleOmega", "scaleAlpha", "scaleLambda2",
                                  "averageInterItemCorrelation", "scaleMean", "scaleSd", "meanSdScoresMethod",
                                  "omegaEstimationMethod", "intervalMethod", "alphaType", "intervalMethodVar",
                                  "ciLevel"))
  scaleTable$addColumnInfo(name = "coefficient", title = gettext("Coefficient"), type = "string")
  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "number")
  scaleTable$addColumnInfo(name = "se", title = gettext("Std. error"), type = "number")
  ci <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)
  scaleTable$addColumnInfo(name = "lower", title = "Lower", type = "number", overtitle = gettextf("%s%% CI", ci))
  scaleTable$addColumnInfo(name = "upper", title = "Upper", type = "number", overtitle = gettextf("%s%% CI", ci))


  scaleTable$position <- 1
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  if (!is.null(model[["scaleResults"]][["error"]][["scaleOmega"]])) {
    model[["footnote"]] <- paste(model[["footnote"]], model[["scaleResults"]][["error"]][["scaleOmega"]])
  }

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

  if (!.is.empty(model) && any(selected)) {
    dt$estimate <- unlist(model[["scaleResults"]][["est"]], use.names = FALSE)
    dt$se <- unlist(model[["scaleResults"]][["se"]])
    dt$lower <- sapply(model[["scaleResults"]][["conf"]], function(x) x[1])
    dt$upper <- sapply(model[["scaleResults"]][["conf"]], function(x) x[2])

  }

  scaleTable$setData(dt)

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
  itemTable$dependOn(options = c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2",
                                 "itemMean", "itemRestCorrelation", "itemSd", "itemVar",
                                 "scaleOmega", "scaleAlpha", "scaleLambda2"))
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


.frequentistLoadingsTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["loadTable"]]$object) ||
      is.null(model[["scaleResults"]][["loadings"]]) || !options[["standardizedLoadings"]])
    return()

  loadTable <- createJaspTable(gettext("Standardized Loadings of the Single-Factor Model"))

  loadTable$dependOn(options = c("scaleOmega", "standardizedLoadings", "omegaEstimationMethod", "variables",
                                 "naAction", "reverseScaledItems"))

  loadTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")
  loadTable$addColumnInfo(name = "loadings", title = gettext("Standardized loading"), type = "number")

  if (options[["scaleOmega"]] && options[["standardizedLoadings"]] && is.null(model[["empty"]])) {

    df <- data.frame(
      variable = options[["variables"]],
      loadings = c(model[["scaleResults"]][["loadings"]]))
    loadTable$setData(df)

    loadTable$position <- 4
    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["loadTable"]] <- loadTable
  }

  return()
}



# ####  Some common and not so common helper functions ####
# unidim. reliability functions
# put here because there is space


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
    for (i in seq_len(nrow(pos))) {
      footnote <- gettextf("%1$sVariables %2$s and %3$s correlated perfectly. ",
                           footnote, variables[pos[i, 1]], variables[pos[i, 2]])
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
  return(gettext("Please enter at least 2 variables to do an analysis."))
}

.is.empty <- function(model) {
  !is.null(model[["empty"]])
}


.sampleListHelper <- function(ll, llName) {
  samps <- lapply(ll, `[[`, llName)
  samps[lengths(samps) == 0] <- NULL
  return(samps)
}


# have this here instead of the Bayesrel package
applyomegaPFA <- function(m, callback = function(){}, loadings = FALSE){

  f <- try(Bayesrel:::pfaArma(m), silent = TRUE)
  if (inherits(f, "try-error")) {
    om <- NA
    l_fa <- NA
    warning("singular bootstrapped covariance matrices encountered when computing omega")
  } else {
    l_fa <- f$loadings
    er_fa <- f$err_var
    om <- sum(l_fa)^2 / (sum(l_fa)^2 + sum(er_fa))
    if (om < 0 || om > 1 || is.na(om)) om <- NA

    if (loadings) {
      mm <- cov2cor(m)
      ff <- try(Bayesrel:::pfaArma(mm), silent = TRUE)
      l_fa <- ff$loadings
    }
  }

  callback()

  if (loadings) {
    return(list(om = om, loadings = l_fa))
  } else {
    return(om)
  }
}


#### Functions from Andries ####
# SE of sample variance
.seVar <- function(x, eps = 1e-16){
  N <- length(x)
  d <- ((x - mean(x))^2 - var(x)) / (N - 1)
  V <- sum(d^2) - sum(outer(d, d))/ N
  if (max(V) < eps) return(NA) else return(sqrt(V))
}

# SE of sample standard deviation
.seSd <- function(x, eps = 1e-16){
  return(.seVar(x)/ (2 * sd(x)))
}


# SE of sample lambda-2
.seLambda2 <- function(X, VC = NULL, eps = 1e-16){
  J <- ncol(X)
  if(is.null(VC)) VC <- .varCM(X)
  C <- var(X)
  C. <- C; diag(C.) <- 0
  vecC. <- matrix(C., nrow = 1)
  C2plus <- sum(C.^2)
  Cplus <- sum(C.)
  Splus <- sum(C)
  C2plusX <- sqrt(J / (J-1) * C2plus)

  e <-  matrix(diag(rep(1, J)), nrow = 1)
  u <-  matrix(1, nrow = 1, ncol = J^2)
  G <- (((C2plusX / C2plus * vecC. + (u - e)) - .lambda2(X)) / Splus)
  return(sqrt(G %*% VC %*% t(G)))
}


# SE of sample lambda-3 (alpha)
.seLambda3 <- function(X, VC = NULL, eps = 1e-16){
  J <- ncol(X)
  return((J / (J - 1)) *.seLambda1(X, VC))
}

# Sample split-half reliability coefficient: input requires 2 scores (sum scores of test halves)
.splithalf <- function(x, y) {
  return(2 * cor(x, y) / (1 + cor(x, y)))
}

# SE of sample split-half reliability coefficient: input requires 2 scores (sum scores of test halves)
.seSplithalf <- function(x, y) {
  k <- cor(x, y)
  seK <- seCor(x, y)
  sh <- 2 * k / (1 + k)
  return((sh/k - sh/(1 + k)) * seCor(x, y))
}

# SE of sample lambda-1
.seLambda1 <- function(X, VC = NULL, eps = 1e-16){
  D <- function(x) diag(as.numeric(x))
  J <- ncol(X)
  if(is.null(VC)) VC <- .varCM(X)
  g0 <- vecC <- matrix(var(X))
  e <-  matrix(diag(rep(1, J)), nrow = 1)
  u <-  matrix(1, nrow = 1, ncol = J^2)
  A1 <- rbind(e, u, u)
  A2 <- matrix(c(1, 0, -1, 1, 0, -1), 2, 3)
  A3 <- matrix(c(-1, 1), nrow = 1, ncol = 2)
  g1 <- log(A1 %*% g0 + eps)
  g2 <- exp(A2 %*% g1)
  g3 <- A3 %*% g2
  G1 <- D(1/(A1 %*% g0)) %*% A1
  G2 <- D(g2) %*% A2 %*% G1
  G <- A3 %*% G2
  V <- G %*% VC %*% t(G)
  return(sqrt(V))
}

# SE of sample covariance matrix
.varCM <- function(X, marg = marg.default, eps = 1e-16) {

  X2nR <- function(X){
    if (!is.matrix(X)) X <- as.matrix(X)
    n <- as.matrix(table(apply(X, 1, paste, collapse="-")))
    R <- matrix(as.numeric(unlist(strsplit(dimnames(n)[[1]], "-"))), ncol = ncol(X), byrow = TRUE)
    return(list(n = n, R = R))
  }
  D <- function(x) diag(as.numeric(x))
  Rij <- function(i, j, x) {
    if (i == j) {
      tmp <- expand.grid(x[[i]]);
      tmp <- unlist(tmp)
      tmp <- matrix(rep(rep(tmp, each = length(tmp)), 2), ncol = 2)
      return(tmp)
    } else {
      tmp <- expand.grid(x[[i]], x[[j]]);
      return(tmp[, ncol(tmp):1])
    }
  }

  res <- X2nR(X)
  R <- res$R
  n <- res$n
  J <- ncol(R)
  K <- length(n)
  marg.default <- list()
  for (j in 1 : J) marg.default[[j]] <- sort(unique(R[, j]))

  A2 <- matrix(c(1, 0, 0, 0,  1, 0, 0, 0,  0, 1, 0, 0,  -1, 0, 1, 1,  0, 0, 0, -1), 4, 5)
  A3 <- matrix(c(-1, 0,  1, 0,  0, 1,  0, -1), 2, 4)
  A4 <- matrix(c(1, -1), 1, 2)
  # nrowG <- .5 * J * (J+1)
  nrowG <- J^2
  G <- matrix(NA, nrowG, K)
  labG <- rep(NA, nrowG)
  iter <- i <- j <- 0
  for (i in (1:(J))) for (j in ((1):J)){  # all variances and covariances (once)
    iter <- iter + 1
    if (i <= j) {
      M2 <- NULL
      A1 <- NULL
      for (a in marg[[i]]) for (b in marg[[j]]) M2 <- rbind(M2,as.numeric(R[,i]==a & R[,j]==b))
      tmp <- Rij(i, j, marg)
      A1 <- t(cbind(tmp[, 1], tmp[, 2], tmp[, 1] * tmp[, 2], 1, 1))
      g1 <- log(A1 %*% M2 %*% n + eps)
      g2 <- exp(A2 %*% g1)
      A3g2 <- matrix(A3 %*% g2, nrow = 2)
      g4 <- matrix(A3g2[1, ] / A3g2[2, ])
      G[(i - 1) * J + j, ]  <- g4 %*% A4 %*% D(1/(A3 %*% g2 + eps)) %*% A3 %*% D(g2) %*% A2 %*% D(1/(A1 %*% M2 %*% n + eps)) %*% A1 %*% M2
      labG[(i - 1) * J + j] <- paste(i, ",", j, sep = "")
      if (i < j) {
        G[(j - 1) * J + i, ] <- G[(i - 1) * J + j, ]
        labG[(j - 1) * J + i] <- paste(j, ",", i, sep = "")
      }
    }
  }
  dimnames(G) <- list(labG, dimnames(n)[[1]])
  Vn <- D(n) - n %*% t(n) / sum(n)
  VC <- G %*% Vn %*% t(G) - (G %*% n %*% t(n) %*% t(G))/sum(n)
  return(VC)
}

# Sample lambda-2
.lambda2 <- function(X){
  J <- ncol(X)
  C <- var(X)
  return(((sum(C)-sum(diag(C))) + (sqrt((J/(J-1))*(sum(C^2)-sum(diag(C^2))))))/sum(C))
}
