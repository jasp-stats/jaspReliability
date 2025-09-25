
#' @export
unidimensionalReliabilityFrequentist <- function(jaspResults, dataset, options) {

  # check for listwise deletion
  datasetOld <- dataset
  dataset <- .handleData(datasetOld, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }

  .checkErrors(dataset, options)

  model <- .frequentistPreCalc(jaspResults, dataset, options, datasetOld)
  options <- .scaleItemBoxAlign(options)

  model[["derivedOptions"]] <- .frequentistDerivedOptions(options)
  model[["bootCor"]] <- .frequentistStdCov(jaspResults, dataset, options, model)
  model[["scaleOmega"]] <- .frequentistOmegaScale(jaspResults, dataset, options, model)
  model[["itemDeletedOmega"]] <- .frequentistOmegaItem(jaspResults, dataset, options, model)
  model[["scaleAlpha"]] <- .frequentistAlphaScale(jaspResults, dataset, options, model)
  model[["itemDeletedAlpha"]] <- .frequentistAlphaItem(jaspResults, dataset, options, model)
  model[["scaleLambda2"]] <- .frequentistLambda2Scale(jaspResults, dataset, options, model)
  model[["itemDeletedLambda2"]] <- .frequentistLambda2Item(jaspResults, dataset, options, model)
  model[["scaleSplithalf"]] <- .frequentistSplithalfScale(jaspResults, dataset, options, model)
  model[["itemDeletedSplithalf"]] <- .frequentistSplithalfItem(jaspResults, dataset, options, model)
  model[["averageInterItemCorrelation"]] <- .frequentistAverageCor(jaspResults, dataset, options, model)
  model[["scaleMean"]] <- .frequentistMean(jaspResults, dataset, options, model)
  model[["scaleVar"]] <- .frequentistVar(jaspResults, dataset, options, model)
  model[["scaleSd"]] <- .frequentistStdDev(jaspResults, dataset, options, model)
  model[["itemRestCorrelation"]] <- .frequentistItemRestCor(jaspResults, dataset, options, model)
  model[["itemMean"]] <- .frequentistMeanItem(jaspResults, dataset, options, model)
  model[["itemVar"]] <- .frequentistVarItem(jaspResults, dataset, options, model)
  model[["itemSd"]] <- .frequentistSdItem(jaspResults, dataset, options, model)

  model[["scaleResults"]] <- .frequentistComputeScaleResults(jaspResults, dataset, options, model)
  model[["itemResults"]] <- .frequentistComputeItemResults(jaspResults, dataset, options, model)

  .frequentistScaleTable(jaspResults, model, options)
  .frequentistItemTable(jaspResults, model, options)
  .frequentistSingleFactorFitTable(jaspResults, model, options)
  .frequentistLoadingsTable(jaspResults, model, options)

  return()

}

.frequentistDerivedOptions <- function(options) {

  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("scaleOmega", "scaleAlpha", "scaleLambda2", "scaleSplithalf",
                                           "averageInterItemCorrelation", "scaleMean", "scaleVar", "scaleSd")]),
    itemDroppedSelected = unlist(options[c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2", "itemDeletedSplithalf",
                                           "itemRestCorrelation", "itemMean", "itemVar", "itemSd")]),
    namesEstimators     = list(
      tables = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2", gettext("Split-half coefficient"),
                 gettext("Average interitem correlation"), gettext("Mean"), gettext("Variance"), gettext("SD")),
      tables_item = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2", gettext("Split-half coefficient"),
                      gettext("Item-rest correlation"), gettext("Mean"), gettext("Variance"), gettext("SD")),
      coefficientsDeleted = c("Coefficient \u03C9", "Coefficient \u03B1", "Guttman's \u03BB2", gettext("Split-half coefficient")))
  )

  return(derivedOptions)
}


.getStateContainerF <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = c("variables", "reverseScaledItems", "bootstrapSamples",
                                                                          "naAction", "bootstrapType", "setSeed",
                                                                          "seed", "ciLevel", "samplesSavingDisabled", "intervalMethod")
  )

  return(jaspResults[["stateContainer"]])
}



###### Precalculate results #####

.frequentistPreCalc <- function(jaspResults, dataset, options, datasetOld) {

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
    if (options[["naAction"]] == "listwise" && nrow(datasetOld) > nrow(dataset)) { # this indicates listwise deletion
      model[["use.cases"]] <- "complete.obs"
      model[["pairwise"]] <- FALSE
      model[["footnote"]] <- gettextf("%1$s Of the observations, %2$1.f complete cases were used. ",
                                      model[["footnote"]], nrow(dataset))

    } else if (anyNA(dataset)) { # when pairwise deletion
      model[["use.cases"]] <- "pairwise.complete.obs"
      model[["pairwise"]] <- TRUE

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
    model[["footnote"]] <- gettextf("%1$s%2$s", model[["footnote"]], .checkLoadings(dataset, options[["variables"]]))
  }

  jaspBase::.setSeedJASP(options)

  # check if interval is checked and bootstrapped covariance sample has to be generated
  if (is.null(model[["bootSamp"]]) &&
      (options[["intervalMethod"]] == "bootstrapped")) {

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


.frequentistStdCov <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["bootCor"]]$object))
    return(.getStateContainerF(jaspResults)[["bootCor"]]$object)

  if (is.null(model[["bootSamp"]])) {
    return()
  } else {
    if (options[["coefficientType"]] == "unstandardized") {
      return()
    } else { # standardized
      out <- model[["bootSamp"]]
      startProgressbar(options[["bootstrapSamples"]])
      for (i in seq_len(nrow(model[["bootSamp"]]))) {
        out[i, , ] <- .cov2cor.callback(model[["bootSamp"]][i, , ], progressbarTick)
      }

      if (options[["samplesSavingDisabled"]])
        return(out)

      stateContainer <- .getStateContainerF(jaspResults)
      stateContainer[["bootCor"]] <- createJaspState(out, dependencies = "coefficientType")
    }
  }
  return(out)
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
        startProgressbar(options[["bootstrapSamples"]])
        type <- ifelse(options[["coefficientType"]] == "unstandardized", "bootSamp", "bootCor")
        missing <- ifelse(options[["naAction"]] == "pairwise", "ml", "listwise")
        jaspBase::.setSeedJASP(options)
        out[["samp"]] <- apply(model[[type]], 1, .applyOmegaCov, missing = missing, callback = progressbarTick)

      }
    } else { # omega with pfa
      # pfa does not work with analytic interval
      if (options[["intervalMethod"]] == "bootstrapped") {
        if (is.null(out[["samp"]])) {
          type <- ifelse(options[["coefficientType"]] == "unstandardized", "bootSamp", "bootCor")
          startProgressbar(options[["bootstrapSamples"]])
          out[["samp"]] <- apply(model[[type]], 1, .applyOmegaPFA, callback = progressbarTick)
        }
      }

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleOmegaObj"]] <- createJaspState(out, dependencies = c("scaleOmega", "omegaEstimationMethod",
                                                                               "coefficientType", "intervalMethod"))
  }

  return(out)
}


.frequentistOmegaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedOmegaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedOmegaObj"]]$object)

  out <- model[["itemDeletedOmega"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedOmega"]] && is.null(model[["empty"]]) && options[["intervalMethod"]] == "bootstrapped") {

    startProgressbar(options[["bootstrapSamples"]] * ncol(dataset))
    type <- ifelse(options[["coefficientType"]] == "unstandardized", "bootSamp", "bootCor")
    missing <- ifelse(options[["naAction"]] == "pairwise", "ml", "listwise")
    jaspBase::.setSeedJASP(options)

    if (options[["omegaEstimationMethod"]] == "cfa") {

      out[["itemSamp"]] <- .frequentistItemDroppedStats(covSamp = model[[type]],
                                                        f1 = .applyOmegaCov,
                                                        missing = missing,
                                                        callback = progressbarTick)

    } else { # omega with pfa
      out[["itemSamp"]] <- .frequentistItemDroppedStats(covSamp = model[[type]],
                                                        f1 = .applyOmegaPFA,
                                                        callback = progressbarTick)

    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedOmegaObj"]] <- createJaspState(out,
                                                               dependencies = c("itemDeletedOmega",
                                                                                "omegaEstimationMethod",
                                                                                "coefficientType",
                                                                                "intervalMethod"))
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

    # should the interval be bootstrapped
    if (options[["intervalMethod"]] == "bootstrapped") {
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["bootstrapSamples"]])
        out[["samp"]] <- apply(model[[if (options[["coefficientType"]] == "unstandardized") "bootSamp" else "bootCor"]],
                               1, Bayesrel:::applyalpha, callback = progressbarTick)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleAlphaObj"]] <- createJaspState(out, dependencies = c("scaleAlpha",
                                                                               "coefficientType"))
  }
  return(out)
}

.frequentistAlphaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedAlphaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedAlphaObj"]]$object)

  out <- model[["itemDeletedAlpha"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedAlpha"]] && is.null(model[["empty"]]) && options[["intervalMethod"]] == "bootstrapped") {


    startProgressbar(options[["bootstrapSamples"]] * ncol(dataset))
    type <- ifelse(options[["coefficientType"]] == "unstandardized", "bootSamp", "bootCor")
    jaspBase::.setSeedJASP(options)

    out[["itemSamp"]] <- .frequentistItemDroppedStats(covSamp = model[[type]],
                                                      f1 = Bayesrel:::applyalpha,
                                                      callback = progressbarTick)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedAlphaObj"]] <- createJaspState(out,
                                                               dependencies = c("itemDeletedAlpha", "coefficientType",
                                                                                "intervalMethod"))
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

    type <- ifelse(options[["coefficientType"]] == "unstandardized", "bootSamp", "bootCor")

    if (options[["intervalMethod"]] == "bootstrapped") {
      if (is.null(out[["samp"]])) {
          startProgressbar(options[["bootstrapSamples"]])
          out[["samp"]] <- apply(model[[type]], 1, Bayesrel:::applylambda2, callback = progressbarTick)
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleLambda2Obj"]] <- createJaspState(out, dependencies = c("scaleLambda2", "coefficientType"))
  }
  return(out)
}

.frequentistLambda2Item <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedLambda2Obj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedLambda2Obj"]]$object)

  out <- model[["itemDeletedLambda2"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedLambda2"]] && is.null(model[["empty"]]) && options[["intervalMethod"]] == "bootstrapped") {


    startProgressbar(options[["bootstrapSamples"]] * ncol(dataset))
    type <- ifelse(options[["coefficientType"]] == "unstandardized", "bootSamp", "bootCor")
    jaspBase::.setSeedJASP(options)

    out[["itemSamp"]] <- .frequentistItemDroppedStats(covSamp = model[[type]],
                                                      f1 = Bayesrel:::applylambda2,
                                                      callback = progressbarTick)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedLambda2Obj"]] <- createJaspState(out,
                                                               dependencies = c("itemDeletedLambda2", "coefficientType",
                                                                                "intervalMethod"))
  }
  return(out)
}

.frequentistSplithalfScale <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleSplithalfObj"]]$object))
    return(model)

  out <- model[["scaleSplithalf"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["scaleSplithalf"]]  && is.null(model[["empty"]])) {

    if (options[["intervalMethod"]] == "bootstrapped") {

      if (is.null(out[["samp"]])) {
        startProgressbar(options[["bootstrapSamples"]])

        nit <- ncol(dataset)
        splits <- split(seq_len(nit), 1:2)
        if (options[["coefficientType"]] == "unstandardized") {
          for (i in seq_len(options[["bootstrapSamples"]])) {
            out[["samp"]][i] <- .splithalfCor(model[["bootSamp"]][i, , ], splits, progressbarTick)
          }
        } else { # either we have the boostrapped cor samples from the standardized coefficients or we have them through
          # the splithalf method
          out[["samp"]] <- apply(model[["bootCor"]], 1, .splithalfCor, splits = splits)
        }
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["scaleSplithalfObj"]] <- createJaspState(out, dependencies = c("scaleSplithalf", "coefficientType"))
  }
  return(out)
}

.frequentistSplithalfItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemDeletedSplithalfObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDeletedSplithalfObj"]]$object)

  out <- model[["itemDeletedSplithalf"]]
  if (is.null(out))
    out <- list()

  if (options[["itemDeletedSplithalf"]] && is.null(model[["empty"]]) && options[["intervalMethod"]] == "bootstrapped") {

    type <- ifelse(options[["coefficientType"]] == "unstandardized", "bootSamp", "bootCor")

    startProgressbar(options[["bootstrapSamples"]] * ncol(dataset))
    jaspBase::.setSeedJASP(options)

    nit <- ncol(dataset) - 1
    splits <- split(seq_len(nit), 1:2)

    out[["itemSamp"]] <- .frequentistItemDroppedStats(covSamp = model[[type]],
                                                      f1 = .splithalfCor,
                                                      callback = progressbarTick,
                                                      splits = splits)

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDeletedSplithalfObj"]] <- createJaspState(out,
                                                                 dependencies = c("itemDeletedSplithalf",
                                                                                  "intervalMethod",
                                                                                  "coefficientType"))
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

    if (options[["intervalMethod"]] == "bootstrapped") {
      if (is.null(out[["samp"]])) {

        startProgressbar(options[["bootstrapSamples"]])
        out[["samp"]] <- numeric(options[["bootstrapSamples"]])

        if (options[["coefficientType"]] == "unstandardized" && is.null(model[["bootCor"]])) {
          model[["bootCor"]] <- model[["bootSamp"]]
          for (i in seq_len(options[["bootstrapSamples"]])) {
            corm <- .cov2cor.callback(model[["bootSamp"]][i, , ], progressbarTick)
            out[["samp"]][i] <- mean(corm[lower.tri(corm)])
          }
        } else { # either we have the boostrapped cor samples from the standardized coefficients or we have them through
          # the splithalf method
          out[["samp"]] <- apply(model[["bootCor"]], 1, function(x) mean(x[lower.tri(x)]))
        }

      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["avgObj"]] <- createJaspState(out, dependencies = "averageInterItemCorrelation")
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

    if (options[["intervalMethodVar"]] == "chisq") {
      out[["se"]] <- out[["est"]] * sqrt(2 / (model[["n"]] - 1))
      chiValueLow <- stats::qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- stats::qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
      out[["conf"]] <- c(((model[["n"]] - 1) * out[["est"]]) / chiValueLow,
                         ((model[["n"]] - 1) * out[["est"]]) / chiValueHigh)
    } else { # wald interval
      out[["se"]] <- .seVar(x = xx)
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
      chiValueLow <- stats::qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- stats::qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
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

  if (options[["itemRestCorrelation"]] && is.null(model[["empty"]]) && options[["intervalMethod"]] == "bootstrapped") {

    startProgressbar(options[["bootstrapSamples"]] * ncol(dataset))
    jaspBase::.setSeedJASP(options)

    out[["itemSamp"]] <- matrix(0, nrow = options[["bootstrapSamples"]], ncol = ncol(dataset))
    n <- model[["n"]]
    if (options[["bootstrapType"]] == "parametric") {
      cc <- cov(dataset, use = model[["use.cases"]])
      for (j in 1:ncol(dataset)) {
        for (b in seq_len(options[["bootstrapSamples"]])) {
          boot_data <- MASS::mvrnorm(n, colMeans(dataset, na.rm = TRUE), cc)
          out[["itemSamp"]][b, j] <- cor(as.matrix(boot_data[, j]), rowSums(as.matrix(boot_data[, -j]), na.rm = TRUE),
                                         use = model[["use.cases"]])
          progressbarTick()
        }
      }
    } else {
      for (j in 1:ncol(dataset)) {
        for (b in seq_len(options[["bootstrapSamples"]])) {
          boot_data <- as.matrix(dataset[sample.int(n, size = n, replace = TRUE), ])
          out[["itemSamp"]][b, j] <- cor(as.matrix(boot_data[, j]), rowSums(as.matrix(boot_data[, -j]), na.rm = TRUE),
                                         use = model[["use.cases"]])
          progressbarTick()
        }
      }
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = c("itemRestCorrelation", "intervalMethod",
                                                                             "bootstrapType"))
  }
  return(out)
}

.frequentistMeanItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemMeanObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemMeanObj"]]$object)

  out <- model[["itemMean"]]
  if (is.null(out))
    out <- list()

  if (options[["itemMean"]]  && is.null(model[["empty"]])) {
    ciValue <- options[["itemCiLevel"]]

    out[["est"]] <- colMeans(dataset, na.rm = TRUE)
    sds <- apply(dataset, 2, sd, na.rm = TRUE)

    zz <- qnorm(1 - (1 - ciValue) / 2)
    ses <- sds / sqrt(model[["n"]])
    out[["lower"]] <- out[["est"]] - zz * ses
    out[["upper"]] <- out[["est"]] + zz * ses

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemMeanObj"]] <- createJaspState(out, dependencies = c("itemMean", "itemCiLevel"))
  }
  return(out)
}

.frequentistVarItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemVarObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemVarObj"]]$object)

  out <- model[["itemVar"]]
  if (is.null(out))
    out <- list()

  if (options[["itemVar"]]  && is.null(model[["empty"]])) {

    ciValue <- options[["itemCiLevel"]]

    out[["est"]] <- apply(dataset, 2, var, na.rm = TRUE)

    if (options[["intervalMethodVar"]] == "chisq") {
      chiValueLow <- stats::qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- stats::qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
      out[["lower"]] <- ((model[["n"]] - 1) * out[["est"]]) / chiValueLow
      out[["upper"]] <- ((model[["n"]] - 1) * out[["est"]]) / chiValueHigh
    } else { # wald interval
      ses <- apply(dataset, 2, .seVar)
      out[["lower"]] <- out[["est"]] - qnorm(1 - (1 - ciValue) / 2) * ses
      out[["upper"]] <- out[["est"]] + qnorm(1 - (1 - ciValue) / 2) * ses
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemVarObj"]] <- createJaspState(out, dependencies = c("itemVar", "intervalMethodVar", "itemCiLevel"))
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

    ciValue <- options[["itemCiLevel"]]
    out[["est"]] <- apply(dataset, 2, sd, na.rm = TRUE)
    ses <- apply(dataset, 2, .seSd)

    if (options[["intervalMethodVar"]] == "chisq") {
      chiValueLow <- stats::qchisq(1 - (1 - ciValue) / 2, df = model[["n"]] - 1)
      chiValueHigh <- stats::qchisq((1 - ciValue) / 2, df = model[["n"]] - 1)
      out[["lower"]] <- sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueLow)
      out[["upper"]] <- sqrt(((model[["n"]] - 1) * out[["est"]]^2) / chiValueHigh)
    } else {
      out[["lower"]] <- out[["est"]] - qnorm(1 - (1 - ciValue) / 2) * ses
      out[["upper"]] <- out[["est"]] + qnorm(1 - (1 - ciValue) / 2) * ses
    }

    if (options[["samplesSavingDisabled"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemSdObj"]] <- createJaspState(out, dependencies = c("itemSd", "intervalMethodVar", "itemCiLevel"))
  }
  return(out)
}




.frequentistComputeScaleResults <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleResultsObj"]]$object)) {
    return(.getStateContainerF(jaspResults)[["scaleResultsObj"]]$object)

  }

  out <- model[["scaleResults"]]
  if (is.null(out))
    out <- list()

  if (is.null(model[["empty"]])) {

    ciValue <- options[["ciLevel"]]

    if (options[["coefficientType"]] == "unstandardized") {
      cc <- model[["data_cov"]]
      dtUse <- dataset
    } else {
      cc <- model[["data_cor"]]
      dtUse <- scale(dataset, scale = TRUE)
    }

    # go one coefficient at a time, because there are too many special options for a generic solution
    # omega
    if (options[["scaleOmega"]]) {
      if (options[["omegaEstimationMethod"]] == "cfa") {
        datasetOmega <- scale(dataset, scale = FALSE)
        omegaO <- .omegaFreqData(datasetOmega,
                                interval = ciValue,
                                omega.int.analytic = TRUE,
                                pairwise = model[["pairwise"]],
                                standardized = options[["coefficientType"]] == "standardized")
        if (is.na(omegaO[["omega"]])) {
          out[["error"]][["scaleOmega"]] <- gettext("Omega calculation with CFA failed.
                                                    Try changing to PFA in 'Advanced Options'. ")
          out[["est"]][["scaleOmega"]] <- NA
          out[["conf"]][["scaleOmega"]] <- c(NA, NA)
        } else {
          out[["fit"]][["scaleOmega"]] <- omegaO[["indices"]]
          out[["est"]][["scaleOmega"]] <- omegaO[["omega"]]

          if (options[["intervalMethod"]] == "analytic") {
            pp <- lavaan::parameterEstimates(omegaO$fit.object)
            out[["se"]][["scaleOmega"]] <- pp[pp$label == "omega", "se"]
            out[["conf"]][["scaleOmega"]] <- c(omegaO$omega_lower, omegaO$omega_upper)
          } else { # omega interval bootstrapped
            if (sum(!is.na(model[["scaleOmega"]][["samp"]])) >= 2) {
              out[["conf"]][["scaleOmega"]] <- quantile(model[["scaleOmega"]][["samp"]],
                                                        probs = c((1 - ciValue)/2, 1 - (1 - ciValue) / 2),
                                                        na.rm = TRUE)
              out[["se"]][["scaleOmega"]] <- sd(model[["scaleOmega"]][["samp"]], na.rm = TRUE)
            } else {
              out[["error"]][["scaleOmega"]] <- gettext("Omega bootstrapped interval calculation with CFA failed.
                                                        Try changing to PFA in 'Advanced Options'.  ")
              out[["conf"]][["scaleOmega"]] <- c(NA, NA)
              out[["se"]][["scaleOmega"]] <- NA
            }
          }

          if (options[["standardizedLoadings"]]) {

            lavObj <- omegaO$fit.object
            loads <- lavaan::inspect(lavObj, what = "std")$lambda
            out[["loadings"]] <- loads
          }
        }
      } else { # omega method is pfa
        omOut <- .applyOmegaPFA(cc, loadings = TRUE)
        out[["est"]][["scaleOmega"]] <- omOut[["om"]]
        if (is.na(omOut[["om"]])) {
          out[["error"]][["scaleOmega"]] <- gettext("Omega calculation with PFA failed. ")
          out[["est"]][["scaleOmega"]] <- NA
        } else {
          if (options[["intervalMethod"]] == "analytic") {
            out[["conf"]][["scaleOmega"]] <- c(NA, NA)
            out[["se"]][["scaleOmega"]] <- NA
            out[["error"]][["scaleOmega"]] <- gettext("The analytic confidence interval is not available for coefficient omega obtained with PFA. Change to 'Bootstrapped' in Advanced Options.")
          } else {
            if (sum(!is.na(model[["scaleOmega"]][["samp"]])) >= 2) {
              out[["conf"]][["scaleOmega"]] <- quantile(model[["scaleOmega"]][["samp"]],
                                                        probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                        na.rm = TRUE)
              out[["se"]][["scaleOmega"]] <- sd(model[["scaleOmega"]][["samp"]], na.rm = TRUE)
            } else {
              out[["error"]][["scaleOmega"]] <- gettext("Omega interval calculation with PFA failed. ")
              out[["conf"]][["scaleOmega"]] <- c(NA, NA)
              out[["se"]][["scaleOmega"]] <- NA
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

      out[["est"]][["scaleAlpha"]] <- Bayesrel:::applyalpha(cc)

      if (options[["intervalMethod"]] == "bootstrapped") {
        samp <- model[["scaleAlpha"]][["samp"]]
        out[["conf"]][["scaleAlpha"]] <- quantile(samp,
                                                    probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2),
                                                    na.rm = TRUE)
        out[["se"]][["scaleAlpha"]] <- sd(samp, na.rm = TRUE)

      } else { # alpha interval analytic

        if (model[["pairwise"]]) {
          out[["se"]][["scaleAlpha"]] <- NA
          out[["error"]][["scaleAlpha"]] <- gettext("The analytic confidence interval is not available for coefficient alpha/lambda2 when data contain missings and pairwise complete observations are used. Try changing to 'Delete listwise' within 'Advanced Options'.")
        } else {
          out[["se"]][["scaleAlpha"]] <- .seLambda3(dtUse, scaleThreshold = options[["hiddenScaleThreshold"]])
        }
        out[["conf"]][["scaleAlpha"]] <- out[["est"]][["scaleAlpha"]] + c(-1, 1) * out[["se"]][["scaleAlpha"]] * qnorm(1 - (1 - ciValue) / 2)

      }
    }


    # lambda 2
    if (options[["scaleLambda2"]]) {
      out[["est"]][["scaleLambda2"]] <- Bayesrel:::applylambda2(cc)
      if (options[["intervalMethod"]] == "bootstrapped") {
        samp <- model[["scaleLambda2"]][["samp"]]
        out[["conf"]][["scaleLambda2"]] <- quantile(samp, probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
        out[["se"]][["scaleLambda2"]] <- sd(samp, na.rm = TRUE)
      } else { # interval analytic
        if (model[["pairwise"]]) {
          out[["se"]][["scaleLambda2"]] <- NA
          if (is.null(out[["error"]][["scaleAlpha"]]))
            out[["error"]][["scaleLambda2"]] <- gettext("The analytic confidence interval is not available for coefficient alpha/lambda2 when data contain missings and pairwise complete observations are used. Try changing to 'Delete listwise' within 'Advanced Options'.")
        } else {
          out[["se"]][["scaleLambda2"]] <- .seLambda2(dtUse, scaleThreshold = options[["hiddenScaleThreshold"]])
        }
        out[["conf"]][["scaleLambda2"]] <- out[["est"]][["scaleLambda2"]] + c(-1, 1) * out[["se"]][["scaleLambda2"]] * qnorm(1 - (1 - ciValue) / 2)

      }
    }

    # split-half
    if (options[["scaleSplithalf"]]) {
      nit <- ncol(dataset)
      splits <- split(seq_len(nit), 1:2)
      out[["est"]][["scaleSplithalf"]] <- .splithalfData(dtUse, splits = splits, useCase = model[["use.cases"]])
      if (options[["intervalMethod"]] == "bootstrapped") {
        samp <- model[["scaleSplithalf"]][["samp"]]
        out[["conf"]][["scaleSplithalf"]] <- quantile(samp, probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
        out[["se"]][["scaleSplithalf"]] <- sd(samp, na.rm = TRUE)
      } else { # interval analytic
        partSums1 <- rowSums(dtUse[, splits[[1]]])
        partSums2 <- rowSums(dtUse[, splits[[2]]])

        out[["se"]][["scaleSplithalf"]] <- .seSplithalf(partSums1, partSums2, model[["use.cases"]])
        out[["conf"]][["scaleSplithalf"]] <- out[["est"]][["scaleSplithalf"]] + c(-1, 1) * out[["se"]][["scaleSplithalf"]] * qnorm(1 - (1 - ciValue) / 2)
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
        out[["error"]][["averageInterItemCorrelation"]] <- gettext("The standard error of the average interitem correlation is not available. ")
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
                                                        "scaleLambda2", "scaleSplithalf",
                                                        "averageInterItemCorrelation",
                                                        "meanSdScoresMethod",
                                                        "omegaEstimationMethod", "intervalMethod",
                                                        "intervalMethodVar", "coefficientType"))


  }

  return(out)

}


.frequentistComputeItemResults <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["itemResultsObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemResultsObj"]]$object)

  out <- model[["itemResults"]]
  if (is.null(out))
    out <- list()

  if (is.null(model[["empty"]])) {

    ciValue <- options[["itemCiLevel"]]

    if (options[["coefficientType"]] == "unstandardized") {
      cc <- model[["data_cov"]]
      dtUse <- dataset
    } else {
      cc <- model[["data_cor"]]
      dtUse <- scale(dataset, scale = TRUE)
    }

    # go one coefficient at a time, because there are too many special options for a generic solution
    # omega
    if (options[["itemDeletedOmega"]]) {

      if (ncol(dataset) == 2) {
        out[["est"]][["itemDeletedOmega"]] <- c(NA, NA)
        out[["lower"]][["itemDeletedOmega"]] <- c(NA, NA)
        out[["upper"]][["itemDeletedOmega"]] <- c(NA, NA)
      } else {
        if (options[["omegaEstimationMethod"]] == "cfa") {
          datasetOmega <- scale(dtUse, scale = FALSE)
          # do we have to compute item dropped values
          for (i in seq_len(ncol(dataset))) {
            outTmp <- .omegaFreqData(data = datasetOmega[, -i], interval = ciValue,
                                     pairwise = model[["pairwise"]],
                                     omega.int.analytic = TRUE,
                                     standardized = options[["coefficientType"]] == "standardized")
            out[["est"]][["itemDeletedOmega"]][i] <- outTmp$omega
            if (options[["intervalMethod"]] == "analytic") {
              out[["lower"]][["itemDeletedOmega"]][i] <- outTmp$omega_lower
              out[["upper"]][["itemDeletedOmega"]][i] <- outTmp$omega_upper
            } else {
              itemSamp <- model[["itemDeletedOmega"]][["itemSamp"]]
              conf <- quantile(itemSamp[, i], probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
              out[["lower"]][["itemDeletedOmega"]][i] <- conf[1]
              out[["upper"]][["itemDeletedOmega"]][i] <- conf[2]
            }
          }

          if (anyNA(c(out[["est"]][["itemDeletedOmega"]],
                      out[["lower"]][["itemDeletedOmega"]],
                      out[["upper"]][["itemDeletedOmega"]])))
            out[["error"]][["itemDeletedOmega"]] <- gettext("Omega item dropped statistics with CFA failed.")

        } else { # pfa

          for (i in seq_len(ncol(dataset))) {

            out[["est"]][["itemDeletedOmega"]][i] <- .applyOmegaPFA(cc[-i, -i], loadings = FALSE)

            itemSamp <- model[["itemDeletedOmega"]][["itemSamp"]]
            conf <- quantile(itemSamp[, i], probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
            out[["lower"]][["itemDeletedOmega"]][i] <- conf[1]
            out[["upper"]][["itemDeletedOmega"]][i] <- conf[2]

          }
          if (anyNA(c(out[["est"]][["itemDeletedOmega"]],
                      out[["lower"]][["itemDeletedOmega"]],
                      out[["upper"]][["itemDeletedOmega"]])))
            out[["error"]][["itemDeletedOmega"]] <- gettext("Omega item dropped statistics with PFA failed.")
        }
      }
    }

    # alpha
    if (options[["itemDeletedAlpha"]]) {

      if (ncol(dataset) == 2) {
        out[["est"]][["itemDeletedAlpha"]] <- c(NA, NA)
        out[["lower"]][["itemDeletedAlpha"]] <- c(NA, NA)
        out[["upper"]][["itemDeletedAlpha"]] <- c(NA, NA)
      } else {
        for (i in seq_len(ncol(dataset))) {

          est <- Bayesrel:::applyalpha(cc[-i, -i])
          out[["est"]][["itemDeletedAlpha"]][i] <- est

          if (options[["intervalMethod"]] == "analytic") {
            if (model[["pairwise"]]) {
              out[["lower"]][["itemDeletedAlpha"]][i] <- NA
              out[["upper"]][["itemDeletedAlpha"]][i] <- NA
              out[["error"]][["itemDeletedAlpha"]] <- gettext("The analytic confidence interval not available for coefficient alpha/lambda2 when data contain missings and pairwise complete observations are used. Try changing to 'Delete listwise' within 'Advanced Options'.")
            } else {
              se <- .seLambda3(dtUse[, -i, drop = FALSE], scaleThreshold = options[["hiddenScaleThreshold"]])
              conf <- est + c(-1, 1) * se * qnorm(1 - (1 - ciValue) / 2)
              out[["lower"]][["itemDeletedAlpha"]][i] <- conf[1]
              out[["upper"]][["itemDeletedAlpha"]][i] <- conf[2]
            }

          } else {
            itemSamp <- model[["itemDeletedAlpha"]][["itemSamp"]]
            conf <- quantile(itemSamp[, i], probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
            out[["lower"]][["itemDeletedAlpha"]][i] <- conf[1]
            out[["upper"]][["itemDeletedAlpha"]][i] <- conf[2]
          }
        }
      }

    }

    # lambda2
    if (options[["itemDeletedLambda2"]]) {

      if (ncol(dataset) == 2) {
        out[["est"]][["itemDeletedLambda2"]] <- c(NA, NA)
        out[["lower"]][["itemDeletedLambda2"]] <- c(NA, NA)
        out[["upper"]][["itemDeletedLambda2"]] <- c(NA, NA)
      } else {
        for (i in seq_len(ncol(dataset))) {

          est <- Bayesrel:::applylambda2(cc[-i, -i])
          out[["est"]][["itemDeletedLambda2"]][i] <- est

          if (options[["intervalMethod"]] == "analytic") {
            if (model[["pairwise"]]) {
              out[["lower"]][["itemDeletedLambda2"]][i] <- NA
              out[["upper"]][["itemDeletedLambda2"]][i] <- NA
              if (is.null(out[["error"]][["itemDeletedAlpha"]]))
                out[["error"]][["itemDeletedLambda2"]] <- gettext("The analytic confidence interval not available for coefficient alpha/lambda2 when data contain missings and pairwise complete observations are used. Try changing to 'Delete listwise' within 'Advanced Options'.")
            } else {
              se <- .seLambda2(dtUse[, -i, drop = FALSE], scaleThreshold = options[["hiddenScaleThreshold"]])
              conf <- est + c(-1, 1) * se * qnorm(1 - (1 - ciValue) / 2)
              out[["lower"]][["itemDeletedLambda2"]][i] <- conf[1]
              out[["upper"]][["itemDeletedLambda2"]][i] <- conf[2]
            }

          } else {
            itemSamp <- model[["itemDeletedLambda2"]][["itemSamp"]]
            conf <- quantile(itemSamp[, i], probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
            out[["lower"]][["itemDeletedLambda2"]][i] <- conf[1]
            out[["upper"]][["itemDeletedLambda2"]][i] <- conf[2]
          }
        }
      }
    }

    # split-half
    if (options[["itemDeletedSplithalf"]]) {

      if (ncol(dataset) == 2) {
        out[["est"]][["itemDeletedSplithalf"]] <- c(NA, NA)
        out[["lower"]][["itemDeletedSplithalf"]] <- c(NA, NA)
        out[["upper"]][["itemDeletedSplithalf"]] <- c(NA, NA)
      }

      for (i in seq_len(ncol(dtUse))) {
        dtCut <- dtUse[, -i, drop = FALSE]
        nit <- ncol(dtCut)
        splits <- split(seq_len(nit), 1:2)
        est <- .splithalfData(dtCut, splits = splits, useCase = model[["use.cases"]])
        out[["est"]][["itemDeletedSplithalf"]][i] <- est

        if (options[["intervalMethod"]] == "analytic") {

          partSums1 <- rowSums(dtCut[, splits[[1]]])
          partSums2 <- rowSums(dtCut[, splits[[2]]])

          se <- .seSplithalf(partSums1, partSums2, model[["use.cases"]])
          conf <- est + c(-1, 1) * se * qnorm(1 - (1 - ciValue) / 2)
          out[["lower"]][["itemDeletedSplithalf"]][i] <- conf[1]
          out[["upper"]][["itemDeletedSplithalf"]][i] <- conf[2]
        } else {
          itemSamp <- model[["itemDeletedSplithalf"]][["itemSamp"]]
          conf <- quantile(itemSamp[, i], probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
          out[["lower"]][["itemDeletedSplithalf"]][i] <- conf[1]
          out[["upper"]][["itemDeletedSplithalf"]][i] <- conf[2]
        }
      }
    }

    # item-rest correlation
    if (options[["itemRestCorrelation"]]) {

      for (i in seq_len(ncol(dataset))) {
        sum1 <- as.matrix(dataset[, i])
        sum2 <- rowSums(as.matrix(dataset[, -i]), na.rm = TRUE)
        cc <- cor(sum1, sum2, use = model[["use.cases"]])
        out[["est"]][["itemRestCorrelation"]][i] <- cc

        if (options[["intervalMethod"]] == "analytic") {
          conf <- .confCorZ(cr= cc, n = model[["n"]], alpha = ciValue)
          out[["lower"]][["itemRestCorrelation"]][i] <- conf[1]
          out[["upper"]][["itemRestCorrelation"]][i] <- conf[2]
        } else {
          itemSamp <- model[["itemRestCorrelation"]][["itemSamp"]]
          conf <- quantile(itemSamp[, i], probs = c((1 - ciValue) / 2, 1 - (1 - ciValue) / 2), na.rm = TRUE)
          out[["lower"]][["itemRestCorrelation"]][i] <- conf[1]
          out[["upper"]][["itemRestCorrelation"]][i] <- conf[2]
        }
      }
    }

    # just copying for mean and sd
    if (options[["itemMean"]]) {
      out[["est"]][["itemMean"]] <- model[["itemMean"]][["est"]]
      out[["lower"]][["itemMean"]] <- model[["itemMean"]][["lower"]]
      out[["upper"]][["itemMean"]] <- model[["itemMean"]][["upper"]]
    }

    if (options[["itemVar"]]) {
      out[["est"]][["itemVar"]] <- model[["itemVar"]][["est"]]
      out[["lower"]][["itemVar"]] <- model[["itemVar"]][["lower"]]
      out[["upper"]][["itemVar"]] <- model[["itemVar"]][["upper"]]
    }

    if (options[["itemSd"]]) {
      out[["est"]][["itemSd"]] <- model[["itemSd"]][["est"]]
      out[["lower"]][["itemSd"]] <- model[["itemSd"]][["lower"]]
      out[["upper"]][["itemSd"]] <- model[["itemSd"]][["upper"]]
    }


    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemResultsObj"]] <- createJaspState(out,
                                                          dependencies = c("itemDeletedOmega",
                                                                           "itemDeletedAlpha", "itemDeletedLambda2",
                                                                           "itemDeletedSplithalf", "itemMean",
                                                                           "itemRestCorrelation", "itemSd", "itemVar",
                                                                           "itemCiLevel", "coefficientType",
                                                                           "intervalMethod", "omegaEstimationMethod",
                                                                           "intervalMethodVar"))

  }

  return(out)

}



##### Output tables ####

.frequentistScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Frequentist Scale Reliability Statistics"))

  scaleTable$dependOn(options = c("scaleOmega", "scaleAlpha", "scaleLambda2", "scaleSplithalf",
                                  "averageInterItemCorrelation", "scaleMean", "scaleSd", "meanSdScoresMethod",
                                  "omegaEstimationMethod", "intervalMethod", "intervalMethodVar",
                                  "ciLevel"))
  scaleTable$addColumnInfo(name = "coefficient", title = gettext("Coefficient"), type = "string")
  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "number")
  scaleTable$addColumnInfo(name = "se", title = gettext("Std. Error"), type = "number")
  ci <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)
  scaleTable$addColumnInfo(name = "lower", title = "Lower", type = "number", overtitle = gettextf("%s%% CI", ci))
  scaleTable$addColumnInfo(name = "upper", title = "Upper", type = "number", overtitle = gettextf("%s%% CI", ci))


  scaleTable$position <- 1
  stateContainer <- .getStateContainerF(jaspResults)
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

  if (!.is.empty(model) && any(selected)) {
    dt$estimate <- unlist(model[["scaleResults"]][["est"]], use.names = FALSE)
    dt$se <- unlist(model[["scaleResults"]][["se"]])
    dt$lower <- sapply(model[["scaleResults"]][["conf"]], function(x) x[1])
    dt$upper <- sapply(model[["scaleResults"]][["conf"]], function(x) x[2])

    errors <- unlist(model[["scaleResults"]][["error"]])
    if (length(errors) > 0) {
      model[["footnote"]] <- paste(model[["footnote"]], paste0(errors, collapse = ""))
    }

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
  itemTable$dependOn(options = c("itemDeletedOmega", "itemDeletedAlpha", "itemDeletedLambda2", "itemDeletedSplithalf",
                                 "itemMean", "itemRestCorrelation", "itemSd", "itemVar",
                                 "scaleOmega", "scaleAlpha", "scaleLambda2", "scaleSplithalf",
                                 "itemCiLevel", "intervalMethodVar", "intervalMethod"))
  # adding the scale options fixes a bug, where the item table would remain displayed
  # after one had checked a scale coefficient box and the item coefficient box and then unchecked the scale coeff box

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemTable$position <- 2
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  overTitles <- format(derivedOptions[["namesEstimators"]][["tables_item"]], digits = 3, drop0trailing = TRUE)
  overTitles <- gettextf("%s (if item dropped)", overTitles)
  ci <- format(100 * options[["itemCiLevel"]], digits = 3, drop0trailing = TRUE)

  selected <- derivedOptions[["itemDroppedSelected"]]
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]]
  idxSelected <- which(selected)
  coefficientsDeleted <- derivedOptions[["namesEstimators"]][["coefficientsDeleted"]]

  footnote <- ""
  if (length(model[["itemsDropped"]]) > 0) {
    itemTable[["variable"]] <- model[["itemsDropped"]]

    if (!is.null(unlist(options[["reverseScaledItems"]])))
      footnote <- .addFootnoteReverseScaledItems(options)
  }

  fewItemProblem <- FALSE
  for (i in idxSelected) {
    if (estimators[i] %in% coefficientsDeleted) {
      itemTable$addColumnInfo(name = paste0("pointEstimate", i), title = gettext("Estimate"), type = "number",
                              overtitle = overTitles[i])
      itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", ci), type = "number",
                              overtitle = overTitles[i])
      itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", ci), type = "number",
                              overtitle = overTitles[i])
      fewItemProblem <- TRUE
    } else {
      itemTable$addColumnInfo(name = paste0("pointEstimate", i), title = gettext("Estimate"), type = "number",
                              overtitle = estimators[i])
      itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", ci), type = "number",
                              overtitle = estimators[i])
      itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", ci), type = "number",
                              overtitle = estimators[i])
    }
  }

  if (is.null(model[["empty"]])) {
    tb <- data.frame(variable = model[["itemsDropped"]])

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])
      newtb <- cbind(pointEstimate = model[["itemResults"]][["est"]][[nm]],
                     lower = model[["itemResults"]][["lower"]][[nm]],
                     upper = model[["itemResults"]][["upper"]][[nm]])
      colnames(newtb) <- paste0(colnames(newtb), i)
      tb <- cbind(tb, newtb)
    }

    itemTable$setData(tb)

    errors <- unlist(model[["itemResults"]][["error"]])
    if (length(errors) > 0) {
      footnote <- paste(model[["footnote"]], paste0(errors, collapse = ""))
    }

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


.handleData <- function(dataset, options) {

  if (options[["naAction"]] == "listwise")
    dataset <- dataset[complete.cases(dataset), ]
  if (length(options$variables) > 0)
    dataset <- dataset[, options[["variables"]]]

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
  if (ncol(dataset) >= 2) {
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



.frequentistItemDroppedStats <- function(covSamp,
                                         f1 = function(){},
                                         callback = function(){},
                                         missing = NULL,
                                         splits = NULL) {

  dd <- dim(covSamp)
  out <- matrix(0, dd[1], dd[3])
  if (!is.null(splits)) { # split half
    for (i in seq_len(dd[3])) {
      out[, i] <- apply(covSamp[, -i, -i], c(1), f1, callback = callback, splits = splits)
    }
  } else {
    if (!is.null(missing)) { # cfa
      for (i in seq_len(dd[3])) {
        out[, i] <- apply(covSamp[, -i, -i], c(1), f1, callback = callback, missing = missing)
      }
    } else {
      for (i in seq_len(dd[3])) {
        out[, i] <- apply(covSamp[, -i, -i], c(1), f1, callback = callback)
      }
    }

  }

  return(out)
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
  if (!options[["scaleSplithalf"]])
    opts[["itemDeletedSplithalf"]] <- FALSE

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
.applyOmegaPFA <- function(m, callback = function(){}, loadings = FALSE){

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

.omegaFreqData <- function(
    data,
    interval,
    omega.int.analytic,
    pairwise,
    n.boot = 1e3,
    parametric = FALSE,
    callback = function(){},
    standardized = FALSE) {

  p <- ncol(data)
  n <- nrow(data)
  file <- Bayesrel:::lavOneFile(data)
  colnames(data) <- file$names

  lam_names <- paste0("l", 1:p)
  err_names <- paste0("e", 1:p)
  model <- paste0("f1 =~ ")
  loadings <- paste(paste(lam_names, "*", file$names, sep = ""),
                    collapse = " + ")
  errors <- paste(paste(file$names, " ~~ ", err_names, "*",
                        file$names, sep = ""), collapse = "\n")
  sum_loads <- paste("loading :=", paste(lam_names, collapse = " + "),
                     "\n")
  sum_errs <- paste("error :=", paste(err_names, collapse = " + "),
                    "\n")
  omega <- "omega := (loading^2) / ((loading^2) + error) \n"
  mod <- paste(model, loadings, "\n", errors,
               "\n", sum_loads, sum_errs, omega)

  if (pairwise) {
    fit <-  Bayesrel:::fitmodelMis(mod, data)
  } else {
    fit <-  Bayesrel:::fitmodel(mod, data)
  }

  if (is.null(fit)) {
    return(list(omega = NA, omega_lower = NA, omega_upper = NA, fit.object = NULL))
  } else {
    if (standardized) {
      params <- lavaan::standardizedsolution(fit, level = interval)
      omega <- params$est.std[params$lhs == "omega"]
    } else {
      params <- lavaan::parameterestimates(fit, level = interval)
      omega <- params$est[params$lhs == "omega"]
    }

    if (omega.int.analytic) {
      om_low <- params$ci.lower[params$lhs == "omega"]
      om_up <- params$ci.upper[params$lhs == "omega"]
      om_obj <- NA
    } else { # omega cfa with bootstrapping:
      if (parametric) {
        if (pairwise) {
          cc <- cov(data, use = "pairwise.complete.obs")
        } else {
          cc <- cov(data)
        }
        om_obj <- numeric(n.boot)
        for (i in 1:n.boot){
          boot_data <- MASS::mvrnorm(n, colMeans(data, na.rm = TRUE), cc)
          fitTmp <-  Bayesrel:::fitmodel(mod, boot_data)
          callback()
          if (!is.null(fitTmp)) {
            if (standardized) {
              params <- lavaan::standardizedsolution(fitTmp, level = interval)
              om_obj[i] <- params$est.std[params$lhs == "omega"]
            } else {
              params <- lavaan::parameterestimates(fitTmp, level = interval)
              om_obj[i] <- params$est[params$lhs == "omega"]
            }

          } else {
            om_obj[i] <- NA
          }
        }

        if (sum(!is.na(om_obj)) > 1) {
          om_low <- quantile(om_obj, prob = (1 - interval) / 2, na.rm = TRUE)
          om_up <- quantile(om_obj, prob = interval + (1 - interval) / 2, na.rm = TRUE)
        } else {
          om_low <- NA
          om_up <- NA
          om_obj <- NA
        }

      } else { # bootstrap non parametric

        om_obj <- numeric(n.boot)
        for (i in 1:n.boot){
          boot_data <- as.matrix(data[sample.int(n, size = n, replace = TRUE), ])
          if (pairwise) {
            fitTmp <- Bayesrel:::fitmodelMis(mod, boot_data)
          } else {
            fitTmp <-  Bayesrel:::fitmodel(mod, boot_data)
          }
          callback()
          if (!is.null(fitTmp)) {
            if (standardized) {
              params <- lavaan::standardizedsolution(fitTmp, level = interval)
              om_obj[i] <- params$est.std[params$lhs=="omega"]
            } else {
              params <- lavaan::parameterestimates(fitTmp, level = interval)
              om_obj[i] <- params$est[params$lhs=="omega"]
            }

          } else {
            om_obj[i] <- NA
          }
        }

        if (sum(!is.na(om_obj)) > 1) {
          om_low <- quantile(om_obj, prob = (1 - interval) / 2, na.rm = TRUE)
          om_up <- quantile(om_obj, prob = interval + (1 - interval) / 2, na.rm = TRUE)
        } else {
          om_low <- NA
          om_up <- NA
          om_obj <- NA
        }
      }
    }

    fit_tmp <- lavaan::fitMeasures(fit)
    indic <- c(fit_tmp["chisq"], fit_tmp["df"], fit_tmp["pvalue"],
               fit_tmp["rmsea"], fit_tmp["rmsea.ci.lower"], fit_tmp["rmsea.ci.upper"],
               fit_tmp["srmr"])
  }
  return(list(omega = omega, omega_lower = om_low, omega_upper = om_up, indices = indic, fit.object = fit,
              omega_boot = om_obj))
}


.applyOmegaCov <- function(
    cc,
    missing,
    callback = function(){})
  {

  p <- ncol(cc)
  file <- Bayesrel:::lavOneFile(cc)
  colnames(cc) <- file$names

  lam_names <- paste0("l", 1:p)
  err_names <- paste0("e", 1:p)
  model <- paste0("f1 =~ ")
  loadings <- paste(paste(lam_names, "*", file$names, sep = ""),
                    collapse = " + ")
  errors <- paste(paste(file$names, " ~~ ", err_names, "*",
                        file$names, sep = ""), collapse = "\n")
  sum_loads <- paste("loading :=", paste(lam_names, collapse = " + "),
                     "\n")
  sum_errs <- paste("error :=", paste(err_names, collapse = " + "),
                    "\n")
  omega <- "omega := (loading^2) / ((loading^2) + error) \n"
  mod <- paste(model, loadings, "\n", errors,
               "\n", sum_loads, sum_errs, omega)

  fit <-  .fitModelCov(mod, sample.cov = cc, sample.nobs = 500, missing)

  callback()

  if (is.null(fit)) {
    return(NA)
  } else {
    params <- lavaan::parameterestimates(fit)
    omega <- params$est[params$lhs == "omega"]
  }
  return(omega)
}


.fitModelCov <- function (mod, sample.cov, sample.nobs = 500, missing = "listwise") {
  out <- tryCatch({
    lavaan::cfa(mod, sample.cov = sample.cov, sample.nobs = sample.nobs, std.lv = TRUE, missing = missing)
  }, error = function(cond) {
    return(NULL)
  }, warning = function(cond) {
    return(NULL)
  }, finally = {
  })
  return(out)
}

#### Functions from Andries ####
# SE of sample variance
.seVar <- function(x, eps = 1e-16) {
  x <- na.omit(x)
  N <- length(x)
  d <- ((x - mean(x))^2 - var(x)) / (N - 1)
  V <- sum(d^2) - sum(outer(d, d))/ N
  if (max(V) < eps) return(NA) else return(sqrt(V))
}

# SE of sample standard deviation
.seSd <- function(x, eps = 1e-16) {
  return(.seVar(x)/ (2 * sd(x, na.rm = TRUE)))
}

# SE of sample lambda-3 (alpha)
.seLambda3 <- function(X, VC = NULL, eps = 1e-16, scaleThreshold = 10){
  J <- ncol(X)
  return((J / (J - 1)) *.seLambda1(X, VC, scaleThreshold = scaleThreshold))
}


# SE of sample lambda-1
.seLambda1 <- function(X, VC = NULL, eps = 1e-16, scaleThreshold){
  D <- function(x) diag(as.numeric(x))
  J <- ncol(X)
  # if(is.null(VC)) VC <- .varCM(X)
  if (is.null(VC)) {
    levs <- sapply(as.data.frame(X), function(col) length(unique(col[!is.na(col)])))
    if (any(levs > scaleThreshold)) {
      # Continuous/mixed: fast, stable normal-theory VC
      VC <- .varVCwishart(stats::var(X), nrow(X))
    } else {
      # Ordinal: keep your existing categorical VC
      VC <- .varCM(X)
    }
  }

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
    n <- as.matrix(table(apply(X, 1, paste, collapse=",")))                               # MODIFIED on 30-aug-2024
    R <- matrix(as.numeric(unlist(strsplit(dimnames(n)[[1]], ","))), ncol = ncol(X), byrow = TRUE)
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

  X <- X - min(X, na.rm = TRUE) + 1
  res <- X2nR(X)
  R <- res$R
  n <- res$n
  J <- ncol(R)
  K <- length(n)
  marg.default <- list()
  for (j in 1 : J) marg.default[[j]] <- sort(unique(R[, j]))
  marg <- marg.default # remove line when not testing
  eps <- 1e-16         # remove line when not testing

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

# SE of sample lambda-2
.seLambda2 <- function(X, VC = NULL, eps = 1e-16, scaleThreshold = 10) {
  J <- ncol(X)
  if (is.null(VC)) {
    levs <- sapply(as.data.frame(X), function(col) length(unique(col[!is.na(col)])))
    if (any(levs > scaleThreshold)) {
      # Continuous/mixed: fast, stable normal-theory VC
      VC <- .varVCwishart(stats::var(X), nrow(X))
    } else {
      # Ordinal: keep your existing categorical VC
      VC <- .varCM(X)
    }
  }

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

# Sample lambda-2
.lambda2 <- function(X){
  J <- ncol(X)
  C <- var(X)
  return(((sum(C)-sum(diag(C))) + (sqrt((J/(J-1))*(sum(C^2)-sum(diag(C^2))))))/sum(C))
}


.splithalfData <- function(X, splits, useCase) {

  partSums1 <- rowSums(X[, splits[[1]], drop = FALSE])
  partSums2 <- rowSums(X[, splits[[2]], drop = FALSE])

  rsh_uncorrected <- cor(partSums1, partSums2, use = useCase)
  rsh <- (2 * rsh_uncorrected) / (1 + rsh_uncorrected)
  return(rsh)
}

.splithalfCor <- function(R, splits, callback = function(){}) {

  R_AA <- R[splits[[1]], splits[[1]]]
  R_BB <- R[splits[[2]], splits[[2]]]
  R_AB <- R[splits[[1]], splits[[2]]]

  Var_XA <- sum(R_AA)
  Var_XB <- sum(R_BB)
  Cov_XA_XB <- sum(R_AB)

  rsh_uncorrected <- Cov_XA_XB / sqrt(Var_XA * Var_XB)
  rsh <- (2 * rsh_uncorrected) / (1 + rsh_uncorrected)

  callback()
  return(rsh)
}

.seSplithalf <- function(x, y, useObs){
  k <- cor(x, y, use = useObs)
  seK <- .seCor(k, length(x))
  sh <- 2 * k / (1 + k)
  return((sh/k - sh/(1 + k)) * seK)
}

.seCor <- function(r, n) { # Bonett (2008)
  return((1 - r^2) / sqrt(n - 3))
}

# fisher transformed correlation interval
.confCorZ <- function(cr, n, alpha) {
  zcor <- 0.5 * log((1 + cr) / (1  - cr))
  seZcor <- 1 / sqrt(n - 3)
  qq <- qnorm(1 - (1 - alpha) / 2)
  zlow <- zcor - qq * seZcor
  zup <- zcor + qq * seZcor
  low <- (exp(2 * zlow) - 1) / (exp(2 * zlow) + 1)
  up <- (exp(2 * zup) - 1) / (exp(2 * zup) + 1)
  return(c(low, up))
}

.varVCwishart <- function(Sigma, n) {
  # Var(vec(S)) for sample covariance with denominator (n-1), under MVN
  J  <- ncol(Sigma); nu <- n - 1L
  VC <- matrix(0, J*J, J*J)
  idx <- function(a,b) (a-1L)*J + b
  for (i in 1:J) for (j in 1:J) for (k in 1:J) for (l in 1:J) {
    VC[idx(i,j), idx(k,l)] <- (Sigma[i,k]*Sigma[j,l] + Sigma[i,l]*Sigma[j,k]) / nu
  }
  VC
}

# # function by Don
# .varVCwishart <- function(Sigma, nu) {
#   J <- ncol(Sigma)
#   # commutation matrix K
#   idx_row <- as.vector(outer(1:J, 1:J, function(i,j) (i-1L)*J + j))
#   idx_col <- as.vector(outer(1:J, 1:J, function(i,j) (j-1L)*J + i))
#   K <- matrix(0, J*J, J*J)
#   K[cbind(idx_row, idx_col)] <- 1
#   # main expression
#   AA <- kronecker(Sigma, Sigma)
#   ((diag(J*J) + K) %*% AA) / nu
# }
