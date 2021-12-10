

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

    ### Chisqs ###
    LL1 <- sum(.dmultinorm(dataset, cdat)) # loglikelihood saturated model
    LR_obs <- apply(implieds, c(1, 2), .LRblav, data = dataset, basell = LL1) # loglikelihoods tested model

    lsm <- apply(model[["singleFactor"]][["loadings"]], 3, mean)
    psm <- apply(model[["singleFactor"]][["residuals"]], 3, mean)
    implM <- lsm %*% t(lsm) + diag(psm) # mean implied covariance matrix
    Dtm <- .LRblav(dataset, implM, LL1) # deviance of the mean model implied cov matrix
    out[["B-LR"]] <- Dtm

    ### BRMSEA ###
    Dm <- mean(LR_obs) # mean deviance of the model implied cov matrices
    pD <- Dm - Dtm # effective number of parameters (free parameters)
    out[["B-RMSEA"]] <- .BRMSEA(LR_obs, pstar, pD, n)
    # rmsea_m <- sqrt((Dtm - (pstar - pD)) / ((pstar - pD) * n))

    startProgressbar(options[["noSamples"]] * options[["noChains"]])

    ### CFI and TLI need a nullmodel:
    res_null <- try(Bayesrel:::omegaSamplerNull(dataset, options[["noSamples"]], options[["noBurnin"]],
                                                options[["noThin"]], options[["noChains"]],
                                                model[["pairwise"]], progressbarTick,
                                                a0 = options[["igShape"]], b0 = options[["igScale"]],
                                                m0 = options[["loadMean"]]), silent = TRUE)

    if (model[["pairwise"]] && inherits(tmp_out, "try-error")) {
      .quitAnalysis(gettext("Sampling the posterior null model failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
    }

    ps_null <- res_null$psi
    ps_null <- apply(ps_null, 3, as.vector)
    impl_null <- array(0, c(nrow(ps_null), k, k))
    for (i in 1:nrow(ps_null)) {
      impl_null[i, , ] <- diag(ps_null[i, ])
    }

    LR_null <- apply(impl_null, 1, .LRblav, data = dataset, basell = LL1)
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

    out[["B-CFI"]] <- cfis
    out[["B-TLI"]] <- tlis


    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["fitMeasuresObj"]] <- createJaspState(out, dependencies = c("omegaScale", "igShape", "igScale",
                                                                                "loadMean", "fitMeasures"))
  }

  return(out)
}
