
.BayesianOmegaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["omegaScaleObj"]]$object))
    return(.getStateContainerB(jaspResults)[["omegaScaleObj"]]$object)

  out <- model[["omegaScale"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaScale"]] && is.null(model[["empty"]])) {

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]])

      dataset <- scale(dataset, scale = FALSE)

      jaspBase::.setSeedJASP(options)

      tmp_out <- try(Bayesrel:::omegaSampler(dataset, options[["noSamples"]], options[["noBurnin"]],
                                             options[["noThin"]], options[["noChains"]],
                                             model[["pairwise"]], progressbarTick), silent = TRUE)
      if (model[["pairwise"]] && inherits(tmp_out, "try-error")) {
        .quitAnalysis(gettext("Sampling the posterior factor model for omega failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
      }
      out[["samp"]] <- tmp_out$omega
      out[["loadings"]] <- apply(tmp_out$lambda, 3, as.vector)
      out[["residuals"]] <- apply(tmp_out$psi, 3, as.vector)
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["omegaScaleObj"]] <- createJaspState(out, dependencies = "omegaScale")
  }

  return(out)
}

.BayesianOmegaItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["omegaItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["omegaItemObj"]]$object)

  out <- model[["omegaItem"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaItem"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemEst"]] <- c(NaN, NaN)
      out[["itemCred"]] <- matrix(NaN, 2, 2)
      colnames(out[["itemCred"]]) <- c("lower", "upper")
      return(out)
    }

    ciValueItem <- options[["credibleIntervalValueItem"]]

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

      dataset <- scale(dataset, scale = FALSE)
      jaspBase::.setSeedJASP(options)

      out[["itemSamp"]] <- array(0, c(options[["noChains"]],
                                      length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]])),
                                      ncol(dataset)))

      for (i in seq_len(ncol(dataset))) {
        out[["itemSamp"]][, , i] <- Bayesrel:::omegaSampler(dataset[, -i],
                                                            options[["noSamples"]], options[["noBurnin"]], options[["noThin"]],
                                                            options[["noChains"]], model[["pairwise"]], progressbarTick)$omega
      }
      dd <- dim(out[["itemSamp"]])
      out[["itemSamp"]] <- matrix(out[["itemSamp"]], dd[1] * dd[2], ncol(dataset))

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["omegaItemObj"]] <- createJaspState(out, dependencies = "omegaItem")
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

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])
      out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::applyalpha, progressbarTick))
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["alphaScaleObj"]] <- createJaspState(out, dependencies = "alphaScale")
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
    ciValueItem <- options[["credibleIntervalValueItem"]]

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesItemDroppedStats(model[["gibbsSamp"]], Bayesrel:::applyalpha, progressbarTick)

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["alphaItemObj"]] <- createJaspState(out, dependencies = "alphaItem")
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

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])
      out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::applylambda2, progressbarTick))

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda2ScaleObj"]] <- createJaspState(out, dependencies = "lambda2Scale")
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
    ciValueItem <- options[["credibleIntervalValueItem"]]

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))
      out[["itemSamp"]] <- .BayesItemDroppedStats(model[["gibbsSamp"]], Bayesrel:::applylambda2, progressbarTick)

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda2ItemObj"]] <- createJaspState(out, dependencies = "lambda2Item")
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
    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(model[["progressbarLength"]])
      out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::applylambda6, progressbarTick))
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda6ScaleObj"]] <- createJaspState(out, dependencies = "lambda6Scale")
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
    ciValueItem <- options[["credibleIntervalValueItem"]]

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(model[["progressbarLength"]] * ncol(dataset))

      out[["itemSamp"]] <- .BayesItemDroppedStats(model[["gibbsSamp"]], Bayesrel:::applylambda6, progressbarTick)

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["lambda6ItemObj"]] <- createJaspState(out, dependencies = "lambda6Item")
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

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {

      dd <- dim(model[["gibbsSamp"]])
      out[["samp"]] <- matrix(0, dd[1], dd[2])

      startProgressbar(dd[1] * 3)
      for (i in seq_len(dd[1])) {
        out[["samp"]][i, ] <- Bayesrel:::glbOnArrayCustom(model[["gibbsSamp"]][i, , , ], callback = progressbarTick)
      }

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["glbScaleObj"]] <- createJaspState(out, dependencies = "glbScale")
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
    ciValueItem <- options[["credibleIntervalValueItem"]]

    if (is.null(out[["itemSamp"]])) {
      # special case glb, because it works with arrays not only matrices, small speedup...
      dd <- dim(model[["gibbsSamp"]])
      out[["itemSamp"]] <- matrix(0, dd[1] * dd[2], dd[3])
      cov_samp <- array(model[["gibbsSamp"]], c(dd[1] * dd[2], dd[3], dd[3]))

      startProgressbar(3 * ncol(dataset))
      for (i in seq_len(dd[3])) {
        out[["itemSamp"]][, i] <- Bayesrel:::glbOnArrayCustom(cov_samp[, -i, -i], callback = progressbarTick)
      }

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["glbItemObj"]] <- createJaspState(out, dependencies = "glbItem")
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

    ciValue <- options[["credibleIntervalValueScale"]]

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
    stateContainer[["avgCorObj"]] <- createJaspState(out, dependencies = "averageItemItemCor")
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
    ciValueItem <- options[["credibleIntervalValueItem"]]

    if (is.null(out[["itemSamp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

      jaspBase::.setSeedJASP(options)
      out[["itemSamp"]] <- .itemRestCor(dataset, options[["noSamples"]], options[["noBurnin"]],
                                        options[["noThin"]], options[["noChains"]], model[["pairwise"]],
                                        callback = progressbarTick)
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = "itemRestCor")
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


.itemRestCor <- function(dataset, n.iter, n.burnin, thin, n.chains, pairwise, callback) {

  ircor_samp <- matrix(0, n.chains * length(seq(1, n.iter - n.burnin, thin)), ncol(dataset))
  for (i in seq(ncol(dataset))) {
    help_dat <- cbind(as.matrix(dataset[, i]), rowMeans(as.matrix(dataset[, -i]), na.rm = TRUE))
    ircor_samp[, i] <- .WishartCorTransform(help_dat, n.iter = n.iter, n.burnin = n.burnin, thin = thin,
                                            n.chains = n.chains, pairwise = pairwise, callback = callback)
  }

  return(ircor_samp)
}

.WishartCorTransform <- function(x, n.iter, n.burnin, thin, n.chains, pairwise, callback) {
  tmp_cov <- Bayesrel:::covSamp(x, n.iter, n.burnin, thin, n.chains, pairwise, callback)$cov_mat
  dd <- dim(tmp_cov)
  tmp_cov <- array(tmp_cov, c(dd[1] * dd[2], dd[3], dd[4]))
  tmp_cor <- apply(tmp_cov, c(1), cov2cor)
  out <- tmp_cor[2, ]
  callback()
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
    samps <- model[selected]
    samps <- lapply(samps, function(x) x[["samp"]])
    samps[sapply(samps, is.null)] <- NULL

    out[["est"]] <- lapply(samps, mean)
    out[["cred"]] <- lapply(samps, function(x) coda::HPDinterval(coda::mcmc(c(x)), prob = ciValue))

    if (options[["rHat"]]) {
      tmp <- lapply(samps, function(x) {lapply(as.data.frame(t(x)), coda::mcmc)})
      out[["rHat"]] <- lapply(tmp, function(x) {coda::gelman.diag(coda::as.mcmc.list(x))[["psrf"]][, 1]})
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
                                                                                 "meanMethod", "sdMethod"))

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
    samps <- model[selected]
    samps <- lapply(samps, function(x) x[["itemSamp"]])
    samps[sapply(samps, is.null)] <- NULL

    out[["est"]] <- lapply(samps, colMeans)
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
                                                                                "itemRestCor", "meanItem", "sdItem"))

  }

  return(out)
}
