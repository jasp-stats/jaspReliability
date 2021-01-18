


.BayesianOmega <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["omegaObj"]]$object))
    return(.getStateContainerB(jaspResults)[["omegaObj"]]$object)

  out <- model[["omega"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaScale"]] && is.null(model[["empty"]])) {

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]])
      dataset <- scale(dataset, scale = F)
      tmp_out <- Bayesrel:::omegaSampler(dataset, options[["noSamples"]], options[["noBurnin"]],
                                               options[["noThin"]], options[["noChains"]],
                                               model[["pairwise"]], progressbarTick)
      out[["samp"]] <- tmp_out$omega
      out[["loadings"]] <- apply(tmp_out$lambda, 3, mean)
      out[["residuals"]] <- apply(tmp_out$psi, 3, mean)
    }
    out[["est"]] <- mean(out[["samp"]])
    out[["cred"]] <- coda::HPDinterval(coda::mcmc(as.vector(out[["samp"]])), prob = ciValue)

    # do we have to compute item dropped values
    if (options[["omegaItem"]]) {
      ciValueItem <- options[["credibleIntervalValueItem"]]

      if (is.null(out[["itemSamp"]])) {
        startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

        dataset <- scale(dataset, scale = F)
        out[["itemSamp"]] <- array(0,
                                   c(options[["noChains"]],
                                     length(seq(1, options[["noSamples"]]-options[["noBurnin"]], options[["noThin"]])),
                                     ncol(dataset)))
        for (i in 1:ncol(dataset)) {
          out[["itemSamp"]][, , i] <- Bayesrel:::omegaSampler(dataset[-i, -i],
                                                   options[["noSamples"]], options[["noBurnin"]], options[["noThin"]],
                                                   options[["noChains"]], model[["pairwise"]], progressbarTick)$omega
        }
      }
      out[["itemEst"]] <- apply(out[["itemSamp"]], 3, mean)
      out[["itemCred"]] <- coda::HPDinterval(coda::mcmc(apply(out[["itemSamp"]], 3, as.vector)),
                                             prob = ciValueItem)
      }

    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["omegaObj"]] <- createJaspState(out, dependencies = c("omegaScale", "omegaItem"))
  }

  return(out)
}



.BayesianAlpha <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["alphaObj"]]$object))
    return(.getStateContainerB(jaspResults)[["alphaObj"]]$object)

  out <- model[["alpha"]]
  if (is.null(out))
    out <- list()

  if (options[["alphaScale"]] && !is.null(model[["gibbsSamp"]])) {

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]])
      out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::applyalpha, progressbarTick))
    }
    out[["est"]] <- mean(out[["samp"]])
    out[["cred"]] <- coda::HPDinterval(coda::mcmc(as.vector(out[["samp"]])), prob = ciValue)
    # do we have to compute item dropped values
    if (options[["alphaItem"]] && !is.null(model[["itemDroppedCovs"]])) {
      ciValueItem <- options[["credibleIntervalValueItem"]]

      if (is.null(out[["itemSamp"]])) {
        startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

        out[["itemSamp"]] <- apply(model[["itemDroppedCovs"]], c(1, 2, 3), Bayesrel:::applyalpha, progressbarTick)

      }
      out[["itemEst"]] <- apply(out[["itemSamp"]], 3, mean)
      out[["itemCred"]] <- coda::HPDinterval(coda::mcmc(apply(out[["itemSamp"]], 3, as.vector)),
                                             prob = ciValueItem)
    }

    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["alphaObj"]] <- createJaspState(out, dependencies = c("alphaItem"))
  }

  return(out)
}



.BayesianLambda2 <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["lambda2Obj"]]$object))
    return(.getStateContainerB(jaspResults)[["lambda2Obj"]]$object)

  out <- model[["lambda2"]]
  if (is.null(out))
    out <- list()

  if (options[["lambda2Scale"]] && !is.null(model[["gibbsSamp"]])) {

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]])
      out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::applylambda2, progressbarTick))
    }
    out[["est"]] <- mean(out[["samp"]])
    out[["cred"]] <- coda::HPDinterval(coda::mcmc(as.vector(out[["samp"]])), prob = ciValue)
    # do we have to compute item dropped values
    if (options[["lambda2Item"]] && !is.null(model[["itemDroppedCovs"]])) {
      ciValueItem <- options[["credibleIntervalValueItem"]]

      if (is.null(out[["itemSamp"]])) {
        startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

        out[["itemSamp"]] <- apply(model[["itemDroppedCovs"]], c(1, 2, 3), Bayesrel:::applylambda2, progressbarTick)

      }
      out[["itemEst"]] <- apply(out[["itemSamp"]], 3, mean)
      out[["itemCred"]] <- coda::HPDinterval(coda::mcmc(apply(out[["itemSamp"]], 3, as.vector)),
                                             prob = ciValueItem)
    }

    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["lambda2Obj"]] <- createJaspState(out, dependencies = c("lambda2Item"))
  }

  return(out)
}




.BayesianLambda6 <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["lambda6Obj"]]$object))
    return(.getStateContainerB(jaspResults)[["lambda6Obj"]]$object)

  out <- model[["lambda6"]]
  if (is.null(out))
    out <- list()

  if (options[["lambda6Scale"]] && !is.null(model[["gibbsSamp"]])) {

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]])
      out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::applylambda6, progressbarTick))
    }
    out[["est"]] <- mean(out[["samp"]])
    out[["cred"]] <- coda::HPDinterval(coda::mcmc(as.vector(out[["samp"]])), prob = ciValue)
    # do we have to compute item dropped values
    if (options[["lambda6Item"]] && !is.null(model[["itemDroppedCovs"]])) {
      ciValueItem <- options[["credibleIntervalValueItem"]]

      if (is.null(out[["itemSamp"]])) {
        startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

        out[["itemSamp"]] <- apply(model[["itemDroppedCovs"]], c(1, 2, 3), Bayesrel:::applylambda6, progressbarTick)

      }
      out[["itemEst"]] <- apply(out[["itemSamp"]], 3, mean)
      out[["itemCred"]] <- coda::HPDinterval(coda::mcmc(apply(out[["itemSamp"]], 3, as.vector)),
                                             prob = ciValueItem)
    }

    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["lambda6Obj"]] <- createJaspState(out, dependencies = c("lambda6Item"))
  }

  return(out)
}



.BayesianGlb <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["glbObj"]]$object))
    return(.getStateContainerB(jaspResults)[["glbObj"]]$object)

  out <- model[["glb"]]
  if (is.null(out))
    out <- list()

  if (options[["glbScale"]] && !is.null(model[["gibbsSamp"]])) {

    ciValue <- options[["credibleIntervalValueScale"]]

    if (is.null(out[["samp"]])) {
      # startProgressbar(options[["noSamples"]] * options[["noChains"]])
      out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::glbOnArray))
      # out[["samp"]] <- coda::mcmc(apply(model[["gibbsSamp"]], MARGIN = c(1, 2), Bayesrel:::glbOnArray_custom))
    }
    out[["est"]] <- mean(out[["samp"]])
    out[["cred"]] <- coda::HPDinterval(coda::mcmc(as.vector(out[["samp"]])), prob = ciValue)
    # do we have to compute item dropped values
    if (options[["glbItem"]] && !is.null(model[["itemDroppedCovs"]])) {
      ciValueItem <- options[["credibleIntervalValueItem"]]

      if (is.null(out[["itemSamp"]])) {
        # startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

        out[["itemSamp"]] <- apply(model[["itemDroppedCovs"]], c(1, 2, 3), Bayesrel:::glbOnArray)

        # out[["itemSamp"]] <- apply(model[["itemDroppedCovs"]], c(1, 2, 3), Bayesrel:::glbOnArray_custom)

      }
      out[["itemEst"]] <- apply(out[["itemSamp"]], 3, mean)
      out[["itemCred"]] <- coda::HPDinterval(coda::mcmc(apply(out[["itemSamp"]], 3, as.vector)),
                                             prob = ciValueItem)
    }

    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["glbObj"]] <- createJaspState(out, dependencies = c("glbItem"))
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
      startProgressbar(options[["noSamples"]] * options[["noChains"]])
      out[["samp"]] <- matrix(0, nrow(model[["gibbsSamp"]]), ncol(model[["gibbsSamp"]]))
      lowerTriangleIndex = which(lower.tri(model[["gibbsSamp"]][1, 1, , ]))
      for (i in 1:nrow(model[["gibbsSamp"]])) {
        for (j in 1:ncol(model[["gibbsSamp"]])) {

          corm <- .cov2cor.callback(model[["gibbsSamp"]][i, j, , ], progressbarTick)
          out[["samp"]][i, j] <- mean(corm[lowerTriangleIndex])
        }
      }
      out[["samp"]] <- coda::mcmc(out[["samp"]])

    }
    out[["est"]] <- mean(out[["samp"]])
    out[["cred"]] <- coda::HPDinterval(coda::mcmc(as.vector(out[["samp"]])), prob = ciValue)


    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["avgCorObj"]] <- createJaspState(out)
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
    if (options[["meanMethod"]] == "sumScores")
      out[["est"]] <- mean(rowSums(dataset, na.rm = T))
    else
      out[["est"]] <- mean(rowMeans(dataset, na.rm = T))

    out[["cred"]] <- c(NA_real_, NA_real_)
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("meanScale", "meanMethod"))
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
    if (options[["sdMethod"]] == "sumScores")
      out[["est"]] <- sd(rowSums(dataset, na.rm = T))
    else
      out[["est"]] <- sd(rowMeans(dataset, na.rm = T))

    out[["cred"]] <- c(NA_real_, NA_real_)
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["sdObj"]] <- createJaspState(out, dependencies = c("sdScale", "sdMethod"))
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

      out[["itemSamp"]] <- .itemRestCor(dataset, options[["noSamples"]], options[["noBurnin"]],
                              options[["noThin"]], options[["noChains"]], model[["pairwise"]],
                              callback = progressbarTick)
    }

    out[["itemEst"]] <- apply(out[["itemSamp"]], 3, mean)
    out[["itemCred"]] <- coda::HPDinterval(coda::mcmc(apply(out[["itemSamp"]], 3, as.vector)),
                                           prob = ciValueItem)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = c("itemRestCor"))
  }

  return(out)
}

.BayesianItemMean <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemMeanObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemMeanObj"]]$object)

  out <- model[["itemMean"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemMean"]] && is.null(model[["empty"]])) {
    out[["itemEst"]] <- numeric(ncol(dataset))
    for (i in 1:ncol(dataset)) {
      out[["itemEst"]][i] <- mean(dataset[, i], na.rm = T)
    }
    out[["itemCred"]] <- matrix(NA_real_, ncol(dataset), 2)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemMeanObj"]] <- createJaspState(out, dependencies = c("itemMean"))
  }
  return(out)
}

.BayesianItemSd <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemSdObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemSdObj"]]$object)

  out <- model[["itemSd"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemSd"]] && is.null(model[["empty"]])) {

    out[["itemEst"]] <- apply(dataset, 2, sd, na.rm = T)
    out[["itemCred"]] <- matrix(NA_real_, ncol(dataset), 2)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["itemSdObj"]] <- createJaspState(out, dependencies = c("itemSd"))
  }
  return(out)
}



