
.frequentistOmega <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["omegaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["omegaObj"]]$object)

  out <- model[["omega"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaScale"]] && is.null(model[["empty"]])) {

    ciValue <- options[["confidenceIntervalValue"]]

    if (options[["omegaEst"]] == "cfa") {
      omegaO <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = T,
                                         pairwise = model[["pairwise"]])
      if (is.na(omegaO[["omega"]])) {
        # out[["omegaFit"]] <- NA
        .quitAnalysis("Omega calculation with CFA failed. \n Try changing to PFA in 'Advanced Options'") # a proper error needs to appear here
      } else {
        out[["omegaFit"]] <- omegaO[["indices"]]
        out[["est"]] <- omegaO[["omega"]]
      }

      if (options[["intervalOn"]]) {
        if (options[["omegaInterval"]] == "omegaAnalytic") {
          out[["conf"]] <- c(omegaO$omega_lower, omegaO$omega_upper)
        } else {
          omegaboot <- out[["omegaBoot"]]
          if (is.null(omegaboot)) {
            omegaboot <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = F,
                                                  pairwise = model[["pairwise"]], parametric = model[["parametric"]])
          }
          if (is.na(omegaboot[["omega_lower"]]) || is.na(omegaboot[["omega_upper"]])) {
            .quitAnalysis("Omega bootstrapped interval calculation with CFA failed. \n Try changing to PFA in 'Advanced Options'")
          } else {
            out[["conf"]] <- c(omegaboot[["omega_lower"]], omegaboot[["omega_upper"]])
          }
        }
      }
      # do we have to compute item dropped values
      if (options[["omegaItem"]]) {
        if (is.null(out[["itemDropped"]])) {
          out[["itemDropped"]] <- numeric(ncol(dataset))
          for (i in 1:ncol(dataset)) {
            out[["itemDropped"]][i] <- Bayesrel:::applyomega_cfa_data(dataset[, -i], interval = ciValue,
                                                                      pairwise = model[["pairwise"]])
          }
        }
        if (anyNA(out[["itemDropped"]]))
          .quitAnalysis("Omega item dropped statistics with CFA failed")
      }

    } else { # omega with pfa
      out[["est"]] <- Bayesrel:::applyomega_pfa(model[["data_cov"]])
      if (is.na(out[["est"]])) {
        .quitAnalysis("Omega calculation with PFA failed")
      }

      if (options[["intervalOn"]]) {
        if (is.null(out[["sampCov"]])) {
          startProgressbar(options[["noSamples"]])
          out[["sampCov"]] <- apply(model[["bootsamp"]], 1, Bayesrel:::applyomega_pfa, callback = progressbarTick)
        }
        if (anyNA(out[["sampCov"]]))
          .quitAnalysis("Omega interval calculation with PFA failed")
        out[["conf"]] <- quantile(out[["sampCov"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2))
      }

      # do we have to compute item dropped values
      if (options[["omegaItem"]]) {
        if (is.null(out[["itemDropped"]]))
          out[["itemDropped"]] <- apply(model[["itemDroppedCovs"]], 1, Bayesrel:::applyomega_pfa)
        if (anyNA(out[["itemDropped"]]))
          .quitAnalysis("Omega item dropped statistics with PFA failed")
      }
    }

    stateContainerF <- .getStateContainerF(jaspResults)
    stateContainerF[["omegaObj"]] <- createJaspState(out, dependencies = c("omegaScale", "omegaItem", "omegaEst",
                                                                           "omegaInterval"))
  }

  return(out)
}


.frequentistAlpha <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["alphaObj"]]$object))
    return(.getStateContainerF(jaspResults)[["alphaObj"]]$object)

  out <- model[["alpha"]]
  if (is.null(out))
    out <- list()

  if (options[["alphaScale"]] && is.null(model[["empty"]])) {

    # alpha unstandardized
    if (options[["alphaMethod"]] == "alphaUnstand") {

      out[["est"]] <- Bayesrel:::applyalpha(model[["data_cov"]])

      # do we need an interval estimate?
      if (options[["intervalOn"]]) {
        ciValue <- options[["confidenceIntervalValue"]]

        # should the interval be analytic
        if (options[["alphaInterval"]] == "alphaAnalytic") {
          out[["conf"]] <- Bayesrel:::ciAlpha(1 - ciValue, nrow(dataset), model[["data_cov"]])

        } else {
          if (is.null(out[["sampCov"]])) {
            startProgressbar(options[["noSamples"]])
            out[["sampCov"]] <- apply(model[["bootsamp"]], 1, Bayesrel:::applyalpha, callback = progressbarTick)
          }
          out[["conf"]] <- quantile(out[["sampCov"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2))
        }
      }

      # do we have to compute item dropped values
      if (options[["alphaItem"]]) {
        if (is.null(out[["itemDropped"]]))
          out[["itemDropped"]] <- apply(model[["itemDroppedCovs"]], 1, Bayesrel:::applyalpha)
      }


    } else { # alpha standardized
      ccor <- model[["data_cor"]]
      out[["est"]] <- Bayesrel:::applyalpha(ccor)

      # do we need an interval estimate?
      if (options[["intervalOn"]]) {
        ciValue <- options[["confidenceIntervalValue"]]

        # should the interval be analytic
        if (options[["alphaInterval"]] == "alphaAnalytic") {
          out[["conf"]] <- Bayesrel:::ciAlpha(1 - ciValue, nrow(dataset), ccor)

        } else {

          if (is.null(out[["sampCor"]])) {
            out[["sampCor"]] <- numeric(options[["noSamples"]])
            for (i in 1:options[["noSamples"]]) {
              out[["sampCor"]][i] <- Bayesrel:::applyalpha(cov2cor(model[["bootsamp"]][i, ,]))
            }
          }
          out[["conf"]] <- quantile(out[["sampCor"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm=T)

        }
      }

      # do we have to compute item dropped values
      if (options[["alphaItem"]]) {
        if (is.null(out[["itemDropped"]])) {
          out[["itemDropped"]] <- numeric(ncol(dataset))
          for (i in 1:ncol(dataset)){
            out[["itemDropped"]][i] <- Bayesrel:::applyalpha(ccor[-i, -i])
          }
        }
      }
    }

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["alphaObj"]] <- createJaspState(out, dependencies = c("alphaScale", "alphaItem", "alphaMethod"))
  }
  return(out)
}


.frequentistLambda2 <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["lambda2Obj"]]$object))
      return(.getStateContainerF(jaspResults)[["lambda2Obj"]]$object)

  out <- model[["lambda2"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["lambda2Scale"]]  && is.null(model[["empty"]])) {

    out[["est"]] <- Bayesrel:::applylambda2(model[["data_cov"]])

    # do we need an interval estimate?
    if (options[["intervalOn"]]) {

      ciValue <- options[["confidenceIntervalValue"]]

      if (is.null(out[["sampCov"]])) {
        startProgressbar(options[["noSamples"]])
        out[["sampCov"]] <- apply(model[["bootsamp"]], 1, Bayesrel:::applylambda2, callback = progressbarTick)
      }
      out[["conf"]] <- quantile(out[["sampCov"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2))

    }

    # do we have to compute item dropped values
    if (options[["lambda2Item"]]) {
      if (is.null(out[["itemDropped"]]))
        out[["itemDropped"]] <- apply(model[["itemDroppedCovs"]], 1, Bayesrel:::applylambda2)
    }

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda2Obj"]] <- createJaspState(out, dependencies = c("lambda2Scale", "lambda2Item"))
  }
  return(out)
}


.frequentistLambda6 <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["lambda6Obj"]]$object))
    return(.getStateContainerF(jaspResults)[["lambda6Obj"]]$object)

  out <- model[["lambda6"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["lambda6Scale"]]  && is.null(model[["empty"]])) {

    out[["est"]] <- Bayesrel:::applylambda6(model[["data_cov"]])
    if (is.na(out[["est"]]))
      .quitAnalysis("Lambda6 calculation failed because the data covariance matrix is not invertible")

    # do we need an interval estimate?
    if (options[["intervalOn"]]) {

      ciValue <- options[["confidenceIntervalValue"]]

      if (is.null(out[["sampCov"]])) {
        startProgressbar(options[["noSamples"]])
        out[["sampCov"]] <- apply(model[["bootsamp"]], 1, Bayesrel:::applylambda6, callback = progressbarTick)
      }

      if (sum(!is.na(out[["sampCov"]])) < 3)
        .quitAnalysis("Lambda6 interval calculation failed because some bootstrapped covariance matrix are not invertible")

      out[["conf"]] <- quantile(out[["sampCov"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2))

    }

    # do we have to compute item dropped values
    if (options[["lambda6Item"]]) {
      if (is.null(out[["itemDropped"]]))
        out[["itemDropped"]] <- apply(model[["itemDroppedCovs"]], 1, Bayesrel:::applylambda6)
      if (anyNA(out[["itemDropped"]]))
        .quitAnalysis("Lambda6 item dropped statistics failed because some bootstrapped covariance matrix are not invertible")
    }

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda6Obj"]] <- createJaspState(out, dependencies = c("lambda6Scale", "lambda6Item"))
  }
  return(out)
}


# check the error handling of the glb !!!!!!!!
.frequentistGlb <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["glbObj"]]$object))
    return(.getStateContainerF(jaspResults)[["glbObj"]]$object)

  out <- model[["glb"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["glbScale"]]  && is.null(model[["empty"]])) {

    out[["est"]] <- Bayesrel:::glbOnArray(model[["data_cov"]])
    # if (is.na(out[["est"]]))
    #   .quitAnalysis("glb calculation failed because the data covariance matrix is not invertible")

    # do we need an interval estimate?
    if (options[["intervalOn"]]) {

      ciValue <- options[["confidenceIntervalValue"]]

      if (is.null(out[["sampCov"]])) {
        # startProgressbar(options[["noSamples"]])
        out[["sampCov"]] <- Bayesrel:::glbOnArray(model[["bootsamp"]])
      }

      # if (sum(!is.na(out[["sampCov"]])) < 3)
      #   .quitAnalysis("glb interval calculation failed because some bootstrapped covariance matrix are not invertible")
      out[["conf"]] <- quantile(out[["sampCov"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2))
    }

    # do we have to compute item dropped values
    if (options[["glbItem"]]) {
      if (is.null(out[["itemDropped"]]))
        out[["itemDropped"]] <- Bayesrel:::glbOnArray(model[["itemDroppedCovs"]])
      # if (anyNA(out[["itemDropped"]]))
      #   .quitAnalysis("glb item dropped statistics failed because some bootstrapped covariance matrix are not invertible")
    }

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["glbObj"]] <- createJaspState(out, dependencies = c("glbScale", "glbItem"))
  }
  return(out)
}


.frequentistAverageCor <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["avgObj"]]$object))
    return(.getStateContainerF(jaspResults)[["avgObj"]]$object)

  out <- model[["average"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["averageInterItemCor"]] && is.null(model[["empty"]])) {
    ciValue <- options[["confidenceIntervalValue"]]

    out[["est"]] <- mean(model[["data_cor"]][lower.tri(model[["data_cor"]])])

    if (options[["intervalOn"]]) {

      if (is.null(out[["sampCov"]])) {
        startProgressbar(options[["noSamples"]])
        out[["sampCov"]] <- numeric(options[["noSamples"]])
        for (i in 1:options[["noSamples"]]) {
          corm <- .cov2cor.callback(model[["bootsamp"]][i, , ], progressbarTick)
          out[["sampCov"]][i] <- mean(corm[lower.tri(corm)])
        }
      }
      out[["conf"]] <- quantile(out[["sampCov"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2))
    }
    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["avgObj"]] <- createJaspState(out, dependencies = c("averageInterItemCor"))
  }
  return(out)
}

.frequentistMean <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["meanObj"]]$object))
    return(.getStateContainerF(jaspResults)[["meanObj"]]$object)

  out <- model[["mean"]]
  if (is.null(out))
    out <- list()
  if (options[["meanScale"]] && is.null(model[["empty"]])) {
    if (options[["meanMethod"]] == "sumScores")
      out[["est"]] <- mean(rowSums(dataset, na.rm = T))
    else
      out[["est"]] <- mean(rowMeans(dataset, na.rm = T))

    if (options[["intervalOn"]])
      out[["conf"]] <- c(NA_real_, NA_real_)
    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("meanScale", "meanMethod"))
  }
  return(out)
}

.frequentistStdDev <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["sdObj"]]$object))
    return(.getStateContainerF(jaspResults)[["sdObj"]]$object)

  out <- model[["sd"]]
  if (is.null(out))
    out <- list()
  if (options[["sdScale"]] && is.null(model[["empty"]])) {
    if (options[["sdMethod"]] == "sumScores")
      out[["est"]] <- sd(rowSums(dataset, na.rm = T))
    else
      out[["est"]] <- sd(rowMeans(dataset, na.rm = T))

    if (options[["intervalOn"]])
      out[["conf"]] <- c(NA_real_, NA_real_)
    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["sdObj"]] <- createJaspState(out, dependencies = c("sdScale", "sdMethod"))
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
    for (i in 1:ncol(dataset)) {
      out[["itemDropped"]][i] <- cor(dataset[, i], rowMeans(dataset[, -i], na.rm = T), use = model[["use.cases"]])
    }

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = c("itemRestCor"))
  }
  return(out)
}

.frequentistItemMean <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemMeanObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemMeanObj"]]$object)

  out <- model[["itemMean"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemMean"]]  && is.null(model[["empty"]])) {
    out[["itemDropped"]] <- numeric(ncol(dataset))
    for (i in 1:ncol(dataset)) {
      out[["itemDropped"]][i] <- mean(dataset[, i], na.rm = T)
    }
    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemMeanObj"]] <- createJaspState(out, dependencies = c("itemMean"))
  }
  return(out)
}

.frequentistItemSd <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemSdObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemSdObj"]]$object)

  out <- model[["itemSd"]]
  if (is.null(out))
    out <- list()
  # is box even checked?
  if (options[["itemSd"]]  && is.null(model[["empty"]])) {
    out[["itemDropped"]] <- apply(dataset, 2, sd, na.rm = T)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemSdObj"]] <- createJaspState(out, dependencies = c("itemSd"))
  }
  return(out)
}



