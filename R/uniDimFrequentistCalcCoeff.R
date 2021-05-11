
.frequentistOmegaScale <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["omegaScaleObj"]]$object))
    return(.getStateContainerF(jaspResults)[["omegaScaleObj"]]$object)

  out <- model[["omegaScale"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaScale"]] && is.null(model[["empty"]])) {

    ciValue <- options[["confidenceIntervalValue"]]

    if (options[["omegaMethod"]] == "cfa") {

      dataset <- scale(dataset, scale = FALSE)
      omegaO <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = TRUE,
                                         pairwise = model[["pairwise"]])
      if (is.na(omegaO[["omega"]])) {
        out[["error"]] <- gettext("Omega calculation with CFA failed. Try changing to PFA in Advanced Options")
      } else {
        out[["omegaFit"]] <- omegaO[["indices"]]
        out[["est"]] <- omegaO[["omega"]]

        if (options[["intervalOn"]]) {

          if (options[["omegaInterval"]] == "omegaAnalytic") {
            out[["conf"]] <- c(omegaO$omega_lower, omegaO$omega_upper)
          } else {

            parametric <- options[["bootType"]] == "parametric"

            omegaboot <- out[["omegaBoot"]]
            if (is.null(omegaboot)) {
              jaspBase::.setSeedJASP(options)
              omegaboot <- Bayesrel:::omegaFreqData(dataset, interval = ciValue, omega.int.analytic = FALSE,
                                                    pairwise = model[["pairwise"]], parametric = parametric,
                                                    n.boot = options[["noSamples"]])
            }
            if (is.na(omegaboot[["omega_lower"]]) || is.na(omegaboot[["omega_upper"]])) {
              out[["error"]] <- gettext("Omega bootstrapped interval calculation with CFA failed.
                                          \n Try changing to PFA in 'Advanced Options'")
            } else {
              out[["conf"]] <- c(omegaboot[["omega_lower"]], omegaboot[["omega_upper"]])
            }
          }
        }
      }

    } else { # omega with pfa
      out[["est"]] <- Bayesrel:::applyomega_pfa(model[["data_cov"]])
      if (is.na(out[["est"]])) {
        out[["error"]] <- gettext("Omega calculation with PFA failed")
      } else {
        if (options[["intervalOn"]]) {
          if (is.null(out[["samp"]])) {
            startProgressbar(options[["noSamples"]])
            out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applyomega_pfa, callback = progressbarTick)
          }
          if (anyNA(out[["samp"]]))
            out[["error"]] <- gettext("Omega interval calculation with PFA failed")
          else
            out[["conf"]] <- quantile(out[["samp"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm = TRUE)
        }
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["omegaScaleObj"]] <- createJaspState(out,
                                                          dependencies = c("omegaScale", "omegaMethod",
                                                                           "omegaInterval", "confidenceIntervalValue"))
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
      out[["itemDropped"]] <- c(NA_real_, NA_real_)
      return(out)
    }

    if (options[["omegaMethod"]] == "cfa") {

      dataset <- scale(dataset, scale = FALSE)
      # do we have to compute item dropped values
        if (is.null(out[["itemDropped"]])) {
          out[["itemDropped"]] <- numeric(ncol(dataset))
          for (i in 1:ncol(dataset)) {
            out[["itemDropped"]][i] <- Bayesrel:::applyomega_cfa_data(dataset[, -i], interval = .95,
                                                                      pairwise = model[["pairwise"]])
          }
        }
        if (anyNA(out[["itemDropped"]]))
          out[["error"]] <- gettext("Omega item dropped statistics with CFA failed")

      } else { # omega with pfa
      # do we have to compute item dropped values
        if (is.null(out[["itemDropped"]]))
          out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applyomega_pfa)

        if (anyNA(out[["itemDropped"]]))
          out[["error"]] <- gettext("Omega item dropped statistics with PFA failed")

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

      out[["est"]] <- Bayesrel:::applyalpha(model[["data_cov"]])

      # do we need an interval estimate?
      if (options[["intervalOn"]]) {
        ciValue <- options[["confidenceIntervalValue"]]

        # should the interval be analytic
        if (options[["alphaInterval"]] == "alphaAnalytic") {
          out[["conf"]] <- Bayesrel:::ciAlpha(1 - ciValue, model[["n"]], model[["data_cov"]])

        } else {
          if (is.null(out[["samp"]])) {
            startProgressbar(options[["noSamples"]])
            out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applyalpha, callback = progressbarTick)
          }
          out[["conf"]] <- quantile(out[["samp"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2))
        }
      }

    } else { # alpha standardized
      ccor <- model[["data_cor"]]
      out[["est"]] <- Bayesrel:::applyalpha(ccor)

      # do we need an interval estimate?
      if (options[["intervalOn"]]) {
        ciValue <- options[["confidenceIntervalValue"]]

        # should the interval be analytic
        if (options[["alphaInterval"]] == "alphaAnalytic") {
          out[["conf"]] <- Bayesrel:::ciAlpha(1 - ciValue, model[["n"]], ccor)

        } else {

          if (is.null(out[["sampCor"]])) {
            out[["sampCor"]] <- numeric(options[["noSamples"]])
            for (i in 1:options[["noSamples"]]) {
              out[["sampCor"]][i] <- Bayesrel:::applyalpha(cov2cor(model[["bootSamp"]][i, ,]))
            }
          }
          out[["conf"]] <- quantile(out[["sampCor"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm = TRUE)

        }
      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["alphaScaleObj"]] <- createJaspState(out,
                                                         dependencies = c("alphaScale", "alphaMethod",
                                                                          "alphaInterval", "confidenceIntervalValue"))
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
      out[["itemDropped"]] <- c(NA_real_, NA_real_)
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

    out[["est"]] <- Bayesrel:::applylambda2(model[["data_cov"]])

    # do we need an interval estimate?
    if (options[["intervalOn"]]) {

      ciValue <- options[["confidenceIntervalValue"]]

      if (is.null(out[["samp"]])) {
        startProgressbar(options[["noSamples"]])
        out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applylambda2, callback = progressbarTick)
      }
      out[["conf"]] <- quantile(out[["samp"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm = TRUE)

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda2ScaleObj"]] <- createJaspState(out,
                                                           dependencies = c("lambda2Scale", "confidenceIntervalValue"))
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
      out[["itemDropped"]] <- c(NA_real_, NA_real_)
      return(out)
    }

    if (is.null(out[["itemDropped"]]))
      out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applylambda2)

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda2ItemObj"]] <- createJaspState(out, dependencies = c("lambda2Item"))
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

    out[["est"]] <- Bayesrel:::applylambda6(model[["data_cov"]])
    if (is.na(out[["est"]])) {
      out[["error"]] <- gettext("Lambda6 calculation failed")
    } else {
      # do we need an interval estimate?
      if (options[["intervalOn"]]) {

        ciValue <- options[["confidenceIntervalValue"]]

        if (is.null(out[["samp"]])) {
          startProgressbar(options[["noSamples"]])
          out[["samp"]] <- apply(model[["bootSamp"]], 1, Bayesrel:::applylambda6, callback = progressbarTick)
        }

        if (sum(!is.na(out[["samp"]])) < 3)
          out[["error"]] <- gettext("Lambda6 interval calculation failed")
        else
          out[["conf"]] <- quantile(out[["samp"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm = TRUE)

      }
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda6ScaleObj"]] <- createJaspState(out,
                                                           dependencies = c("lambda6Scale", "confidenceIntervalValue"))
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
      out[["itemDropped"]] <- c(NA_real_, NA_real_)
      return(out)
    }

    if (is.null(out[["itemDropped"]]))
      out[["itemDropped"]] <- .freqItemDroppedStats(model[["data_cov"]], Bayesrel:::applylambda6)
    if (anyNA(out[["itemDropped"]]))
      out[["error"]] <- gettext("Lambda6 item dropped statistics failed")

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["lambda6ItemObj"]] <- createJaspState(out, dependencies = c("lambda6Item"))
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

    out[["est"]] <- Bayesrel:::glbOnArray_custom(model[["data_cov"]])


    # do we need an interval estimate?
    if (options[["intervalOn"]]) {
      ciValue <- options[["confidenceIntervalValue"]]
      if (is.null(out[["samp"]])) {
        startProgressbar(options[["noSamples"]] %/% 500 + 1)
        out[["samp"]] <- Bayesrel:::glbOnArray_custom(model[["bootSamp"]], callback = progressbarTick)

      }
      out[["conf"]] <- quantile(out[["samp"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm = TRUE)
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["glbScaleObj"]] <- createJaspState(out,
                                                       dependencies = c("glbScale", "confidenceIntervalValue"))
  }
  return(out)
}

# check the error handling of the glb !!!!!!!!
.frequentistGlbItem <- function(jaspResults, dataset, options, model) {

  if (!is.null(.getStateContainerF(jaspResults)[["glbItemObj"]]$object))
    return(.getStateContainerF(jaspResults)[["glbItemObj"]]$object)

  out <- model[["glbItem"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["glbItem"]]  && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      out[["itemDropped"]] <- c(NA_real_, NA_real_)
      return(out)
    }

    # do we have to compute item dropped values
    if (is.null(out[["itemDropped"]])) {
      # special case glb since it has build in array functionality, but it might be only slightly faster
      itemDroppedCovs <- array(0, c(model[["k"]], model[["k"]]-1, model[["k"]]-1))
      for (i in 1:model[["k"]]) {
        itemDroppedCovs[i, , ] <- model[["data_cov"]][-i, -i]
      }
      out[["itemDropped"]] <- c(Bayesrel:::glbOnArray_custom(itemDroppedCovs, callback = progressbarTick))
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["glbItemObj"]] <- createJaspState(out, dependencies = c("glbItem"))
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
    ciValue <- options[["confidenceIntervalValue"]]

    out[["est"]] <- mean(model[["data_cor"]][lower.tri(model[["data_cor"]])])

    if (options[["intervalOn"]]) {

      if (is.null(out[["samp"]])) {
        startProgressbar(options[["noSamples"]])
        out[["samp"]] <- numeric(options[["noSamples"]])
        for (i in 1:options[["noSamples"]]) {
          corm <- .cov2cor.callback(model[["bootSamp"]][i, , ], progressbarTick)
          out[["samp"]][i] <- mean(corm[lower.tri(corm)])
        }
      }
      out[["conf"]] <- quantile(out[["samp"]], probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm = TRUE)
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["avgObj"]] <- createJaspState(out,
                                                  dependencies = c("averageInterItemCor", "confidenceIntervalValue"))
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

    if (options[["meanMethod"]] == "sumScores") {
      out[["est"]] <- mean(rowSums(dataset, na.rm = TRUE))
      sdmean <- sd(rowSums(dataset, na.rm = TRUE))
    } else {
      out[["est"]] <- mean(rowMeans(dataset, na.rm = TRUE))
      sdmean <- sd(rowMeans(dataset, na.rm = TRUE))
    }

    if (options[["intervalOn"]]) {
      zz <- qnorm(1-(1-ciValue)/2)
      out[["conf"]] <- c(out[["est"]] - zz*(sdmean/sqrt(model[["n"]])),
                         out[["est"]] + zz*(sdmean/sqrt(model[["n"]])))
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["meanObj"]] <- createJaspState(out, dependencies = c("meanScale", "meanMethod"))
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

    out[["est"]] <- if (options[["sdMethod"]] == "sumScores")
      sd(rowSums(dataset, na.rm = TRUE))
    else
      sd(rowMeans(dataset, na.rm = TRUE))

    if (options[["intervalOn"]]) {
      chiValueLow <- qchisq(1-(1-ciValue)/2, df = model[["n"]]-1)
      chiValueHigh <- qchisq((1-ciValue)/2, df = model[["n"]]-1)
      out[["conf"]] <- c(sqrt(((model[["n"]]-1) * out[["est"]]^2) / chiValueLow),
                         sqrt(((model[["n"]]-1) * out[["est"]]^2) / chiValueHigh))
    }

    if (options[["disableSampleSave"]])
      return(out)

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
      out[["itemDropped"]][i] <- cor(as.matrix(dataset[, i]), rowMeans(as.matrix(dataset[, -i]), na.rm = TRUE), use = model[["use.cases"]])
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemRestObj"]] <- createJaspState(out, dependencies = c("itemRestCor"))
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
    stateContainer[["meanItemObj"]] <- createJaspState(out, dependencies = c("meanItem"))
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
    stateContainer[["sdItemObj"]] <- createJaspState(out, dependencies = c("sdItem"))
  }
  return(out)
}



