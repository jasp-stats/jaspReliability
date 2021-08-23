
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
                                                                                 "meanMethod", "sdMethod",
                                                                                 "omegaMethod", "omegaInterval",
                                                                                 "alphaInterval", "alphaMethod"))


  }

  return(out)

}
