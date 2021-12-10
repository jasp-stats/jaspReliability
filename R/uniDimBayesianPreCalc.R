

.BayesianPreCalc <- function(jaspResults, dataset, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object)) {
    if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object$gibbsSamp))
      return(.getStateContainerB(jaspResults)[["modelObj"]]$object)
  }


  derivedOptions <- .BayesianDerivedOptions(options)

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
    model[["itemsDropped"]] <- .unv(colnames(dataset))
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

    cc <- cov(dataset, use = model[["use.cases"]])
    model[["data_cov"]] <- cc
    model[["k"]] <- ncol(dataset)
    model[["n"]] <- nrow(dataset)

    # check if any items correlate negatively with the scale
    model[["footnote"]] <- sprintf("%s%s", model[["footnote"]], .checkLoadings(dataset, options[["variables"]]))
  }

  jaspBase::.setSeedJASP(options)

  # check if posterior cov sample already exists and any of the relevant coefficients are checked
  if (is.null(model[["gibbsSamp"]]) &&
      (options[["alphaScale"]] || options[["lambda2Scale"]] || options[["lambda6Scale"]] || options[["glbScale"]] ||
       options[["averageInterItemCor"]])
  ) {

    startProgressbar(options[["noSamples"]] * options[["noChains"]])
    dataset <- scale(dataset, scale = FALSE)
    c_out <- try(Bayesrel:::covSamp(dataset, options[["noSamples"]], options[["noBurnin"]],
                                    options[["noThin"]], options[["noChains"]],
                                    model[["pairwise"]], progressbarTick, k0 = options[["iwScale"]],
                                    df0 = options[["iwDf"]]), silent = TRUE)
    if (model[["pairwise"]] && inherits(c_out, "try-error")) {
      .quitAnalysis(gettext("Sampling the posterior covariance matrix for either one of [alpha, lambda2, lambda6, glb] failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
    }

    model[["gibbsSamp"]] <- c_out$cov_mat

  }

  model[["progressbarLength"]] <- options[["noChains"]] *
    length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]]))

  model[["itemsDropped"]] <- .unv(colnames(dataset))

  if (options[["disableSampleSave"]])
    return(model)

  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model, dependencies = c("iwScale", "iwDf"))

  return(model)
}


# standardize covariance matrix sample
.BayesianStdCov <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["gibbsCor"]]$object))
    return(.getStateContainerB(jaspResults)[["gibbsCor"]]$object)

  if (is.null(model[["gibbsSamp"]])) {
    return()
  } else {
    if (options[["stdCoeffs"]] == "unstand") {
      return()
    } else { # standardized
      out <- model[["gibbsSamp"]]
      startProgressbar(model[["progressbarLength"]])
      for (i in seq_len(nrow(model[["gibbsSamp"]]))) {
        for (j in seq_len(ncol(model[["gibbsSamp"]]))) {
          out[i, j, , ] <- .cov2cor.callback(model[["gibbsSamp"]][i, j, , ], progressbarTick)
        }
      }
      stateContainer <- .getStateContainerB(jaspResults)
      stateContainer[["gibbsCorObj"]] <- createJaspState(out, dependencies = "stdCoeffs")
    }
  }
  return(out)
}


.BayesianSingleFactorModel <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["singleFactorObj"]]$object))
    return(.getStateContainerB(jaspResults)[["singleFactorObj"]]$object)

  out <- model[["singleFactor"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaScale"]] && is.null(model[["empty"]])) {

    if (is.null(out[["loadings"]])) {
      startProgressbar(options[["noSamples"]] * options[["noChains"]])

      dataset <- scale(dataset, scale = FALSE)

      jaspBase::.setSeedJASP(options)

      tmp_out <- try(Bayesrel:::omegaSampler(dataset, options[["noSamples"]], options[["noBurnin"]],
                                             options[["noThin"]], options[["noChains"]],
                                             model[["pairwise"]], progressbarTick,
                                             a0 = options[["igShape"]], b0 = options[["igScale"]],
                                             m0 = options[["loadMean"]]), silent = TRUE)
      if (model[["pairwise"]] && inherits(tmp_out, "try-error")) {
        .quitAnalysis(gettext("Sampling the posterior factor model for omega failed. Try changing to 'Exclude cases listwise' in 'Advanced Options'"))
      }

      out[["loadings"]] <- tmp_out$lambda
      out[["residuals"]] <- tmp_out$psi
      out[["factor_var"]] <- tmp_out$phi
    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["singleFactorObj"]] <- createJaspState(out, dependencies = c("omegaScale", "igShape", "igScale",
                                                                                 "loadMean"))
  }

  return(out)
}


.BayesianSingleFactorModelItem <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["singleFactorItemObj"]]$object))
    return(.getStateContainerB(jaspResults)[["singleFactorItemObj"]]$object)

  out <- model[["singleFactorItem"]]
  if (is.null(out))
    out <- list()

  if (options[["omegaItem"]] && is.null(model[["empty"]])) {

    if (ncol(dataset) == 2) {
      return()
    }

    if (is.null(out[["itemLoadings"]])) {

      startProgressbar(options[["noSamples"]] * options[["noChains"]] * ncol(dataset))

      dataset <- scale(dataset, scale = FALSE)

      jaspBase::.setSeedJASP(options)

      out[["itemLoadings"]] <- array(0, c(options[["noChains"]],
                                          length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]])),
                                          ncol(dataset), ncol(dataset) - 1))
      out[["itemResiduals"]] <- array(0, c(options[["noChains"]],
                                           length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]])),
                                           ncol(dataset), ncol(dataset) - 1))

      for (i in seq_len(ncol(dataset))) {
        tmp_out <- Bayesrel:::omegaSampler(dataset[, -i],
                                           options[["noSamples"]], options[["noBurnin"]], options[["noThin"]],
                                           options[["noChains"]], model[["pairwise"]], progressbarTick,
                                           a0 = options[["igShape"]], b0 = options[["igScale"]],
                                           m0 = options[["loadMean"]])
        out[["itemLoadings"]][, , i, ] <- tmp_out$lambda
        out[["itemResiduals"]][, , i, ] <- tmp_out$psi
      }

      dd <- dim(out[["itemLoadings"]])
      out[["itemLoadings"]] <- array(out[["itemLoadings"]], c(dd[1] * dd[2], ncol(dataset), ncol(dataset) - 1))
      out[["itemResiduals"]] <- array(out[["itemResiduals"]], c(dd[1] * dd[2], ncol(dataset), ncol(dataset) - 1))

    }

    if (options[["disableSampleSave"]])
      return(out)

    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["singleFactorItemObj"]] <- createJaspState(out, dependencies = c("omegaItem", "igShape", "igScale",
                                                                                     "loadMean"))
  }

  return(out)
}

