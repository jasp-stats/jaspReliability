

.BayesianPreCalc <- function(jaspResults, dataset, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object)) {
    if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object$gibbsSamp))
      return(.getStateContainerB(jaspResults)[["modelObj"]]$object)
  }


  derivedOptions <- .BayesianDerivedOptions(options)

  # what if no coefficient boxes are checked?
  if(!any(derivedOptions[["selectedEstimators"]]) && !any(derivedOptions[["itemDroppedSelected"]])) {
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

    # what if exactly two variables are entered:
    if (length(options[["variables"]]) == 2)
      model[["twoItems"]] <- TRUE

    # check for missings and determine the missing handling
    if (anyNA(dataset)) {
      if (options[["missingValues"]] == "excludeCasesPairwise") {
        model[["use.cases"]] <- "pairwise.complete.obs"
        model[["pairwise"]] <- TRUE
        model[["footnote"]] <- gettextf("%s Of the observations, pairwise complete cases were used. ",
                                        model[["footnote"]])
      } else {
        model[["use.cases"]] <- "complete.obs"
        model[["pairwise"]] <- FALSE
        model[["footnote"]] <- gettextf("%s Of the observations, %1.f complete cases were used. ",
                                        model[["footnote"]], nrow(dataset))
      }
    } else {
      model[["use.cases"]] <- "everything"
      model[["pairwise"]] <- FALSE
    }

    cc <- cov(dataset, use = model[["use.cases"]])
    model[["data_cov"]] <- cc
    model[["k"]] <- ncol(dataset)
    model[["n"]] <- nrow(dataset)

    # check if any items correlate negatively with the scale
    model[["footnote"]] <- .checkLoadings(dataset, options[["variables"]])
  }

  jaspBase::.setSeedJASP(options)

  # check if posterior cov sample already exists and any of the relevant coefficients are checked
  if (is.null(model[["gibbsSamp"]]) &&
      (options[["alphaScale"]] || options[["lambda2Scale"]] || options[["lambda6Scale"]] || options[["glbScale"]] ||
       options[["averageInterItemCor"]])
      ) {

    startProgressbar(options[["noSamples"]] * options[["noChains"]])
    dataset <- scale(dataset, scale = FALSE)
    model[["gibbsSamp"]] <- Bayesrel:::covSamp(dataset, options[["noSamples"]], options[["noBurnin"]],
                                               options[["noThin"]], options[["noChains"]],
                                               model[["pairwise"]], progressbarTick)$cov_mat

  }

  model[["progressbarLength"]] <- options[["noChains"]] *
    length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]]))

  model[["itemsDropped"]] <- .unv(colnames(dataset))

  if (options[["disableSampleSave"]])
    return(model)

  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model)

  return(model)
}

