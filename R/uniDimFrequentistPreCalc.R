
.frequentistPreCalc <- function(jaspResults, dataset, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["modelObj"]]$object)) {
    if (!is.null(.getStateContainerF(jaspResults)[["modelObj"]]$object$bootSamp))
      return(.getStateContainerF(jaspResults)[["modelObj"]]$object)
  }



  derivedOptions <- .frequentistDerivedOptions(options)

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

    k <- ncol(dataset); n <- nrow(dataset)
    model[["k"]] <- ncol(dataset)
    model[["n"]] <- nrow(dataset)
    cc <- cov(dataset, use = model[["use.cases"]])
    model[["data_cov"]] <- cc
    model[["data_cor"]] <- cor(dataset, use = model[["use.cases"]])
    model[["itemsDropped"]] <- .unv(colnames(dataset))


    # check if any items correlate negatively with the scale
    model[["footnote"]] <- .checkLoadings(dataset, options[["variables"]])
  }

  jaspBase::.setSeedJASP(options)

  # check if interval is checked and bootstrapped covariance sample has to be generated
  if (is.null(model[["bootSamp"]]) &&
      options[["intervalOn"]] &&
      (
        (options[["omegaScale"]] && options[["omegaMethod"]] == "pfa") ||
        (options[["alphaScale"]] && !(options[["alphaInterval"]] == "alphaAnalytic")) ||
        options[["lambda2Scale"]] ||
        options[["lambda6Scale"]] ||
        options[["glbScale"]] ||
        options[["averageInterItemCor"]]
      )
  ) {


    boot_cov <- array(0, c(options[["noSamples"]], k, k))

    startProgressbar(options[["noSamples"]])

    if (options[["bootType"]] == "parametric") {
      model[["parametric"]] <- TRUE
      for (i in 1:options[["noSamples"]]) {
        boot_data <- MASS::mvrnorm(n, colMeans(dataset, na.rm = TRUE), cc)
        boot_cov[i, , ] <- cov(boot_data)
        progressbarTick()
      }
    } else {
      model[["parametric"]] <- FALSE
      for (i in 1:options[["noSamples"]]){
        boot_data <- as.matrix(dataset[sample.int(n, size = n, replace = TRUE), ])
        boot_cov[i, , ] <- cov(boot_data, use = model[["use.cases"]])
        progressbarTick()
      }
    }
    model[["bootSamp"]] <- boot_cov
  }

  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model)

  return(model)
}

