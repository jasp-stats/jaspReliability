
.frequentistPreCalc <- function(jaspResults, dataset, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["modelObj"]]$object))
    return(.getStateContainerF(jaspResults)[["modelObj"]]$object)

  derivedOptions <- .frequentistDerivedOptions(options)

  # what if no coefficient boxes are checked?
  if(!any(derivedOptions[["selectedEstimators"]]) & !any(derivedOptions[["itemDroppedSelected"]])) {
    variables <- options[["variables"]]
    # if (length(options[["reverseScaledItems"]]) > 0L) {
    #   dataset <- .reverseScoreItems(dataset, options)
    # }
    empty <-  T
    model <- list(empty = empty)
    model[["footnote"]] <- .checkLoadings(dataset, variables)
    return(model)
  }

  # what if too few variables are entered:
  if (length(options[["variables"]]) < 3) {
    empty <-  T
    model <- list(empty = empty)
    model[["footnote"]] <- gettext("Please enter at least 3 Variables to do an analysis")
    return(model)
  }

  model <- jaspResults[["modelObj"]]$object

  if (is.null(model)) {
    model <- list()
    # check for missings and determine the missing handling
    if (anyNA(dataset)) {
      if (options[["missingValues"]] == "excludeCasesPairwise") {
        model[["use.cases"]] <- "pairwise.complete.obs"
        model[["pairwise"]] <- TRUE
        model[["footnote"]] <- gettextf("%s Of the observations, pairwise complete cases were used. ",
                                        model[["footnote"]])
      } else {
        pos <- which(is.na(dataset), arr.ind = TRUE)[, 1]
        dataset <- dataset[-pos, ]
        model[["use.cases"]] <- "complete.obs"
        model[["pairwise"]] <- FALSE
        model[["footnote"]] <- gettextf("%s Of the observations, %1.f complete cases were used. ",
                                        model[["footnote"]], nrow(dataset))
      }
    } else {
      model[["use.cases"]] <- "everything"
      model[["pairwise"]] <- FALSE
    }

    # # check for inverse scored items
    # if (length(options[["reverseScaledItems"]]) > 0L) {
    #   dataset <- .reverseScoreItems(dataset, options)
    # }
    p <- ncol(dataset); n <- nrow(dataset)
    cc <- cov(dataset, use = model[["use.cases"]])
    model[["data_cov"]] <- cc
    model[["data_cor"]] <- cor(dataset, use = model[["use.cases"]])

    # check if any items correlate negatively with the scale
    model[["footnote"]] <- .checkLoadings(dataset, options[["variables"]])
  }

  # check if interval is checked and bootstrapped covariance sample has to be generated
  if (is.null(model[["bootsamp"]]) &&
      options[["intervalOn"]] &&
      (
        (options[["omegaScale"]] && options[["omegaEst"]] == "pfa") ||
        (options[["omegaScale"]] && options[["omegaEst"]] == "cfa" && options[["omegaInterval"]] == "omegaBoot") ||
        (options[["alphaScale"]] && !(options[["alphaInterval"]] == "alphaAnalytic")) ||
        options[["lambda2Scale"]] ||
        options[["lambda6Scale"]] ||
        options[["glbScale"]] ||
        options[["averageInterItemCor"]]
      )
  ) {


    boot_cov <- array(0, c(options[["noSamples"]], p, p))

    startProgressbar(options[["noSamples"]])

    if (options[["setSeed"]])
      set.seed(options[["seedValue"]])

    if (options[["bootType"]] == "para") {
      model[["parametric"]] <- TRUE
      for (i in 1:options[["noSamples"]]) {
        boot_data <- MASS::mvrnorm(n, colMeans(dataset, na.rm = T), cc)
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
    model[["bootsamp"]] <- boot_cov
  }
  model[["itemsDropped"]] <- .unv(colnames(dataset))

  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model)

  return(model)
}

.frequentistItemDroppedMats <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerF(jaspResults)[["itemDroppedObj"]]$object))
    return(.getStateContainerF(jaspResults)[["itemDroppedObj"]]$object)

  out <- model[["itemsDroppedCovs"]]
  if (is.null(out) && is.null(model[["empty"]])) {
    cc <- model[["data_cov"]]
    p <- ncol(dataset)
    if (options[["omegaItem"]] || options[["alphaItem"]] || options[["lambda2Item"]] ||
         options[["lambda6Item"]] || options[["glbItem"]]) {
      Ctmp <- array(0, c(p, p - 1, p - 1))
      for (i in 1:p){
        Ctmp[i, , ] <- cc[-i, -i]
      }
      out <- Ctmp
    }

    stateContainer <- .getStateContainerF(jaspResults)
    stateContainer[["itemDroppedObj"]] <- createJaspState(out, dependencies = c("omegaItem", "alphaItem", "lambda2Item",
                                                                                "lambda6Item", "glbItem"))
  }
  return(out)
}

