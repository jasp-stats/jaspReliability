

.BayesianPreCalc <- function(jaspResults, dataset, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object)) {
    if (!is.null(.getStateContainerB(jaspResults)[["modelObj"]]$object$gibbsSamp))
      return(.getStateContainerB(jaspResults)[["modelObj"]]$object)
  }


  derivedOptions <- .BayesianDerivedOptions(options)

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

    cc <- cov(dataset, use = model[["use.cases"]])
    model[["data_cov"]] <- cc
    model[["k"]] <- ncol(dataset)
    model[["n"]] <- nrow(dataset)

    # check if any items correlate negatively with the scale
    model[["footnote"]] <- .checkLoadings(dataset, options[["variables"]])
  }

  # check if posterior cov sample already exists and any of the relevant coefficients are checked
  if (is.null(model[["gibbsSamp"]]) &&
      (options[["alphaScale"]] || options[["lambda2Scale"]] || options[["lambda6Scale"]] || options[["glbScale"]] ||
       options[["averageInterItemCor"]])
      ) {

    if (options[["setSeed"]])
      set.seed(options[["seedValue"]])

    startProgressbar(options[["noSamples"]] * options[["noChains"]])
    dataset <- scale(dataset, scale = F)
    model[["gibbsSamp"]] <- Bayesrel:::covSamp(dataset, options[["noSamples"]], options[["noBurnin"]],
                                               options[["noThin"]], options[["noChains"]],
                                               model[["pairwise"]], progressbarTick)$cov_mat

  }

  model[["itemsDropped"]] <- .unv(colnames(dataset))

  stateContainerB <- .getStateContainerB(jaspResults)
  stateContainerB[["modelObj"]] <- createJaspState(model)

  return(model)
}

# change options when scale box is unchecked
.scaleItemBoxAlign <- function(options) {
  opts <- options
  if(!options[["omegaScale"]])
    opts[["omegaItem"]] <- F
  if(!options[["alphaScale"]])
    opts[["alphaItem"]] <- F
  if(!options[["lambda2Scale"]])
    opts[["lambda2Item"]] <- F
  if(!options[["lambda6Scale"]])
    opts[["lambda6Item"]] <- F
  if(!options[["glbScale"]])
    opts[["glbItem"]] <- F

  return(opts)

}


.BayesianItemDroppedMats <- function(jaspResults, dataset, options, model) {
  if (!is.null(.getStateContainerB(jaspResults)[["itemDroppedObj"]]$object))
    return(.getStateContainerB(jaspResults)[["itemDroppedObj"]]$object)

  out <- model[["itemDroppedCovs"]]
  if (is.null(out)) {
    if (is.null(model[["empty"]]) &&
        (options[["alphaItem"]] || options[["lambda2Item"]] || options[["lambda6Item"]] || options[["glbItem"]])) {
      cc <- model[["gibbsSamp"]]
      p <- ncol(dataset)
      out <- array(0, c(options[["noChains"]],
                         length(seq(1, options[["noSamples"]] - options[["noBurnin"]], options[["noThin"]])),
                         p, p - 1, p - 1))
      for (i in 1:p){
        out[, , i, , ] <- cc[ , , -i, -i]
      }
    }

    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["itemDroppedObj"]] <- createJaspState(out, dependencies = c("alphaItem", "lambda2Item",
                                                                                "lambda6Item", "glbItem"))
  }
  return(out)
}

