

.BayesianScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Bayesian Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("credibleIntervalValueScale","meanScale", "sdScale", "rHat",
                                  "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                  "averageInterItemCor", "meanMethod", "sdMethod"))

  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")
  interval <- gettextf("%s%% CI",
                       format(100*options[["credibleIntervalValueScale"]], digits = 3, drop0trailing = TRUE))
  intervalLow <- gettextf("%s lower bound", interval)
  intervalUp <- gettextf("%s upper bound", interval)

  if (options[["rHat"]]) {
    allData <- data.frame(estimate = c("Posterior mean", intervalLow, intervalUp, "R-hat"))
  } else {
    allData <- data.frame(estimate = c("Posterior mean", intervalLow, intervalUp))
  }

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idxSelected <- which(selected)

  if (.is.empty(model)) {
    scaleTable$setData(allData)
    nvar <- length(options[["variables"]])
    for (i in idxSelected) {
      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
    }

    scaleTable$addFootnote(model[["footnote"]])
    scaleTable$position <- 1
    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["scaleTable"]] <- scaleTable
    return()
  }

  for (j in seq_along(idxSelected)) {
    i <- idxSelected[j]
    nm <- names(idxSelected[j])

    scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
    newData <- data.frame(est = c(unlist(model[[nm]][["est"]], use.names = FALSE),
                                  unlist(model[[nm]][["cred"]], use.names = FALSE)))

    if (options[["rHat"]]) {
      if (opts[i] == "mean" || opts[i] == "sd") {
        rhat <- NA_real_
      } else {
        # browser()
        tmp <- lapply(as.data.frame(t(model[[nm]][["samp"]])), coda::mcmc)
        rhat <- coda::gelman.diag(coda::as.mcmc.list(tmp))[["psrf"]][, 1]
      }
      newData <- rbind(newData, rhat)
    }
    colnames(newData) <- paste0(colnames(newData), i)
    allData <- cbind(allData, newData)
  }

  scaleTable$setData(allData)

  if (!is.null(model[["footnote"]]))
    scaleTable$addFootnote(model[["footnote"]])

  scaleTable$position <- 1
  stateContainerB <- .getStateContainerB(jaspResults)
  stateContainerB[["scaleTable"]] <- scaleTable

  return()
}


.BayesianItemTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["itemTable"]]$object) ||
      !any(model[["derivedOptions"]][["itemDroppedSelected"]]))
    return()

  derivedOptions <- model[["derivedOptions"]]

  overTitles <- format(derivedOptions[["namesEstimators"]][["tables_item"]], digits = 3, drop0trailing = TRUE)
  overTitles <- gettextf("%s (if item dropped)", overTitles)

  cred <- format(100*options[["credibleIntervalValueItem"]], digits = 3, drop0trailing = TRUE)
  itemTable <- createJaspTable(gettext("Bayesian Individual Item Reliability Statistics"))

  itemTable$dependOn(options = c("omegaItem",  "alphaItem",  "lambda2Item",  "lambda6Item", "glbItem",
                                 "credibleIntervalValueItem", "itemRestCor", "itemMean", "itemSd",
                                 "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale"))

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  selected <- derivedOptions[["itemDroppedSelected"]]
  idxSelected <- which(selected)
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]]
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]

  for (i in idxSelected) {
    if (estimators[i] %in% coefficients) {
      if (estimators[i] == "Item-rest correlation") { # no item deleted for item rest cor
        itemTable$addColumnInfo(name = paste0("postMean", i), title = gettext("Posterior Mean"), type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
      } else {
        itemTable$addColumnInfo(name = paste0("postMean", i), title = gettext("Posterior Mean"), type = "number",
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", cred), type = "number",
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", cred), type = "number",
                                overtitle = overTitles[i])
      }
    } else {
      itemTable$addColumnInfo(name = paste0("postMean", i), title = estimators[i], type = "number")
    }
  }

  if (is.null(model[["empty"]])) {
    tb <- data.frame(variable = model[["itemsDropped"]])
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (i %in% c(1:6)) { # check this when more estimators are included !!!!!!!!!!!!!!!!!!!!!
        newtb <- cbind(postMean = model[[nm]][["itemEst"]], model[[nm]][["itemCred"]])
      } else {
        newtb <- cbind(postMean = model[[nm]][["itemEst"]])
      }
      colnames(newtb) <- paste0(colnames(newtb), i)
      tb <- cbind(tb, newtb)
    }
    itemTable$setData(tb)

    if (!is.null(unlist(options[["reverseScaledItems"]]))) {
      itemTable$addFootnote(.addFootnoteReverseScaledItems(options))
    }

  } else if (length(model[["itemsDropped"]]) > 0) {
    itemTable[["variable"]] <- model[["itemsDropped"]]

    if (!is.null(unlist(options[["reverseScaledItems"]]))) {
      itemTable$addFootnote(.addFootnoteReverseScaledItems(options))
    }
  }


  itemTable$position <- 2
  stateContainerB <- .getStateContainerB(jaspResults)
  stateContainerB[["itemTable"]] <- itemTable

  return()
}



.BayesianProbTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["probTable"]]$object))
    return()

  probTable <- createJaspTable(
    gettextf("Probability that Reliability Statistic is Larger than %.2f and Smaller than %.2f",
             options[["probTableValueLow"]], options[["probTableValueHigh"]]))
  probTable$dependOn(options = c("probTableValueLow", "probTable", "probTableValueHigh"))

  overTitle <- gettext("Probability")
  probTable$addColumnInfo(name = "statistic", title = gettext("Statistic"),   type = "string")
  probTable$addColumnInfo(name = "prior",     title = gettext("Prior"), type = "number", overtitle = overTitle )
  probTable$addColumnInfo(name = "posterior", title = gettext("Posterior"), type = "number", overtitle = overTitle )

  derivedOptions <- model[["derivedOptions"]]
  order_end    <- derivedOptions[["order_end"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  idxSelected  <- which(selected)

  if (options[["probTable"]] && is.null(model[["empty"]])) {

    n.item <- model[["k"]]
    priors <- Bayesrel:::priors[[as.character(n.item)]]
    end <- length(priors[[1]][["x"]])
    poslow <- end - sum(priors[[1]][["x"]] > options[["probTableValueLow"]])
    poshigh <- end - sum(priors[[1]][["x"]] > options[["probTableValueHigh"]])
    # since the priors are only available in density form, the prior probability for the estimator being larger than
    # a cutoff is given by calculating the relative probability of the density from the cutoff to 1.
    # check this with R package

    probsPost <- numeric(sum(selected))
    probsPrior <- numeric(sum(selected))

    # the prior names dont match the model names, thus rename the priors in their original order
    names(priors) <- c("alphaScale", "lambda2Scale", "lambda6Scale", "glbScale", "omegaScale")

    for (i in 1:length(idxSelected)) {
        nm <- names(idxSelected[i])
      samp_tmp <- as.vector(model[[nm]][["samp"]])
      probsPost[i] <- mean(samp_tmp > options[["probTableValueLow"]]) -
        mean(samp_tmp > options[["probTableValueHigh"]])
      probsPrior[i] <- sum(priors[[nm]][["y"]][poslow:end]) / sum(priors[[nm]][["y"]]) -
        sum(priors[[nm]][["y"]][poshigh:end]) / sum(priors[[nm]][["y"]])

    }

    df <- data.frame(statistic = opts[idxSelected], prior = probsPrior, posterior = probsPost)
    probTable$setData(df)

    probTable$position <- 3
    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["probTable"]] <- probTable
  }

  return()
}
