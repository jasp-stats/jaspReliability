

.BayesianScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Bayesian Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("credibleIntervalValueScale", "meanScale", "sdScale", "rHat",
                                  "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                  "averageInterItemCor", "scoresMethod"))

  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  interval <- gettextf("%s%% CI",
                       format(100 * options[["credibleIntervalValueScale"]], digits = 3, drop0trailing = TRUE))
  intervalLow <- gettextf("%s lower bound", interval)
  intervalUp <- gettextf("%s upper bound", interval)

  pointEst <- gettextf("Posterior %s", options[["pointEst"]])

  if (options[["rHat"]]) {
    allData <- data.frame(estimate = c(pointEst, intervalLow, intervalUp, "R-hat"))
  } else {
    allData <- data.frame(estimate = c(pointEst, intervalLow, intervalUp))
  }

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idxSelected <- which(selected)

  for (i in idxSelected) {
    scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
  }

  if (.is.empty(model)) {
    scaleTable$setData(allData)
    if (model[["footnote"]] != "") {
      scaleTable$addFootnote(model[["footnote"]])
    }
    scaleTable$position <- 1
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["scaleTable"]] <- scaleTable
    return()
  }

  for (j in seq_along(idxSelected)) {
    i <- idxSelected[j]
    nm <- names(idxSelected[j])

    newData <- data.frame(est = c(unlist(model[["scaleResults"]][["est"]][[nm]], use.names = FALSE),
                                  unlist(model[["scaleResults"]][["cred"]][[nm]], use.names = FALSE)))

    if (options[["rHat"]]) {
      newData <- rbind(newData, model[["scaleResults"]][["rHat"]][[nm]])
    }
    colnames(newData) <- paste0(colnames(newData), i)
    allData <- cbind(allData, newData)
  }

  scaleTable$setData(allData)

  if (model[["footnote"]] != "") {
    scaleTable$addFootnote(model[["footnote"]])
  }

  scaleTable$position <- 1
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  return()
}


.BayesianItemTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["itemTable"]]$object) ||
      !any(model[["derivedOptions"]][["itemDroppedSelected"]]))
    return()

  derivedOptions <- model[["derivedOptions"]]

  itemTable <- createJaspTable(gettext("Bayesian Individual Item Reliability Statistics"))

  itemTable$dependOn(options = c("omegaItem",  "alphaItem",  "lambda2Item",  "lambda6Item", "glbItem",
                                 "credibleIntervalValueItem", "itemRestCor", "meanItem", "sdItem",
                                 "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale"))

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemTable$position <- 2
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  overTitles <- format(derivedOptions[["namesEstimators"]][["tables_item"]], digits = 3, drop0trailing = TRUE)
  overTitles <- gettextf("%s (if item dropped)", overTitles)
  cred <- format(100 * options[["credibleIntervalValueItem"]], digits = 3, drop0trailing = TRUE)

  selected <- derivedOptions[["itemDroppedSelected"]]
  idxSelected <- which(selected)
  estimators <- derivedOptions[["namesEstimators"]][["tables_item"]]
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]


  itemTable[["variable"]] <- model[["itemsDropped"]]

  footnote <- ""
  if (!is.null(unlist(options[["reverseScaledItems"]])))
    footnote <- .addFootnoteReverseScaledItems(options)

  pointEst <- gettextf("Posterior %s", options[["pointEst"]])

  fewItemProblem <- FALSE
  for (i in idxSelected) {
    if (estimators[i] %in% coefficients) {
      if (estimators[i] == "Item-rest correlation") { # no item deleted for item rest cor
        itemTable$addColumnInfo(name = paste0("postMean", i), title = pointEst, type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
        itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", cred), type = "number",
                                overtitle = gettext("Item-rest correlation"))
      } else {
        itemTable$addColumnInfo(name = paste0("postMean", i), title = pointEst, type = "number",
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("lower", i), title = gettextf("Lower %s%% CI", cred), type = "number",
                                overtitle = overTitles[i])
        itemTable$addColumnInfo(name = paste0("upper", i), title = gettextf("Upper %s%% CI", cred), type = "number",
                                overtitle = overTitles[i])
        fewItemProblem <- TRUE
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

      if (i %in% 1:6) {
        rows <- length(options[["variables"]])
        if (rows < 3 && nm != "itemRestCor") {
          newtb <- cbind(postMean = rep(NaN, rows), matrix(NaN, rows, 2, dimnames = list(NULL, c("lower", "upper"))))
        } else {
          newtb <- cbind(postMean = model[["itemResults"]][["est"]][[nm]], model[["itemResults"]][["cred"]][[nm]])
        }
      } else {
        newtb <- cbind(postMean = model[["itemResults"]][["est"]][[nm]])
      }
      colnames(newtb) <- paste0(colnames(newtb), i)
      tb <- cbind(tb, newtb)
    }
    itemTable$setData(tb)

    if (length(options[["variables"]]) < 3 && fewItemProblem) {
      footnote <- gettextf("%s Please enter at least 3 variables for the if item dropped statistics.", footnote)
    }
  }

  if (footnote != "") {
    itemTable$addFootnote(footnote)
  }

  itemTable$position <- 2
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  return()
}



.BayesianProbTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["probTable"]]$object))
    return()

  # check if the values are in the proper order, meaning the lower field value is smaller than the upper field
  if (options[["probTableValueLow"]] > options[["probTableValueHigh"]]) {
    low <- options[["probTableValueHigh"]]
    high <- options[["probTableValueLow"]]
    footnote <- gettext("The bounds you entered have been rearranged in increasing order to provide meaningful results.")
  } else {
    low <- options[["probTableValueLow"]]
    high <- options[["probTableValueHigh"]]
    footnote <- ""
  }
  probTable <- createJaspTable(
    gettextf("Probability that Reliability Statistic is Larger than %.2f and Smaller than %.2f", low, high))
  probTable$dependOn(options = c("probTableValueLow", "probTable", "probTableValueHigh"))

  overTitle <- gettext("Probability")
  probTable$addColumnInfo(name = "statistic", title = gettext("Statistic"), type = "string")
  probTable$addColumnInfo(name = "prior", title = gettext("Prior"), type = "number", overtitle = overTitle)
  probTable$addColumnInfo(name = "posterior", title = gettext("Posterior"), type = "number", overtitle = overTitle)

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  idxSelected  <- which(selected)

  if (options[["probTable"]] && is.null(model[["empty"]])) {

    n.item <- model[["k"]]

    end <- 512
    xx <- seq(0, 1, length.out = 512)
    poslow <- end - sum(xx > low)
    poshigh <- end - sum(xx > high)

    probsPost <- numeric(sum(selected))
    probsPrior <- numeric(sum(selected))

    for (i in seq_len(length(idxSelected))) {
      nm <- names(idxSelected[i])

      samp_tmp <- as.vector(model[[nm]][["samp"]])
      probsPost[i] <- mean(samp_tmp > low) - mean(samp_tmp > high)

      if (nm == "omegaScale") {
        startProgressbar(2e3)
      } else {
        startProgressbar(4e3)
      }

      prior <- .samplePrior(n.item, nm, progressbarTick, options[["iwScale"]], options[["iwDf"]],
                            options[["igShape"]], options[["igScale"]])

      probsPrior[i] <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]]) -
        sum(prior[["y"]][poshigh:end]) / sum(prior[["y"]])

    }

    df <- data.frame(statistic = opts[idxSelected], prior = probsPrior, posterior = probsPost)
    probTable$setData(df)

    if (footnote != "") {
      probTable$addFootnote(footnote)
    }

    probTable$position <- 3
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["probTable"]] <- probTable
  }

  return()
}


.BayesianLoadingsTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["loadTable"]]$object))
    return()

  loadTable <- createJaspTable(gettextf("Omega Single-Factor Model"))

  loadTable$dependOn(options = c("omegaScale", "dispLoadings"))

  loadTable$addColumnInfo(name = "loadings", title = gettext("Standardized loadings"), type = "number")

  derivedOptions <- model[["derivedOptions"]]

  if (options[["omegaScale"]] && options[["dispLoadings"]] && is.null(model[["empty"]])) {

    df <- data.frame(loadings = model[["omegaScale"]][["loadingsStd"]])
    loadTable$setData(df)

    loadTable$position <- 4
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["loadTable"]] <- loadTable
  }

  return()
}

