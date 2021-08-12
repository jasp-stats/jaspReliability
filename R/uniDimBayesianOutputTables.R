

.BayesianScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Bayesian Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("credibleIntervalValueScale", "meanScale", "sdScale", "rHat",
                                  "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                  "averageInterItemCor", "meanMethod", "sdMethod"))

  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  scaleTable$position <- 1
  stateContainer <- .getStateContainerB(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  interval <- gettextf("%s%% CI",
                       format(100 * options[["credibleIntervalValueScale"]], digits = 3, drop0trailing = TRUE))
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

  for (i in idxSelected) {
    scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
  }

  if (.is.empty(model)) {
    scaleTable$setData(allData)
    scaleTable$addFootnote(model[["footnote"]])
    return()
  }

  for (j in seq_along(idxSelected)) {
    i <- idxSelected[j]
    nm <- names(idxSelected[j])

    newData <- data.frame(est = c(unlist(model[[nm]][["est"]], use.names = FALSE),
                                  unlist(model[[nm]][["cred"]], use.names = FALSE)))

    if (options[["rHat"]]) {
      if (opts[i] == "mean" || opts[i] == "sd") {
        rhat <- NA_real_
      } else {
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


  twoItemProblem <- FALSE
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
        twoItemProblem <- TRUE
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

    if (!is.null(model[["twoItems"]]) && twoItemProblem)
      footnote <- gettextf("%s Please enter at least 3 variables for the if item dropped statistics.", footnote)
  }

  itemTable$addFootnote(footnote)

  return()
}



.BayesianProbTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["probTable"]]$object))
    return()

  # check if the values are in the proper order, meaning the lower field value is smaller than the upper field
  if (options[["probTableValueLow"]] > options[["probTableValueHigh"]]) {
    low <- options[["probTableValueHigh"]]
    high <- options[["probTableValueLow"]]
    note <- gettext("The bounds you entered have been rearranged in increasing order to provide meaningful results.")
  } else {
    low <- options[["probTableValueLow"]]
    high <- options[["probTableValueHigh"]]
    note <- ""
  }
  probTable <- createJaspTable(
    gettextf("Probability that Reliability Statistic is Larger than %.2f and Smaller than %.2f",
             low, high))
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
    if (n.item < 3 || n.item > 50) {
      noPriorSaved <- TRUE
    } else {
      priors <- Bayesrel:::priors[[as.character(n.item)]]
      # the prior names dont match the model names, thus rename the priors in their original order
      names(priors) <- c("alphaScale", "lambda2Scale", "lambda6Scale", "glbScale", "omegaScale")
      noPriorSaved <- FALSE
    }

    end <- 512
    xx <- seq(0, 1, length.out = 512)
    poslow <- end - sum(xx > low)
    poshigh <- end - sum(xx > high)
    # since the priors are only available in density form, the prior probability for the estimator being larger than
    # a cutoff is given by calculating the relative probability of the density from the cutoff to 1.
    # check this with R package
    probsPost <- numeric(sum(selected))
    probsPrior <- numeric(sum(selected))

    for (i in seq_len(length(idxSelected))) {
        nm <- names(idxSelected[i])
      samp_tmp <- as.vector(model[[nm]][["samp"]])
      probsPost[i] <- mean(samp_tmp > low) -
        mean(samp_tmp > high)

      if (noPriorSaved) {
        startProgressbar(4e3)
        prior <- .samplePrior(n.item, nm, progressbarTick)
      } else {
        prior <- priors[[nm]]
      }
      probsPrior[i] <- sum(prior[["y"]][poslow:end]) / sum(prior[["y"]]) -
        sum(prior[["y"]][poshigh:end]) / sum(prior[["y"]])

    }

    df <- data.frame(statistic = opts[idxSelected], prior = probsPrior, posterior = probsPost)
    probTable$setData(df)
    probTable$addFootnote(note)

    probTable$position <- 3
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["probTable"]] <- probTable
  }

  return()
}
