

.frequentistScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Frequentist Scale Reliability Statistics"))

  scaleTable$dependOn(options = c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                  "averageInterItemCor", "meanScale", "sdScale", "scoresMethod",
                                  "omegaMethod", "omegaInterval", "alphaMethod", "alphaInterval"))
  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  if (options[["intervalOn"]]) {
    interval <- gettextf("%s%% CI",
                         format(100 * options[["confidenceIntervalValue"]], digits = 3, drop0trailing = TRUE))
    intervalLow <- gettextf("%s lower bound", interval)
    intervalUp <- gettextf("%s upper bound", interval)
    allData <- data.frame(estimate = c(gettext("Point estimate"), intervalLow, intervalUp))
  } else {
    allData <- data.frame(estimate = c(gettext("Point estimate")))
  }

  scaleTable$position <- 1
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  if (!is.null(model[["scaleResults"]][["error"]][["omegaScale"]])) {
    model[["footnote"]] <- paste(model[["footnote"]], model[["scaleResults"]][["error"]][["omegaScale"]])
  }
  if (!is.null(model[["scaleResults"]][["error"]][["lambda6Scale"]])) {
    model[["footnote"]] <- paste(model[["footnote"]], model[["scaleResults"]][["error"]][["lambda6Scale"]])
  }

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idxSelected <- which(selected)


  # if no coefficients selected:
  if (.is.empty(model)) {
    scaleTable$setData(allData)


    for (i in idxSelected) {
      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
    }

    if (model[["footnote"]] != "") {
      scaleTable$addFootnote(model[["footnote"]])
    }

    return()
  }

  if (options[["intervalOn"]]) {
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
      newData <- data.frame(est = c(unlist(model[["scaleResults"]][["est"]][[nm]], use.names = FALSE),
                                    unlist(model[["scaleResults"]][["conf"]][[nm]], use.names = FALSE)))
      colnames(newData) <- paste0(colnames(newData), i)
      allData <- cbind(allData, newData)
    }
  } else {
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
      newData <- data.frame(est = c(unlist(model[["scaleResults"]][["est"]][[nm]], use.names = FALSE)))
      colnames(newData) <- paste0(colnames(newData), i)
      allData <- cbind(allData, newData)
    }
  }

  scaleTable$setData(allData)

  if (model[["footnote"]] != "") {
    scaleTable$addFootnote(model[["footnote"]])
  }
  return()
}



.frequentistItemTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["itemTable"]]$object) ||
      !any(model[["derivedOptions"]][["itemDroppedSelected"]])) {
    return()
  }

  derivedOptions <- model[["derivedOptions"]]

  itemTable <- createJaspTable(gettext("Frequentist Individual Item Reliability Statistics"))
  itemTable$dependOn(options = c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item", "glbItem",
                                 "meanItem", "itemRestCor", "sdItem",
                                 "omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale"))
  # adding the scale options fixes a bug, where the item table would remain displayed
  # after one had checked a scale coefficient box and the item coefficient box and then unchecked the scale coeff box

  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemTable$position <- 2
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  footnote <- ""

  if (!is.null(model[["omegaItem"]][["error"]])) {
    footnote <- paste(footnote, model[["omegaItem"]][["error"]])
  }
  if (!is.null(model[["lambda6Item"]][["error"]])) {
    footnote <- paste(footnote, model[["omegaItem"]][["error"]])
  }

  selected <- derivedOptions[["itemDroppedSelected"]]
  coefficientsTable <- derivedOptions[["namesEstimators"]][["tables_item"]]
  overTitle <- gettext("If item dropped")
  idxSelected <- which(selected)
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]


  if (length(model[["itemsDropped"]]) > 0) {
    itemTable[["variable"]] <- model[["itemsDropped"]]

    if (!is.null(unlist(options[["reverseScaledItems"]])))
      footnote <- .addFootnoteReverseScaledItems(options)
  }

  fewItemProblem <- FALSE
  for (i in idxSelected) {
    if (coefficientsTable[i] %in% coefficients) {
      itemTable$addColumnInfo(name = paste0("pointEst", i), title = coefficientsTable[i], type = "number",
                              overtitle = overTitle)
      fewItemProblem <- TRUE
    } else {
      itemTable$addColumnInfo(name = paste0("pointEst", i), title = coefficientsTable[i], type = "number")
    }
  }

  if (is.null(model[["empty"]])) {
    tb <- data.frame(variable = model[["itemsDropped"]])

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])
      newtb <- cbind(pointEst = model[[nm]][["itemDropped"]])
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

  return()
}



# once the package is updated check this again and apply:
.frequentistSingleFactorFitTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["fitTable"]]$object) ||
      is.null(model[["scaleResults"]][["fit"]][["omegaScale"]]) || !options[["fitMeasures"]])
    return()

  fitTable <- createJaspTable(gettextf("Fit Measures of Single Factor Model Fit"))

  fitTable$addColumnInfo(name = "measure", title = gettext("Fit measure"),   type = "string")
  fitTable$addColumnInfo(name = "value",  title = gettext("Value"), type = "number")

  opts <- c("Chi-Square", "df", "p.value", "RMSEA", "Lower CI RMSEA", "Upper CI RMSEA", "SRMR")
  allData <- data.frame(
    measure = opts,
    value = as.vector(model[["scaleResults"]][["fit"]][["omegaScale"]])
  )
  fitTable$setData(allData)


  fitTable$dependOn(options = c("variables", "omegaScale", "reverseScaledItems", "fitMeasures", "missingValues",
                                "omegaMethod"))
  fitTable$position <- 3
  stateContainer <- .getStateContainerF(jaspResults)
  stateContainer[["fitTable"]] <- fitTable

  return()

}
