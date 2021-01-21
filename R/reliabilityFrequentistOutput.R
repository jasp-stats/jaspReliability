

.frequentistScaleTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["scaleTable"]]$object))
      return()

  scaleTable <- createJaspTable(gettext("Frequentist Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                  "averageInterItemCor", "meanScale", "sdScale", "meanMethod", "sdMethod"))
  scaleTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")

  if (options[["intervalOn"]]) {
    intervalLow <- gettextf("%s%% CI",
                            format(100*options[["confidenceIntervalValue"]], digits = 3, drop0trailing = TRUE))
    intervalUp <- gettextf("%s%% CI",
                           format(100*options[["confidenceIntervalValue"]], digits = 3, drop0trailing = TRUE))
    intervalLow <- gettextf("%s lower bound", intervalLow)
    intervalUp <- gettextf("%s upper bound", intervalUp)
    allData <- data.frame(estimate = c(gettext("Point estimate"), intervalLow, intervalUp))
  } else {
    allData <- data.frame(estimate = c(gettext("Point estimate")))
  }

  derivedOptions <- model[["derivedOptions"]]
  opts     <- derivedOptions[["namesEstimators"]][["tables"]]
  selected <- derivedOptions[["selectedEstimators"]]
  idxSelected <- which(selected)
  # match names from model object with the names from selected coefficients to find their index:
  idMatchedNames <- which(!is.na(charmatch(names(model), names(selected[idxSelected]))))

  # if no coefficients selected:
  if (!is.null(model[["empty"]])) {
    scaleTable$setData(allData)
    nvar <- length(options[["variables"]])
    if (sum(selected) > 0L && nvar < 3) {
      for (i in idxSelected) {
        scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
      }
    }
    scaleTable$addFootnote(model[["footnote"]])
    scaleTable$position <- 1
    stateContainerF <- .getStateContainerF(jaspResults)
    stateContainerF[["scaleTable"]] <- scaleTable
    return()
  }

  if (options[["intervalOn"]]) {
    z <- 1
    for (i in idxSelected) {
      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
      newData <- data.frame(est = c(unlist(model[[idMatchedNames[z]]][["est"]], use.names = F),
                                    unlist(model[[idMatchedNames[z]]][["conf"]], use.names = F)))
      colnames(newData) <- paste0(colnames(newData), i)
      allData <- cbind(allData, newData)
      z <- z+1
    }
  } else {
    z <- 1
    for (i in idxSelected) {
      scaleTable$addColumnInfo(name = paste0("est", i), title = opts[i], type = "number")
      newData <- data.frame(est = c(unlist(model[[idMatchedNames[z]]][["est"]], use.names = F)))
      colnames(newData) <- paste0(colnames(newData), i)
      allData <- cbind(allData, newData)
      z <- z+1
    }
  }

  scaleTable$setData(allData)

  # if (!is.null(model[["error"]]))
  #   scaleTable$setError(model[["error"]])

  if (!is.null(model[["footnote"]]))
    scaleTable$addFootnote(model[["footnote"]])

  scaleTable$position <- 1
  stateContainerF <- .getStateContainerF(jaspResults)
  stateContainerF[["scaleTable"]] <- scaleTable

  return()
}



.frequentistItemTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["itemTable"]]$object) ||
      !any(model[["derivedOptions"]][["itemDroppedSelected"]]))
    return()

  derivedOptions <- model[["derivedOptions"]]

  # fixes issue that unchecking the scale coefficient box, does not uncheck the item-dropped coefficient box:
  for (i in 1:5) {
    if (!derivedOptions[["selectedEstimators"]][i]) {
      derivedOptions[["itemDroppedSelected"]][i] <- derivedOptions[["selectedEstimators"]][i]
    }
  }

  itemTable <- createJaspTable(gettext("Frequentist Individual Item Reliability Statistics"))
  itemTable$dependOn(options = c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item", "glbItem",
                                 "itemMean", "itemRestCor", "itemSd"))
  itemTable$addColumnInfo(name = "variable", title = gettext("Item"), type = "string")

  itemDroppedSelected <- derivedOptions[["itemDroppedSelected"]]
  coefficientsTable <- derivedOptions[["namesEstimators"]][["tables_item"]]
  overTitle <- gettext("If item dropped")

  idxSelected <- which(itemDroppedSelected)
  coefficients <- derivedOptions[["namesEstimators"]][["coefficients"]]
  idMatchedNames <- which(!is.na(charmatch(names(model), names(itemDroppedSelected[idxSelected]))))

  for (i in idxSelected) {
    if (coefficientsTable[i] %in% coefficients) {
      itemTable$addColumnInfo(name = paste0("pointEst", i), title = coefficientsTable[i], type = "number",
                              overtitle = overTitle)
    } else {
      itemTable$addColumnInfo(name = paste0("pointEst", i), title = coefficientsTable[i], type = "number")
    }
  }

  if (is.null(model[["empty"]])) {
    tb <- data.frame(variable = model[["itemsDropped"]])
    z <- 1
    for (i in idxSelected) {
      newtb <- cbind(pointEst = model[[idMatchedNames[z]]][["itemDropped"]])
      colnames(newtb) <- paste0(colnames(newtb), i)
      tb <- cbind(tb, newtb)
      z <- z+1
    }
    itemTable$setData(tb)

    if (!is.null(unlist(options[["reverseScaledItems"]]))) {
      itemTable$addFootnote(sprintf(ngettext(length(options[["reverseScaledItems"]]),
                                             "The following item was reverse scaled: %s. ",
                                             "The following items were reverse scaled: %s. "),
                                    paste(options[["reverseScaledItems"]], collapse = ", ")))
    }

  } else if (length(model[["itemsDropped"]]) > 0) {
    itemTable[["variables"]] <- model[["itemsDropped"]]

    if (!is.null(unlist(options[["reverseScaledItems"]]))) {
      itemTable$addFootnote(sprintf(ngettext(length(options[["reverseScaledItems"]]),
                                             "The following item was reverse scaled: %s. ",
                                             "The following items were reverse scaled: %s. "),
                                    paste(options[["reverseScaledItems"]], collapse = ", ")))
    }
  }

  itemTable$position <- 2
  stateContainerF <- .getStateContainerF(jaspResults)
  stateContainerF[["itemTable"]] <- itemTable

  return()
}



# once the package is updated check this again and apply:
.frequentistSingleFactorFitTable <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerF(jaspResults)[["fitTable"]]$object) || is.null(model[["omega"]][["omegaFit"]]) ||
      !options[["fitMeasures"]])
    return()

  fitTable <- createJaspTable(gettextf("Fit Measures of Single Factor Model Fit"))

  fitTable$addColumnInfo(name = "measure", title = gettext("Fit Measure"),   type = "string")
  fitTable$addColumnInfo(name = "value",  title = gettext("Value"), type = "number")

  opts <- c("Chi-Square", "df", "p.value", "RMSEA", "Lower CI RMSEA", "Upper CI RMSEA", "SRMR")
  allData <- data.frame(
    measure = opts,
    value = as.vector(model[["omega"]][["omegaFit"]])
  )
  fitTable$setData(allData)

  # if (!is.null(model[["error"]]))
  #   fitTable$setError(model[["error"]])

  fitTable$dependOn(options = c("variables", "omegaScale", "reverseScaledItems", "fitMeasures", "missingValues",
                                "omegaMethod"))
  fitTable$position <- 3
  stateContainerF <- .getStateContainerF(jaspResults)
  stateContainerF[["fitTable"]] <- fitTable

  return()

}
