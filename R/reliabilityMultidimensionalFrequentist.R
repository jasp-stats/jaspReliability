#' @importFrom jaspBase createJaspContainer createJaspHtml createJaspPlot createJaspState createJaspTable progressbarTick startProgressbar
#' @importFrom stats cor cov sd complete.cases quantile qnorm

# Dependencies that invalidate the fitted model (and therefore everything downstream).
# ciLevel is included because the Wald interval is produced by the fit itself; refitting is cheap.
.multiDimFreqBaseDependencies <- c(
  "factors", "reverseScaledItems", "modelType", "naAction", "ciLevel",
  "intervalMethod", "bootstrapSamples", "setSeed", "seed", "samplesSavingDisabled"
)

# the bootstrap samples are reused across confidence levels, so they are kept outside the state
# container and depend on everything above except ciLevel
.multiDimFreqBootDependencies <- setdiff(.multiDimFreqBaseDependencies, "ciLevel")

#' @export
reliabilityMultidimensionalFrequentistInternal <- function(jaspResults, dataset, options) {

  allItems <- .multiDimGetItems(options)
  ready    <- .multiDimReady(options)

  # option names match the other reliability analyses, so the shared helpers only need the item list
  optTmp <- options
  optTmp[["variables"]] <- allItems

  if (ready) {
    dataset <- dataset[, allItems, drop = FALSE]
    if (length(options[["reverseScaledItems"]]) > 0L)
      dataset <- .reverseScoreItems(dataset, optTmp)
    # under listwise deletion the model sees only the complete cases, so the observation and
    # variance checks must be evaluated on those rows rather than on the full dataset
    .checkErrors(.multiDimFreqAnalysisData(dataset, options), optTmp)
  }

  model <- .multiDimFreqComputeModel(jaspResults, dataset, options, ready, allItems)

  .multiDimFreqScaleTable(jaspResults, dataset, model, options, ready)
  .multiDimFreqItemTable(jaspResults, dataset, model, options, ready, allItems)
  .multiDimFreqFitTable(jaspResults, model, options, ready)
  .multiDimFreqLoadingsTable(jaspResults, model, options, ready, allItems)

  return()
}


#### Helpers ####

.getStateContainerMDF <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = .multiDimFreqBaseDependencies)
  return(jaspResults[["stateContainer"]])
}

# descriptive statistics computed outside the factor model (scale scores, item-rest correlations)
# must use the same rows as the fit: complete cases under listwise deletion, everything under FIML
.multiDimFreqAnalysisData <- function(dataset, options) {
  if (options[["naAction"]] != "listwise")
    return(dataset)
  return(dataset[complete.cases(dataset), , drop = FALSE])
}

.multiDimFreqMissing <- function(options) {
  if (options[["naAction"]] == "listwise") "listwise" else "fiml"
}

# the names shown as column headers of the loadings table, taken from the factors form so the
# table matches what the user typed rather than showing f1, f2, ...
.multiDimFreqFactorTitles <- function(options) {
  facs  <- options[["factors"]]
  keep  <- vapply(facs, function(f) length(unlist(f[["indicators"]])) > 0L, logical(1))
  vapply(facs[keep], function(f) as.character(f[["title"]]), character(1))
}

# one fit of the factor model. lavaan reports convergence and identification problems through
# warnings; these are suppressed here and reported through the diagnostics instead, so the user
# gets a translated footnote rather than raw lavaan output
.multiDimFreqFit <- function(dataset, options, modelSyntax, fitMeasures) {
  try(suppressWarnings(Bayesrel::omegasCFA(
    data         = as.data.frame(dataset),
    model        = modelSyntax,
    model.type   = options[["modelType"]],
    interval     = options[["ciLevel"]],
    missing      = .multiDimFreqMissing(options),
    fit.measures = fitMeasures
  )), silent = TRUE)
}

# builds the model syntax for the reduced item set used by the if-item-dropped statistics
.multiDimFreqReducedSyntax <- function(factors, item) {
  reduced <- lapply(factors, function(f) setdiff(unlist(f[["indicators"]]), item))
  sizes   <- vapply(reduced, length, integer(1))
  if (any(sizes[sizes > 0L] < 2L))        # would under-identify a factor
    return(NULL)
  lines <- character(0)
  fi    <- 0L
  for (f in reduced) {
    if (length(f) == 0L) next
    fi    <- fi + 1L
    lines <- c(lines, sprintf("f%d =~ %s", fi, paste(f, collapse = " + ")))
  }
  if (fi < 2L)                            # fewer than two factors left
    return(NULL)
  return(paste(lines, collapse = "\n"))
}

.multiDimFreqOmegaTLabel <- function() gettext("McDonald's ωₜ")
.multiDimFreqOmegaHLabel <- function() gettext("McDonald's ωₕ")

# shown alongside omega_h when the general factor is not identified. The likelihood is flat along
# the direction that separates the two structural loadings, so the estimate is one arbitrary point
# on that ridge and no standard error exists for it.
.multiDimFreqOmegaHNotIdentifiedNote <- function() {
  return(gettext("McDonald's ωₕ is not identified with only two group factors: the loadings of the general factor are determined only up to their product, so the value shown is arbitrary and has no confidence interval. Assign the items to at least three group factors, or choose the bi-factor model. McDonald's ωₜ is unaffected."))
}

# footnotes derived from the state of the fitted model; these are the conditions under which
# lavaan still returns numbers that should not be read as ordinary estimates
.multiDimFreqDiagnosticsFootnotes <- function(model, options) {
  diagnostics <- model[["fit"]][["diagnostics"]]
  if (is.null(diagnostics))
    return(character(0))

  notes <- character(0)
  if (!diagnostics[["se.available"]] && model[["intervalMethod"]] == "analytic") {
    # when the general factor is known to be unidentified, that footnote already names the cause,
    # so this one only reports the consequence instead of guessing at it a second time
    notes <- c(notes, if (.multiDimGeneralFactorUnidentified(options))
      gettext("Standard errors could not be computed for this model, so no confidence intervals are reported.")
    else
      gettext("Standard errors could not be computed, so no confidence intervals are reported. The factor model is probably not identified."))
  }
  if (!diagnostics[["admissible"]])
    notes <- c(notes, gettext("The factor model solution is inadmissible, for example because of a negative variance estimate. The coefficients may fall outside the interval [0, 1]."))

  return(notes)
}


#### Compute ####

.multiDimFreqComputeModel <- function(jaspResults, dataset, options, ready, allItems) {

  if (!is.null(.getStateContainerMDF(jaspResults)[["modelObj"]]$object))
    return(.getStateContainerMDF(jaspResults)[["modelObj"]]$object)

  model <- list(modelType      = options[["modelType"]],
                k              = length(allItems),
                items          = allItems,
                intervalMethod = options[["intervalMethod"]])

  if (!ready)
    return(model)

  # caught here rather than by Bayesrel's own stop(), so the message is translated and actionable
  if (options[["modelType"]] == "biFactor" && .multiDimHasCrossLoadings(options)) {
    model[["error"]] <- gettext("The bi-factor model does not support items that load on more than one factor. Assign each item to a single factor, or select the second-order or correlated-factors model.")
    return(model)
  }

  fit <- .multiDimFreqFit(dataset, options, .multiDimBuildModel(options), options[["fitMeasures"]])

  if (inherits(fit, "try-error")) {
    model[["error"]] <- jaspBase::.extractErrorMessage(fit)
    return(model)
  }

  # a model the optimizer never reached produces coefficients that must not be interpreted at all,
  # so this is an error rather than a footnote
  if (!fit[["diagnostics"]][["converged"]]) {
    model[["error"]] <- gettext("The factor model did not converge, so the coefficients cannot be interpreted. Consider a different model, or check the items assigned to each factor.")
    return(model)
  }

  model[["fit"]] <- fit

  stateContainer <- .getStateContainerMDF(jaspResults)
  stateContainer[["modelObj"]] <- createJaspState(model)

  return(model)
}


# Percentile bootstrap for the omega coefficients: the Wald interval of a ratio such as omega is
# symmetric and can run outside [0, 1], so resampling is offered as an alternative. The samples
# are cached independently of ciLevel, so changing the confidence level only recomputes quantiles.
.multiDimFreqBootstrap <- function(jaspResults, dataset, options, model, ready) {

  if (options[["intervalMethod"]] != "bootstrapped" || !ready || !is.null(model[["error"]]))
    return(NULL)

  if (!is.null(jaspResults[["bootSamples"]]$object))
    return(jaspResults[["bootSamples"]]$object)

  nBoot       <- options[["bootstrapSamples"]]
  modelSyntax <- .multiDimBuildModel(options)
  correlated  <- options[["modelType"]] == "correlated"
  bootData    <- .multiDimFreqAnalysisData(dataset, options)
  n           <- nrow(bootData)

  samples <- matrix(NA_real_, nBoot, 2L, dimnames = list(NULL, c("omegaT", "omegaH")))

  jaspBase::.setSeedJASP(options)
  startProgressbar(nBoot, label = gettext("Bootstrapping the omega coefficients"))

  for (i in seq_len(nBoot)) {
    resample <- bootData[sample.int(n, size = n, replace = TRUE), , drop = FALSE]
    fitBoot  <- .multiDimFreqFit(resample, options, modelSyntax, FALSE)
    # a resample that fails to converge contributes nothing rather than a misleading value
    if (!inherits(fitBoot, "try-error") && isTRUE(fitBoot[["diagnostics"]][["converged"]])) {
      samples[i, "omegaT"] <- fitBoot[["omega_t"]][["est"]]
      if (!correlated)
        samples[i, "omegaH"] <- fitBoot[["omega_h"]][["est"]]
    }
    progressbarTick()
  }

  if (options[["samplesSavingDisabled"]])
    return(samples)

  bootState <- createJaspState(samples)
  bootState$dependOn(options = .multiDimFreqBootDependencies)
  jaspResults[["bootSamples"]] <- bootState

  return(samples)
}


# omega if an item is dropped: the model is refit once per item with that item removed from its
# factor(s). Items whose removal would leave a factor with fewer than two indicators are skipped.
.multiDimFreqItemDeleted <- function(jaspResults, dataset, options, allItems) {

  stateContainer <- .getStateContainerMDF(jaspResults)
  if (!is.null(stateContainer[["itemDeletedObj"]]$object))
    return(stateContainer[["itemDeletedObj"]]$object)

  k          <- length(allItems)
  correlated <- options[["modelType"]] == "correlated"
  omegaT     <- rep(NA_real_, k)
  omegaH     <- rep(NA_real_, k)

  startProgressbar(k, label = gettext("Computing item statistics (refitting the model once per item)"))

  for (i in seq_len(k)) {
    syntax <- .multiDimFreqReducedSyntax(options[["factors"]], allItems[i])
    if (is.null(syntax)) {
      progressbarTick()
      next
    }
    reducedData <- dataset[, setdiff(allItems, allItems[i]), drop = FALSE]
    fitReduced  <- .multiDimFreqFit(reducedData, options, syntax, FALSE)

    if (!inherits(fitReduced, "try-error") && isTRUE(fitReduced[["diagnostics"]][["converged"]])) {
      omegaT[i] <- fitReduced[["omega_t"]][["est"]]
      if (!correlated)
        omegaH[i] <- fitReduced[["omega_h"]][["est"]]
    }
    progressbarTick()
  }

  out <- list(omegaT = omegaT, omegaH = omegaH)
  stateContainer[["itemDeletedObj"]] <- createJaspState(out, dependencies = c("itemDeletedOmegaT", "itemDeletedOmegaH"))

  return(out)
}


#### Tables ####

.multiDimFreqScaleTable <- function(jaspResults, dataset, model, options, ready) {

  if (!is.null(.getStateContainerMDF(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Frequentist Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("meanSdScoresMethod"))
  scaleTable$position <- 1

  ci <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)
  scaleTable$addColumnInfo(name = "coefficient", title = gettext("Coefficient"), type = "string")
  scaleTable$addColumnInfo(name = "estimate",    title = gettext("Estimate"),    type = "number")
  scaleTable$addColumnInfo(name = "lower", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI", ci))
  scaleTable$addColumnInfo(name = "upper", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI", ci))

  stateContainer <- .getStateContainerMDF(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  if (!ready) {
    scaleTable$addFootnote(gettext("Please assign at least two factors with at least two items each."))
    return()
  }
  if (!is.null(model[["error"]])) {
    scaleTable$setError(model[["error"]])
    return()
  }

  fit        <- model[["fit"]]
  correlated <- options[["modelType"]] == "correlated"
  bootSamples <- .multiDimFreqBootstrap(jaspResults, dataset, options, model, ready)

  # the interval is taken from the resampling distribution when the user asked for it, and from
  # the Wald interval of the fit otherwise
  interval <- function(key, waldConf) {
    if (is.null(bootSamples))
      return(if (is.null(waldConf)) c(NA_real_, NA_real_) else waldConf)
    samp <- bootSamples[, key]
    if (all(is.na(samp)))
      return(c(NA_real_, NA_real_))
    probs <- c((1 - options[["ciLevel"]]) / 2, 1 - (1 - options[["ciLevel"]]) / 2)
    unname(quantile(samp, probs = probs, na.rm = TRUE))
  }

  rows <- list()
  addCoefRow <- function(label, key, est, waldConf) {
    ci <- interval(key, waldConf)
    list(coefficient = label, estimate = est, lower = ci[1], upper = ci[2])
  }
  addStatRow <- function(label, value) {
    list(coefficient = label, estimate = value, lower = NA_real_, upper = NA_real_)
  }

  rows[[length(rows) + 1L]] <- addCoefRow(.multiDimFreqOmegaTLabel(), "omegaT",
                                          fit[["omega_t"]][["est"]], fit[["omega_t"]][["conf"]])

  footnotes <- .multiDimFreqDiagnosticsFootnotes(model, options)
  if (.multiDimHasOmegaH(options)) {
    rows[[length(rows) + 1L]] <- addCoefRow(.multiDimFreqOmegaHLabel(), "omegaH",
                                            fit[["omega_h"]][["est"]], fit[["omega_h"]][["conf"]])
    if (.multiDimGeneralFactorUnidentified(options))
      footnotes <- c(footnotes, .multiDimFreqOmegaHNotIdentifiedNote())
  } else {
    footnotes <- c(footnotes, .multiDimOmegaHFootnote(options))
  }

  scoreData <- .multiDimFreqAnalysisData(dataset, options)
  pairwise  <- options[["naAction"]] != "listwise"
  cc <- cor(dataset, use = if (pairwise) "pairwise.complete.obs" else "complete.obs")
  rows[[length(rows) + 1L]] <- addStatRow(gettext("Average interitem correlation"), mean(cc[lower.tri(cc)]))

  scores <- if (options[["meanSdScoresMethod"]] == "sumScores")
    rowSums(scoreData, na.rm = TRUE) else rowMeans(scoreData, na.rm = TRUE)
  rows[[length(rows) + 1L]] <- addStatRow(gettext("Mean"), mean(scores))
  rows[[length(rows) + 1L]] <- addStatRow(gettext("SD"), sd(scores))

  scaleTable$setData(do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE)))

  if (options[["intervalMethod"]] == "bootstrapped" && !is.null(bootSamples)) {
    failed <- sum(is.na(bootSamples[, "omegaT"]))
    if (failed > 0L)
      footnotes <- c(footnotes, gettextf("%1$i of %2$i bootstrap samples did not converge and were excluded from the interval.",
                                         failed, nrow(bootSamples)))
  }
  for (note in footnotes)
    scaleTable$addFootnote(note)

  return()
}


.multiDimFreqItemTable <- function(jaspResults, dataset, model, options, ready, allItems) {

  showOmegaT <- options[["itemDeletedOmegaT"]]
  showOmegaH <- options[["itemDeletedOmegaH"]] && .multiDimHasOmegaH(options)
  showRest   <- options[["itemRestCorrelation"]]
  if (!(showOmegaT || showOmegaH || showRest) ||
      !is.null(.getStateContainerMDF(jaspResults)[["itemTable"]]$object))
    return()

  itemTable <- createJaspTable(gettext("Frequentist Individual Item Reliability Statistics"))
  itemTable$dependOn(options = c("itemDeletedOmegaT", "itemDeletedOmegaH", "itemRestCorrelation"))
  itemTable$position <- 2

  itemTable$addColumnInfo(name = "item", title = gettext("Item"), type = "string")
  if (showOmegaT)
    itemTable$addColumnInfo(name = "omegaT", title = gettext("McDonald's ωₜ (if item dropped)"), type = "number")
  if (showOmegaH)
    itemTable$addColumnInfo(name = "omegaH", title = gettext("McDonald's ωₕ (if item dropped)"), type = "number")
  if (showRest)
    itemTable$addColumnInfo(name = "itemRestCorrelation", title = gettext("Item-rest correlation"), type = "number")

  stateContainer <- .getStateContainerMDF(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  if (!ready)
    return()
  if (!is.null(model[["error"]])) {
    itemTable$setError(model[["error"]])
    return()
  }

  df <- data.frame(item = jaspBase::decodeColNames(allItems), stringsAsFactors = FALSE)

  footnotes <- character(0)
  if (showOmegaT || showOmegaH) {
    deleted <- .multiDimFreqItemDeleted(jaspResults, dataset, options, allItems)
    if (showOmegaT) df[["omegaT"]] <- deleted[["omegaT"]]
    if (showOmegaH) df[["omegaH"]] <- deleted[["omegaH"]]
    if (anyNA(c(if (showOmegaT) df[["omegaT"]], if (showOmegaH) df[["omegaH"]])))
      footnotes <- c(footnotes, gettext("Empty cells indicate that the model could not be estimated without this item, for instance because a factor would be left with fewer than two items."))
    footnotes <- c(footnotes, gettext("No confidence intervals are reported for the if-item-dropped coefficients."))
  }
  if (showRest) {
    restData <- .multiDimFreqAnalysisData(dataset, options)
    df[["itemRestCorrelation"]] <- vapply(seq_along(allItems), function(i) {
      rest <- rowSums(restData[, -i, drop = FALSE], na.rm = TRUE)
      cor(restData[, i], rest, use = "pairwise.complete.obs")
    }, numeric(1))
  }

  itemTable$setData(df)

  for (note in footnotes)
    itemTable$addFootnote(note)
  if (length(options[["reverseScaledItems"]]) > 0L)
    itemTable$addFootnote(.addFootnoteReverseScaledItems(options))

  return()
}


.multiDimFreqFitTable <- function(jaspResults, model, options, ready) {

  if (!options[["fitMeasures"]] || !is.null(.getStateContainerMDF(jaspResults)[["fitTable"]]$object))
    return()

  fitTable <- createJaspTable(gettext("Fit Measures of the Factor Model"))
  fitTable$dependOn(options = "fitMeasures")
  fitTable$position <- 3

  fitTable$addColumnInfo(name = "measure", title = gettext("Fit measure"), type = "string")
  fitTable$addColumnInfo(name = "value",   title = gettext("Value"),       type = "number")

  stateContainer <- .getStateContainerMDF(jaspResults)
  stateContainer[["fitTable"]] <- fitTable

  if (!ready)
    return()
  if (!is.null(model[["error"]])) {
    fitTable$setError(model[["error"]])
    return()
  }

  measures <- model[["fit"]][["model"]][["fit.measures"]]
  srmr     <- model[["fit"]][["model"]][["srmr.summary"]]
  if (is.null(measures) || is.null(srmr)) {
    fitTable$setError(gettext("The fit measures could not be computed for this factor model."))
    return()
  }

  chisqNames <- c("chisq", "df", "pvalue", "cfi", "tli",
                  "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "aic", "bic")
  srmrNames  <- c("usrmr", "usrmr.ci.lower", "usrmr.ci.upper", "usrmr.closefit.pvalue")

  labels <- c(gettext("χ²"), gettext("df"), gettext("p"), gettext("CFI"), gettext("TLI"),
              gettext("RMSEA"), gettext("RMSEA 90% CI lower"), gettext("RMSEA 90% CI upper"),
              gettext("RMSEA p-value"), gettext("AIC"), gettext("BIC"),
              gettext("SRMR"), gettext("SRMR 90% CI lower"), gettext("SRMR 90% CI upper"),
              gettext("SRMR p-value"))

  fitTable$setData(data.frame(
    measure = labels,
    value   = c(unname(measures[chisqNames]), unname(unlist(srmr[srmrNames, ]))),
    stringsAsFactors = FALSE
  ))

  for (note in .multiDimFreqDiagnosticsFootnotes(model, options))
    fitTable$addFootnote(note)

  return()
}


.multiDimFreqLoadingsTable <- function(jaspResults, model, options, ready, allItems) {

  if (!options[["standardizedLoadings"]] ||
      !is.null(.getStateContainerMDF(jaspResults)[["loadingsContainer"]]$object))
    return()

  loadingsContainer <- createJaspContainer(gettext("Standardized Factor Loadings"))
  loadingsContainer$dependOn(options = "standardizedLoadings")
  loadingsContainer$position <- 4

  stateContainer <- .getStateContainerMDF(jaspResults)
  stateContainer[["loadingsContainer"]] <- loadingsContainer

  biFactor    <- options[["modelType"]] == "biFactor"
  secondOrder <- options[["modelType"]] == "secondOrder"
  titles      <- .multiDimFreqFactorTitles(options)

  itemTable <- createJaspTable(gettext("Item Loadings"))
  itemTable$position <- 1
  itemTable$addColumnInfo(name = "item", title = gettext("Item"), type = "string")
  if (biFactor)
    itemTable$addColumnInfo(name = "general", title = gettext("General"), type = "number")
  for (i in seq_along(titles))
    itemTable$addColumnInfo(name = paste0("factor", i), title = titles[i], type = "number")
  loadingsContainer[["itemLoadings"]] <- itemTable

  # in the second-order model the general factor loads on the group factors rather than on the
  # items, so those loadings need a table of their own
  if (secondOrder) {
    factorTable <- createJaspTable(gettext("Group Factor Loadings on the General Factor"))
    factorTable$position <- 2
    factorTable$addColumnInfo(name = "factor",  title = gettext("Factor"),  type = "string")
    factorTable$addColumnInfo(name = "loading", title = gettext("Loading"), type = "number")
    loadingsContainer[["factorLoadings"]] <- factorTable
  }

  if (!ready)
    return()
  if (!is.null(model[["error"]])) {
    itemTable$setError(model[["error"]])
    return()
  }

  loadings <- model[["fit"]][["loadings"]]
  specific <- loadings[["specific"]]

  df <- data.frame(item = jaspBase::decodeColNames(allItems), stringsAsFactors = FALSE)
  if (biFactor)
    df[["general"]] <- loadings[["general"]]
  for (i in seq_along(titles)) {
    # an item that does not load on a factor gets an empty cell rather than a zero, which would
    # read as an estimated loading of zero
    column <- specific[, i]
    column[column == 0] <- NA_real_
    df[[paste0("factor", i)]] <- column
  }
  itemTable$setData(df)

  if (secondOrder)
    loadingsContainer[["factorLoadings"]]$setData(data.frame(
      factor  = titles,
      loading = loadings[["general"]],
      stringsAsFactors = FALSE
    ))

  for (note in .multiDimFreqDiagnosticsFootnotes(model, options))
    itemTable$addFootnote(note)

  return()
}
