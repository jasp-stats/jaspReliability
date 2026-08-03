#' @importFrom jaspBase createJaspContainer createJaspHtml createJaspPlot createJaspState createJaspTable progressbarTick startProgressbar
#' @importFrom stats cor cov median sd complete.cases

# Dependencies that invalidate the fitted model (and therefore everything downstream).
.multiDimBaseDependencies <- c(
  "factors", "reverseScaledItems", "modelType",
  "samples", "burnin", "thinning", "chains",
  "setSeed", "seed", "naAction", "samplesSavingDisabled",
  "igShapeManifest", "igScaleManifest", "loadMeanManifest", "loadScaleManifest",
  "igShapeLatent", "igScaleLatent", "loadMeanLatent", "loadScaleLatent",
  "igShapeGFactor", "igScaleGFactor", "latentCorDf"
)

# translated coefficient labels shared by tables and plots
.multiDimOmegaTLabel <- function() gettext("McDonald's ωₜ")
.multiDimOmegaHLabel <- function() gettext("McDonald's ωₕ")

#' @export
reliabilityMultidimensionalBayesian <- function(jaspResults, dataset, options) {

  allItems <- .multiDimGetItems(options)
  ready    <- .multiDimReady(options)

  # option names match the unidimensional analyses, so the shared helpers only need the item list
  optTmp <- options
  optTmp[["variables"]] <- allItems

  if (ready) {
    dataset <- dataset[, allItems, drop = FALSE]
    if (length(options[["reverseScaledItems"]]) > 0L)
      dataset <- .reverseScoreItems(dataset, optTmp)
    .checkErrors(dataset, optTmp, Bayes = TRUE)
  }

  model <- .multiDimComputeModel(jaspResults, dataset, options, ready, allItems)

  .multiDimScaleTable(jaspResults, dataset, model, options, ready)
  .multiDimItemTable(jaspResults, dataset, model, options, ready, allItems)
  .multiDimProbabilityTable(jaspResults, model, options, ready)
  .multiDimFitTable(jaspResults, dataset, model, options, ready)
  .multiDimPosteriorPlot(jaspResults, model, options, ready)
  .multiDimTracePlot(jaspResults, model, options, ready)
  .multiDimPPCPlot(jaspResults, dataset, model, options, ready)

  return()
}


#### Helpers ####

# all unique items across the factors form
.multiDimGetItems <- function(options) {
  facs <- options[["factors"]]
  items <- unlist(lapply(facs, function(f) unlist(f[["indicators"]])), use.names = FALSE)
  return(unique(items))
}

# a multidimensional model needs at least two factors with at least two items each
# an item assigned to more than one factor is a cross-loading; Bayesrel samples these for the
# second-order and correlated models, but stops for the bi-factor model
.multiDimHasCrossLoadings <- function(options) {
  items <- unlist(lapply(options[["factors"]], function(f) unlist(f[["indicators"]])), use.names = FALSE)
  return(anyDuplicated(items) > 0L)
}


.multiDimReady <- function(options) {
  facs <- options[["factors"]]
  if (length(facs) < 2)
    return(FALSE)
  nWithTwo <- sum(vapply(facs, function(f) length(unlist(f[["indicators"]])) >= 2L, logical(1)))
  return(nWithTwo >= 2L)
}

# build the lavaan-style group-factor model syntax used by Bayesrel::bomegas
.multiDimBuildModel <- function(options) {
  facs  <- options[["factors"]]
  lines <- character(0)
  fi    <- 0L
  for (i in seq_along(facs)) {
    inds <- unlist(facs[[i]][["indicators"]])
    if (length(inds) == 0L)
      next
    fi <- fi + 1L
    lines <- c(lines, sprintf("f%d =~ %s", fi, paste(inds, collapse = " + ")))
  }
  return(paste(lines, collapse = "\n"))
}

.getStateContainerMD <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = .multiDimBaseDependencies)
  return(jaspResults[["stateContainer"]])
}

# descriptive statistics computed outside the sampler (scale scores, item-rest correlations) must
# use the same rows as the fit: complete cases under listwise deletion, everything otherwise
.multiDimAnalysisData <- function(dataset, options) {
  if (options[["naAction"]] != "listwise")
    return(dataset)
  return(dataset[complete.cases(dataset), , drop = FALSE])
}

.multiDimPointEst <- function(chains, pointEstimate) {
  .getPointEstFun(pointEstimate)(as.vector(chains))
}

.multiDimRhat <- function(chains) {
  # chains is a [n.chains x n.samples] matrix
  mcmcList <- coda::as.mcmc.list(lapply(seq_len(nrow(chains)), function(i) coda::mcmc(chains[i, ])))
  coda::gelman.diag(mcmcList, autoburnin = FALSE)[["psrf"]][1, 1]
}


#### Compute ####

.multiDimComputeModel <- function(jaspResults, dataset, options, ready, allItems) {

  if (!is.null(.getStateContainerMD(jaspResults)[["modelObj"]]$object))
    return(.getStateContainerMD(jaspResults)[["modelObj"]]$object)

  model <- list(modelType = options[["modelType"]], k = length(allItems), items = allItems)

  if (!ready)
    return(model)

  # caught here rather than by Bayesrel's own stop(), so the message is translated and actionable
  if (options[["modelType"]] == "biFactor" && .multiDimHasCrossLoadings(options)) {
    model[["error"]] <- gettext("The bi-factor model does not support items that load on more than one factor. Assign each item to a single factor, or select the second-order or correlated-factors model.")
    return(model)
  }

  modelSyntax <- .multiDimBuildModel(options)
  missing     <- if (options[["naAction"]] == "listwise") "listwise" else "impute"

  # the g-factor priors differ between models: the second-order/bi-factor models use an inverse-gamma on
  # the g-factor variance (p0 = shape, R0 = scale), whereas the correlated model uses an inverse-Wishart
  # on the latent correlation matrix (p0 = degrees of freedom, R0 = scale matrix). For the correlated
  # case we expose only the df (latentCorDf, clamped to >= number of factors so the prior stays proper)
  # and let Bayesrel build its default scale matrix.
  correlated <- options[["modelType"]] == "correlated"
  nFactors   <- sum(vapply(options[["factors"]], function(f) length(unlist(f[["indicators"]])) > 0L, logical(1)))
  p0 <- if (correlated) max(options[["latentCorDf"]], nFactors) else options[["igShapeGFactor"]]
  R0 <- if (correlated) NA else options[["igScaleGFactor"]]

  jaspBase::.setSeedJASP(options)
  startProgressbar(options[["chains"]] * options[["samples"]],
                   label = gettext("Sampling the posterior distribution of the omega coefficients"))

  fit <- try(Bayesrel::bomegas(
    data       = as.matrix(dataset),
    model      = modelSyntax,
    model.type = options[["modelType"]],
    n.iter     = options[["samples"]],
    n.burnin   = options[["burnin"]],
    n.chains   = options[["chains"]],
    thin       = options[["thinning"]],
    interval   = options[["ciLevel"]],
    missing    = missing,
    a0         = options[["igShapeManifest"]],
    b0         = options[["igScaleManifest"]],
    l0         = options[["loadMeanManifest"]],
    A0         = options[["loadScaleManifest"]],
    c0         = options[["igShapeLatent"]],
    d0         = options[["igScaleLatent"]],
    beta0      = options[["loadMeanLatent"]],
    B0         = options[["loadScaleLatent"]],
    p0         = p0,
    R0         = R0,
    param.out  = FALSE,   # only omega_t/omega_h/implCovs are used; TRUE also crashes bi-factor in Bayesrel <= 0.7.8
    callback   = progressbarTick,
    disableMcmcCheck = TRUE   # MCMC bounds are validated in .checkErrors
  ), silent = TRUE)

  if (inherits(fit, "try-error")) {
    model[["error"]] <- jaspBase::.extractErrorMessage(fit)
    return(model)
  }

  model[["fit"]] <- fit

  if (!options[["samplesSavingDisabled"]]) {
    stateContainer <- .getStateContainerMD(jaspResults)
    stateContainer[["modelObj"]] <- createJaspState(model)
  }

  return(model)
}


#### Tables ####

.multiDimScaleTable <- function(jaspResults, dataset, model, options, ready) {

  if (!is.null(.getStateContainerMD(jaspResults)[["scaleTable"]]$object))
    return()

  scaleTable <- createJaspTable(gettext("Bayesian Scale Reliability Statistics"))
  scaleTable$dependOn(options = c("meanSdScoresMethod", "ciLevel", "rHat", "pointEstimate"))
  scaleTable$position <- 1

  pointEstimate <- if (options[["pointEstimate"]] == "mean") gettext("Posterior mean") else gettext("Posterior median")
  ci  <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)

  scaleTable$addColumnInfo(name = "coefficient", title = gettext("Coefficient"), type = "string")
  scaleTable$addColumnInfo(name = "estimate",    title = pointEstimate,          type = "number")
  scaleTable$addColumnInfo(name = "lower", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI", ci))
  scaleTable$addColumnInfo(name = "upper", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI", ci))
  if (options[["rHat"]])
    scaleTable$addColumnInfo(name = "rHat", title = gettext("R-hat"), type = "number")

  stateContainer <- .getStateContainerMD(jaspResults)
  stateContainer[["scaleTable"]] <- scaleTable

  if (!ready) {
    scaleTable$addFootnote(gettext("Please assign at least two factors with at least two items each."))
    return()
  }
  if (!is.null(model[["error"]])) {
    scaleTable$setError(model[["error"]])
    return()
  }

  fit       <- model[["fit"]]
  ciValue   <- options[["ciLevel"]]
  pointFun  <- options[["pointEstimate"]]
  correlated <- options[["modelType"]] == "correlated"

  rows <- list()
  footnote <- ""

  addCoefRow <- function(label, chains) {
    est  <- .multiDimPointEst(chains, pointFun)
    cred <- coda::HPDinterval(coda::mcmc(as.vector(chains)), prob = ciValue)
    row  <- list(coefficient = label, estimate = est, lower = cred[1], upper = cred[2])
    if (options[["rHat"]])
      row[["rHat"]] <- .multiDimRhat(chains)
    row
  }
  addStatRow <- function(label, value) {
    row <- list(coefficient = label, estimate = value, lower = NA_real_, upper = NA_real_)
    if (options[["rHat"]])
      row[["rHat"]] <- NA_real_
    row
  }

  # always display all coefficients; McDonald's omega_h exists only for the non-correlated models
  rows[[length(rows) + 1L]] <- addCoefRow(.multiDimOmegaTLabel(), fit[["omega_t"]][["chains"]])

  if (correlated) {
    footnote <- gettext("McDonald's ωₕ is shown only for the second-order and bi-factor models, not the correlated-factors model.")
  } else {
    rows[[length(rows) + 1L]] <- addCoefRow(.multiDimOmegaHLabel(), fit[["omega_h"]][["chains"]])
  }

  pairwise <- options[["naAction"]] != "listwise"
  cc <- cor(dataset, use = if (pairwise) "pairwise.complete.obs" else "complete.obs")
  rows[[length(rows) + 1L]] <- addStatRow(gettext("Average interitem correlation"), mean(cc[lower.tri(cc)]))

  scoreData <- .multiDimAnalysisData(dataset, options)
  scores    <- if (options[["meanSdScoresMethod"]] == "sumScores")
    rowSums(scoreData, na.rm = TRUE) else rowMeans(scoreData, na.rm = TRUE)
  rows[[length(rows) + 1L]] <- addStatRow(gettext("Mean"), mean(scores))
  rows[[length(rows) + 1L]] <- addStatRow(gettext("SD"), sd(scores))

  if (length(rows) > 0L)
    scaleTable$setData(do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE)))

  if (footnote != "")
    scaleTable$addFootnote(footnote)

  return()
}


.multiDimItemTable <- function(jaspResults, dataset, model, options, ready, allItems) {

  showOmegaT <- options[["itemDeletedOmegaT"]]
  showOmegaH <- options[["itemDeletedOmegaH"]] && options[["modelType"]] != "correlated"
  showRest   <- options[["itemRestCorrelation"]]
  if (!(showOmegaT || showOmegaH || showRest) ||
      !is.null(.getStateContainerMD(jaspResults)[["itemTable"]]$object))
    return()

  itemTable <- createJaspTable(gettext("Bayesian Individual Item Reliability Statistics"))
  itemTable$dependOn(options = c("itemDeletedOmegaT", "itemDeletedOmegaH", "itemRestCorrelation", "pointEstimate",
                                 "ciLevel"))
  itemTable$position <- 2

  pointEstimate <- if (options[["pointEstimate"]] == "mean") gettext("Posterior mean") else gettext("Posterior median")
  cred          <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)

  # point estimate + credible interval per omega, mirroring the unidimensional Bayesian item table;
  # the item-rest correlation is a plug-in statistic, not a posterior, so it gets no interval
  addOmegaColumns <- function(prefix, overtitle) {
    itemTable$addColumnInfo(name = prefix,                 title = pointEstimate, type = "number", overtitle = overtitle)
    itemTable$addColumnInfo(name = paste0(prefix, "Lower"),
                            title = gettextf("Lower %s%% CI", cred), type = "number", overtitle = overtitle)
    itemTable$addColumnInfo(name = paste0(prefix, "Upper"),
                            title = gettextf("Upper %s%% CI", cred), type = "number", overtitle = overtitle)
  }

  itemTable$addColumnInfo(name = "item", title = gettext("Item"), type = "string")
  if (showOmegaT)
    addOmegaColumns("omegaT", gettext("McDonald's ωₜ (if item dropped)"))
  if (showOmegaH)
    addOmegaColumns("omegaH", gettext("McDonald's ωₕ (if item dropped)"))
  if (showRest)
    itemTable$addColumnInfo(name = "itemRestCorrelation", title = gettext("Item-rest correlation"), type = "number")

  stateContainer <- .getStateContainerMD(jaspResults)
  stateContainer[["itemTable"]] <- itemTable

  if (!ready)
    return()
  if (!is.null(model[["error"]])) {
    itemTable$setError(model[["error"]])
    return()
  }

  df <- data.frame(item = jaspBase::decodeColNames(allItems), stringsAsFactors = FALSE)

  footnote <- ""
  if (showOmegaT || showOmegaH) {
    del      <- .multiDimItemDeletedOmega(jaspResults, dataset, options, ready, allItems)
    pointFun <- .getPointEstFun(options[["pointEstimate"]])
    ciValue  <- options[["ciLevel"]]

    # summarise the cached chains at display time, so switching point estimate or CI level costs no refit
    summarizeChains <- function(chains) {
      if (is.null(chains))
        return(c(NA_real_, NA_real_, NA_real_))
      cred <- coda::HPDinterval(coda::mcmc(chains), prob = ciValue)
      c(pointFun(chains), cred[1], cred[2])
    }
    addOmegaData <- function(df, prefix, chainList) {
      summaries <- vapply(chainList, summarizeChains, numeric(3))
      df[[prefix]]                  <- summaries[1, ]
      df[[paste0(prefix, "Lower")]] <- summaries[2, ]
      df[[paste0(prefix, "Upper")]] <- summaries[3, ]
      df
    }

    if (showOmegaT) df <- addOmegaData(df, "omegaT", del[["omtChains"]])
    if (showOmegaH) df <- addOmegaData(df, "omegaH", del[["omhChains"]])
    if (anyNA(c(if (showOmegaT) df[["omegaT"]], if (showOmegaH) df[["omegaH"]])))
      footnote <- gettext("Empty cells indicate that dropping the item would leave a factor with fewer than two items.")
  }
  if (showRest) {
    restData <- .multiDimAnalysisData(dataset, options)
    df$itemRestCorrelation <- vapply(seq_along(allItems), function(i) {
      rest <- rowSums(restData[, -i, drop = FALSE], na.rm = TRUE)
      cor(restData[, i], rest, use = "pairwise.complete.obs")
    }, numeric(1))
  }

  itemTable$setData(df)

  if (footnote != "")
    itemTable$addFootnote(footnote)
  if (length(options[["reverseScaledItems"]]) > 0L)
    itemTable$addFootnote(.addFootnoteReverseScaledItems(options))

  return()
}


# omega-if-item-deleted for multidimensional models: Bayesrel has no built-in support, so refit the
# full Gibbs sampler once per item with that item removed from its factor(s). Items whose removal would
# leave a factor with < 2 indicators are skipped (NA), since the reduced model would be under-identified.
.multiDimItemDeletedOmega <- function(jaspResults, dataset, options, ready, allItems) {

  sc <- .getStateContainerMD(jaspResults)
  if (!is.null(sc[["itemDeletedObj"]]$object))
    return(sc[["itemDeletedObj"]]$object)

  k          <- length(allItems)
  omtChains  <- vector("list", k)
  omhChains  <- vector("list", k)
  correlated <- options[["modelType"]] == "correlated"
  facs       <- options[["factors"]]
  missing    <- if (options[["naAction"]] == "listwise") "listwise" else "impute"

  startProgressbar(k, label = gettext("Computing item statistics (refitting the model once per item)"))
  for (i in seq_len(k)) {
    item    <- allItems[i]
    reduced <- lapply(facs, function(f) setdiff(unlist(f[["indicators"]]), item))
    sizes   <- vapply(reduced, length, integer(1))
    if (any(sizes[sizes > 0L] < 2L)) {       # would under-identify a factor
      progressbarTick()
      next
    }
    lines <- character(0); fi <- 0L
    for (f in reduced) {
      if (length(f) == 0L) next
      fi    <- fi + 1L
      lines <- c(lines, sprintf("f%d =~ %s", fi, paste(f, collapse = " + ")))
    }
    p0     <- if (correlated) max(options[["latentCorDf"]], fi) else options[["igShapeGFactor"]]
    R0     <- if (correlated) NA else options[["igScaleGFactor"]]
    datRed <- as.matrix(dataset[, setdiff(allItems, item), drop = FALSE])

    jaspBase::.setSeedJASP(options)
    fitRed <- try(Bayesrel::bomegas(
      data     = datRed, model = paste(lines, collapse = "\n"), model.type = options[["modelType"]],
      n.iter   = options[["samples"]], n.burnin = options[["burnin"]], n.chains = options[["chains"]],
      thin     = options[["thinning"]], interval = options[["ciLevel"]], missing = missing,
      a0       = options[["igShapeManifest"]], b0 = options[["igScaleManifest"]], l0 = options[["loadMeanManifest"]],
      A0       = options[["loadScaleManifest"]],
      c0       = options[["igShapeLatent"]], d0 = options[["igScaleLatent"]], beta0 = options[["loadMeanLatent"]],
      B0       = options[["loadScaleLatent"]],
      p0       = p0, R0 = R0, param.out = FALSE, callback = function(){}, disableMcmcCheck = TRUE
    ), silent = TRUE)

    if (!inherits(fitRed, "try-error")) {
      omtChains[[i]] <- as.vector(fitRed[["omega_t"]][["chains"]])
      if (!correlated)
        omhChains[[i]] <- as.vector(fitRed[["omega_h"]][["chains"]])
    }
    progressbarTick()
  }

  # the chains are cached rather than a point estimate, so the point-estimate type and CI level are
  # display-time choices that do not trigger k additional model refits
  out <- list(omtChains = omtChains, omhChains = omhChains)
  sc[["itemDeletedObj"]] <- createJaspState(out, dependencies = c("itemDeletedOmegaT", "itemDeletedOmegaH"))
  return(out)
}


.multiDimProbabilityTable <- function(jaspResults, model, options, ready) {

  if (!options[["probabilityTable"]] || !is.null(.getStateContainerMD(jaspResults)[["probabilityTable"]]$object))
    return()

  if (options[["probabilityTableLowerBound"]] > options[["probabilityTableUpperBound"]]) {
    low  <- options[["probabilityTableUpperBound"]]
    high <- options[["probabilityTableLowerBound"]]
    footnote <- gettext("The bounds you entered have been rearranged in increasing order to provide meaningful results.")
  } else {
    low  <- options[["probabilityTableLowerBound"]]
    high <- options[["probabilityTableUpperBound"]]
    footnote <- ""
  }

  probabilityTable <- createJaspTable(
    gettextf("Probability that Reliability Coefficient is Larger than %1$.2f and Smaller than %2$.2f", low, high))
  probabilityTable$dependOn(options = c("probabilityTable", "probabilityTableLowerBound", "probabilityTableUpperBound"))
  probabilityTable$position <- 3
  probabilityTable$addColumnInfo(name = "coefficient", title = gettext("Coefficient"), type = "string")
  probabilityTable$addColumnInfo(name = "posterior",   title = gettext("Posterior"),   type = "number",
                                 overtitle = gettext("Probability"))

  stateContainer <- .getStateContainerMD(jaspResults)
  stateContainer[["probabilityTable"]] <- probabilityTable

  if (!ready)
    return()
  if (!is.null(model[["error"]])) {
    probabilityTable$setError(model[["error"]])
    return()
  }

  fit  <- model[["fit"]]
  rows <- list()
  probInRange <- function(chains) {
    samp <- as.vector(chains)
    mean(samp > low) - mean(samp > high)
  }
  rows[[length(rows) + 1L]] <- list(coefficient = .multiDimOmegaTLabel(), posterior = probInRange(fit[["omega_t"]][["chains"]]))
  if (options[["modelType"]] != "correlated")
    rows[[length(rows) + 1L]] <- list(coefficient = .multiDimOmegaHLabel(), posterior = probInRange(fit[["omega_h"]][["chains"]]))

  if (length(rows) > 0L)
    probabilityTable$setData(do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE)))

  if (footnote != "")
    probabilityTable$addFootnote(footnote)

  return()
}


.multiDimFitTable <- function(jaspResults, dataset, model, options, ready) {

  if (!options[["fitMeasures"]] || !is.null(.getStateContainerMD(jaspResults)[["fitTable"]]$object))
    return()

  fitTable <- createJaspTable(gettext("Fit Measures of the Factor Model"))
  fitTable$dependOn(options = c("fitMeasures", "fitMeasuresCiLevel", "fitMeasuresCutoffRmsea"))
  fitTable$position <- 4

  fitTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")
  fitTable$addColumnInfo(name = "lr",    title = "B-LR",    type = "number")
  fitTable$addColumnInfo(name = "srmr",  title = "B-SRMR",  type = "number")
  fitTable$addColumnInfo(name = "rmsea", title = "B-RMSEA", type = "number")

  stateContainer <- .getStateContainerMD(jaspResults)
  stateContainer[["fitTable"]] <- fitTable

  if (!ready)
    return()
  if (!is.null(model[["error"]])) {
    fitTable$setError(model[["error"]])
    return()
  }

  ciFit  <- options[["fitMeasuresCiLevel"]]
  fitOut <- try(Bayesrel::multiFit(model[["fit"]], data = as.matrix(dataset), ppc = FALSE,
                                   cutoff = options[["fitMeasuresCutoffRmsea"]], ci = ciFit), silent = TRUE)
  if (inherits(fitOut, "try-error")) {
    fitTable$setError(jaspBase::.extractErrorMessage(fitOut))
    return()
  }

  srmrCi <- as.numeric(coda::HPDinterval(coda::mcmc(fitOut[["srmr_samp"]]), prob = ciFit))
  cred   <- gettextf("%s%% CI", format(100 * ciFit, digits = 3, drop0trailing = TRUE))

  df <- data.frame(
    estimate = c(gettext("Point estimate"), gettextf("%s lower bound", cred),
                 gettextf("%s upper bound", cred), gettext("Relative to cutoff")),
    lr    = c(fitOut[["LR"]],            NA_real_,        NA_real_,         NA_real_),
    srmr  = c(fitOut[["srmr_pointEst"]], srmrCi[1],       srmrCi[2],        NA_real_),
    rmsea = c(fitOut[["rmsea_pointEst"]], fitOut[["rmsea_ci"]][1], fitOut[["rmsea_ci"]][2], fitOut[["p_rmsea"]]),
    stringsAsFactors = FALSE
  )
  fitTable$setData(df)
  fitTable$addFootnote(gettextf("'Relative to cutoff' denotes the posterior probability that the B-RMSEA is smaller than the cutoff of %.2f.",
                                options[["fitMeasuresCutoffRmsea"]]))

  return()
}


#### Plots ####

.multiDimPosteriorPlot <- function(jaspResults, model, options, ready) {

  if (!options[["posteriorPlot"]] || !is.null(.getStateContainerMD(jaspResults)[["posteriorPlots"]]$object))
    return()

  plotContainer <- createJaspContainer(gettext("Posterior Plots"))
  plotContainer$dependOn(options = c("posteriorPlot", "posteriorPlotFixedRange", "posteriorPlotPriorDisplayed", "posteriorPlotShaded",
                                     "probabilityTable", "probabilityTableLowerBound", "probabilityTableUpperBound",
                                     "ciLevel"))
  plotContainer$position <- 6

  if (!ready || !is.null(model[["error"]])) {
    stateContainer <- .getStateContainerMD(jaspResults)
    stateContainer[["posteriorPlots"]] <- plotContainer
    return()
  }

  fit     <- model[["fit"]]
  ciValue <- options[["ciLevel"]]

  if (options[["posteriorPlotShaded"]] && options[["probabilityTable"]]) {
    low  <- min(options[["probabilityTableLowerBound"]], options[["probabilityTableUpperBound"]])
    high <- max(options[["probabilityTableLowerBound"]], options[["probabilityTableUpperBound"]])
    shade <- c(low, high)
  } else {
    shade <- NULL
  }

  coefs <- list()
  coefs[["omegaT"]] <- list(chains = fit[["omega_t"]][["chains"]], label = .multiDimOmegaTLabel())
  if (options[["modelType"]] != "correlated")
    coefs[["omegaH"]] <- list(chains = fit[["omega_h"]][["chains"]], label = .multiDimOmegaHLabel())

  # prior samples of the omegas, drawn with the same prior parameterization as the Gibbs sampler
  priorSamples <- NULL
  if (options[["posteriorPlotPriorDisplayed"]]) {
    priorFun <- switch(fit[["model.type"]],
      "second-order" = Bayesrel:::omegasSecoPrior,
      "bi-factor"    = Bayesrel:::omegasBifPrior,
      "correlated"   = Bayesrel:::omegasCorrPrior)
    priorSamples <- try(priorFun(fit, nsamp = 2e3), silent = TRUE)
    if (inherits(priorSamples, "try-error"))
      priorSamples <- NULL
  }
  priorKeys <- c(omegaT = "omt_prior", omegaH = "omh_prior")

  for (nm in names(coefs)) {
    cred <- coda::HPDinterval(coda::mcmc(as.vector(coefs[[nm]][["chains"]])), prob = ciValue)
    priorDens <- if (!is.null(priorSamples))
      stats::density(priorSamples[[priorKeys[[nm]]]], from = 0, to = 1, n = 512)
    else
      NULL
    p <- .makeSinglePosteriorPlot(list(samp = coefs[[nm]][["chains"]]), cred, coefs[[nm]][["label"]],
                                  options[["posteriorPlotFixedRange"]], shade,
                                  priorTrue = !is.null(priorDens), priorSample = priorDens)
    plotObj <- createJaspPlot(plot = p, title = coefs[[nm]][["label"]])
    plotContainer[[nm]] <- plotObj
  }

  stateContainer <- .getStateContainerMD(jaspResults)
  stateContainer[["posteriorPlots"]] <- plotContainer

  return()
}


.multiDimTracePlot <- function(jaspResults, model, options, ready) {

  if (!options[["tracePlot"]] || !is.null(.getStateContainerMD(jaspResults)[["tracePlots"]]$object))
    return()

  plotContainer <- createJaspContainer(gettext("Convergence Traceplot"))
  plotContainer$dependOn(options = c("tracePlot"))
  plotContainer$position <- 7

  if (!ready || !is.null(model[["error"]])) {
    stateContainer <- .getStateContainerMD(jaspResults)
    stateContainer[["tracePlots"]] <- plotContainer
    return()
  }

  fit   <- model[["fit"]]
  coefs <- list()
  coefs[["omegaT"]] <- list(chains = fit[["omega_t"]][["chains"]], label = .multiDimOmegaTLabel())
  if (options[["modelType"]] != "correlated")
    coefs[["omegaH"]] <- list(chains = fit[["omega_h"]][["chains"]], label = .multiDimOmegaHLabel())

  for (nm in names(coefs)) {
    p <- .makeTracePlot(list(samp = coefs[[nm]][["chains"]]), coefs[[nm]][["label"]])
    plotObj <- createJaspPlot(plot = p, title = coefs[[nm]][["label"]], width = 400)
    plotContainer[[nm]] <- plotObj
  }

  stateContainer <- .getStateContainerMD(jaspResults)
  stateContainer[["tracePlots"]] <- plotContainer

  return()
}


.multiDimPPCPlot <- function(jaspResults, dataset, model, options, ready) {

  if (!options[["posteriorPredictiveCheck"]] || !is.null(.getStateContainerMD(jaspResults)[["ppcPlot"]]$object))
    return()

  plot <- createJaspPlot(title = gettext("Posterior Predictive Check"), width = 350)
  plot$dependOn(options = "posteriorPredictiveCheck")
  plot$position <- 8

  stateContainer <- .getStateContainerMD(jaspResults)
  stateContainer[["ppcPlot"]] <- plot

  if (!ready || !is.null(model[["error"]]))
    return()

  implCovs <- model[["fit"]][["implCovs"]]
  pairwise <- options[["naAction"]] != "listwise"
  cdat     <- cov(dataset, use = if (pairwise) "pairwise.complete.obs" else "complete.obs")
  k        <- ncol(cdat)
  n        <- nrow(dataset)
  nsamp    <- dim(implCovs)[1]

  eeImpl <- matrix(0, nsamp, k)
  for (i in seq_len(nsamp)) {
    dtmp <- MASS::mvrnorm(n, rep(0, k), implCovs[i, , ])
    eeImpl[i, ] <- eigen(cov(dtmp), only.values = TRUE)$values
  }

  eframe <- data.frame(number = seq_len(k), eigen_value = eigen(cdat, only.values = TRUE)$values)
  eframe$eigen_sim_low <- apply(eeImpl, 2, quantile, prob = .025)
  eframe$eigen_sim_up  <- apply(eeImpl, 2, quantile, prob = .975)

  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(eframe$eigen_sim_up)))

  g <- ggplot2::ggplot(eframe, mapping = ggplot2::aes(x = number, y = eigen_value)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = eigen_sim_low, ymax = eigen_sim_up), color = "grey55",
                           width = 0.2, linewidth = 1) +
    ggplot2::geom_point(size = 2.25) +
    ggplot2::scale_y_continuous(name = gettext("Eigenvalue"), breaks = yBreaks, limits = range(yBreaks)) +
    ggplot2::scale_x_continuous(name = gettext("Eigenvalue No."), breaks = seq_len(k),
                                expand = ggplot2::expansion(mult = c(.1, .1))) +
    jaspGraphs::geom_rangeframe(sides = "bl") +
    jaspGraphs::themeJaspRaw()

  plot$plotObject <- g

  return()
}
