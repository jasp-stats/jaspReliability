#
# Copyright (C) 2026 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Bayesian kappa for every rater pair: the cells of the pair's agreement table (categories used by the pair) get a
# Dirichlet prior with the same concentration per cell; the posterior is Dirichlet(prior + counts), and kappa is
# computed from every posterior draw of the cell probabilities. Chance agreement uses each rater's own marginals
# (Cohen) or the marginals pooled over both raters (Fleiss, raters treated as interchangeable).

#' @export
raterAgreementBayesianInternal <- function(jaspResults, dataset, options) {

  ready <- length(options[["variables"]]) > 1

  dataset <- .raterAgreementHandleData(dataset, options)

  if (!options[["cohensKappa"]] && !options[["fleissKappa"]])
    .raterAgreementBayesianPlaceholderTable(jaspResults, options, ready)

  for (coefficient in c("cohensKappa", "fleissKappa")) {
    if (!options[[coefficient]])
      next
    .raterAgreementBayesianComputeSamples(jaspResults, dataset, options, ready, coefficient)
    .raterAgreementBayesianTable(jaspResults, options, ready, coefficient)
    if (options[["posteriorPlot"]])
      .raterAgreementBayesianPosteriorPlots(jaspResults, options, ready, coefficient)
  }

  return()
}

.raterAgreementBayesianPlaceholderTable <- function(jaspResults, options, ready) {
  if (!is.null(jaspResults[["placeholder"]]))
    return()

  jaspTable <- createJaspTable(title = gettext("Agreement Coefficient"))
  jaspTable$info <- gettext("Empty placeholder table shown while no agreement coefficient is selected; check a coefficient to obtain results.")
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"),        type = "string")
  jaspTable$addColumnInfo(name = "mean",    title = gettext("Posterior mean"), type = "number")
  jaspTable$addColumnInfo(name = "lower",   title = gettext("Lower"),          type = "number")
  jaspTable$addColumnInfo(name = "upper",   title = gettext("Upper"),          type = "number")
  if (ready)
    jaspTable$addFootnote(gettext("Check one of the coefficients to start the analysis."))
  jaspTable$dependOn(options = c("cohensKappa", "fleissKappa", "variables"))
  jaspResults[["placeholder"]] <- jaspTable
}

.raterAgreementBayesianWeighted <- function(options, coefficient) {
  return(coefficient == "cohensKappa" && options[["cohensKappaType"]] == "weighted")
}

.raterAgreementBayesianName <- function(options, coefficient) {
  if (coefficient == "fleissKappa")
    return(gettext("Fleiss' kappa"))
  if (.raterAgreementBayesianWeighted(options, coefficient))
    return(gettext("Weighted Cohen's kappa"))
  return(gettext("Cohen's kappa"))
}

.raterAgreementBayesianSampleDependencies <- function(coefficient) {
  dependencies <- c("variables", "variables.types", "dataStructure", coefficient,
                    "dirichletPriorConcentration", "samples", "setSeed", "seed")
  if (coefficient == "cohensKappa")
    dependencies <- c(dependencies, "cohensKappaType", "weightType")
  return(dependencies)
}

.raterAgreementBayesianComputeSamples <- function(jaspResults, dataset, options, ready, coefficient) {
  stateName <- paste0(coefficient, "Samples")
  if (!ready || !is.null(jaspResults[[stateName]]))
    return()

  samplesState <- createJaspState()
  samplesState$dependOn(.raterAgreementBayesianSampleDependencies(coefficient))
  jaspResults[[stateName]] <- samplesState

  samplesState$object <- .raterAgreementBayesianFit(dataset, options, coefficient)
  return()
}

.raterAgreementBayesianFit <- function(dataset, options, coefficient) {

  weighted <- .raterAgreementBayesianWeighted(options, coefficient)

  if (ncol(dataset) < 2)
    return(list(error = gettextf("%s requires at least 2 raters/measurements.", .raterAgreementBayesianName(options, coefficient))))

  if (weighted && any(options[["variables.types"]] == "nominal"))
    return(list(error = gettext("Weighted Cohen's kappa requires ordinal variables. Remove nominal variables or change their type.")))

  # all raters share one set of category codes; for weighted kappa the codes are the positions on the common
  # ordinal scale, so a pair that never uses an interior category still gets the right distances
  unionLevels <- .raterAgreementUnionLevels(dataset)
  if (weighted && !unionLevels[["ordered"]])
    return(list(error = .raterAgreementAmbiguousOrderMessage(gettext("Weighted Cohen's kappa"))))

  codes         <- .raterAgreementUnionCodes(dataset, unionLevels[["levels"]])
  possiblePairs <- utils::combn(ncol(dataset), 2)

  pairs <- lapply(seq_len(ncol(possiblePairs)), function(j) {
    pair    <- possiblePairs[, j]
    ratings <- codes[stats::complete.cases(codes[, pair]), pair, drop = FALSE]
    # the seed is set per pair, so a pair's posterior does not depend on which other raters are in the analysis
    jaspBase::.setSeedJASP(options)
    pairFit <- .raterAgreementBayesianPairSamples(
      ratings        = ratings,
      nLevels        = length(unionLevels[["levels"]]),
      weightExponent = if (weighted) ifelse(options[["weightType"]] == "quadratic", 2, 1) else NULL,
      chance         = if (coefficient == "fleissKappa") "pooled" else "cohen",
      prior          = options[["dirichletPriorConcentration"]],
      samples        = options[["samples"]]
    )
    pairFit[["label"]] <- paste(colnames(dataset)[pair], collapse = " - ")
    return(pairFit)
  })

  return(list(
    pairs      = pairs,
    nSubjects  = sum(rowSums(!is.na(dataset)) >= 2), # subjects contribute pairwise
    nRaters    = ncol(dataset),
    anyMissing = anyNA(dataset),
    nLevels    = length(unique(stats::na.omit(as.vector(codes))))
  ))
}

# same algorithm as MCMCpack::rdirichlet, so results can be reproduced with it
.raterAgreementBayesianRdirichlet <- function(n, alpha) {
  l     <- length(alpha)
  draws <- matrix(stats::rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  return(draws / as.vector(draws %*% rep(1, l)))
}

.raterAgreementBayesianPairSamples <- function(ratings, nLevels, weightExponent, chance, prior, samples) {

  n <- nrow(ratings)
  if (n < 3)
    return(list(n = n, failedFor = gettext("fewer than 3 jointly rated subjects/items")))

  categories <- sort(unique(as.vector(ratings)))
  k          <- length(categories)
  if (k < 2)
    return(list(n = n, failedFor = gettext("the ratings do not vary")))

  counts <- table(factor(ratings[, 1], levels = categories), factor(ratings[, 2], levels = categories))

  # agreement weights: identity for unweighted kappa, 1 - |i - j|^e / (nLevels - 1)^e on the common scale otherwise
  weights <- if (is.null(weightExponent)) {
    diag(k)
  } else {
    1 - abs(outer(categories, categories, "-"))^weightExponent / (nLevels - 1)^weightExponent
  }

  # every row of draws holds one posterior draw of the k x k cell probabilities (column-major, as as.vector(counts))
  draws        <- .raterAgreementBayesianRdirichlet(samples, as.vector(counts) + prior)
  cellRow      <- rep(seq_len(k), times = k)
  cellColumn   <- rep(seq_len(k), each = k)
  rowMarginals <- draws %*% outer(cellRow, seq_len(k), "==")
  colMarginals <- draws %*% outer(cellColumn, seq_len(k), "==")

  observed <- as.vector(draws %*% as.vector(weights))
  if (chance == "cohen") {
    expected <- rowSums((rowMarginals %*% weights) * colMarginals)
  } else {
    pooled   <- 0.5 * (rowMarginals + colMarginals)
    expected <- rowSums((pooled %*% weights) * pooled)
  }
  kappa <- (observed - expected) / (1 - expected)

  finite <- is.finite(kappa)
  if (!any(finite))
    return(list(n = n, failedFor = gettext("the coefficient is not estimable")))

  return(list(
    n         = n,
    kappa     = kappa[finite],
    observed  = observed[finite],
    expected  = expected[finite],
    nDropped  = sum(!finite)
  ))
}

.raterAgreementBayesianTable <- function(jaspResults, options, ready, coefficient) {
  if (!is.null(jaspResults[[coefficient]]))
    return()

  name     <- .raterAgreementBayesianName(options, coefficient)
  weighted <- .raterAgreementBayesianWeighted(options, coefficient)

  jaspTable <- createJaspTable(title = name)
  jaspTable$info <- if (coefficient == "fleissKappa") {
    gettext("Fleiss' kappa: chance-corrected agreement with the category proportions pooled over both raters of a pair, for raters that are interchangeable (e.g., different raters per subject). Summarized by the posterior mean and credible interval.")
  } else {
    gettext("Cohen's kappa: chance-corrected agreement between two raters, with each rater's own category proportions. Summarized by the posterior mean and credible interval.")
  }
  jaspTable$position <- if (coefficient == "cohensKappa") 1 else 2
  jaspTable$dependOn(c(.raterAgreementBayesianSampleDependencies(coefficient), "ci", "ciLevel", "observedAndChanceAgreement"))
  jaspResults[[coefficient]] <- jaspTable

  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")

  if (!ready)
    return()

  fit <- jaspResults[[paste0(coefficient, "Samples")]]$object
  if (!is.null(fit[["error"]])) {
    jaspTable$setError(fit[["error"]])
    return()
  }

  valid <- vapply(fit[["pairs"]], function(pair) is.null(pair[["failedFor"]]), logical(1L))
  if (!any(valid)) {
    failedFor <- vapply(fit[["pairs"]], `[[`, character(1L), "failedFor")
    jaspTable$setError(gettextf("%1$s could not be computed for any rater pair: %2$s.", name, paste(unique(failedFor), collapse = "; ")))
    return()
  }

  if (fit[["anyMissing"]])
    jaspTable$addColumnInfo(name = "n", title = gettext("n"), type = "integer")
  jaspTable$addColumnInfo(name = "mean", title = gettext("Posterior mean"), type = "number")
  if (options[["ci"]]) {
    ciPercent <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)
    overtitle <- gettextf("%s%% Credible Interval", ciPercent)
    jaspTable$addColumnInfo(name = "lower", title = gettext("Lower"), type = "number", overtitle = overtitle)
    jaspTable$addColumnInfo(name = "upper", title = gettext("Upper"), type = "number", overtitle = overtitle)
  }
  if (options[["observedAndChanceAgreement"]]) {
    jaspTable$addColumnInfo(name = "observed", title = gettext("Observed agreement"), type = "number")
    jaspTable$addColumnInfo(name = "expected", title = gettext("Chance agreement"),   type = "number")
  }

  # failed pairs keep their row with empty cells; they are explained in a footnote
  rows <- lapply(fit[["pairs"]], function(pair) {
    row <- data.frame(ratings = pair[["label"]], n = pair[["n"]], mean = NA_real_, lower = NA_real_,
                      upper = NA_real_, observed = NA_real_, expected = NA_real_)
    if (!is.null(pair[["failedFor"]]))
      return(row)
    interval         <- coda::HPDinterval(coda::mcmc(pair[["kappa"]]), prob = options[["ciLevel"]])
    row[["mean"]]     <- mean(pair[["kappa"]])
    row[["lower"]]    <- interval[1L, "lower"]
    row[["upper"]]    <- interval[1L, "upper"]
    row[["observed"]] <- mean(pair[["observed"]])
    row[["expected"]] <- mean(pair[["expected"]])
    return(row)
  })
  jaspTable$showSpecifiedColumnsOnly <- TRUE
  jaspTable$setData(do.call(rbind, rows))

  .raterAgreementBayesianFootnotes(jaspTable, fit, options, coefficient, weighted, valid)
  return()
}

.raterAgreementBayesianFootnotes <- function(jaspTable, fit, options, coefficient, weighted, valid) {

  labels <- vapply(fit[["pairs"]], `[[`, character(1L), "label")
  if (any(!valid)) {
    failedFor <- vapply(fit[["pairs"]][!valid], `[[`, character(1L), "failedFor")
    jaspTable$addFootnote(paste0(gettext("Some rater pairs could not be computed: "),
                                 paste0(labels[!valid], " (", failedFor, ")", collapse = "; "), "."),
                          symbol = gettext("Note:"))
  }

  footnote <- gettextf("%1$i subjects/items and %2$i raters/measurements.", fit[["nSubjects"]], fit[["nRaters"]])
  if (fit[["anyMissing"]])
    footnote <- gettextf("%1$s Based on pairwise complete cases.", footnote)
  footnote <- paste(footnote, gettextf(
    "Dirichlet prior with concentration %1$s on each cell of a pair's agreement table; %2$i posterior samples.",
    format(options[["dirichletPriorConcentration"]], drop0trailing = TRUE), options[["samples"]]))
  if (options[["ci"]])
    footnote <- paste(footnote, gettext("Credible intervals are highest posterior density intervals."))
  if (coefficient == "fleissKappa")
    footnote <- paste(footnote, gettext("Chance agreement uses the category proportions pooled over both raters; for two raters this equals Fleiss' kappa and Scott's pi."))
  if (weighted && fit[["nLevels"]] < 3)
    footnote <- paste(footnote, gettext("If there are only 2 levels, weighted kappa is equal to unweighted kappa."))
  jaspTable$addFootnote(footnote)

  nDropped <- vapply(fit[["pairs"]][valid], `[[`, numeric(1L), "nDropped")
  if (any(nDropped > 0))
    jaspTable$addFootnote(gettextf("%i posterior samples with undefined kappa (all probability on a single category) were excluded.", sum(nDropped)),
                          symbol = gettext("Note:"))
  return()
}

.raterAgreementBayesianPosteriorPlots <- function(jaspResults, options, ready, coefficient) {
  plotsName <- paste0(coefficient, "PosteriorPlots")
  if (!ready || !is.null(jaspResults[[plotsName]]))
    return()

  fit <- jaspResults[[paste0(coefficient, "Samples")]]$object
  if (is.null(fit) || !is.null(fit[["error"]]))
    return() # the table shows the error

  name      <- .raterAgreementBayesianName(options, coefficient)
  container <- createJaspContainer(title = gettextf("Posterior Distribution: %s", name))
  container$position <- if (coefficient == "cohensKappa") 3 else 4
  container$dependOn(c(.raterAgreementBayesianSampleDependencies(coefficient), "posteriorPlot", "ci", "ciLevel"))
  jaspResults[[plotsName]] <- container

  for (j in seq_along(fit[["pairs"]])) {
    pair <- fit[["pairs"]][[j]]
    if (!is.null(pair[["failedFor"]]))
      next
    posteriorPlot <- createJaspPlot(title = pair[["label"]], width = 480, height = 320)
    posteriorPlot$position <- j
    container[[paste0("pair", j)]] <- posteriorPlot

    plotObject <- try(.raterAgreementBayesianPosteriorPlot(pair[["kappa"]], options, name))
    if (jaspBase::isTryError(plotObject))
      posteriorPlot$setError(jaspBase::.extractErrorMessage(plotObject))
    else
      posteriorPlot$plotObject <- plotObject
  }
  return()
}

.raterAgreementBayesianPosteriorPlot <- function(kappa, options, name) {

  if (stats::sd(kappa) == 0)
    stop(gettext("The posterior samples do not vary, so their density cannot be plotted."))

  density    <- stats::density(kappa)
  plotData   <- data.frame(x = density[["x"]], y = density[["y"]])
  xBreaks    <- jaspGraphs::getPrettyAxisBreaks(plotData[["x"]])
  yBreaks    <- jaspGraphs::getPrettyAxisBreaks(c(0, plotData[["y"]]))

  plotObject <- ggplot2::ggplot(plotData, ggplot2::aes(x = x, y = y))

  if (options[["ci"]]) {
    interval   <- coda::HPDinterval(coda::mcmc(kappa), prob = options[["ciLevel"]])
    inInterval <- plotData[plotData[["x"]] >= interval[1L, "lower"] & plotData[["x"]] <= interval[1L, "upper"], ]
    plotObject <- plotObject + ggplot2::geom_area(data = inInterval, fill = "grey80")
  }

  plotObject <- plotObject +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_x_continuous(name = name, breaks = xBreaks, limits = range(xBreaks)) +
    ggplot2::scale_y_continuous(name = gettext("Density"), breaks = yBreaks, limits = range(yBreaks)) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()

  return(plotObject)
}
