# Property-based (invariance) tests for Rater Agreement.
#
# Every agreement coefficient must be unchanged by operations that are mathematically
# no-ops on the ratings:
#
#   * permuting the raters
#   * permuting the subjects/items
#   * relabelling the categories (order-preserving for ordinal, any bijection for nominal)
#   * expressing the same data with raters in rows instead of in columns
#   * a strictly increasing per-rater transform (Kendall's W ranks within raters)
#   * an affine transform of the scores (Krippendorff's interval alpha)
#
# These properties -- rather than hand-picked examples -- are what the historical bugs of
# this analysis violated: a cross-rater level merge that depended on rater order, an
# encoder that recoded numeric columns when any column was discrete, and a transpose that
# dropped factor levels. Several fixtures therefore let the raters declare only SUBSETS of
# the scale (no single rater observes every category), because that is the case in which
# the cross-rater level merge becomes observable.

raOptions <- function(variables, types, ...) {
  opts                 <- analysisOptions("raterAgreement")
  opts$dataStructure   <- "ratersInColumns"
  opts$ci              <- FALSE
  opts$variables       <- variables
  opts$variables.types <- types
  extra <- list(...)
  for (name in names(extra)) opts[[name]] <- extra[[name]]
  return(opts)
}

raRun <- function(dataset, opts) runAnalysis("raterAgreement", dataset, opts)

# a coefficient that never got computed would make every comparison vacuously true
raCoefficient <- function(results, table, field) {
  rows <- results[["results"]][[table]][["data"]]
  expect_gt(length(rows), 0)
  return(rows[[1]][[field]])
}

raAllCoefficients <- function(results, kendall = TRUE) {
  out <- list(
    cohen  = raCoefficient(results, "cohensKappa", "cKappa"),
    fleiss = raCoefficient(results, "fleissKappa", "fKappa"),
    alpha  = raCoefficient(results, "krippendorffsAlpha", "kAlpha")
  )
  if (kendall) out$kendall <- raCoefficient(results, "kendallW", "W")
  return(out)
}

raAllOptions <- function(variables, types, kendall = TRUE, weighted = TRUE, method = "ordinal") {
  raOptions(variables, types,
            cohensKappa              = TRUE,
            cohensKappaType          = if (weighted) "weighted" else "unweighted",
            fleissKappa              = TRUE,
            krippendorffsAlpha       = TRUE,
            krippendorffsAlphaMethod = method,
            kendallW                 = kendall)
}

raOrdinalLevels <- c("low", "medium", "high")

raOrdinalData <- function(seed, nSubjects = 20, nRaters = 4, pMissing = 0) {
  set.seed(seed)
  dataset <- as.data.frame(lapply(seq_len(nRaters), function(j)
    factor(sample(raOrdinalLevels, nSubjects, replace = TRUE), levels = raOrdinalLevels, ordered = TRUE)))
  names(dataset) <- paste0("r", seq_len(nRaters))
  if (pMissing > 0)
    for (j in seq_len(nRaters)) {
      missingRows <- which(runif(nSubjects) < pMissing)
      if (length(missingRows) > 0) dataset[missingRows, j] <- NA
    }
  return(dataset)
}

raScaleData <- function(seed, nSubjects = 20, nRaters = 4) {
  set.seed(seed)
  dataset <- as.data.frame(matrix(round(rnorm(nSubjects * nRaters), 4), nSubjects, nRaters))
  names(dataset) <- paste0("r", seq_len(nRaters))
  return(dataset)
}

# low < medium < high, but no rater declares the whole scale; the pairwise constraints
# (low < high, medium < high, low < medium) still determine one total order
raPartialLevelData <- function(columns) {
  dataset <- data.frame(
    r1 = factor(c("low","high","low","high","low","high","high","low"),
                levels = c("low", "high"),   ordered = TRUE),
    r2 = factor(c("medium","high","medium","high","high","medium","high","medium"),
                levels = c("medium", "high"), ordered = TRUE),
    r3 = factor(c("low","medium","low","medium","low","medium","medium","low"),
                levels = c("low", "medium"),  ordered = TRUE)
  )
  return(dataset[, columns, drop = FALSE])
}

raExpectSame <- function(actual, expected, label) {
  expect_equal(actual, expected, tolerance = 1e-10, label = label)
}


# ==== Permuting the raters ====
test_that("Coefficients are invariant to the order of the raters", {
  permutation <- c(3, 1, 4, 2)
  for (seed in 1:3) {
    dataset <- raOrdinalData(seed)
    original <- raRun(dataset, raAllOptions(names(dataset), rep("ordinal", 4)))
    permuted <- raRun(dataset[, permutation],
                      raAllOptions(names(dataset)[permutation], rep("ordinal", 4)))

    a <- raAllCoefficients(original)
    b <- raAllCoefficients(permuted)
    raExpectSame(b$cohen,   a$cohen,   "weighted Cohen's kappa (average)")
    raExpectSame(b$fleiss,  a$fleiss,  "Fleiss' kappa")
    raExpectSame(b$alpha,   a$alpha,   "ordinal Krippendorff's alpha")
    raExpectSame(b$kendall, a$kendall, "Kendall's W")
  }
})

test_that("Coefficients are invariant to rater order when raters declare only part of the scale", {
  # the case the cross-rater level merge exists for: a greedy merge makes these differ
  opts <- function(cols) raOptions(cols, rep("ordinal", 3),
                                   cohensKappa = TRUE, cohensKappaType = "weighted",
                                   krippendorffsAlpha = TRUE, krippendorffsAlphaMethod = "ordinal")
  reference <- raRun(raPartialLevelData(c("r1", "r2", "r3")), opts(c("r1", "r2", "r3")))
  for (columns in list(c("r2", "r3", "r1"), c("r3", "r1", "r2"), c("r3", "r2", "r1"))) {
    permuted <- raRun(raPartialLevelData(columns), opts(columns))
    raExpectSame(raCoefficient(permuted, "cohensKappa", "cKappa"),
                 raCoefficient(reference, "cohensKappa", "cKappa"),
                 paste("weighted Cohen's kappa for", paste(columns, collapse = ",")))
    raExpectSame(raCoefficient(permuted, "krippendorffsAlpha", "kAlpha"),
                 raCoefficient(reference, "krippendorffsAlpha", "kAlpha"),
                 paste("ordinal alpha for", paste(columns, collapse = ",")))
  }
})


# ==== Permuting the subjects/items ====
test_that("Coefficients are invariant to the order of the subjects, also with missing data", {
  for (seed in 1:3) {
    dataset <- raOrdinalData(seed, pMissing = 0.15)
    set.seed(seed + 100)
    shuffled <- dataset[sample(nrow(dataset)), ]

    opts <- raAllOptions(names(dataset), rep("ordinal", 4))
    a <- raAllCoefficients(raRun(dataset, opts))
    b <- raAllCoefficients(raRun(shuffled, opts))
    raExpectSame(b$cohen,   a$cohen,   "weighted Cohen's kappa (average)")
    raExpectSame(b$fleiss,  a$fleiss,  "Fleiss' kappa")
    raExpectSame(b$alpha,   a$alpha,   "ordinal Krippendorff's alpha")
    raExpectSame(b$kendall, a$kendall, "Kendall's W")
  }
})

test_that("Kendall's W is invariant to subject order when raters mix ordinal and scale types", {
  # a shared cross-rater encoding would rank the numeric columns by first appearance here
  dataset <- data.frame(
    ord = factor(c("low","medium","high","medium","low","high","high","low","medium"),
                 levels = raOrdinalLevels, ordered = TRUE),
    sc1 = c(2.5, 7.1, 9.9, 4.0, 1.2, 8.8, 9.1, 0.5, 5.5),
    sc2 = c(1.0, 6.0, 9.0, 5.0, 2.0, 8.0, 9.5, 0.1, 4.4)
  )
  opts <- raOptions(names(dataset), c("ordinal", "scale", "scale"), kendallW = TRUE)
  reference <- raCoefficient(raRun(dataset, opts), "kendallW", "W")
  for (seed in 1:3) {
    set.seed(seed)
    shuffled <- raCoefficient(raRun(dataset[sample(nrow(dataset)), ], opts), "kendallW", "W")
    raExpectSame(shuffled, reference, paste("Kendall's W, subject permutation", seed))
  }
})


# ==== Relabelling the categories ====
test_that("Ordinal coefficients are invariant to order-preserving relabelling", {
  for (seed in 1:2) {
    dataset  <- raOrdinalData(seed)
    relabels <- c(low = "1", medium = "5", high = "9") # order preserved, spacing irrelevant
    relabelled <- as.data.frame(lapply(dataset, function(column)
      factor(relabels[as.character(column)], levels = unname(relabels), ordered = TRUE)))
    names(relabelled) <- names(dataset)

    opts <- raOptions(names(dataset), rep("ordinal", 4),
                      cohensKappa = TRUE, cohensKappaType = "weighted",
                      krippendorffsAlpha = TRUE, krippendorffsAlphaMethod = "ordinal",
                      kendallW = TRUE)
    a <- raRun(dataset, opts)
    b <- raRun(relabelled, opts)
    raExpectSame(raCoefficient(b, "cohensKappa", "cKappa"),
                 raCoefficient(a, "cohensKappa", "cKappa"), "weighted Cohen's kappa")
    raExpectSame(raCoefficient(b, "krippendorffsAlpha", "kAlpha"),
                 raCoefficient(a, "krippendorffsAlpha", "kAlpha"), "ordinal alpha")
    raExpectSame(raCoefficient(b, "kendallW", "W"),
                 raCoefficient(a, "kendallW", "W"), "Kendall's W")
  }
})

test_that("Nominal coefficients are invariant to any bijective relabelling of the categories", {
  for (seed in 1:2) {
    dataset <- as.data.frame(lapply(raOrdinalData(seed), function(column)
      factor(as.character(column), levels = raOrdinalLevels)))
    names(dataset) <- paste0("r", seq_len(ncol(dataset)))
    relabels <- c(low = "zeta", medium = "alpha", high = "mu") # scrambles alphabetical order too
    relabelled <- as.data.frame(lapply(dataset, function(column)
      factor(relabels[as.character(column)], levels = unname(relabels))))
    names(relabelled) <- names(dataset)

    opts <- raOptions(names(dataset), rep("nominal", 4),
                      cohensKappa = TRUE, fleissKappa = TRUE,
                      krippendorffsAlpha = TRUE, krippendorffsAlphaMethod = "nominal")
    a <- raRun(dataset, opts)
    b <- raRun(relabelled, opts)
    raExpectSame(raCoefficient(b, "cohensKappa", "cKappa"),
                 raCoefficient(a, "cohensKappa", "cKappa"), "unweighted Cohen's kappa")
    raExpectSame(raCoefficient(b, "fleissKappa", "fKappa"),
                 raCoefficient(a, "fleissKappa", "fKappa"), "Fleiss' kappa")
    raExpectSame(raCoefficient(b, "krippendorffsAlpha", "kAlpha"),
                 raCoefficient(a, "krippendorffsAlpha", "kAlpha"), "nominal alpha")
  }
})


# ==== Raters in rows vs raters in columns ====
test_that("Raters in rows reproduces raters in columns for labelled ordered factors", {
  for (seed in 1:2) {
    dataset <- raOrdinalData(seed, nSubjects = 10, nRaters = 3)
    transposed <- as.data.frame(lapply(seq_len(nrow(dataset)), function(i)
      factor(as.character(unlist(dataset[i, ])), levels = raOrdinalLevels, ordered = TRUE)))
    names(transposed) <- paste0("s", seq_len(nrow(dataset)))

    columnOpts <- raAllOptions(names(dataset), rep("ordinal", 3))
    rowOpts    <- raAllOptions(names(transposed), rep("ordinal", ncol(transposed)))
    rowOpts$dataStructure <- "ratersInRows"

    a <- raAllCoefficients(raRun(dataset, columnOpts))
    b <- raAllCoefficients(raRun(transposed, rowOpts))
    raExpectSame(b$cohen,   a$cohen,   "weighted Cohen's kappa (average)")
    raExpectSame(b$fleiss,  a$fleiss,  "Fleiss' kappa")
    raExpectSame(b$alpha,   a$alpha,   "ordinal Krippendorff's alpha")
    raExpectSame(b$kendall, a$kendall, "Kendall's W")
  }
})


# ==== Transforms that preserve the quantity each coefficient measures ====
test_that("Kendall's W is invariant to strictly increasing per-rater transforms", {
  for (seed in 1:3) {
    dataset <- raScaleData(seed)
    transformed <- dataset
    transformed$r1 <- 3 * dataset$r1 + 7   # affine, increasing
    transformed$r2 <- exp(dataset$r2)      # increasing, nonlinear
    transformed$r3 <- dataset$r3^3         # increasing (odd power)
    transformed$r4 <- dataset$r4 - 100     # shift

    opts <- raOptions(names(dataset), rep("scale", 4), kendallW = TRUE)
    raExpectSame(raCoefficient(raRun(transformed, opts), "kendallW", "W"),
                 raCoefficient(raRun(dataset, opts), "kendallW", "W"),
                 paste("Kendall's W under monotone transforms, seed", seed))
  }
})

test_that("Krippendorff's alpha respects the invariances of its measurement level", {
  dataset <- raScaleData(4)
  intervalOpts <- raOptions(names(dataset), rep("scale", 4),
                            krippendorffsAlpha = TRUE, krippendorffsAlphaMethod = "interval")
  reference <- raCoefficient(raRun(dataset, intervalOpts), "krippendorffsAlpha", "kAlpha")

  # interval data are invariant to any affine rescaling
  rescaled <- as.data.frame(lapply(dataset, function(column) 2.5 * column + 13))
  raExpectSame(raCoefficient(raRun(rescaled, intervalOpts), "krippendorffsAlpha", "kAlpha"),
               reference, "interval alpha under an affine transform")

  # ratio data have a meaningful zero: invariant to scaling, but NOT to a shift
  positive  <- as.data.frame(lapply(dataset, function(column) abs(column) + 0.5))
  ratioOpts <- raOptions(names(dataset), rep("scale", 4),
                         krippendorffsAlpha = TRUE, krippendorffsAlphaMethod = "ratio")
  ratioReference <- raCoefficient(raRun(positive, ratioOpts), "krippendorffsAlpha", "kAlpha")
  raExpectSame(raCoefficient(raRun(as.data.frame(lapply(positive, function(x) 4 * x)), ratioOpts),
                             "krippendorffsAlpha", "kAlpha"),
               ratioReference, "ratio alpha under scaling")

  shifted <- raCoefficient(raRun(as.data.frame(lapply(positive, function(x) x + 10)), ratioOpts),
                           "krippendorffsAlpha", "kAlpha")
  expect_gt(abs(shifted - ratioReference), 1e-6) # negative control: the test can detect a change
})
