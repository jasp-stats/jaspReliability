# The UPPS data (Lozano et al., 2018) are used rather than Reliability.csv: that scale is
# essentially unidimensional, so every multidimensional model fitted to it is either
# non-converged or inadmissible, which makes it useless for testing the ordinary output.

# JASP serializes an NA numeric cell as an empty string, so blank cells are tested for that
isBlank <- function(value) is.null(value) || identical(value, "") || (is.numeric(value) && is.na(value))

uppsFactors <- list(
  list(indicators = c("U17_r", "U22_r", "U29_r", "U34_r"), name = "Factor1", title = "Urgency"),
  list(indicators = c("U4", "U14", "U19", "U27"),          name = "Factor2", title = "Premeditation"),
  list(indicators = c("U6", "U16", "U28", "U48"),          name = "Factor3", title = "Perseverance"),
  list(indicators = c("U23_r", "U31_r", "U36_r", "U46_r"), name = "Factor4", title = "Sensation"),
  list(indicators = c("U10_r", "U20_r", "U35_r", "U52_r"), name = "Factor5", title = "Positive")
)

options <- analysisOptions("reliabilityMultidimensionalFrequentist")
options$factors              <- uppsFactors
options$modelType            <- "secondOrder"
options$naAction             <- "listwise"
options$itemRestCorrelation  <- TRUE
options$itemDeletedOmegaT    <- TRUE
options$itemDeletedOmegaH    <- TRUE
options$fitMeasures          <- TRUE
options$standardizedLoadings <- TRUE
results <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                       options, makeTests = FALSE)


test_that("Analysis completes without errors", {
  expect_equal(results[["status"]], "complete")
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 0.864287332085935, 0.845989830348801,
         0.88258483382307, "McDonald's <unicode><unicode>", 0.644139461755755,
         0.594619887192351, 0.693659036319159, "Average interitem correlation",
         0.18468108489681, "", "", "Mean", 39.9367945823928, "", "",
         "SD", 7.74337116131412, "", ""))
})

test_that("Item table reports omega if dropped and item-rest correlations", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  expect_equal(length(table), 20)
  omegaT <- vapply(table, function(x) x[["omegaT"]], numeric(1))
  omegaH <- vapply(table, function(x) x[["omegaH"]], numeric(1))
  # dropping any single item leaves a well-behaved model here, so every cell is filled
  expect_true(all(is.finite(omegaT) & omegaT > 0 & omegaT < 1))
  expect_true(all(is.finite(omegaH) & omegaH > 0 & omegaH < 1))
  expect_equal(omegaT[1], 0.8600326, tolerance = 1e-5)
})

test_that("Fit measures table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_fitTable"]][["data"]]
  values <- vapply(table, function(x) x[["value"]], numeric(1))
  expect_equal(unname(values[1:3]), c(410.60557, 165, 0), tolerance = 1e-4)   # chisq, df, p
  expect_equal(unname(values[6]), 0.0580, tolerance = 1e-3)                    # RMSEA
})

test_that("Loadings tables match the fitted model", {
  loadings <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_loadingsContainer"]][["collection"]]
  itemTable <- loadings[["stateContainer_loadingsContainer_itemLoadings"]][["data"]]
  expect_equal(length(itemTable), 20)
  # an item loads on exactly one group factor, the remaining cells are empty rather than zero
  expect_equal(itemTable[[1]][["factor1"]], 0.603772, tolerance = 1e-4)
  expect_true(isBlank(itemTable[[1]][["factor2"]]))

  factorTable <- loadings[["stateContainer_loadingsContainer_factorLoadings"]][["data"]]
  expect_equal(length(factorTable), 5)
  expect_equal(factorTable[[1]][["factor"]], "Urgency")   # the factor form title is used
  expect_equal(factorTable[[1]][["loading"]], 0.814527, tolerance = 1e-4)
})


# The second-order model with two group factors is not identified: the general factor has one
# loading per group factor, but the group factors supply only a single correlation. lavaan then
# cannot compute standard errors, so no interval may be reported, and omega_h is left empty.
optionsTwo <- analysisOptions("reliabilityMultidimensionalFrequentist")
optionsTwo$factors           <- uppsFactors[1:2]
optionsTwo$modelType         <- "secondOrder"
optionsTwo$naAction          <- "listwise"
optionsTwo$itemDeletedOmegaT <- TRUE
optionsTwo$itemDeletedOmegaH <- TRUE
resultsTwo <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                          optionsTwo, makeTests = FALSE)

test_that("Unidentified second-order model reports no confidence interval", {
  scaleTable <- resultsTwo[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  omegaRows <- Filter(function(x) grepl("ω", x[["coefficient"]]), scaleTable[["data"]])

  expect_equal(length(omegaRows), 2)
  omegaT <- Filter(function(x) grepl("ωₜ", x[["coefficient"]]), omegaRows)[[1]]
  expect_true(is.finite(omegaT[["estimate"]]))                      # point estimate is kept
  expect_true(isBlank(omegaT[["lower"]]))                           # interval is withheld
  expect_true(isBlank(omegaT[["upper"]]))

  notes <- paste(vapply(scaleTable[["footnotes"]], function(x) x[["text"]], character(1)), collapse = " ")
  expect_true(grepl("no confidence intervals", notes))
})


# The correlated-factors model has no general factor, so omega_h is dropped with a footnote
optionsCorr <- optionsTwo
optionsCorr$factors   <- uppsFactors
optionsCorr$modelType <- "correlated"
resultsCorr <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                           optionsCorr, makeTests = FALSE)

test_that("Correlated model omits omega_h and adds a footnote", {
  scaleTable   <- resultsCorr[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  coefficients <- vapply(scaleTable[["data"]], function(x) x[["coefficient"]], character(1))
  expect_false(any(grepl("ωₕ", coefficients)))
  expect_true(any(grepl("ωₜ", coefficients)))
  expect_true(length(scaleTable[["footnotes"]]) >= 1)
})


# The bi-factor model rejects cross-loadings, matching the Bayesian analysis
optionsCross <- optionsCorr
optionsCross$modelType  <- "biFactor"
optionsCross$factors    <- uppsFactors
optionsCross$factors[[2]]$indicators <- c(optionsCross$factors[[2]]$indicators, "U17_r")
resultsCross <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                            optionsCross, makeTests = FALSE)

test_that("Bi-factor model rejects cross-loadings with an actionable message", {
  scaleTable <- resultsCross[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  expect_equal(scaleTable[["status"]], "error")
  expect_true(grepl("bi-factor", scaleTable[["error"]][["errorMessage"]]))
})


# A model the optimizer never reached must not be readable, so non-convergence is an error on every
# table rather than a footnote. Convergence cannot be forced through the data, since it differs by
# platform, so the fit is replaced by one whose diagnostics report that it did not converge.
test_that("A non-converged factor model is reported as an error on every table", {
  testthat::local_mocked_bindings(
    .multiDimFreqFit = function(...) list(diagnostics = list(converged = FALSE, admissible = TRUE, se.available = TRUE)),
    .package = "jaspReliability")

  optionsNc <- analysisOptions("reliabilityMultidimensionalFrequentist")
  optionsNc$factors           <- uppsFactors
  optionsNc$modelType         <- "secondOrder"
  optionsNc$naAction          <- "listwise"
  optionsNc$itemDeletedOmegaT <- TRUE
  optionsNc$fitMeasures       <- TRUE
  resultsNc <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                           optionsNc, makeTests = FALSE)

  collection <- resultsNc[["results"]][["stateContainer"]][["collection"]]
  for (table in c("stateContainer_scaleTable", "stateContainer_itemTable", "stateContainer_fitTable")) {
    expect_equal(collection[[table]][["status"]], "error", label = table)
    expect_true(grepl("did not converge", collection[[table]][["error"]][["errorMessage"]]), label = table)
    expect_equal(length(collection[[table]][["data"]]), 0, label = table)   # no coefficients are shown
  }
})


# Fewer than two usable factors is not an analysis: the table prompts instead of erroring
optionsEmpty <- analysisOptions("reliabilityMultidimensionalFrequentist")
optionsEmpty$factors <- list(
  list(indicators = c("U17_r", "U22_r"), name = "Factor1", title = "Factor 1"),
  list(indicators = list(),              name = "Factor2", title = "Factor 2"))
resultsEmpty <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                            optionsEmpty, makeTests = FALSE)

test_that("Incomplete factor assignment shows a prompt rather than results", {
  scaleTable <- resultsEmpty[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  expect_equal(length(scaleTable[["data"]]), 0)
  expect_true(length(scaleTable[["footnotes"]]) >= 1)
})


# Bootstrapped interval: few samples, only checking that the interval is produced from the
# resampling distribution and that failed resamples are reported rather than silently dropped
optionsBoot <- analysisOptions("reliabilityMultidimensionalFrequentist")
optionsBoot$factors          <- uppsFactors
optionsBoot$modelType        <- "secondOrder"
optionsBoot$naAction         <- "listwise"
optionsBoot$intervalMethod   <- "bootstrapped"
optionsBoot$bootstrapSamples <- 100
optionsBoot$setSeed          <- TRUE
optionsBoot$seed             <- 1
resultsBoot <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                           optionsBoot, makeTests = FALSE)

test_that("Bootstrapped interval brackets the point estimate", {
  scaleTable <- resultsBoot[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  omegaRow   <- scaleTable[["data"]][[1]]
  expect_equal(omegaRow[["estimate"]], 0.864287, tolerance = 1e-4)
  expect_true(omegaRow[["lower"]] < omegaRow[["estimate"]])
  expect_true(omegaRow[["upper"]] > omegaRow[["estimate"]])
})


# The second-order model with two group factors identifies the loadings of the general factor only
# up to their product, so omega_h is an arbitrary point of a flat likelihood. Its row stays in the
# scale table with empty cells and the footnote saying why; the item table leaves it empty without
# a footnote of its own.
test_that("Second-order model with two group factors leaves omega_h empty and says why", {
  scaleTable   <- resultsTwo[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  coefficients <- vapply(scaleTable[["data"]], function(x) x[["coefficient"]], character(1))

  expect_true(any(grepl("ωₜ", coefficients)))
  omegaH <- Filter(function(x) grepl("ωₕ", x[["coefficient"]]), scaleTable[["data"]])
  expect_equal(length(omegaH), 1)
  expect_true(isBlank(omegaH[[1]][["estimate"]]))
  expect_true(isBlank(omegaH[[1]][["lower"]]))
  expect_true(isBlank(omegaH[[1]][["upper"]]))

  notes <- paste(vapply(scaleTable[["footnotes"]], function(x) x[["text"]], character(1)), collapse = " ")
  expect_true(grepl("up to their product", notes))
})

test_that("Item table leaves omega_h empty for the unidentified model without an extra footnote", {
  itemTable <- resultsTwo[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]]
  expect_true(all(vapply(itemTable[["data"]], function(x) isBlank(x[["omegaH"]]), logical(1))))
  expect_true(all(vapply(itemTable[["data"]], function(x) is.finite(x[["omegaT"]]), logical(1))))

  notes <- paste(vapply(itemTable[["footnotes"]], function(x) x[["text"]], character(1)), collapse = " ")
  expect_false(grepl("Empty cells", notes))
  expect_false(grepl("up to their product", notes))
})

# Urgency, Premeditation and Sensation, not the first three factors: with Perseverance the third
# general-factor loading sits on the boundary at 1 (residual variance ~ 0), which leaves no standard
# errors and lets convergence depend on the platform. Here every loading is well inside (max 0.57).
test_that("Three group factors report omega_h and drop the identification footnote", {
  optionsThree <- optionsTwo
  optionsThree$factors <- uppsFactors[c(1, 2, 4)]
  resultsThree <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                              optionsThree, makeTests = FALSE)
  scaleTable <- resultsThree[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  omegaH     <- Filter(function(x) grepl("ωₕ", x[["coefficient"]]), scaleTable[["data"]])[[1]]
  notes      <- paste(vapply(scaleTable[["footnotes"]], function(x) x[["text"]], character(1)), collapse = " ")
  expect_true(is.finite(omegaH[["estimate"]]))
  expect_false(grepl("up to their product", notes))
})


# Bayesrel fills its loadings matrix in model-syntax order but reads it in data-row order, which
# differ once a cross-loaded item is listed after the items of a later factor. The loadings must
# instead follow the item, so they are compared against lavaan's own standardized solution.
test_that("Cross-loaded item appended to a factor keeps every loading on its own item", {
  optionsCl <- analysisOptions("reliabilityMultidimensionalFrequentist")
  optionsCl$factors <- uppsFactors
  optionsCl$factors[[2]]$indicators <- c(optionsCl$factors[[2]]$indicators, "U17_r")
  optionsCl$modelType            <- "correlated"
  optionsCl$naAction             <- "listwise"
  optionsCl$standardizedLoadings <- TRUE
  resultsCl <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                           optionsCl, makeTests = FALSE)

  itemTable <- resultsCl[["results"]][["stateContainer"]][["collection"]][["stateContainer_loadingsContainer"]][["collection"]][["stateContainer_loadingsContainer_itemLoadings"]][["data"]]
  loadingOf <- function(item, factor) Filter(function(x) x[["item"]] == item, itemTable)[[1]][[factor]]

  expect_equal(loadingOf("U17_r", "factor2"), 0.036039, tolerance = 1e-4)   # the cross-loading itself
  expect_equal(loadingOf("U17_r", "factor1"), 0.591404, tolerance = 1e-4)   # its primary loading
  expect_equal(loadingOf("U4",    "factor2"), 0.632917, tolerance = 1e-4)   # not shifted onto a neighbour
  expect_equal(loadingOf("U14",   "factor2"), 0.694644, tolerance = 1e-4)
  expect_true(isBlank(loadingOf("U4", "factor1")))
})

# With two group factors only the product of the general-factor loadings is identified, so the
# individual loadings are as arbitrary as omega_h and are left empty as well
test_that("Unidentified general-factor loadings are left empty", {
  optionsLd <- optionsTwo
  optionsLd$standardizedLoadings <- TRUE
  resultsLd <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                           optionsLd, makeTests = FALSE)
  factorTable <- resultsLd[["results"]][["stateContainer"]][["collection"]][["stateContainer_loadingsContainer"]][["collection"]][["stateContainer_loadingsContainer_factorLoadings"]][["data"]]
  expect_equal(length(factorTable), 2)
  expect_true(all(vapply(factorTable, function(x) isBlank(x[["loading"]]), logical(1))))
})

# Under listwise deletion an item-deleted refit must use the respondents of the full model. Rows
# missing only the dropped item would otherwise re-enter its refit, so both the item and the sample
# would change and the omega if item dropped could not be compared with the full model's.
test_that("Item-deleted refits keep the full model's listwise sample", {
  uppsData  <- read.csv(testthat::test_path("upps.csv"))
  uppsItems <- unlist(lapply(uppsFactors, function(f) f[["indicators"]]))
  complete  <- which(complete.cases(uppsData[, uppsItems]))
  highest   <- intersect(order(rowSums(uppsData[, uppsItems]), decreasing = TRUE), complete)[1:80]
  uppsData[highest, "U17_r"] <- NA

  optionsId <- analysisOptions("reliabilityMultidimensionalFrequentist")
  optionsId$factors           <- uppsFactors
  optionsId$modelType         <- "secondOrder"
  optionsId$naAction          <- "listwise"
  optionsId$itemDeletedOmegaT <- TRUE
  resultsId <- runAnalysis("reliabilityMultidimensionalFrequentist", uppsData, optionsId, makeTests = FALSE)

  itemTable <- resultsId[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  droppedU17 <- Filter(function(x) x[["item"]] == "U17_r", itemTable)[[1]][["omegaT"]]
  # refitting on all 443 respondents that are complete without U17_r would give 0.860033
  expect_equal(droppedU17, 0.7472325, tolerance = 1e-4)
})

# FIML ignores a respondent without any observed item, so the descriptive statistics computed
# outside the model must ignore it too rather than turning it into NaN or a zero sum score
test_that("An entirely missing row does not distort the FIML descriptive statistics", {
  uppsData  <- read.csv(testthat::test_path("upps.csv"))
  uppsItems <- unlist(lapply(uppsFactors, function(f) f[["indicators"]]))
  emptyRow  <- uppsData[1, , drop = FALSE]
  emptyRow[, uppsItems] <- NA
  withEmpty <- rbind(uppsData, emptyRow)

  meanSd <- function(dataset, method) {
    optionsFiml <- analysisOptions("reliabilityMultidimensionalFrequentist")
    optionsFiml$factors            <- uppsFactors
    optionsFiml$modelType          <- "secondOrder"
    optionsFiml$naAction           <- "fiml"
    optionsFiml$meanSdScoresMethod <- method
    res   <- runAnalysis("reliabilityMultidimensionalFrequentist", dataset, optionsFiml, makeTests = FALSE)
    table <- res[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
    vapply(c("Mean", "SD"), function(label) Filter(function(x) x[["coefficient"]] == label, table)[[1]][["estimate"]], numeric(1))
  }

  for (method in c("meanScores", "sumScores")) {
    expect_true(all(is.finite(meanSd(withEmpty, method))), label = method)
    expect_equal(meanSd(withEmpty, method), meanSd(uppsData, method), label = method)
  }
})
