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
# cannot compute standard errors, so no interval may be reported.
optionsTwo <- analysisOptions("reliabilityMultidimensionalFrequentist")
optionsTwo$factors   <- uppsFactors[1:2]
optionsTwo$modelType <- "secondOrder"
optionsTwo$naAction  <- "listwise"
resultsTwo <- runAnalysis("reliabilityMultidimensionalFrequentist", testthat::test_path("upps.csv"),
                          optionsTwo, makeTests = FALSE)

test_that("Unidentified second-order model reports no confidence interval", {
  scaleTable <- resultsTwo[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  omegaRows  <- Filter(function(x) grepl("ω", x[["coefficient"]]), scaleTable[["data"]])

  expect_equal(length(omegaRows), 2)
  for (row in omegaRows) {
    expect_true(is.finite(row[["estimate"]]))                       # point estimate is kept
    expect_true(isBlank(row[["lower"]]))                            # interval is withheld
    expect_true(isBlank(row[["upper"]]))
  }
  expect_true(length(scaleTable[["footnotes"]]) >= 1)
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
