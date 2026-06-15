
# 2-factor second-order model on the Reliability example data, with one crossloading
# item (Question_12 loads on both factors)
f1 <- paste0("Question_", sprintf("%02d", 1:12))
f2 <- paste0("Question_", sprintf("%02d", 12:23))

options <- analysisOptions("reliabilityMultidimensionalBayesian")
options$factors <- list(
  list(indicators = f1, name = "Factor1", title = "Factor 1"),
  list(indicators = f2, name = "Factor2", title = "Factor 2")
)
options$modelType <- "secondOrder"
# scale coefficients are always displayed; item statistics are opt-in
options$itemRestCor <- TRUE
options$probTable <- TRUE
options$probTableValueLow <- 0.7
options$probTableValueHigh <- 0.9
options$noSamples <- 200
options$noBurnin <- 50
options$noChains <- 2
options$rHat <- TRUE
options$setSeed <- TRUE
options$seed <- 1
options$fitMeasures <- TRUE
set.seed(1)
results <- runAnalysis("reliabilityMultidimensionalBayesian", testthat::test_path("Reliability.csv"), options, makeTests = FALSE)

test_that("Analysis completes without errors", {
  expect_equal(results[["status"]], "complete")
})

test_that("Bayesian Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 0.7788594113243, 0.7655848422395,
         0.998727899094338, 0.791322733994566, "McDonald's <unicode><unicode>",
         0.757721516474659, 0.741581961183294, 1.00323445815622, 0.771598841259063,
         "Average interitem correlation", 0.130074006470925, "", "",
         "", "Mean", 61.437183975107, "", "", "", "SD", 9.01983558246063,
         "", "", ""))
})

test_that("Bayesian Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("Question_01", 0.491798910912644, "Question_02", -0.105371122635483,
         "Question_03", -0.434761937472932, "Question_04", 0.533243866223406,
         "Question_05", 0.453605309151933, "Question_06", 0.477809305171136,
         "Question_07", 0.559884251471989, "Question_08", 0.493112167091624,
         "Question_09", -0.0811768317813592, "Question_10", 0.345790483737811,
         "Question_11", 0.540473158732128, "Question_12", 0.518535898276985,
         "Question_13", 0.558606070855461, "Question_14", 0.528139241251555,
         "Question_15", 0.457047628630301, "Question_16", 0.524877464899673,
         "Question_17", 0.568384535026814, "Question_18", 0.578263548919627,
         "Question_19", -0.247956846021939, "Question_20", 0.265420749483938,
         "Question_21", 0.514201155527712, "Question_22", -0.120800284171491,
         "Question_23", -0.0127987366827629))
})

test_that("Probability table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_probabilityTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 1, "McDonald's <unicode><unicode>", 1))
})

test_that("Fit measures table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_fitTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("Point estimate", 4494.12650123064, 0.089187044699341, 0.063017234802653,
         "90% CI lower bound", "", 0.0888221787024989, 0.0641069635806036,
         "90% CI upper bound", "", 0.0895624113228447, 0.0688465705668997,
         "Relative to cutoff", "", 0, ""))
})


# correlated-factors model: McDonald's omega_h is undefined and must be dropped with a footnote
optionsCorr <- options
optionsCorr$modelType <- "correlated"
set.seed(1)
resultsCorr <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsCorr, makeTests = FALSE)

test_that("Correlated model omits omega_h and adds a footnote", {
  scaleTable   <- resultsCorr[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  coefficients <- vapply(scaleTable[["data"]], function(x) x[["coefficient"]], character(1))
  expect_false(any(grepl("ωₕ", coefficients)))   # no McDonald's omega_h row
  expect_true(length(scaleTable[["footnotes"]]) >= 1)
})


# omega-if-item-deleted: per-item refit. Use 2 factors x 3 items so dropping an item leaves valid
# (2-item) factors and every item yields a refit value.
optionsDel <- analysisOptions("reliabilityMultidimensionalBayesian")
optionsDel$factors <- list(
  list(indicators = paste0("Question_", sprintf("%02d", 1:3)), name = "Factor1", title = "Factor 1"),
  list(indicators = paste0("Question_", sprintf("%02d", 4:6)), name = "Factor2", title = "Factor 2")
)
optionsDel$modelType        <- "secondOrder"
optionsDel$itemDeletedOmegaT <- TRUE
optionsDel$itemDeletedOmegaH <- TRUE
optionsDel$noSamples <- 100
optionsDel$noBurnin  <- 30
optionsDel$noChains  <- 2
optionsDel$setSeed   <- TRUE
optionsDel$seed      <- 1
set.seed(1)
resultsDel <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsDel, makeTests = FALSE)

test_that("Omega-if-item-deleted produces a populated item table", {
  itemTable <- resultsDel[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  expect_equal(length(itemTable), 6L)                                   # one row per item
  omtDropped <- vapply(itemTable, function(x) x[["omegaT"]], numeric(1))
  expect_true(all(is.finite(omtDropped)))                              # every item refit succeeded
})
