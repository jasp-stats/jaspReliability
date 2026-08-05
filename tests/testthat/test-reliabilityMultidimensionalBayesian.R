
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
options$itemRestCorrelation <- TRUE
options$probabilityTable <- TRUE
options$probabilityTableLowerBound <- 0.7
options$probabilityTableUpperBound <- 0.9
options$samples <- 200
options$burnin <- 50
options$chains <- 2
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


# bi-factor model: regression test, previously crashed in Bayesrel when param.out = TRUE
# (psis array sized without the g-factor). No crossloading: Bayesrel rejects those for bi-factor models.
optionsBif <- options
optionsBif$modelType <- "biFactor"
optionsBif$factors <- list(
  list(indicators = paste0("Question_", sprintf("%02d", 1:12)),  name = "Factor1", title = "Factor 1"),
  list(indicators = paste0("Question_", sprintf("%02d", 13:23)), name = "Factor2", title = "Factor 2")
)
set.seed(1)
resultsBif <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsBif, makeTests = FALSE)

test_that("Bi-factor model runs and reports both omegas", {
  expect_equal(resultsBif[["status"]], "complete")
  scaleTable   <- resultsBif[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  coefficients <- vapply(scaleTable[["data"]], function(x) x[["coefficient"]], character(1))
  expect_true(any(grepl("ωₜ", coefficients)))
  expect_true(any(grepl("ωₕ", coefficients)))
  ests <- vapply(scaleTable[["data"]][grepl("ω", coefficients)], function(x) x[["estimate"]], numeric(1))
  expect_true(all(is.finite(ests) & ests > 0 & ests <= 1))
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
optionsDel$samples <- 100
optionsDel$burnin  <- 30
optionsDel$chains  <- 2
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

test_that("Omega-if-item-deleted reports credible intervals bracketing the point estimate", {
  itemTable <- resultsDel[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  for (coefficient in c("omegaT", "omegaH")) {
    est   <- vapply(itemTable, function(x) x[[coefficient]],                    numeric(1))
    lower <- vapply(itemTable, function(x) x[[paste0(coefficient, "Lower")]],   numeric(1))
    upper <- vapply(itemTable, function(x) x[[paste0(coefficient, "Upper")]],   numeric(1))
    expect_true(all(is.finite(c(lower, upper))), label = coefficient)
    expect_true(all(lower <= est & est <= upper), label = coefficient)
    expect_true(all(lower >= 0 & upper <= 1), label = coefficient)
  }
})


# plots: posterior densities (prior displayed, shaded probability region, fixed x-range),
# traceplots, and the posterior predictive check. The scale table anchors the chains numerically,
# so a plot snapshot failure with a passing table points at rendering, not sampling.
optionsPlots <- analysisOptions("reliabilityMultidimensionalBayesian")
optionsPlots$factors <- list(
  list(indicators = paste0("Question_", sprintf("%02d", 1:4)), name = "Factor1", title = "Factor 1"),
  list(indicators = paste0("Question_", sprintf("%02d", 5:8)), name = "Factor2", title = "Factor 2")
)
optionsPlots$samples <- 100
optionsPlots$burnin  <- 30
optionsPlots$chains  <- 2
optionsPlots$setSeed <- TRUE
optionsPlots$seed    <- 1
optionsPlots$posteriorPlot               <- TRUE
optionsPlots$posteriorPlotFixedRange     <- TRUE
optionsPlots$posteriorPlotPriorDisplayed <- TRUE
optionsPlots$posteriorPlotShaded         <- TRUE
optionsPlots$probabilityTable            <- TRUE
optionsPlots$tracePlot                   <- TRUE
optionsPlots$posteriorPredictiveCheck    <- TRUE
set.seed(1)
resultsPlots <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsPlots, makeTests = FALSE)

test_that("Bayesian Scale Reliability Statistics table anchors the plot run", {
  table <- resultsPlots[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 0.578143408799999, 0.556397832369372,
         0.59545593378841, "McDonald's <unicode><unicode>", 0.483847464320034,
         0.456430592167667, 0.505980528973344, "Average interitem correlation",
         0.101362519780686, "", "", "Mean", 19.4791909762738, "", "",
         "SD", 3.57246477368688, "", ""))
})

test_that("Posterior plot omega_t matches", {
  plotName <- resultsPlots[["results"]][["stateContainer"]][["collection"]][["stateContainer_posteriorPlots"]][["collection"]][["stateContainer_posteriorPlots_omegaT"]][["data"]]
  testPlot <- resultsPlots[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "posterior-omega-t")
})

test_that("Posterior plot omega_h matches", {
  plotName <- resultsPlots[["results"]][["stateContainer"]][["collection"]][["stateContainer_posteriorPlots"]][["collection"]][["stateContainer_posteriorPlots_omegaH"]][["data"]]
  testPlot <- resultsPlots[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "posterior-omega-h")
})

test_that("Traceplot omega_t matches", {
  plotName <- resultsPlots[["results"]][["stateContainer"]][["collection"]][["stateContainer_tracePlots"]][["collection"]][["stateContainer_tracePlots_omegaT"]][["data"]]
  testPlot <- resultsPlots[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "trace-omega-t")
})

test_that("Traceplot omega_h matches", {
  plotName <- resultsPlots[["results"]][["stateContainer"]][["collection"]][["stateContainer_tracePlots"]][["collection"]][["stateContainer_tracePlots_omegaH"]][["data"]]
  testPlot <- resultsPlots[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "trace-omega-h")
})

test_that("Posterior predictive check plot matches", {
  plotName <- resultsPlots[["results"]][["stateContainer"]][["collection"]][["stateContainer_ppcPlot"]][["data"]]
  testPlot <- resultsPlots[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "posterior-predictive-check")
})


# missing values: debMiss30 has 30% missing observations; run both Bayesian imputation and
# listwise deletion
optionsMiss <- analysisOptions("reliabilityMultidimensionalBayesian")
optionsMiss$factors <- list(
  list(indicators = c("contNormal", "contcor1"), name = "Factor1", title = "Factor 1"),
  list(indicators = c("contcor2", "debMiss30"), name = "Factor2", title = "Factor 2")
)
optionsMiss$samples <- 100
optionsMiss$burnin  <- 30
optionsMiss$chains  <- 2
optionsMiss$setSeed <- TRUE
optionsMiss$seed    <- 1
optionsMiss$itemRestCorrelation <- TRUE
optionsMiss$naAction <- "imputation"
set.seed(1)
resultsMissImp <- runAnalysis("reliabilityMultidimensionalBayesian", "test.csv", optionsMiss, makeTests = FALSE)

test_that("Missing data with Bayesian imputation: scale table matches", {
  table <- resultsMissImp[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 0.0233447741856189, 0.00049167696902388,
         0.0910937592021331, "McDonald's <unicode><unicode>", 0.0166544859440133,
         6.12745270505382e-06, 0.0622094079174373, "Average interitem correlation",
         0.0991999010928576, "", "", "Mean", 5.54108024245, "", "", "SD",
         20.3923746243839, "", ""))
})

# imputation keeps all rows, so item-rest correlations use pairwise complete observations
test_that("Missing data with Bayesian imputation: item-rest correlations match", {
  table <- resultsMissImp[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("contNormal", -0.109440767266664, "contcor1", 0.00472160500708882,
         "contcor2", 0.0422940370884689, "debMiss30", -0.122450817493202))
})

optionsMiss$naAction <- "listwise"
optionsMiss$posteriorPredictiveCheck <- TRUE
set.seed(1)
resultsMissLw <- runAnalysis("reliabilityMultidimensionalBayesian", "test.csv", optionsMiss, makeTests = FALSE)

test_that("Missing data with listwise deletion: scale table matches", {
  table <- resultsMissLw[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 0.563750589461712, 0.000541705738093422,
         0.919261583745341, "McDonald's <unicode><unicode>", 0.27518219902584,
         3.25887297714518e-05, 0.694773241637106, "Average interitem correlation",
         0.115180274099673, "", "", "Mean", 8.21151657631428, "", "", "SD",
         23.8940400508469, "", ""))
})

# listwise deletion must restrict the item-rest correlations to the complete cases the fit used,
# so these differ from the pairwise values above (70 of 100 rows are complete)
test_that("Missing data with listwise deletion: item-rest correlations use complete cases", {
  table <- resultsMissLw[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("contNormal", -0.153547377551663, "contcor1", -0.0553604988574586,
         "contcor2", 0.033074200678642, "debMiss30", -0.122450817493202))
})

# the synthetic datasets in the PPC must have as many rows as the fit had complete cases (70, not
# 100); simulating with the full row count shrinks the eigenvalue bands by roughly 15%
test_that("Missing data with listwise deletion: posterior predictive check uses the complete-case n", {
  plotName <- resultsMissLw[["results"]][["stateContainer"]][["collection"]][["stateContainer_ppcPlot"]][["data"]]
  ppcFrame <- resultsMissLw[["state"]][["figures"]][[plotName]][["obj"]][["data"]]
  expect_equal(ppcFrame[["eigen_value"]],
               c(579.201872609291, 1.71880320289344, 1.20706696531347, 0.293018153611386))
  expect_equal(ppcFrame[["eigen_sim_low"]],
               c(371.578866898854, 1.05942659277699, 0.741965436972952, 0.486848523306179))
  expect_equal(ppcFrame[["eigen_sim_up"]],
               c(836.927650710798, 2.08658774538687, 1.35100742373142, 1.00129862123668))
})

# listwise deletion can leave too few rows to analyse even when every column on its own has enough
# observations, so the data checks must run on the complete cases rather than on all rows
test_that("Listwise deletion with too few complete cases is rejected", {
  optionsFew <- analysisOptions("reliabilityMultidimensionalBayesian")
  optionsFew$factors <- list(
    list(indicators = c("i1", "i2"), name = "Factor1", title = "Factor 1"),
    list(indicators = c("i3", "i4"), name = "Factor2", title = "Factor 2")
  )
  optionsFew$samples  <- 60
  optionsFew$burnin   <- 20
  optionsFew$chains   <- 2
  optionsFew$setSeed  <- TRUE
  optionsFew$seed     <- 1
  optionsFew$naAction <- "listwise"

  set.seed(7)
  n  <- 30
  dt <- data.frame(i1 = rnorm(n), i2 = rnorm(n), i3 = rnorm(n), i4 = rnorm(n))
  dt$i1[1:10]  <- NA   # every column keeps at least 20 observations,
  dt$i2[11:18] <- NA   # but only two rows are complete
  dt$i3[19:24] <- NA
  dt$i4[25:28] <- NA
  expect_equal(sum(complete.cases(dt)), 2)

  resultsFew <- runAnalysis("reliabilityMultidimensionalBayesian", dt, optionsFew, makeTests = FALSE)
  expect_equal(resultsFew[["status"]], "validationError")
  expect_match(resultsFew[["results"]][["errorMessage"]], "Number of observations")
})

# a factor with a single indicator is under-identified; the model must not be handed to Bayesrel,
# which would stop with an opaque "invalid 'n' argument"
test_that("A factor with only one indicator leaves the analysis not ready", {
  optionsOne <- analysisOptions("reliabilityMultidimensionalBayesian")
  optionsOne$factors <- list(
    list(indicators = c("Question_01", "Question_02"), name = "Factor1", title = "Factor 1"),
    list(indicators = c("Question_03", "Question_04"), name = "Factor2", title = "Factor 2"),
    list(indicators = c("Question_05"),                name = "Factor3", title = "Factor 3")
  )
  optionsOne$samples <- 60
  optionsOne$burnin  <- 20
  optionsOne$chains  <- 2
  optionsOne$setSeed <- TRUE
  optionsOne$seed    <- 1

  resultsOne <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsOne,
                            makeTests = FALSE)
  expect_equal(resultsOne[["status"]], "complete")

  scaleTable <- resultsOne[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  expect_null(scaleTable[["error"]])
  footnotes <- vapply(scaleTable[["footnotes"]], function(f) f[["text"]], character(1))
  expect_true(any(grepl("at least two factors with at least two items", footnotes)))
})

# with pairwise deletion the covariance matrix need not be positive semi-definite, so an observed
# eigenvalue can be negative; the y-axis must still show it instead of clipping it away
test_that("Posterior predictive check keeps observed eigenvalues inside the plot range", {
  optionsNeg <- analysisOptions("reliabilityMultidimensionalBayesian")
  optionsNeg$factors <- list(
    list(indicators = c("i1", "i2"), name = "Factor1", title = "Factor 1"),
    list(indicators = c("i3", "i4"), name = "Factor2", title = "Factor 2")
  )
  optionsNeg$samples  <- 200
  optionsNeg$burnin   <- 50
  optionsNeg$chains   <- 2
  optionsNeg$setSeed  <- TRUE
  optionsNeg$seed     <- 1
  optionsNeg$naAction <- "imputation"
  optionsNeg$posteriorPredictiveCheck <- TRUE

  set.seed(2)
  n  <- 80
  a  <- rnorm(n)
  dt <- data.frame(i1 =  a + rnorm(n, sd = .3), i2 =  a + rnorm(n, sd = .3),
                   i3 = -a + rnorm(n, sd = .3), i4 =  a + rnorm(n, sd = .3))
  dt$i1[1:20]  <- NA   # staggered missingness, so no two columns share the same respondents
  dt$i2[21:40] <- NA
  dt$i3[41:55] <- NA
  dt$i4[56:70] <- NA

  resultsNeg <- runAnalysis("reliabilityMultidimensionalBayesian", dt, optionsNeg, makeTests = FALSE)
  plotName   <- resultsNeg[["results"]][["stateContainer"]][["collection"]][["stateContainer_ppcPlot"]][["data"]]
  ppcPlot    <- resultsNeg[["state"]][["figures"]][[plotName]][["obj"]]
  ppcFrame   <- ppcPlot[["data"]]

  expect_lt(min(ppcFrame[["eigen_value"]]), 0)   # the case that used to be clipped at zero
  yLimits <- ggplot2::layer_scales(ppcPlot)$y$get_limits()
  expect_true(all(ppcFrame[["eigen_value"]] >= yLimits[1] & ppcFrame[["eigen_value"]] <= yLimits[2]))
})


# reverse-scaled items: Question_02 is recoded before the analysis and flagged in a footnote
optionsRev <- analysisOptions("reliabilityMultidimensionalBayesian")
optionsRev$factors <- list(
  list(indicators = paste0("Question_", sprintf("%02d", 1:4)), name = "Factor1", title = "Factor 1"),
  list(indicators = paste0("Question_", sprintf("%02d", 5:8)), name = "Factor2", title = "Factor 2")
)
optionsRev$samples <- 100
optionsRev$burnin  <- 30
optionsRev$chains  <- 2
optionsRev$setSeed <- TRUE
optionsRev$seed    <- 1
optionsRev$reverseScaledItems  <- "Question_02"
optionsRev$itemRestCorrelation <- TRUE
set.seed(1)
resultsRev <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsRev, makeTests = FALSE)

test_that("Reverse-scaled item: scale table matches", {
  table <- resultsRev[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 0.641395804019822, 0.619352924110949,
         0.655345811411494, "McDonald's <unicode><unicode>", 0.557421766815874,
         0.530671217737979, 0.578844398457801, "Average interitem correlation",
         0.122401576829532, "", "", "Mean", 22.2322053675613, "", "",
         "SD", 3.69553171121895, "", ""))
})

test_that("Reverse-scaled item: item-rest correlations match and footnote is shown", {
  itemTable <- resultsRev[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]]
  jaspTools::expect_equal_tables(itemTable[["data"]],
    list("Question_01", 0.426281726411376, "Question_02", 0.0743575770421569,
         "Question_03", -0.511411265341137, "Question_04", 0.47685808495137,
         "Question_05", 0.43453575979988, "Question_06", 0.411639313079386,
         "Question_07", 0.507494397783249, "Question_08", 0.367333855613921))
  footnotes <- vapply(itemTable[["footnotes"]], function(f) f[["text"]], character(1))
  expect_true(any(grepl("reverse", footnotes, ignore.case = TRUE)))
})


# cross-loadings with the bi-factor model. Sampling of cross-loaded items is already covered by the
# main fixture at the top of this file (Question_12 loads on both factors); Bayesrel refuses them for
# the bi-factor model, and the analysis reports that before the sampler is started.
test_that("Cross-loaded item is rejected by the bi-factor model with a clear error", {
  opts <- options
  opts$modelType <- "biFactor"
  set.seed(1)
  res <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", opts, makeTests = FALSE)
  scaleTable <- res[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  expect_match(scaleTable[["error"]][["errorMessage"]], "bi-factor model does not support")
})
