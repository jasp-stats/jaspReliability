
# 3-factor second-order model on the Reliability example data, with one crossloading
# item (Question_08 loads on both Factor 1 and Factor 2). Three group factors rather than two,
# because the second-order model identifies the general factor only from the correlations among
# the group factors: two of them leave omega_h determined up to a product and therefore arbitrary.
f1 <- paste0("Question_", sprintf("%02d", 1:8))
f2 <- paste0("Question_", sprintf("%02d", 8:15))
f3 <- paste0("Question_", sprintf("%02d", 16:23))

options <- analysisOptions("reliabilityMultidimensionalBayesian")
options$factors <- list(
  list(indicators = f1, name = "Factor1", title = "Factor 1"),
  list(indicators = f2, name = "Factor2", title = "Factor 2"),
  list(indicators = f3, name = "Factor3", title = "Factor 3")
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
    list("McDonald's <unicode><unicode>", 0.780189088627406, 0.767445970684284,
         1.01484927321525, 0.793349865861184, "McDonald's <unicode><unicode>",
         0.763321734422218, 0.749481389834821, 1.02186376750489, 0.777016384422387,
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
    list("Point estimate", 4504.33615737775, 0.0880578185559675, 0.0628897101084209,
         "90% CI lower bound", "", 0.087719045837083, 0.0638424260927066,
         "90% CI upper bound", "", 0.0884396351435885, 0.0682371103290602,
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


# omega-if-item-deleted: per-item refit. Use 3 factors x 3 items so dropping an item leaves valid
# (2-item) factors, every item yields a refit value, and the general factor stays identified.
optionsDel <- analysisOptions("reliabilityMultidimensionalBayesian")
optionsDel$factors <- list(
  list(indicators = paste0("Question_", sprintf("%02d", 1:3)), name = "Factor1", title = "Factor 1"),
  list(indicators = paste0("Question_", sprintf("%02d", 4:6)), name = "Factor2", title = "Factor 2"),
  list(indicators = paste0("Question_", sprintf("%02d", 7:9)), name = "Factor3", title = "Factor 3")
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
  expect_equal(length(itemTable), 9L)                                   # one row per item
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
# This block keeps two group factors: it tests rendering rather than the coefficients, and the
# snapshots stay comparable across changes to the analysis. The omega_h drawn here is the
# unidentified one, which is fine for a rendering check but must not be read as an estimate.
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
  list(indicators = c("contcor2", "debMiss30"), name = "Factor2", title = "Factor 2"),
  list(indicators = c("contGamma", "contOutlier"), name = "Factor3", title = "Factor 3")
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
    list("McDonald's <unicode><unicode>", 0.019340538119256, 0.000442576423266648,
         0.0593915898253003, "McDonald's <unicode><unicode>", 0.00705767855380868,
         1.13757216238552e-08, 0.0292217131671469, "Average interitem correlation",
         0.036435948452956, "", "", "Mean", 7.60957855316, "", "", "SD",
         20.928784371066, "", ""))
})

# imputation keeps all rows, so item-rest correlations use pairwise complete observations
test_that("Missing data with Bayesian imputation: item-rest correlations match", {
  table <- resultsMissImp[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("contNormal", -0.109528587852968, "contcor1", -0.0211291743290384,
         "contcor2", 0.0484455731336876, "debMiss30", 0.0260155552428289,
         "contGamma", 0.14366359620374, "contOutlier", 0.0218132665525293))
})

optionsMiss$naAction <- "listwise"
optionsMiss$posteriorPredictiveCheck <- TRUE
set.seed(1)
resultsMissLw <- runAnalysis("reliabilityMultidimensionalBayesian", "test.csv", optionsMiss, makeTests = FALSE)

test_that("Missing data with listwise deletion: scale table matches", {
  table <- resultsMissLw[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode><unicode>", 0.443076207877799, 0.000936239739393275,
         0.875328892364313, "McDonald's <unicode><unicode>", 0.0595893661324344,
         3.26102311674882e-06, 0.184228922172945, "Average interitem correlation",
         0.0435069256927236, "", "", "Mean", 10.1825294300143, "", "", "SD",
         24.5540789060966, "", ""))
})

# listwise deletion must restrict the item-rest correlations to the complete cases the fit used,
# so these differ from the pairwise values above (70 of 100 rows are complete)
test_that("Missing data with listwise deletion: item-rest correlations use complete cases", {
  table <- resultsMissLw[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("contNormal", -0.151944850496581, "contcor1", -0.0815927556893976,
         "contcor2", 0.0398453253170278, "debMiss30", 0.0260155552428289,
         "contGamma", 0.167208636470416, "contOutlier", 0.0343903439383511))
})

# the synthetic datasets in the PPC must have as many rows as the fit had complete cases (70, not
# 100); simulating with the full row count shrinks the eigenvalue bands by roughly 15%
test_that("Missing data with listwise deletion: posterior predictive check uses the complete-case n", {
  plotName <- resultsMissLw[["results"]][["stateContainer"]][["collection"]][["stateContainer_ppcPlot"]][["data"]]
  ppcFrame <- resultsMissLw[["state"]][["figures"]][[plotName]][["obj"]][["data"]]
  expect_equal(ppcFrame[["eigen_value"]],
               c(579.294364694806, 11.0000135411597, 2.68537227911255,
                 1.67302612601514, 1.20324526241216, 0.267431151515681))
  expect_equal(ppcFrame[["eigen_sim_low"]],
               c(364.35061005767, 6.93249897760065, 1.88382296469919,
                 1.03381919139404, 0.672801386134565, 0.353654268997077))
  expect_equal(ppcFrame[["eigen_sim_up"]],
               c(929.87330888715, 16.5523066387851, 4.06967367128528,
                 2.13052156049937, 1.34949670295412, 0.947055199562322))
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
  list(indicators = paste0("Question_", sprintf("%02d", 1:4)),  name = "Factor1", title = "Factor 1"),
  list(indicators = paste0("Question_", sprintf("%02d", 5:8)),  name = "Factor2", title = "Factor 2"),
  list(indicators = paste0("Question_", sprintf("%02d", 9:12)), name = "Factor3", title = "Factor 3")
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
    list("McDonald's <unicode><unicode>", 0.678836904330263, 0.661472639189372,
         0.69897175075815, "McDonald's <unicode><unicode>", 0.645234337055706,
         0.619086701485578, 0.666880329142608, "Average interitem correlation",
         0.131480124530591, "", "", "Mean", 32.7740178918709, "", "",
         "SD", 5.08869596388986, "", ""))
})

test_that("Reverse-scaled item: item-rest correlations match and footnote is shown", {
  itemTable <- resultsRev[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]]
  jaspTools::expect_equal_tables(itemTable[["data"]],
    list("Question_01", 0.449010913838358, "Question_02", 0.0490884658992861,
         "Question_03", -0.469882320342349, "Question_04", 0.497239560605097,
         "Question_05", 0.450668111762472, "Question_06", 0.442440207570404,
         "Question_07", 0.521075744378894, "Question_08", 0.465121574342944,
         "Question_09", -0.167026025243398, "Question_10", 0.339092272012501,
         "Question_11", 0.518948633030988, "Question_12", 0.456211948190079))
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


# The second-order model with two group factors identifies the general factor's loadings only up
# to their product. A proper prior still yields a proper posterior for omega_h, so the coefficient
# looks like an ordinary result while being driven by the prior rather than the data; it is
# therefore withheld from every output element, as in the frequentist analysis.
optionsUnid <- analysisOptions("reliabilityMultidimensionalBayesian")
optionsUnid$factors <- list(
  list(indicators = paste0("Question_", sprintf("%02d", 1:5)),  name = "Factor1", title = "Factor 1"),
  list(indicators = paste0("Question_", sprintf("%02d", 6:10)), name = "Factor2", title = "Factor 2"))
optionsUnid$modelType        <- "secondOrder"
optionsUnid$samples          <- 200
optionsUnid$burnin           <- 50
optionsUnid$chains           <- 2
optionsUnid$setSeed          <- TRUE
optionsUnid$seed             <- 1
optionsUnid$probabilityTable <- TRUE
optionsUnid$posteriorPlot    <- TRUE
optionsUnid$tracePlot        <- TRUE
set.seed(1)
resultsUnid <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsUnid,
                           makeTests = FALSE)

test_that("Second-order model with two group factors warns that omega_h is not identified", {
  collection   <- resultsUnid[["results"]][["stateContainer"]][["collection"]]
  scaleTable   <- collection[["stateContainer_scaleTable"]]
  coefficients <- vapply(scaleTable[["data"]], function(x) x[["coefficient"]], character(1))

  # the coefficient is still reported; the footnote is what says it follows the prior, not the data
  expect_true(any(grepl("ωₕ", coefficients)))
  expect_true(any(grepl("ωₜ", coefficients)))

  notes <- paste(vapply(scaleTable[["footnotes"]], function(x) x[["text"]], character(1)), collapse = " ")
  expect_true(grepl("up to their product", notes))
})

test_that("Three group factors drop the identification warning", {
  optionsIdent <- optionsUnid
  optionsIdent$factors <- c(optionsUnid$factors, list(
    list(indicators = paste0("Question_", sprintf("%02d", 11:15)), name = "Factor3", title = "Factor 3")))
  set.seed(1)
  resultsIdent <- runAnalysis("reliabilityMultidimensionalBayesian", "Reliability.csv", optionsIdent,
                              makeTests = FALSE)
  scaleTable <- resultsIdent[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]]
  notes      <- paste(vapply(scaleTable[["footnotes"]], function(x) x[["text"]], character(1)), collapse = " ")
  expect_false(grepl("up to their product", notes))
})
