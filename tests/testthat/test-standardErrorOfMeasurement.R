

# test the sem functions for ordinal data

options <- list(
  .meta = list(
    variables = list(
      hasTypes = TRUE,
      shouldEncode = TRUE
    )
  ),
  anova = TRUE,
  combinedPointPlot = TRUE,
  feldt = TRUE,
  feldtNumberOfSplits = "5",
  hideTable = FALSE,
  histogramCounts = FALSE,
  irt = TRUE,
  keats = FALSE,
  lord = FALSE,
  lord2 = FALSE,
  lord2NumberOfSplits = "5",
  minimumGroupSize = 10,
  mollenkopfFeldt = TRUE,
  mollenkopfFeldtNumberOfSplits = "5",
  mollenkopfFeldtPolyDegree = 3,
  plotHeight = 320,
  plotWidth = 480,
  pointPlots = TRUE,
  reliabilityValue = 0.5,
  thorndike = TRUE,
  userReliability = FALSE,
  variables = list(
    types = c("ordinal", "ordinal", "ordinal", "ordinal", "ordinal"),
    value = c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
  )
)

results <- runAnalysis("standardErrorOfMeasurement", data = "asrm.csv", options = options)
# results <- runAnalysis("standardErrorOfMeasurement", dataset = Bayesrel::asrm, options = options, makeTests = TRUE)



# binary data
options <- analysisOptions("standardErrorOfMeasurement")
options$anova <- FALSE
options$combinedPointPlot <- TRUE
options$feldt <- FALSE
options$feldtNumberOfSplits <- 2
options$hideTable <- FALSE
options$histogramCounts <- FALSE
options$irt <- TRUE
options$keats <- TRUE
options$lord <- TRUE
options$lord2 <- TRUE
options$lord2NumberOfSplits <- 2
options$minimumGroup <- 10
options$mollenkopfFeldt <- FALSE
options$mollenkopfFeldtNumberOfSplits <- 2
options$mollenkopfFeldtPolyDegree <- 3
options$plotHeight <- 320
options$plotWidth <- NULL
options$pointPlots <- TRUE
options$reliabilityValue <- 0.5
options$thorndike <- FALSE
options$userReliability <- FALSE
options$variables <- c("V1", "V2", "V3", "V4", "V5")

results <- runAnalysis("standardErrorOfMeasurement", data = "binaryTestDt.csv", options = options)
# results <- runAnalysis("standardErrorOfMeasurement", data = "tests/testthat/binaryTestDt.csv", options = options, makeTests = TRUE)

test_that("Standard error of measurement table results match", {
  table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(7, "", 0, 0, 0.781735959970572, 0, 11, 0.79134461970981, 0.88737278046041,
                                      1, 0.781735959970572, 1, 53, 0.944934004161848, 1.08680526188138,
                                      1.22474487139159, 1.31751624529027, 2, 55, 0.980178580758181,
                                      1.08680526188138, 1.22474487139159, 1.32801971507819, 3, 48,
                                      0.826797499898133, 0.88737278046041, 1, 1, 4, 26, "", 0, 0,
                                      0, 5))
})

test_that("Combined binary plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_combinedPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "combined-plot-binary")
})

test_that("IRT-2PL plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_irtPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "irt-2pl")
})

test_that("Keats plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_keatsPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "keats")
})

test_that("Lord's compound plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_lordPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lord-s-compound")
})


# and some extra plots and stuff

options <- list(.meta = list(variables = list(hasTypes = TRUE,
                                              shouldEncode = TRUE)), anova = TRUE, ciLevelPlots = 0.99,
                ciLevelTable = 0.1, combinedPointPlot = FALSE, feldt = TRUE,
                feldtNumberOfSplits = "5", hideTable = FALSE, histogramCounts = FALSE,
                irt = FALSE, keats = TRUE, lord = FALSE, lord2 = FALSE, lord2NumberOfSplits = "5",
                minimumGroupSize = 10, mollenkopfFeldt = FALSE, mollenkopfFeldtNumberOfSplits = "5",
                mollenkopfFeldtPolyDegree = 3, plotHeight = 320, plotWidth = 480,
                pointPlots = FALSE, reliabilityValue = 0.5, sumScoreCiPlots = TRUE,
                sumScoreCiPlotsCutoff = TRUE, sumScoreCiPlotsCutoffValue = 3,
                sumScoreCiTable = TRUE, thorndike = FALSE, userReliability = FALSE,
                variables = c("V1", "V2", "V3", "V4", "V5"), variables.types = c("nominal",
                                                                                 "nominal", "nominal", "nominal", "nominal"))

results <- runAnalysis("standardErrorOfMeasurement", data = "binaryTestDt.csv", options = options)
# results <- runAnalysis("standardErrorOfMeasurement", data = "tests/testthat/binaryTestDt.csv", options = options, makeTests = TRUE)

