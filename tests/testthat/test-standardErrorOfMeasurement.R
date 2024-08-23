

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
  sumScoreCiTable = FALSE,
  sumScoreCiPlots = FALSE,
  variables.types = c("ordinal", "ordinal", "ordinal", "ordinal", "ordinal"),
  variables = c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
)

results <- runAnalysis("standardErrorOfMeasurement", dataset = testthat::test_path("asrm.csv"), options = options, makeTests = FALSE)


# binary data
options <- list(.meta = list(variables = list(hasTypes = TRUE,
                                              shouldEncode = TRUE)), anova = FALSE, ciLevelPlots = 0.95,
                ciLevelTable = 0.95, combinedPointPlot = FALSE, feldt = FALSE,
                feldtNumberOfSplits = "2", hideTable = FALSE, histogramCounts = FALSE,
                irt = TRUE, keats = TRUE, lord = TRUE, lord2 = TRUE, lord2NumberOfSplits = "2",
                minimumGroupSize = 10, mollenkopfFeldt = FALSE, mollenkopfFeldtNumberOfSplits = "2",
                mollenkopfFeldtPolyDegree = 3, plotHeight = 320, plotWidth = 480,
                pointPlots = TRUE, reliabilityValue = 0.5, sumScoreCiPlots = FALSE,
                sumScoreCiPlotsCutoff = FALSE, sumScoreCiPlotsCutoffValue = 0,
                sumScoreCiTable = FALSE, thorndike = FALSE, userReliability = FALSE,
                variables = c("V1", "V2", "V3", "V4"), variables.types = c("nominal",
                                                                           "nominal", "nominal", "nominal"))

results <- runAnalysis("standardErrorOfMeasurement", dataset = "binaryTestDt.csv", options = options)

test_that("Standard error of measurement table results match", {
  table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(
                                     7.0, 0.49433301750863673, 0.0, 0.0, 0.84852813742385702, 0,
                                     18.0, 0.76104647303067507, 0.9151894370278959, 1.0, 0.84852813742385702, 1,
                                     56.0, 0.84058341505800005, 1.0567697356551156, 1.1547005383792515, 1.2955969390869324, 2,
                                     69.0, 0.81285834219580733, 0.9151894370278959, 1.0, 1.0, 3,
                                     50.0, 0.49433301750863673, 0.0, 0.0, 0.0, 4))
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

results <- runAnalysis("standardErrorOfMeasurement", data = testthat::test_path("binaryTestDt.csv"), options = options)

