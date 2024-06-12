

# test the sem functions for ordinal data


options <- analysisOptions("standardErrorOfMeasurement")

options$anova <- TRUE
options$combinedPointPlot <- TRUE
options$feldt <- TRUE
options$feldtNumberOfSplits <- 3
options$hideTable <- FALSE
options$histogramCounts <- TRUE
options$irt <- TRUE
options$keats <- FALSE
options$lord <- FALSE
options$lord2 <- FALSE
options$lord2NumberOfSplits <- 2
options$minimumGroupSize <- 10
options$mollenkopfFeldt <- TRUE
options$mollenkopfFeldtNumberOfSplits <- 3
options$mollenkopfFeldtPolyDegree <- 3
options$plotHeight <- 320
options$plotWidth <- 480
options$pointPlots <- TRUE
options$reliabilityValue <- 0.5
options$thorndike <- TRUE
options$userReliability <- TRUE
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")

results <- runAnalysis("standardErrorOfMeasurement", data = "asrm.csv", options = options)
# results <- runAnalysis("standardErrorOfMeasurement", data = Bayesrel::asrm, options = options, makeTests = TRUE)



test_that("Standard error of measurement table results match", {
  table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(1.50791562224467, 3, 1.1115554667022, 1.4287178895478, 2.00941503560864,
                                      5, 0.925820099772551, 1.50791562224467, 0, 1.1115554667022,
                                      1.49703106451365, 2.01714592845412, 6, 0.925820099772551, 1.50791562224467,
                                      2, 1.1115554667022, 1.54510894064586, 1.97233761068898, 7, 0.925820099772551,
                                      1.50791562224467, 4, 1.1115554667022, 1.53986972817318, 1.88705837675172,
                                      8, 0.925820099772551, 1.50791562224467, 6, 1.1115554667022,
                                      1.51867216206786, 1.77229381515355, 9, 0.925820099772551, 1.1878488323372,
                                      7, 1.23883906227654, 1.52424619129076, 1.64004587459843, 10,
                                      0.866025403784439, 1.1878488323372, 5, 1.23883906227654, 1.5235523748015,
                                      1.50555348569793, 11, 0.866025403784439, 1.72451684054791, 7,
                                      1.67355347987449, 1.48885399722157, 1.38973901151331, 12, 1.64189930669738,
                                      1.72451684054791, 9, 1.67355347987449, 1.44988865277158, 1.32035101644324,
                                      13, 1.64189930669738, 1.62117399689396, 13, 1.68530017693897,
                                      1.46250310836144, 1.32743166609819, 14, 1.53589529557661, 1.43811745632331,
                                      11, 1.24979337135159, 1.47211242437243, 1.43078483885436, 15,
                                      1.57249078621379, 1.91287503750007, 6, 1.73681574314025, 1.40088808515813,
                                      1.63035724851625, 16, 1.84390889145858, 1.91287503750007, 1,
                                      1.73681574314025, 1.22842056164903, 1.91163804277982, 17, 1.84390889145858,
                                      1.91287503750007, 2, 1.73681574314025, 0.961775139235799, 2.25733419588177,
                                      18, 1.84390889145858, 1.91287503750007, 1, 1.73681574314025,
                                      "", 2.65352302119895, 19, 1.84390889145858, 1.91287503750007,
                                      0, 1.73681574314025, "", 3.0904315124867, 20, 1.84390889145858,
                                      1.91287503750007, 0, 1.73681574314025, "", 3.56143677050039,
                                      21, 1.84390889145858, 1.91287503750007, 1, 1.73681574314025,
                                      "", 4.06200762001967, 22, 1.84390889145858))
})

test_that("Combined plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_combinedPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "combined-plot")
})

test_that("Histogram of counts per sum score group plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_histPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "histogram-of-counts-per-sum-score-group")
})

test_that("ANOVA plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_anovaPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "anova")
})

test_that("Feldt plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_feldtPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "feldt")
})

test_that("IRT-GRM plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_irtPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "irt-grm")
})

test_that("Mollenkopf-Feldt plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_mfPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "mollenkopf-feldt")
})

test_that("Thorndike plot matches", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_pointPlotsContainer"]][["collection"]][["semMainContainer_pointPlotsContainer_thorndikePlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "thorndike")
})



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
