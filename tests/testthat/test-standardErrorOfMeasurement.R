

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

test_that("Standard error of measurement table results match", {

table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
jaspTools::expect_equal_tables(table,
                               list(1.50791562224467, 3, 1.25166555703457, 1.4287178895478,
                                      0.810972506211266, 5, 0.925820099772551, 1.50791562224467, 0,
                                      1.25166555703457, 1.49703106451365, 1.07055114141815, 6, 0.925820099772551,
                                      1.50791562224467, 2, 1.25166555703457, 1.54510894064586, 1.23298708781764,
                                      7, 0.925820099772551, 1.50791562224467, 4, 1.25166555703457,
                                      1.53986972817318, 1.34041735804632, 8, 0.925820099772551, 1.50791562224467,
                                      6, 1.25166555703457, 1.51867216206786, 1.41137323284577, 9, 0.925820099772551,
                                      1.1878488323372, 7, 1.11024271620619, 1.52424619129076, 1.45689201232324,
                                      10, 0.866025403784439, 1.1878488323372, 5, 1.11024271620619,
                                      1.5235523748015, 1.48489344540488, 11, 0.866025403784439, 1.72451684054791,
                                      7, 1.65123060170286, 1.48885399722157, 1.50186786343546, 12,
                                      1.64189930669738, 1.72451684054791, 9, 1.65123060170286, 1.44988865277158,
                                      1.51365121198181, 13, 1.64189930669738, 1.62117399689396, 13,
                                      1.55757359471666, 1.46250310836144, 1.52578416955645, 14, 1.53589529557661,
                                      1.43811745632331, 11, 1.37119119354531, 1.47211242437243, 1.54361561010075,
                                      15, 1.57249078621379, 1.91287503750007, 6, 1.58244501686629,
                                      1.40088808515813, 1.57221038715823, 16, 1.84390889145858, 1.91287503750007,
                                      1, 1.58244501686629, 1.22842056164903, 1.61611151602076, 17,
                                      1.84390889145858, 1.91287503750007, 2, 1.58244501686629, 0.961775139235799,
                                      1.67903900932887, 18, 1.84390889145858, 1.91287503750007, 1,
                                      1.58244501686629, 0.22985441247086, 1.76363892004825, 19, 1.84390889145858,
                                      1.91287503750007, 0, 1.58244501686629, 0.22985441247086, 1.87138377130906,
                                      20, 1.84390889145858, 1.91287503750007, 0, 1.58244501686629,
                                      0.22985441247086, 2.00265997698042, 21, 1.84390889145858, 1.91287503750007,
                                      1, 1.58244501686629, 0.22985441247086, 2.1569982806834, 22, 1.84390889145858))
})


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

results <- runAnalysis("standardErrorOfMeasurement", data = testthat::test_path("binaryTestDt.csv"), options = options, makeTests = F)


test_that("Sum score CI table results match", {
  table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(-0.0853758863090898, -0.0771097189955713, 0, 0, 0.0853758863090898,
                                        0.0771097189955713, 0, 0.91462411369091, 0.922890281004429, 0.888491541244813,
                                        1, 1.08537588630909, 1.07710971899557, 1.11150845875519, 1.87104202115792,
                                        1.87226439925586, 1.8634305870228, 2, 2.12895797884208, 2.12773560074414,
                                        2.1365694129772, 2.87297101344401, 2.87413111981302, 2.8634305870228,
                                        3, 3.12702898655599, 3.12586888018698, 3.1365694129772, 3.88445146704636,
                                        3.88566143267109, 3.88849154124481, 4, 4.11554853295364, 4.11433856732891,
                                        4.11150845875519, 5, 5, 5, 5, 5, 5, 5))
})

test_that("Sum score CI plots matches for ", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciPlotsContainer"]][["collection"]][["semMainContainer_ciPlotsContainer_anovaPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "anova-ci-plot")
})

test_that("Sum score CI plots matches for ", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciPlotsContainer"]][["collection"]][["semMainContainer_ciPlotsContainer_feldtPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "feldt-ci-plot")
})

test_that("Sum score CI plots matches for ", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciPlotsContainer"]][["collection"]][["semMainContainer_ciPlotsContainer_keatsPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "keats-ci-plot")
})

