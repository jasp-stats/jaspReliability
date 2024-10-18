context("Unidimensional Reliability Frequentist -- Verification project")

options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$scaleOmega <- TRUE
options$scaleAlpha <- TRUE
options$scaleLambda2 <- TRUE
options$scaleMean <- TRUE
options$itemDeletedOmega <- TRUE
options$itemDeletedAlpha <- TRUE
options$itemDeletedLambda2 <- TRUE
options$itemRestCorrelation <- TRUE
options$itemMean <- TRUE
options$itemSd <- TRUE
options$omegaIntervalMethod <- "analytic"
options$meanSdScoresMethod <- "meanScores"

options$variables <- c(paste("Question", c(1, 4:8), sep="_0"),
                       paste("Question", 10, sep="_"))

set.seed(1)

results <- jaspTools::runAnalysis("unidimensionalReliabilityFrequentist", testthat::test_path("Reliability.csv"), options, makeTests = F)

test_that("Main (scale) table results match R, SPSS, SAS and MiniTab", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.762830469875068, 0.748803468173598,
                                      0.00715676502839475, 0.776857471576537, "Coefficient <unicode>",
                                      0.757381263545318, 0.742450705275325, 0.00761777174874793, 0.77231182181531,
                                      "Guttman's <unicode>2", 0.764552908800386, 0.750193337649084,
                                      0.0073264464370614, 0.778912479951688, "Mean", 2.50730677335111,
                                      2.48349947134054, 0.012146805858866, 2.53111407536169))
})

test_that("Main (item) table results match R, SPSS, SAS and MiniTab", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.73831872363898, 0.727701679603248, 0.736546525283195, 0.482407830305306,
                                      2.37417347335667, 0.828022142655586, "Question_01", 0.725416778506327,
                                      0.714226029693282, 0.724004514548664, 0.537847564692381, 2.78607545702061,
                                      0.94854824592574, "Question_04", 0.735857963823041, 0.724655278339561,
                                      0.734980705918487, 0.48904564866847, 2.72228704784131, 0.964690420727732,
                                      "Question_05", 0.735979687155174, 0.730026221541046, 0.735364197126949,
                                      0.473324346566251, 2.22714896927266, 1.12200227531637, "Question_06",
                                      0.705632108445892, 0.704421830534527, 0.708784316995843, 0.572648610122148,
                                      2.92376507195644, 1.10236003533831, "Question_07", 0.747791822403911,
                                      0.741268949281306, 0.749193530651702, 0.407453391591576, 2.23687281213536,
                                      0.872570436301406, "Question_08", 0.752113195750575, 0.748493326739936,
                                      0.75495617610098, 0.369190204284382, 2.28082458187476, 0.877129349455149,
                                      "Question_10"))
})
