# analytic confidence interval
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$itemDeletedAlpha <- TRUE
options$scaleAlpha <- TRUE
options$averageInterItemCorrelation <- TRUE
options$ciLevel <- 0.9
options$omegafitMeasures <- TRUE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$scaleSplithalf <- TRUE
options$itemSplithalf <- TRUE
options$itemRestCorrelation <- TRUE
options$itemDeletedOmega <- TRUE
options$omegaEstimationMethod <- "pfa"
options$itemMean <- TRUE
options$scaleMean <- TRUE
options$scaleVar <- TRUE
options$bootstrapSamples <- 300
options$itemSd <- TRUE
options$scaleSd <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options, makeTests = F)


test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.55175372583945, 0.484048636403245, 0.525547257736439, 0.133621852799609,
                                      -0.18874858754, 1.05841360919316, "contNormal", 0.230256036016288,
                                      0.189422705634723, 0.196049691047581, 0.458747451806099, 0.05254867287,
                                      1.01183864387684, "contcor1", 0.28356284370461, 0.282664399460191,
                                      0.283087797840773, 0.363943642284291, 0.06968807084, 1.0041493380131,
                                      "contcor2", 0.671633441261487, 0.535041083185576, 0.60068713700511,
                                      0.139002685382132, 3, 1.4213381090374, "facFive"))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.525637056655307, "", "", "", "Coefficient <unicode>",
                                      0.448585203881531, "", "", "", "Guttman's <unicode>2", 0.490572327059451,
                                      "", "", "", "Split-half coefficient", 0.576970187390583, 0.435670122838468,
                                      0.0859043395940329, 0.7182702519427, "Average interitem correlation",
                                      0.191748160936288, "", "", "", "Mean", 2.93348815617, 2.4742476519306,
                                      0.279198402042951, 3.3927286604094, "Variance", 7.79517477033372,
                                      6.26269762927648, 1.16186834179699, 10.0163406043478, "SD",
                                      2.79198402042951, 2.50253823732555, 0.208072168983663, 3.16486028196314
                                 ))
})

# special options test
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$intervalMethod <- "bootstrapped"
options$itemDeletedAlpha <- TRUE
options$scaleAlpha <- TRUE
options$averageInterItemCorrelation <- TRUE
options$bootstrapType <- "parametric"
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$itemRestCorrelation <- TRUE
options$itemDeletedOmega <- TRUE
options$omegaEstimationMethod <- "pfa"
options$itemMean <- TRUE
options$scaleMean <- TRUE
options$meanSdScoresMethod <- "meanScores"
options$bootstrapSamples <- 300
options$reverseScaledItems <- "debMiss30"
options$itemSd <- TRUE
options$scaleVar <- TRUE
options$itemVar <- TRUE
options$scaleSd <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "debMiss30")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options, makeTests = F)

test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0192598705478673, 0.0144681758066313, 0.0162871794638591, 0.197671423943988,
                                      -0.18874858754, 1.1202393681253, 1.05841360919316, "contNormal",
                                      0.030930851258905, 0.0233203074068998, 0.02887162264849, 0.180747669167931,
                                      0.05254867287, 1.02381744124251, 1.01183864387684, "contcor1",
                                      0.0466790689222428, 0.0346056750405157, 0.0377062864873144,
                                      0.0513962438424752, 0.06968807084, 1.00831589303215, 1.0041493380131,
                                      "contcor2", 0.671633441261487, 0.535041083185576, 0.60068713700511,
                                      0.122450817493202, 15.9882068024571, 579.15817042274, 24.0657052758223,
                                      "debMiss30"))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.0471510241880039, "", "", "", "Coefficient <unicode>",
                                      0.0338270034216142, -0.0094111158069592, 0.0232798957064179,
                                      0.0786531222994487, "Guttman's <unicode>2", 0.0392336793613145,
                                      0.0047666061688859, 0.0220147009193286, 0.0878218264568904,
                                      "Average interitem correlation", 0.184127369413486, 0.0939902180995728,
                                      0.0503714335248153, 0.286242462411683, "Mean", 2.764059782725,
                                      1.68976265760682, 0.548120849970767, 3.83835690784319, "Variance",
                                      30.0436466172676, 23.1605276208449, 4.3735671591732, 40.5435826102092,
                                      "SD", 5.48120849970767, 4.81253858382921, 0.398960116131913,
                                      6.36738428322095))
})


# omega test
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$omegaFitMeasures <- TRUE
options$bootstrapSamples <- 100
options$intervalMethod <- "bootstrapped"
options$naAction <- "listwise"
options$bootstrapType <- "parametric"
options$omegaEstimationMethod <- "cfa"
options$standardizedLoadings <- TRUE
options$setSeed <- TRUE
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
options$setSeed <- TRUE
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", testthat::test_path("asrm_mis.csv"), options, makeTests = F)

test_that("Fit Measures of Single Factor Model Fit table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_fitTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Chi-Square", 12.788508304247, "df", 5, "p.value", 0.0254433712709828,
                                      "RMSEA", 0.15724319758923, "Lower 90% CI RMSEA", 0.0507316506074521,
                                      "Upper 90% CI RMSEA", 0.266560548199575, "SRMR", 0.0708026289801857
                                 ))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  if (jaspBase::getOS() == "linux") {
    jaspTools::expect_equal_tables(table,
                                   list("Coefficient <unicode>", 0.7917101, 0.694692780642854,
                                        0.0428578933018329, 0.85215986924174))
  } else if (jaspBase::getOS() == "osx") {
    jaspTools::expect_equal_tables(table,
                                   list("Coefficient <unicode>", 0.791710063361508, 0.688063057831179,
                                        0.0428578933018329, 0.854457118896833))
  }
})

test_that("Standardized Loadings of the Single-Factor Model table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_loadTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.844967739355576, "asrm_1", 0.754980853055071, "asrm_2", 0.443805410137925,
                                      "asrm_3", 0.551535460946667, "asrm_4", 0.654599358675392, "asrm_5"
                                 ))
})


# standardized coefficients work
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$itemDeletedAlpha <- TRUE
options$scaleAlpha <- TRUE
options$ciLevel <- 0.9
options$itemDeletedGreatestLowerBound <- TRUE
options$scaleGreatestLowerBound <- TRUE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$itemDeletedLambda6 <- TRUE
options$scaleLambda6 <- TRUE
options$itemDeletedOmega <- TRUE
options$omegaEstimationMethod <- "cfa"
options$coefficientType <- "standardized"
options$setSeed <- TRUE
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", testthat::test_path("asrm.csv"), options, makeTests = F)

test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("", 0.609462461310165, 0.617228009973447, "asrm_1", "", 0.639830176986499,
                                      0.639938573049249, "asrm_2", "", 0.751857116024424, 0.759115178994372,
                                      "asrm_3", "", 0.682331846089402, 0.704706497328503, "asrm_4"
                                 ))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.743458875550109, 0.672822439743217,
                                      0.0464758999280928, 0.814095311357002, "Coefficient <unicode>",
                                      0.734183438809739, "", "", "", "Guttman's <unicode>2", 0.744181125723515,
                                      "", "", ""))
})
