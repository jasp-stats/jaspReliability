options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$itemDeletedAlpha <- TRUE
options$scaleAlpha <- TRUE
options$averageInterItemCorrelation <- TRUE
options$ciLevel <- 0.9
options$omegafitMeasures <- TRUE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$itemRestCorrelation <- TRUE
options$itemDeletedOmega <- TRUE
options$omegaEstimationMethod <- "pfa"
options$itemMean <- TRUE
options$scaleMean <- TRUE
options$bootstrapSamples <- 300
options$itemSd <- TRUE
options$scaleSd <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options, makeTests = FALSE)


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
                                 list("Coefficient <unicode>", 0.525637056655307, 0.403771687499453,
                                      "", 0.626298737852686, "Coefficient <unicode>", 0.448585203881531,
                                      0.279137827611669, "", 0.586595705707241, "Guttman's <unicode>2",
                                      0.490572327059451, 0.359901964004142, "", 0.60076050080959,
                                      "Average interitem correlation", 0.191748160936288, 0.0950770018748302,
                                      "", 0.280026170127809, "mean", 2.93348815617, 2.4742476519306,
                                      "", 3.3927286604094, "sd", 2.79198402042951, 2.50253823732555,
                                      "", 3.16486028196314))
})

# special options test
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$alphaIntervalMethod <- "bootstrapped"
options$itemDeletedAlpha <- TRUE
options$alphaType <- "standardized"
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
options$scaleSd <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "debMiss30")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options, makeTests = FALSE)

test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0192598705478673, 0.499115899922193, 0.0162871794638591, 0.197671423943988,
                                      -0.18874858754, 1.05841360919316, "contNormal", 0.030930851258905,
                                      0.173583399426978, 0.02887162264849, 0.180747669167931, 0.05254867287,
                                      1.01183864387684, "contcor1", 0.0466790689222428, 0.325209102320569,
                                      0.0377062864873144, 0.0513962438424752, 0.06968807084, 1.0041493380131,
                                      "contcor2", 0.671633441261487, 0.542545781005174, 0.60068713700511,
                                      0.122450817493202, 15.9882068024571, 24.0657052758223, "debMiss30"
                                 ))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.0471510241880039, 0.00702485789714753,
                                      "", 0.115136402927551, "Coefficient <unicode>", 0.474438267321141,
                                      0.293267132435947, "", 0.615996053946591, "Guttman's <unicode>2",
                                      0.0392336793613145, 0.0047666061688859, "", 0.0878218264568904,
                                      "Average interitem correlation", 0.184127369413486, 0.0939902180995728,
                                      "", 0.286242462411683, "mean", 2.764059782725, 1.68976265760682,
                                      "", 3.83835690784319, "sd", 5.48120849970767, 4.81253858382921,
                                      "", 6.36738428322095))
})


# omega test
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$omegaFitMeasures <- TRUE
options$bootstrapSamples <- 100
options$omegaIntervalMethod <- "bootstrapped"
options$naAction <- "listwise"
options$bootstrapType <- "parametric"
options$omegaEstimationMethod <- "cfa"
options$standardizedLoadings <- TRUE
options$setSeed <- TRUE
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
options$setSeed <- TRUE
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "asrm_mis.csv", options)

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
                                   list("McDonald's <unicode>", 0.7917101, 0.694692780642854,
                                        "", 0.85215986924174))
  } else if (jaspBase::getOS() == "osx") {
    jaspTools::expect_equal_tables(table,
                                   list("McDonald's <unicode>", 0.791710063361508, 0.688063057831179,
                                        "", 0.854457118896833))
  }
})

test_that("Standardized Loadings of the Single-Factor Model table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_loadTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.844967739355576, "asrm_1", 0.754980853055071, "asrm_2", 0.443805410137925,
                                      "asrm_3", 0.551535460946667, "asrm_4", 0.654599358675392, "asrm_5"
                                 ))
})



# disabled sample saving test
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$itemDeletedAlpha <- TRUE
options$scaleAlpha <- TRUE
options$averageInterItemCorrelation <- TRUE
options$ciLevel <- 0.9
options$omegafitMeasures <- TRUE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$itemRestCorrelation <- TRUE
options$itemDeletedOmega <- TRUE
options$omegaEstimationMethod <- "pfa"
options$itemMean <- TRUE
options$scaleMean <- TRUE
options$bootstrapSamples <- 300
options$itemSd <- TRUE
options$samplesSavingDisabled <- TRUE
options$scaleSd <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options, makeTests = FALSE)


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
                                 list("Coefficient <unicode>", 0.525637056655307, 0.403771687499453,
                                      "", 0.626298737852686, "Coefficient <unicode>", 0.448585203881531,
                                      0.279137827611669, "", 0.586595705707241, "Guttman's <unicode>2",
                                      0.490572327059451, 0.359901964004142, "", 0.60076050080959,
                                      "Average interitem correlation", 0.191748160936288, 0.0950770018748302,
                                      "", 0.280026170127809, "mean", 2.93348815617, 2.4742476519306,
                                      "", 3.3927286604094, "sd", 2.79198402042951, 2.50253823732555,
                                      "", 3.16486028196314))
})
