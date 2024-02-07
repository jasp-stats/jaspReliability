options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$itemDeletedAlpha <- TRUE
options$scaleAlpha <- TRUE
options$averageInterItemCorrelation <- TRUE
options$ciLevel <- 0.9
options$omegafitMeasures <- TRUE
options$itemDeletedGreatestLowerBound <- TRUE
options$scaleGreatestLowerBound <- TRUE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$itemDeletedLambda6 <- TRUE
options$scaleLambda6 <- TRUE
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
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options)


test_that("Frequentist Individual Item Reliability Statistics table results match for main options", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.55175372583945, 0.484048636403245, 0.525547257736439, 0.552088397696874,
                                      0.560933950469576, 0.133621852799609, -0.18874858754, 1.05841360919316,
                                      "contNormal", 0.230256036016286, 0.189422705634724, 0.196049691047581,
                                      0.134230003173562, 0.223967986152454, 0.458747451806099, 0.05254867287,
                                      1.01183864387684, "contcor1", 0.28356284370461, 0.282664399460191,
                                      0.283087797840773, 0.222670519008231, 0.284134273672075, 0.363943642284291,
                                      0.06968807084, 1.0041493380131, "contcor2", 0.671633441261486,
                                      0.535041083185576, 0.600687137005109, 0.558313196445623, 0.685678770613101,
                                      0.139002685382132, 3, 1.4213381090374, "facFive"))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.525637056655306, 0.448585203881531, 0.490572327059451, 0.516365401424283,
                                      0.567466299403832, 0.191748160936288, 2.93348815617, 2.79198402042951,
                                      "Point estimate", 0.403771687499453, 0.279137827611668, 0.359901964004142,
                                      0.396402831936408, 0.459081486688799, 0.0950770018748302, 2.4742476519306,
                                      2.50253823732555, "90% CI lower bound", 0.626298737852686, 0.586595705707241,
                                      0.60076050080959, 0.628395738756313, 0.685340010046143, 0.28002617012781,
                                      3.3927286604094, 3.16486028196314, "90% CI upper bound"))
})

options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$alphaIntervalMethod <- "bootstrapped"
options$itemDeletedAlpha <- TRUE
options$alphaType <- "standardized"
options$scaleAlpha <- TRUE
options$averageInterItemCorrelation <- TRUE
options$bootstrapType <- "parametric"
options$itemDeletedGreatestLowerBound <- TRUE
options$scaleGreatestLowerBound <- TRUE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$itemDeletedLambda6 <- TRUE
options$scaleLambda6 <- TRUE
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
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options)


test_that("Frequentist Individual Item Reliability Statistics table results match for special options", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0192598705478672, 0.499115899922193, 0.0162871794638591, 0.529785814609134,
                                      0.0265192931578479, 0.197671423943988, -0.18874858754, 1.05841360919316,
                                      "contNormal", 0.0309308512589049, 0.173583399426978, 0.02887162264849,
                                      0.142434129645294, 0.0490584701419857, 0.180747669167931, 0.05254867287,
                                      1.01183864387684, "contcor1", 0.0466790689222427, 0.325209102320569,
                                      0.0377062864873144, 0.252472824196733, 0.0611690101495298, 0.0513962438424752,
                                      0.06968807084, 1.0041493380131, "contcor2", 0.671633441261486,
                                      0.542545781005174, 0.600687137005109, 0.558313196445623, 0.685678770613101,
                                      0.122450817493202, 15.9882068024571, 24.0657052758223, "debMiss30"
                                 ))
})

test_that("Frequentist Scale Reliability Statistics table results match for special options", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]

  if (jaspBase::getOS() == "linux") {
    jaspTools::expect_equal_tables(table,
                                   list(0.0471510241880039, 0.474438267321141, 0.0392336793613145, 0.517855238417142,
                                        0.0684849632061198, 0.184127369413486, 2.764059782725, 5.48120849970767,
                                        "Point estimate", 0.00596163495689328, 0.289407838775642, 0.00364474587371632,
                                        0.407256236850348, 0.0170876010177963, 0.0924170105791601, 1.68976265760681,
                                        4.81253858382921, "95% CI lower bound", 0.115485828937115, 0.615613625191656,
                                        0.0900905447143638, 0.636260988220761, 0.220873107497465, 0.285913270464004,
                                        3.83835690784319, 6.36738428322095, "95% CI upper bound"))
  } else {
    jaspTools::expect_equal_tables(table,
                                   list(0.0471510241880039, 0.474438267321141, 0.0392336793613145, 0.517855238417142,
                                        0.0684849632061198, 0.184127369413486, 2.764059782725, 5.48120849970767,
                                        "Point estimate", 0.00702485789714729, 0.293267132435948, 0.00476660616888555,
                                        0.394302611317106, 0.0149501819485384, 0.0939902180995729, 1.68976265760681,
                                        4.81253858382921, "95% CI lower bound", 0.11513640292755, 0.61599605394659,
                                        0.0878218264568898, 0.635557616289631, 0.216185154128616, 0.286242462411683,
                                        3.83835690784319, 6.36738428322095, "95% CI upper bound"))
  }

})


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
                                   list(0.7917101, "Point estimate", 0.682155240979747, "95% CI lower bound",
                                        0.844326391786941, "95% CI upper bound"))
  } else if (jaspBase::getOS() == "osx") {
    jaspTools::expect_equal_tables(table,
                                   list(0.791710063361508, "Point estimate", 0.679125749864971, "95% CI lower bound",
                                        0.849422292711498, "95% CI upper bound"))
  }
})

test_that("Standardized Loadings of the Single-Factor Model table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_loadTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.844967739355576, "asrm_1", 0.754980853055071, "asrm_2", 0.443805410137925,
                                      "asrm_3", 0.551535460946667, "asrm_4", 0.654599358675392, "asrm_5"
                                 ))
})


options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$itemDeletedAlpha <- TRUE
options$scaleAlpha <- TRUE
options$averageInterItemCorrelation <- TRUE
options$ciLevel <- 0.9
options$omegafitMeasures <- TRUE
options$itemDeletedGreatestLowerBound <- TRUE
options$scaleGreatestLowerBound <- TRUE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$itemDeletedLambda6 <- TRUE
options$scaleLambda6 <- TRUE
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
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options)


test_that("Frequentist Individual Item Reliability Statistics table results match for main options with disabled sample saving", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.55175372583945, 0.484048636403245, 0.525547257736439, 0.552088397696874,
                                      0.560933950469576, 0.133621852799609, -0.18874858754, 1.05841360919316,
                                      "contNormal", 0.230256036016286, 0.189422705634724, 0.196049691047581,
                                      0.134230003173562, 0.223967986152454, 0.458747451806099, 0.05254867287,
                                      1.01183864387684, "contcor1", 0.28356284370461, 0.282664399460191,
                                      0.283087797840773, 0.222670519008231, 0.284134273672075, 0.363943642284291,
                                      0.06968807084, 1.0041493380131, "contcor2", 0.671633441261486,
                                      0.535041083185576, 0.600687137005109, 0.558313196445623, 0.685678770613101,
                                      0.139002685382132, 3, 1.4213381090374, "facFive"))
})

test_that("Frequentist Scale Reliability Statistics table results match with disabled sample saving", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.525637056655306, 0.448585203881531, 0.490572327059451, 0.516365401424283,
                                      0.567466299403832, 0.191748160936288, 2.93348815617, 2.79198402042951,
                                      "Point estimate", 0.403771687499453, 0.279137827611668, 0.359901964004142,
                                      0.396402831936408, 0.459081486688799, 0.0950770018748302, 2.4742476519306,
                                      2.50253823732555, "90% CI lower bound", 0.626298737852686, 0.586595705707241,
                                      0.60076050080959, 0.628395738756313, 0.685340010046143, 0.28002617012781,
                                      3.3927286604094, 3.16486028196314, "90% CI upper bound"))
})
