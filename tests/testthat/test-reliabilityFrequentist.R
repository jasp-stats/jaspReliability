options <- analysisOptions("reliabilityFrequentist")
options$alphaItem <- TRUE
options$alphaScale <- TRUE
options$averageInterItemCor <- TRUE
options$confidenceIntervalValue <- 0.9
options$fitMeasures <- TRUE
options$glbItem <- TRUE
options$glbScale <- TRUE
options$guttman2Item <- TRUE
options$guttman2Scale <- TRUE
options$guttman6Item <- TRUE
options$guttman6Scale <- TRUE
options$itemRestCor <- TRUE
options$mcDonaldItem <- TRUE
options$meanItem <- TRUE
options$meanScale <- TRUE
options$noSamples <- 300
options$sdItem <- TRUE
options$sdScale <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive")
set.seed(1)
results <- runAnalysis("reliabilityFrequentist", "test.csv", options)


test_that("Fit Measures of Single Factor Model Fit table results match for main options", {
  table <- results[["results"]][["fitTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("", ""))
})

test_that("Frequentist Individual Item Reliability Statistics table results match for main options", {
  table <- results[["results"]][["itemTable"]][["data"]]
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

test_that("Frequentist Scale Reliability Statistics table results match for main options", {
  table <- results[["results"]][["scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.525637056655304, 0.448585203881531, 0.490572327059451, 0.516365401424283,
                                      0.567466299403832, 0.191748160936288, 0.7333720390425, 1.51568528273943,
                                      "Point estimate", 0.403771687499446, 0.279137827611668, 0.359901964004142,
                                      0.396402831936408, 0.459081486688799, 0.0950770018748302, "",
                                      "", "90% CI lower bound", 0.626298737852679, 0.586595705707241,
                                      0.60076050080959, 0.628395738756313, 0.685340010046143, 0.28002617012781,
                                      "", "", "90% CI upper bound"))
})

options <- analysisOptions("reliabilityFrequentist")
options$alphaInterval <- "alphaBoot"
options$alphaItem <- TRUE
options$alphaMethod <- "alphaStand"
options$alphaScale <- TRUE
options$averageInterItemCor <- TRUE
options$bootType <- "bootPara"
options$glbItem <- TRUE
options$glbScale <- TRUE
options$guttman2Item <- TRUE
options$guttman2Scale <- TRUE
options$guttman6Item <- TRUE
options$guttman6Scale <- TRUE
options$itemRestCor <- TRUE
options$mcDonaldItem <- TRUE
options$meanItem <- TRUE
options$meanScale <- TRUE
options$noSamples <- 300
options$reverseScaledItems <- "debMiss30"
options$sdItem <- TRUE
options$sdScale <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "debMiss30")
set.seed(1)
results <- runAnalysis("reliabilityFrequentist", "test.csv", options)


test_that("Frequentist Individual Item Reliability Statistics table results match for special options", {
  table <- results[["results"]][["itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0192598705478672, 0.498142461253454, 0.0162871794638591, 0.529785814609134,
                                      0.0265192931578479, 0.197671423943988, -0.18874858754, 1.05841360919316,
                                      "contNormal", 0.0309308512589049, 0.184853444016974, 0.02887162264849,
                                      0.142434129645294, 0.0490584701419857, 0.180747669167931, 0.05254867287,
                                      1.01183864387684, "contcor1", 0.0466790689222427, 0.332754536370897,
                                      0.0377062864873144, 0.252472824196733, 0.0611690101495298, 0.0513962438424752,
                                      0.06968807084, 1.0041493380131, "contcor2", 0.671633441261486,
                                      0.542545781005174, 0.600687137005109, 0.558313196445623, 0.685678770613101,
                                      0.122450817493202, 15.9882068024571, 24.0657052758223, "debMiss30"
                                 ))
})

test_that("Frequentist Scale Reliability Statistics table results match for special options", {
  table <- results[["results"]][["scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0471510241880036, 0.477845698575834, 0.0392336793613145, 0.517855238417142,
                                      0.0684849632061072, 0.184127369413486, 2.764059782725, 8.00605828917035,
                                      "Point estimate", 0.00702485789714698, 0.293267132435948, 0.00476660616888555,
                                      0.394302611317106, 0.0149501819485526, 0.0939902180995729, "",
                                      "", "95% CI lower bound", 0.115136402927549, 0.615996053946591,
                                      0.0878218264568898, 0.635557616289631, 0.216185154128633, 0.286242462411683,
                                      "", "", "95% CI upper bound"))
})
