options <- analysisOptions("reliabilityUniDimBayesian")
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive", "debMiss30")
options$alphaScale <- TRUE
options$lambda2Scale <- TRUE
options$lambda6Scale <- TRUE
options$glbScale <- TRUE
options$averageInterItemCor <- TRUE
options$meanScale <- TRUE
options$sdScale <- TRUE
options$scoresMethod <- "meanScores"
options$credibleIntervalValueItem <- 0.9
options$omegaItem <- TRUE
options$alphaItem <- TRUE
options$lambda2Item <- TRUE
options$lambda6Item <- TRUE
options$glbItem <- TRUE
options$plotItem <- TRUE
options$orderItem <- TRUE
options$itemRestCor <- TRUE
options$meanItem <- TRUE
options$sdItem <- TRUE
options$plotPosterior <- TRUE
options$fixXRange <- TRUE
options$dispPrior <- TRUE
options$probTable <- TRUE
options$probTableValueLow <- 0.1
options$probTableValueHigh <- 0.3
options$shadePlots <- TRUE
options$noSamples <- 200
options$rHat <- TRUE
options$tracePlot <- TRUE
options$setSeed <- TRUE
options$reverseScaledItems <- "debMiss30"
options$orderType <- "orderItemKL"
set.seed(1)
results <- runAnalysis("reliabilityUniDimBayesian", "test.csv", options)

test_that("Bayesian Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(4.9262740258841e-06, -0.0302604620110933, -0.00730189174965356,
                                      0.381529087630189, 0.00468437245991029, 0.0403066583678939,
                                      0.0322661722934692, 0.0212125066975763, 0.0325973854277802,
                                      0.494218519800715, 0.0831973460090512, 0.198533853285123, -0.18874858754,
                                      1.05841360919316, 0.0908764111289075, 0.0716108498281969, 0.0783885886016254,
                                      0.606959485289979, 0.158707927788854, 0.343183173277219, "contNormal",
                                      9.40024640626554e-07, -0.0157617947104415, -0.00111243856322425,
                                      0.0307771184754301, 0.00866952978122558, 0.0233738672906997,
                                      0.0372377719573929, 0.0263018982813707, 0.0385176798582972,
                                      0.199614543969512, 0.0901500218057822, 0.181125240163059, 0.05254867287,
                                      1.01183864387684, 0.101904192043363, 0.0684025775131465, 0.0730991379330716,
                                      0.391472022288837, 0.162670671684959, 0.314266959992938, "contcor1",
                                      3.87289002980396e-08, -0.00463647527725275, 0.00875019171674774,
                                      0.129610929479583, 0.0225485982857244, -0.0833355720929837,
                                      0.0472986075910507, 0.0386674316516181, 0.0486615947429129,
                                      0.293085737284156, 0.102499450157703, 0.0586411881810892, 0.06968807084,
                                      1.0041493380131, 0.109551905583045, 0.0807291934771201, 0.0876232547094318,
                                      0.440233848836187, 0.194147501919935, 0.230464183178138, "contcor2",
                                      3.06076329375633e-08, -0.0148054356280015, 0.0021066241823039,
                                      0.420553558691172, 0.0159597213180006, -0.0799600621489579,
                                      0.0315765590294507, 0.0308430722235896, 0.0408228091905911,
                                      0.520482384062624, 0.102023432380012, 0.0706298883887492, 3,
                                      1.4213381090374, 0.0878435386563505, 0.0718948260426243, 0.0786323358702756,
                                      0.628165183999919, 0.190816990258425, 0.207264946462154, "facFive",
                                      0.358369568041053, 0.285237073953381, 0.38010239162854, 0.423372619911139,
                                      0.47791304140055, 0.013013921988528, 0.479869257412775, 0.434004318276764,
                                      0.491435548543508, 0.519513315848713, 0.585673590137766, 0.183357847129134,
                                      15.9882068024571, 24.0657052758223, 0.60016655999309, 0.570258196480784,
                                      0.610285976634783, 0.629418255267226, 0.701113938910745, 0.345204213665747,
                                      "debMiss30"))
})

test_that("Cronbach's alpha scale plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "alpha-scale", dir="uniDimBayesian")
})

test_that("Greatest Lower Bound scale plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_glb"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "greatest-lower-bound-scale", dir="uniDimBayesian")
})

test_that("Guttman's lambda2 scale plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lambda2-scale", dir="uniDimBayesian")
})

test_that("Guttman's lambda6 scale plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_lambda6"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lambda6-scale", dir="uniDimBayesian")
})

test_that("McDonald's omega scale plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "omega-scale", dir="uniDimBayesian")
})

test_that("Cronbach's alpha item plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "alpha-item", dir="uniDimBayesian")
})

test_that("Greatest Lower Bound item plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_glb"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "greatest-lower-bound-item", dir="uniDimBayesian")
})

test_that("Guttman's lambda2 item plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lambda2-item", dir="uniDimBayesian")
})

test_that("Guttman's lambda6 item plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_lambda6"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lambda6-item", dir="uniDimBayesian")
})

test_that("McDonald's omega item plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "omega-item", dir="uniDimBayesian")
})

test_that("Cronbach's alpha traceplot plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "alpha-tp", dir="uniDimBayesian")
})

test_that("Greatest Lower Bound traceplot plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_glb"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "greatest-lower-bound-tp", dir="uniDimBayesian")
})

test_that("Guttman's lambda2 traceplot plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lambda2-tp", dir="uniDimBayesian")
})

test_that("Guttman's lambda6 traceplot plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_lambda6"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lambda6-tp", dir="uniDimBayesian")
})

test_that("McDonald's omega traceplot plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "omega-tp", dir="uniDimBayesian")
})

test_that("Probability that Reliability Statistic is Larger than 0.10 and Smaller than 0.30 table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_probTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0577777777777778, 0.26017225362508, "McDonald's <unicode>",
                                      0.02, 0.220163381782161, "Cronbach's <unicode>", 0.0466666666666667,
                                      0.115986652096294, "Guttman's <unicode>2", 0, 0.0617103799372504,
                                      "Guttman's <unicode>6", 0.6, 0.00475933170216669, "Greatest Lower Bound"
                                 ))
})

test_that("Bayesian Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0335378719395441, 0.0379993015242684, 0.0509632083185619, 0.507213118664045,
                                      0.125134713254144, 0.142232495514273, 2.8581975155295, 4.39188432578635,
                                      "Posterior mean", 4.88588161069134e-07, -0.0247593284030748,
                                      0.00211400244586846, 0.358233671861917, 0.0351173024406283,
                                      0.0612110191516504, "", "", "95% CI lower bound", 0.107247682851598,
                                      0.0925517739486723, 0.109789052574058, 0.618692688177314, 0.26368113845796,
                                      0.221618551302267, "", "", "95% CI upper bound", 1.03110267965238,
                                      1.02512927209581, 1.04150159018568, 1.02543279542889, 1.07445147891652,
                                      1.01625374834437, "", "", "R-hat"))
})



options <- analysisOptions("reliabilityUniDimBayesian")
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive", "debMiss30")
options$alphaScale <- TRUE
options$lambda2Scale <- TRUE
options$lambda6Scale <- TRUE
options$glbScale <- TRUE
options$averageInterItemCor <- TRUE
options$meanScale <- TRUE
options$sdScale <- TRUE
options$scoresMethod <- "meanScores"
options$credibleIntervalValueItem <- 0.9
options$omegaItem <- TRUE
options$alphaItem <- TRUE
options$lambda2Item <- TRUE
options$lambda6Item <- TRUE
options$glbItem <- TRUE
options$orderItem <- TRUE
options$itemRestCor <- TRUE
options$meanItem <- TRUE
options$sdItem <- TRUE
options$probTable <- TRUE
options$probTableValueLow <- 0.1
options$probTableValueHigh <- 0.3
options$shadePlots <- TRUE
options$noSamples <- 200
options$rHat <- TRUE
options$disableSampleSave <- TRUE
options$setSeed <- TRUE
options$reverseScaledItems <- "debMiss30"
set.seed(1)
results <- runAnalysis("reliabilityUniDimBayesian", "test.csv", options)

test_that("Bayesian Individual Item Reliability Statistics table results match with disabled sample saving", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(4.9262740258841e-06, -0.0302604620110933, -0.00730189174965356,
                                      0.381529087630189, 0.00468437245991029, 0.0403066583678939,
                                      0.0322661722934692, 0.0212125066975763, 0.0325973854277802,
                                      0.494218519800715, 0.0831973460090512, 0.198533853285123, -0.18874858754,
                                      1.05841360919316, 0.0908764111289075, 0.0716108498281969, 0.0783885886016254,
                                      0.606959485289979, 0.158707927788854, 0.343183173277219, "contNormal",
                                      9.40024640626554e-07, -0.0157617947104415, -0.00111243856322425,
                                      0.0307771184754301, 0.00866952978122558, 0.0233738672906997,
                                      0.0372377719573929, 0.0263018982813707, 0.0385176798582972,
                                      0.199614543969512, 0.0901500218057822, 0.181125240163059, 0.05254867287,
                                      1.01183864387684, 0.101904192043363, 0.0684025775131465, 0.0730991379330716,
                                      0.391472022288837, 0.162670671684959, 0.314266959992938, "contcor1",
                                      3.87289002980396e-08, -0.00463647527725275, 0.00875019171674774,
                                      0.129610929479583, 0.0225485982857244, -0.0833355720929837,
                                      0.0472986075910507, 0.0386674316516181, 0.0486615947429129,
                                      0.293085737284156, 0.102499450157703, 0.0586411881810892, 0.06968807084,
                                      1.0041493380131, 0.109551905583045, 0.0807291934771201, 0.0876232547094318,
                                      0.440233848836187, 0.194147501919935, 0.230464183178138, "contcor2",
                                      3.06076329375633e-08, -0.0148054356280015, 0.0021066241823039,
                                      0.420553558691172, 0.0159597213180006, -0.0799600621489579,
                                      0.0315765590294507, 0.0308430722235896, 0.0408228091905911,
                                      0.520482384062624, 0.102023432380012, 0.0706298883887492, 3,
                                      1.4213381090374, 0.0878435386563505, 0.0718948260426243, 0.0786323358702756,
                                      0.628165183999919, 0.190816990258425, 0.207264946462154, "facFive",
                                      0.358369568041053, 0.285237073953381, 0.38010239162854, 0.423372619911139,
                                      0.47791304140055, 0.013013921988528, 0.479869257412775, 0.434004318276764,
                                      0.491435548543508, 0.519513315848713, 0.585673590137766, 0.183357847129134,
                                      15.9882068024571, 24.0657052758223, 0.60016655999309, 0.570258196480784,
                                      0.610285976634783, 0.629418255267226, 0.701113938910745, 0.345204213665747,
                                      "debMiss30"))
})

test_that("Probability that Reliability Statistic is Larger than 0.10 and Smaller than 0.30 table results match with disabled sample saving", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_probTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0577777777777778, 0.26017225362508, "McDonald's <unicode>",
                                      0.02, 0.220163381782161, "Cronbach's <unicode>", 0.0466666666666667,
                                      0.115986652096294, "Guttman's <unicode>2", 0, 0.0617103799372504,
                                      "Guttman's <unicode>6", 0.6, 0.00475933170216669, "Greatest Lower Bound"
                                 ))
})

test_that("Bayesian Scale Reliability Statistics table results match with disabled sample saving", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0335378719395441, 0.0379993015242684, 0.0509632083185619, 0.507213118664045,
                                      0.125134713254144, 0.142232495514273, 2.8581975155295, 4.39188432578635,
                                      "Posterior mean", 4.88588161069134e-07, -0.0247593284030748,
                                      0.00211400244586846, 0.358233671861917, 0.0351173024406283,
                                      0.0612110191516504, "", "", "95% CI lower bound", 0.107247682851598,
                                      0.0925517739486723, 0.109789052574058, 0.618692688177314, 0.26368113845796,
                                      0.221618551302267, "", "", "95% CI upper bound", 1.03110267965238,
                                      1.02512927209581, 1.04150159018568, 1.02543279542889, 1.07445147891652,
                                      1.01625374834437, "", "", "R-hat"))
})
