options <- analysisOptions("reliabilityBayesian")
options$alphaItem <- TRUE
options$alphaScale <- TRUE
options$averageInterItemCor <- TRUE
options$credibleIntervalValueItem <- 0.9
options$dispPrior <- TRUE
options$fixXRange <- TRUE
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
options$noSamples <- 200
options$orderItem <- TRUE
options$orderType <- "orderItemKS"
options$plotItem <- TRUE
options$plotPosterior <- TRUE
options$probTable <- TRUE
options$probTableValueHigh <- 0.3
options$probTableValueLow <- 0.1
options$rHat <- TRUE
options$reverseScaledItems <- "debMiss30"
options$sdItem <- TRUE
options$sdScale <- TRUE
options$setSeed <- TRUE
options$shadePlots <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive", "debMiss30")
set.seed(1)
results <- runAnalysis("reliabilityBayesian", "debug.csv", options)


test_that("Bayesian Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(1.21695646331964e-07, -0.0302604620110933, -0.00730189174965356,
                                      0.381529087630189, 0.00468437245991029, 0.0452950595494772,
                                      0.0252487392158195, 0.0212125066975763, 0.0325973854277802,
                                      0.494218519800715, 0.0831973460090512, 0.200786958388514, -0.18874858754,
                                      1.05841360919316, 0.0669845350224338, 0.0716108498281969, 0.0783885886016254,
                                      0.606959485289979, 0.158707927788854, 0.353534950462945, "contNormal",
                                      1.67680115627664e-08, -0.0157617947104415, -0.00111243856322425,
                                      0.03077711847543, 0.00866952978122558, 0.0104533283359072, 0.0412635419994771,
                                      0.0263018982813707, 0.0385176798582972, 0.199614543969512, 0.0901500218057822,
                                      0.173831322632861, 0.05254867287, 1.01183864387684, 0.10970303862964,
                                      0.0684025775131465, 0.0730991379330716, 0.391472022288837, 0.162670671684959,
                                      0.315647290004467, "contcor1", 8.60724513641752e-07, -0.00463647527725275,
                                      0.00875019171674774, 0.129610929479583, 0.0225485982857244,
                                      -0.112744886614557, 0.0537890577165174, 0.0386674316516181,
                                      0.0486615947429129, 0.293085737284156, 0.102499450157703, 0.0568995666036369,
                                      0.06968807084, 1.0041493380131, 0.13453163652543, 0.0807291934771201,
                                      0.0876232547094318, 0.440233848836187, 0.194147501919935, 0.203466216010952,
                                      "contcor2", 8.83114909806779e-09, -0.0148054356280018, 0.0021066241823039,
                                      0.420553558691172, 0.0159597213180006, -0.0709816138759204,
                                      0.0312186169316103, 0.0308430722235896, 0.0408228091905911,
                                      0.520482384062624, 0.102023432380012, 0.0844325768836279, 3,
                                      1.4213381090374, 0.0807083859087435, 0.0718948260426243, 0.0786323358702756,
                                      0.628165183999919, 0.190816990258425, 0.255178425193561, "facFive",
                                      0.35942193964239, 0.285237073953381, 0.38010239162854, 0.423372619911139,
                                      0.47791304140055, -0.00395440079343337, 0.484237373873063, 0.434004318276764,
                                      0.491435548543508, 0.519513315848713, 0.585673590137766, 0.170405246301504,
                                      15.9882068024571, 24.0657052758223, 0.589811814300273, 0.570258196480784,
                                      0.610285976634783, 0.629418255267226, 0.701113938910745, 0.338540231176324,
                                      "debMiss30"))
})

test_that("Cronbach's α plot matches", {
  plotName <- results[["results"]][["plotContainer"]][["collection"]][["plotContainer_Cronbach's α"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "cronbach-s-α-scale", dir="reliabilityBayesian")
})

test_that("Greatest Lower Bound plot matches", {
  plotName <- results[["results"]][["plotContainer"]][["collection"]][["plotContainer_Greatest Lower Bound"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "greatest-lower-bound-scale", dir="reliabilityBayesian")
})

test_that("Guttman's λ2 plot matches", {
  plotName <- results[["results"]][["plotContainer"]][["collection"]][["plotContainer_Guttman's λ2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-λ2-scale", dir="reliabilityBayesian")
})

test_that("Guttman's λ6 plot matches", {
  plotName <- results[["results"]][["plotContainer"]][["collection"]][["plotContainer_Guttman's λ6"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-λ6-scale", dir="reliabilityBayesian")
})

test_that("McDonald's ω plot matches", {
  plotName <- results[["results"]][["plotContainer"]][["collection"]][["plotContainer_McDonald's ω"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "mcdonald-s-ω-scale", dir="reliabilityBayesian")
})

test_that("Cronbach's α plot matches", {
  plotName <- results[["results"]][["plotContainerItem"]][["collection"]][["plotContainerItem_Cronbach's α"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "cronbach-s-α-item", dir="reliabilityBayesian")
})

test_that("Greatest Lower Bound plot matches", {
  plotName <- results[["results"]][["plotContainerItem"]][["collection"]][["plotContainerItem_Greatest Lower Bound"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "greatest-lower-bound-item", dir="reliabilityBayesian")
})

test_that("Guttman's λ2 plot matches", {
  plotName <- results[["results"]][["plotContainerItem"]][["collection"]][["plotContainerItem_Guttman's λ2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-λ2-item", dir="reliabilityBayesian")
})

test_that("Guttman's λ6 plot matches", {
  plotName <- results[["results"]][["plotContainerItem"]][["collection"]][["plotContainerItem_Guttman's λ6"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-λ6-item", dir="reliabilityBayesian")
})

test_that("McDonald's ω plot matches", {
  plotName <- results[["results"]][["plotContainerItem"]][["collection"]][["plotContainerItem_McDonald's ω"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "mcdonald-s-ω-item", dir="reliabilityBayesian")
})

test_that("Probability that Reliability Statistic is Larger than 0.10 and Smaller than 0.30 table results match", {
  table <- results[["results"]][["probTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0688888888888889, 0.26017225362508, "McDonald's <unicode>",
                                      0.02, 0.220163381782161, "Cronbach's <unicode>", 0.0466666666666667,
                                      0.115986652096294, "Guttman's <unicode>2", 0, 0.0617103799372504,
                                      "Guttman's <unicode>6", 0.6, 0.00475933170216669, "Greatest Lower Bound"
                                 ))
})

test_that("Bayesian Scale Reliability Statistics table results match", {
  table <- results[["results"]][["scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0335139258203281, 0.0379993015242684, 0.0509632083185619, 0.507213118664045,
                                      0.125134713254144, 0.142232495514273, 2.8581975155295, 6.94729970125454,
                                      "Posterior mean", 7.26490380601197e-08, -0.0247593284030748,
                                      0.00211400244586846, 0.358233671861917, 0.0351173024406283,
                                      0.0612110191516504, "", "", "95% CI lower bound", 0.119049805227077,
                                      0.0925517739486723, 0.109789052574058, 0.618692688177314, 0.26368113845796,
                                      0.221618551302267, "", "", "95% CI upper bound", 1.02785492867472,
                                      1.02512927209581, 1.04150159018568, 1.02543279542889, 1.07445147891652,
                                      1.01625374834437, "", "", "R-hat"))
})