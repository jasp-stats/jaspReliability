
# common options
options <- analysisOptions("unidimensionalReliabilityBayesian")
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive", "debMiss30")
options$scaleAlpha <- TRUE
options$scaleLambda2 <- TRUE
options$averageInterItemCorrelation <- TRUE
options$scaleSplithalf <- TRUE
options$scaleMean <- TRUE
options$scaleVar <- TRUE
options$scaleSd <- TRUE
options$meanSdScoresMethod <- "meanScores"
options$itemCiLevel <- 0.9
options$itemDeletedOmega <- TRUE
options$itemDeletedAlpha <- TRUE
options$itemDeletedLambda2 <- TRUE
options$itemDeletedPlot <- TRUE
options$itemDeletedPlotOrdered <- TRUE
options$itemRestCorrelation <- TRUE
options$itemDeletedSplithalf <- TRUE
options$itemMean <- TRUE
options$itemVar <- TRUE
options$itemSd <- TRUE
options$posteriorPlot <- TRUE
options$posteriorPlotFixedRange <- TRUE
options$posteriorPlotPriorDisplayed <- TRUE
options$probabilityTable <- TRUE
options$probabilityTableLowerBound <- 0.1
options$probabilityTableUpperBound <- 0.3
options$posteriorPlotShaded <- TRUE
options$samples <- 200
options$rHat <- TRUE
options$effectiveSampleSize <- TRUE
options$tracePlot <- TRUE
options$setSeed <- TRUE
options$reverseScaledItems <- "debMiss30"
options$itemDeletedPlotOrderedType <- "kullbackLeibler"
options$inverseWishartPriorDf <- length(options$variables)
options$inverseWishartPriorScale <- 0.0000000001
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityBayesian", "test.csv", options, makeTests = F)

test_that("Bayesian Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0000049262740258844, -0.0302604620110933, -0.00730189174965356,
                                      0.44201519518913, 0.0494291189911074, 0.0322661722934692, 0.0212125066975763,
                                      0.0325973854277802, 0.586113047921276, 0.201506484849513, -0.18874858754,
                                      1.1202393681253, 1.05841360919316, 0.0908764111289078, 0.071610849828197,
                                      0.0783885886016254, 0.7160906715891, 0.347636993674243, "contNormal",
                                      0.000000940024640628503, -0.0157617947104415, -0.00111243856322423,
                                      0.0184364938164071, 0.0250257087054518, 0.0372377719573931,
                                      0.0263018982813707, 0.0385176798582973, 0.252887867184234, 0.18043738105305,
                                      0.05254867287, 1.02381744124251, 1.01183864387684, 0.101904192043364,
                                      0.0684025775131465, 0.0730991379330716, 0.532343590284936, 0.311403152449827,
                                      "contcor1", 0.0000000387289002980385, -0.00463647527725275,
                                      0.00875019171674773, 0.119568735823148, -0.0805010870340148,
                                      0.0472986075910507, 0.038667431651618, 0.0486615947429129, 0.354628587010249,
                                      0.0552126591050449, 0.06968807084, 1.00831589303215, 1.0041493380131,
                                      0.109551905583045, 0.0807291934771201, 0.0876232547094318, 0.538305160532101,
                                      0.228317875033916, "contcor2", 0.0000000306076329374839, -0.0148054356280015,
                                      0.00210662418230391, 0.509396711247221, -0.0786484764858476,
                                      0.0315765590294508, 0.0308430722235896, 0.040822809190591, 0.618684356999201,
                                      0.0698072970377176, 3, 2.02020202020202, 1.4213381090374, 0.0878435386563505,
                                      0.0718948260426243, 0.0786323358702756, 0.741107989525121, 0.204256118691893,
                                      "facFive", 0.358369568041053, 0.285237073953381, 0.38010239162854,
                                      0.506580145691714, 0.0164539097397541, 0.479869257412775, 0.434004318276763,
                                      0.491435548543508, 0.626447969042049, 0.184196159162177, 15.9882068024571,
                                      579.15817042274, 24.0657052758223, 0.600166559993091, 0.570258196480785,
                                      0.610285976634783, 0.718853492788103, 0.34352987256172, "debMiss30"
                                 ))
})

test_that("Coefficient alpha plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "coefficient-alpha-scale")
})

test_that("Guttman's lambda2 scale plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-lambda2-scale")
})

test_that("Coefficient omega plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "coefficient-omega-scale")
})

test_that("Split-half coefficient plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_splithalf"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "split-half-coefficient-tp")
})

test_that("Coefficient alpha plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "coefficient-alpha-item")
})

test_that("Guttman's lambda2 plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-lambda2-item")
})

test_that("Coefficient omega plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "coefficient-omega-item")
})

test_that("Split-half coefficient plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_splithalf"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "split-half-coefficient-item")
})

test_that("Coefficient alpha plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "coefficient-alpha-tp")
})

test_that("Guttman's lambda2 plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-lambda2-tp")
})

test_that("Coefficient omega plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerTP"]][["collection"]][["stateContainer_plotContainerTP_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "coefficient-omega-tp")
})

test_that("Split-half coefficient plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_splithalf"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "split-half-coefficient-scale")
})

test_that("Probability that Reliability Coefficient is Larger than 0.10 and Smaller than 0.30 table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_probabilityTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0577777777777778, 0.257579474218955, "Coefficient <unicode>",
                                      0.02, 0.215729634402404, "Coefficient <unicode>", 0.0466666666666667,
                                      0.117470293201849, "Guttman's <unicode>2", 0.00444444444444447,
                                      0.192241844381343, "Split-half coefficient"))
})

test_that("Bayesian Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 261.905719143696, 0.0335378719395442,
                                      0.000000488588161068121, 1.03110267965238, 0.107247682851598,
                                      "Coefficient <unicode>", 278.717853339024, 0.0379993015242684,
                                      -0.0247593284030748, 1.02512927209581, 0.0925517739486723, "Guttman's <unicode>2",
                                      271.566230782654, 0.0509632083185619, 0.00211400244586848, 1.04150159018568,
                                      0.109789052574058, "Split-half coefficient", 411.30231862117,
                                      0.580365168313184, 0.437096067658011, 1.01119370237162, 0.740706655357837,
                                      "Average interitem correlation", 414.536145329942, 0.142232495514273,
                                      0.0612110191516504, 1.01625374834437, 0.221618551302267, "Mean",
                                      "", 2.8581975155295, "", "", "", "Variance", "", 19.2886479310878,
                                      "", "", "", "SD", "", 4.39188432578635, "", "", ""))
})


# disabled sample saving
options <- analysisOptions("unidimensionalReliabilityBayesian")
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive", "debMiss30")
options$scaleAlpha <- TRUE
options$scaleLambda2 <- TRUE
options$samples <- 200
options$samplesSavingDisabled <- TRUE
options$setSeed <- TRUE
options$reverseScaledItems <- "debMiss30"
options$inverseWishartPriorDf <- length(options$variables)
options$inverseWishartPriorScale <- 0.0000000001

set.seed(1)
results <- runAnalysis("unidimensionalReliabilityBayesian", "test.csv", options, makeTests = F)
test_that("Bayesian Scale Reliability Statistics table results match", {
table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
jaspTools::expect_equal_tables(table,
                               list("Coefficient <unicode>", 0.0335378719395442, 4.88588161068121e-07,
                                    0.107247682851598, "Coefficient <unicode>", 0.0379993015242684,
                                    -0.0247593284030748, 0.0925517739486723, "Guttman's <unicode>2",
                                    0.0509632083185619, 0.00211400244586848, 0.109789052574058
                               ))
})


# adjusted priors
options <- analysisOptions("unidimensionalReliabilityBayesian")
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
options$scaleOmega <- TRUE
options$scaleLambda2 <- TRUE
options$averageInterItemCorrelation <- TRUE
options$itemCiLevel <- 0.95
options$itemDeletedOmega <- TRUE
options$itemDeletedLambda2 <- TRUE
options$itemRestCorrelation <- TRUE
options$posteriorPlot <- TRUE
options$posteriorPlotFixedRange <- TRUE
options$posteriorPlotPriorDisplayed <- TRUE
options$probabilityTable <- TRUE
options$probabilityTableLowerBound <- 0.7
options$probabilityTableUpperBound <- 1
options$samples <- 100
options$chains <- 2
options$rHat <- TRUE
options$setSeed <- TRUE
options$naAction <- "listwise"
options$inverseWishartPriorScale <- 1
options$inverseWishartPriorDf <- 10
options$inverseGammaPriorShape <- 6
options$inverseGammaPriorScale <- 10
options$normalPriorMean <- 1

set.seed(1)
results <- runAnalysis("unidimensionalReliabilityBayesian", testthat::test_path("asrm_mis.csv"), options, makeTests =F)

test_that("Bayesian Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.365317706663268, 0.629664746841127, 0.569474629223085, 0.547833261943788,
                                      0.713418617285292, 0.715574360243638, 0.719358502963356, 0.800896152983193,
                                      0.808079872344535, "asrm_1", 0.425180238104353, 0.647546515358842,
                                      0.506244667189157, 0.590585184223742, 0.737906643927935, 0.626413941162412,
                                      0.720617721441406, 0.821130476744635, 0.754537273415984, "asrm_2",
                                      0.511558346241041, 0.710291628124685, 0.326897431951586, 0.647299797106082,
                                      0.784623165070472, 0.492987161374292, 0.800023243314042, 0.850735328744899,
                                      0.646935311378515, "asrm_3", 0.419641220251401, 0.687663875830053,
                                      0.431363154771682, 0.615730018474419, 0.767034290352854, 0.580923015620642,
                                      0.72516307294578, 0.849444065720693, 0.716515367604935, "asrm_4",
                                      0.380908166253114, 0.624933985182541, 0.553082199033924, 0.559120671647617,
                                      0.733553360969622, 0.658218071150869, 0.717678239451181, 0.81305056069536,
                                      0.779528963734584, "asrm_5"))
})

test_that("Guttman's lambda2 plot matches with adjusted priors", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "lambda2-prior-adjusted")
})

test_that("McDonald's omega plot matches with adjusted priors", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "omega-prior-adjusted")
})

test_that("Probability that Reliability Statistic is Larger than 0.70 and Smaller than 1.00 table results match with adjusted priors", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_probabilityTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.3, 0.521665798540962, "McDonald's <unicode>", 0.99, 0.103837082022117,
                                      "Guttman's <unicode>2"))
})


test_that("Bayesian Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.658799634973158, 0.526011816625903,
                                      1.02778087840674, 0.795894536283031, "Guttman's <unicode>2",
                                      0.789270505536385, 0.71521555858729, 0.990343839523617, 0.859660650371794,
                                      "Average interitem correlation", 0.419218581777162, 0.316743509774531,
                                      0.997596520647971, 0.507230342509293))
})



# standardization and median
options <- analysisOptions("unidimensionalReliabilityBayesian")
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive")
options$scaleAlpha <- TRUE
options$scaleLambda2 <- TRUE
options$averageInterItemCorrelation <- TRUE
options$itemDeletedOmega <- TRUE
options$itemDeletedAlpha <- TRUE
options$itemDeletedLambda2 <- TRUE
options$posteriorPlot <- TRUE
options$posteriorPlotPriorDisplayed <- TRUE
options$probabilityTable <- TRUE
options$probabilityTableLowerBound <- 0.2
options$probabilityTableUpperBound <- 0.5
options$posteriorPlotShaded <- TRUE
options$samples <- 200
options$rHat <- TRUE
options$setSeed <- TRUE
options$inverseWishartPriorDf <- length(options$variables)
options$inverseWishartPriorScale <- 0.0000000001
options$standardizedLoadings <- TRUE
options$coefficientType <- "standardized"
options$pointEstimate <- "median"
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityBayesian", "test.csv", options, makeTests = FALSE)

test_that("Bayesian Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.499864243624905, 0.371987972410805, 0.519686705794261, 0.620670635416665,
                                      0.553348700422337, 0.620992269207978, 0.714612686739472, 0.674415156445438,
                                      0.710080098725572, "contNormal", 5.73312727316968e-06, -0.10711429087776,
                                      -0.00316547939066864, 0.161564148144393, 0.194804130030221,
                                      0.220807457791194, 0.372255405615577, 0.445366740296484, 0.439504615186215,
                                      "contcor1", 0.0393067987346507, 0.0189905734187845, 0.0792132799731725,
                                      0.280634498409812, 0.299360676376119, 0.314435040862574, 0.456896489326536,
                                      0.49846757129663, 0.494379320017436, "contcor2", 0.518994671551128,
                                      0.381555713653778, 0.523617983089419, 0.627501070070124, 0.542402992030407,
                                      0.614551863581106, 0.714884024479502, 0.692929487440639, 0.707153922148802,
                                      "facFive"))
})

test_that("Standardized Loadings of the Single-Factor Model table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_loadTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.118960690052957, "contNormal", 0.808474151944779, "contcor1",
                                      0.768612766587743, "contcor2", 0.144763227189446, "facFive"
                                 ))
})

test_that("Coefficient alpha plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
<<<<<<< HEAD
  jaspTools::expect_equal_plots(testPlot, "alpha-std")
})

=======
  jaspTools::expect_equal_plots(testPlot, "coefficient-alpha-std")
})
>>>>>>> b809c4f (remove lambda6 and glb)

test_that("Guttman's lambda2 plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_lambda2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "guttman-s-lambda2-std")
})

test_that("Coefficient omega plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainer"]][["collection"]][["stateContainer_plotContainer_omega"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "coefficient-omega-std")
})

test_that("Probability that Reliability Coefficient is Larger than 0.20 and Smaller than 0.50 table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_probabilityTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.211111111111111, 0.322257519253291, "Coefficient <unicode>",
                                      0.546666666666667, 0.327288686204119, "Coefficient <unicode>",
                                      0.206666666666667, 0.233840603713808, "Guttman's <unicode>2"
                                 ))
})

test_that("Bayesian Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.549475923084279, 0.435706048978678,
                                      1.00879071329335, 0.664754350432929, "Coefficient <unicode>",
                                      0.485978025107501, 0.32252569270118, 1.00688014001917, 0.639694091225906,
                                      "Guttman's <unicode>2", 0.553414961378168, 0.443001980608279,
                                      1.00846765229577, 0.669397876270056, "Average interitem correlation",
                                      0.191174460794261, 0.106359073588414, 1.01104129919821, 0.307409602568224
                                 ))
})


# fit indices
# results were compared to blavFitIndices and lavaan omegaFitMeasures with the same data but 2000 obs
options <- analysisOptions("unidimensionalReliabilityBayesian")
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
options$scaleOmega <- TRUE
options$samples <- 200
options$chains <- 3
options$setSeed <- TRUE
options$omegaFitMeasures <- TRUE
options$omegaFitMeasuresCutoffRmsea <- .1
options$omegaFitMeasuresCutoffCfiTli <- .85
options$omegaPosteriorPredictiveCheck <- TRUE

set.seed(1)
results <- runAnalysis("unidimensionalReliabilityBayesian", testthat::test_path("asrm.csv"), options)

test_that("Fit Measures for the Single-Factor Model table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_fitTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.929544182303969, "Posterior mean", 13.0900880096089, 0.1319747009136,
                                      0.061860112775025, 0.860237253334787, 0.845602177064211, "90% CI lower bound",
                                      "", 0.0233139770705626, "", 0.693722044278387, 1, "90% CI upper bound",
                                      "", 0.222796365776379, "", 1, 0.946666666666667, "Relative to cutoff",
                                      "", 0.215555555555556, "", 0.595555555555556))
})

test_that("Posterior Predictive Check Omega plot matches", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_omegaPPC"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "posterior-predictive-check-omega")
})


# item plot ordered by mean
options <- analysisOptions("unidimensionalReliabilityBayesian")
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive", "debMiss30")
options$scaleAlpha <- TRUE
options$itemDeletedAlpha <- TRUE
options$itemDeletedPlot <- TRUE
options$itemDeletedPlotOrdered <- TRUE
options$itemDeletedPlotOrderedType <- "mean"
options$inverseWishartPriorDf <- length(options$variables)
options$inverseWishartPriorScale <- 0.0000000001
options$samples <- 200
options$chains <- 2
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityBayesian", "test.csv", options)

test_that("Cronbach's alpha plot item deleted plot matches ordered by mean", {
  plotName <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_plotContainerItem"]][["collection"]][["stateContainer_plotContainerItem_alpha"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "cronbach-s-alpha-item-mean-ordered")
})
