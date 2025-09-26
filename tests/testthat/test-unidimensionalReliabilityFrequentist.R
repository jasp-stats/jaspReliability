# analytic confidence interval
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$scaleOmega <- FALSE
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
options$itemMean <- TRUE
options$scaleMean <- TRUE
options$scaleVar <- TRUE
options$itemSd <- TRUE
options$scaleSd <- TRUE
options$setSeed <- TRUE
options$hiddenScaleThreshold <- 10
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", testthat::test_path("asrm.csv"), options, makeTests = F)


test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.590410772412584, 0.601267685347127, 0.499746429932828, 2.47877140395024,
                                      0.781401765957935, 0.722847443345244, 0.729317086131165, 0.649992231778526,
                                      2.67948717948718, 0.904441992620488, 0.855284114277903, 0.857366486915203,
                                      0.762266302693691, 2.88020295502412, 1.07383333365897, "asrm_1",
                                      0.607066121073447, 0.611182594693247, 0.475168884559224, 2.40584963440731,
                                      0.715913141676863, 0.731880726916123, 0.734730750315596, 0.630980009291163,
                                      2.58974358974359, 0.828641470508619, 0.856695332758799, 0.858278905937945,
                                      0.748422870116717, 2.77363754507987, 0.983836265835256, "asrm_2",
                                      0.697511924639158, 0.705239284978459, 0.233284698172959, 1.74991096107246,
                                      0.773970709653438, 0.793001084660615, 0.797208661700744, 0.433320479366847,
                                      1.94871794871795, 0.895840835489655, 0.888490244682072, 0.889178038423029,
                                      0.598171292118922, 2.14752493636344, 1.06362128102266, "asrm_3",
                                      0.644062273806467, 0.663149641165271, 0.336262234815853, 2.27675483504211,
                                      0.719376699106489, 0.765425655300394, 0.77443685066795, 0.519890282401405,
                                      2.46153846153846, 0.832650402814226, 0.886789036794321, 0.885724060170629,
                                      0.665435729129832, 2.64632208803482, 0.988596023981458, "asrm_4",
                                      0.591511112232894, 0.611205509225008, 0.46192077618739, 2.64223013295526,
                                      0.843801741639071, 0.73282849740501, 0.742662751599646, 0.620653571262307,
                                      2.85897435897436, 0.97666752473882, 0.874145882577127, 0.874119993974284,
                                      0.74086071394996, 3.07571858499345, 1.15958585793662, "asrm_5"
                                 ))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.789985434486178, 0.707653368279996,
                                      0.0500543421354605, 0.87231750069236, "Guttman's <unicode>2",
                                      0.796033564448564, 0.716973502479691, 0.048065104805344, 0.875093626417437,
                                      "Split-half coefficient", 0.78330036202316, 0.700984258147298,
                                      0.0500446377276892, 0.865616465899022, "Average interitem correlation",
                                      0.430268419211834, "", "", "", "Mean", 12.5384615384615, 11.9279578200137,
                                      0.371159906537889, 13.1489652569094, "Variance", 10.7452547452547,
                                      8.401175763325, 1.73175460560579, 14.3179701120402, "SD", 3.27799553771123,
                                      2.8984781805846, 0.273114903855726, 3.78390936889881))
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
options$hiddenScaleThreshold <- 10
options$variables <- c("contNormal", "contcor1", "contcor2", "debMiss30")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options, makeTests = F)

test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(3.71246492618525e-05, -0.032393035453925, -0.00527770593544932,
                                      -0.0093520914049856, -0.396193843016565, 0.863588070980427,
                                      0.929294394140214, 0.0192598705478673, 0.0144681758066313, 0.0162871794638591,
                                      0.198462015554176, -0.18874858754, 1.1202393681253, 1.05841360919316,
                                      0.075346314048224, 0.0573146999248539, 0.0608997002452966, 0.375215704411214,
                                      0.018696667936565, 1.51175115136297, 1.22953289966677, "contNormal",
                                      0.00154980754347149, -0.0124418579566558, -0.000420287856866242,
                                      -0.0492644454996841, -0.145768057146445, 0.789256791250214,
                                      0.888401255768031, 0.030930851258905, 0.0233203074068998, 0.02887162264849,
                                      0.175373623695169, 0.05254867287, 1.02381744124251, 1.01183864387684,
                                      0.0795418575844302, 0.0581030751697449, 0.0652065770036478,
                                      0.309692874988231, 0.250865402886445, 1.38163078322628, 1.17542791494259,
                                      "contcor1", 0.00518404408870834, -0.00124256598622769, 0.00365304316618711,
                                      -0.143295963248798, -0.127121582920542, 0.777306709422078, 0.881649992583269,
                                      0.0466790689222428, 0.0346056750405157, 0.0377062864873144,
                                      0.0481023174485831, 0.06968807084, 1.00831589303215, 1.0041493380131,
                                      0.0984160395516825, 0.0708166876811237, 0.0765991255472786,
                                      0.212626162715547, 0.266497724600542, 1.36071160825196, 1.16649543859029,
                                      "contcor2", 0.533395306905801, 0.316676004777833, 0.470440814190599,
                                      -0.061581356597092, 11.2714152421404, 446.470728863, 21.129853971644,
                                      0.671633441261487, 0.535041083185576, 0.60068713700511, 0.122450817493202,
                                      15.9882068024571, 579.15817042274, 24.0657052758223, 0.77012646283331,
                                      0.671487308870688, 0.696679032007608, 0.30241630495653, 20.7049983627739,
                                      781.567811192939, 27.9565343201359, "debMiss30"))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.0471510241880039, 0.00702485789714753,
                                      0.0281148884501927, 0.115136402927551, "Coefficient <unicode>",
                                      0.0338270034216142, -0.0094111158069592, 0.0232798957064179,
                                      0.0786531222994487, "Guttman's <unicode>2", 0.0392336793613145,
                                      0.0047666061688859, 0.0220147009193286, 0.0878218264568904,
                                      "Average interitem correlation", 0.184127369413486, 0.0939902180995728,
                                      0.0503714335248153, 0.286242462411683, "Mean", 2.764059782725,
                                      1.68976265760682, 0.548120849970767, 3.83835690784319, "Variance",
                                      30.0436466172676, 23.1605276208449, 4.27021798715751, 40.5435826102092,
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
options$hiddenScaleThreshold <- 10
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
options$setSeed <- TRUE
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", testthat::test_path("asrm_mis.csv"), options,
                       makeTests = F)

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
                                   list("Coefficient <unicode>", 0.791710063361508, 0.694692772995628,
                                        0.0405022504460226, 0.852159859484347))
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
options$hiddenScaleThreshold <- 10
options$variables <- c("asrm_1", "asrm_2", "asrm_3", "asrm_4")
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", testthat::test_path("asrm.csv"), options,
                       makeTests = F)

test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.512767157451441, 0.395912376152677, 0.416710118574818, 0.642387810117694,
                                      0.609462461310165, 0.617228009973447, 0.772008462783947, 0.823012546467653,
                                      0.817745901372076, "asrm_1", 0.501955170018076, 0.438382254987982,
                                      0.439733875765864, 0.640141520920559, 0.639830176986499, 0.639938573049249,
                                      0.778327871823042, 0.841278098985017, 0.840143270332634, "asrm_2",
                                      0.688686832630021, 0.629654786651993, 0.645198187231801, 0.769749240010078,
                                      0.751857116024424, 0.759115178994372, 0.850811647390136, 0.874059445396854,
                                      0.873032170756944, "asrm_3", "", 0.457779216045379, 0.514928488572061,
                                      "", 0.682331846089402, 0.704706497328503, "", 0.906884476133425,
                                      0.894484506084945, "asrm_4"))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Coefficient <unicode>", 0.743458875550109, 0.672822439743217,
                                      0.0464758999280928, 0.814095311357002, "Coefficient <unicode>",
                                      0.734183438809739, 0.618202085150078, 0.0705116563317664, 0.850164792469401,
                                      "Guttman's <unicode>2", 0.744181125723515, 0.636954210756928,
                                      0.0651893355187589, 0.851408040690102))
})



# check lambda2 analytic confidence interval
options <- analysisOptions("unidimensionalReliabilityFrequentist")
options$intervalMethod <- "analytic"
options$scaleOmega <- FALSE
options$itemDeletedLambda2 <- TRUE
options$scaleLambda2 <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2")
options$hiddenScaleThreshold <- 10
set.seed(1)
results <- runAnalysis("unidimensionalReliabilityFrequentist", "test.csv", options, makeTests = F)

test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.711442852097357, 0.79299280264282, 0.874542753188283, "contNormal",
                                      -0.307281370858232, 0.0618194975467093, 0.430920365951651, "contcor1",
                                      -0.00732988641121907, 0.277152727398941, 0.561635341209102,
                                      "contcor2"))
})

test_that("Frequentist Scale Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("Guttman's <unicode>2", 0.60068713700511, 0.483142517715655, 0.0599728465505652,
                                      0.718231756294564))
})
