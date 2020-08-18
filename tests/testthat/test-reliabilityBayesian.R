options <- analysisOptions("reliabilityBayesian")
options$alphaItem <- TRUE
options$alphaScale <- TRUE
options$averageInterItemCor <- TRUE
options$credibleIntervalValueItem <- 0.9
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
options$probTable <- TRUE
options$probTableValueHigh <- 0.3
options$probTableValueLow <- 0.1
options$rHat <- TRUE
options$sdItem <- TRUE
options$sdScale <- TRUE
options$setSeed <- TRUE
options$variables <- c("contNormal", "contcor1", "contcor2", "facFive")
set.seed(1)
results <- runAnalysis("reliabilityBayesian", "test.csv", options)


test_that("Bayesian Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.427497950838586, 0.322839035775697, 0.415560596274979, 0.460126468167324,
                                      0.462217995813286, -0.00667393596077111, 0.536821057757968,
                                      0.478840471900155, 0.530636782588017, 0.557650604644619, 0.579809425144588,
                                      0.133506113650053, -0.18874858754, 1.05841360919316, 0.645317651287762,
                                      0.623910211428705, 0.641397904428437, 0.647355364352853, 0.674014454076118,
                                      0.29672201942606, "contNormal", 6.35226389505136e-06, -0.0323260135507699,
                                      0.0144431047172635, -0.0344804943993657, 0.0552617033353108,
                                      0.329095209105205, 0.162115516084883, 0.177735751878795, 0.218708835034426,
                                      0.142196330094696, 0.257988822303677, 0.460268175585853, 0.05254867287,
                                      1.01183864387684, 0.34119418887803, 0.419352791268417, 0.399845316350987,
                                      0.313449059076203, 0.455179333809157, 0.587613246058499, "contcor1",
                                      6.85535342629442e-10, 0.0719061425848915, 0.115551736246703,
                                      0.0716873154472548, 0.144873119849732, 0.215557278637012, 0.218405612580144,
                                      0.271006291731917, 0.293634196111037, 0.228048704209436, 0.326628062147507,
                                      0.363545921738566, 0.06968807084, 1.0041493380131, 0.387088136395753,
                                      0.489896555650716, 0.48120824685717, 0.402266814062928, 0.491063656585353,
                                      0.490416177139362, "contcor2", 0.513820677991841, 0.391046110981327,
                                      0.508002079092401, 0.468462071868989, 0.59672910974336, -0.015499152577849,
                                      0.623096869189241, 0.531550029641048, 0.605588492293822, 0.564545971550885,
                                      0.686909358231204, 0.135200755784571, 3, 1.4213381090374, 0.718234159722534,
                                      0.656828083588266, 0.698243447514041, 0.658511069529606, 0.778923851563923,
                                      0.297086393097092, "facFive"))
})

test_that("Probability that Reliability Statistic is Larger than 0.10 and Smaller than 0.30 table results match", {
  table <- results[["results"]][["probTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.00444444444444447, 0.263179999504362, "McDonald's <unicode>",
                                      0.0844444444444444, 0.212923338143927, "Cronbach's <unicode>",
                                      0.00444444444444447, 0.119459139776366, "Guttman's <unicode>2",
                                      0.00222222222222224, 0.107044445506913, "Guttman's <unicode>6",
                                      0, 0.0131635877156484, "Greatest Lower Bound"))
})

test_that("Bayesian Scale Reliability Statistics table results match", {
  table <- results[["results"]][["scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.479979214722957, 0.442328719603318, 0.498219292597923, 0.525924332046715,
                                      0.590636686593891, 0.19347485623366, 0.7333720390425, 1.51568528273943,
                                      "Posterior mean", 0.342879651825553, 0.254618334986826, 0.353936473332957,
                                      0.391309056840046, 0.461847741296067, 0.106359073588414, "",
                                      "", "95% CI lower bound", 0.61114368191617, 0.615128875906764,
                                      0.636878063449655, 0.636673491920215, 0.719984361008782, 0.307409602568225,
                                      "", "", "95% CI upper bound", 1.00710690801089, 1.0048425336041,
                                      1.00475109953701, 1.00873819260755, 1.00811400872571, 1.01104129919821,
                                      "", "", "R-hat"))
})