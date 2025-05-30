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
results <- jaspTools::runAnalysis("unidimensionalReliabilityFrequentist",
                                  testthat::test_path("Reliability.csv"), options, makeTests = F)

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

test_that("Frequentist Individual Item Reliability Statistics table results match", {
  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.722669473705523, 0.710623104470516, 0.720142584684343, 0.452182921613983,
                                      2.34216691307561, 0.805993418667595, 0.73831872363898, 0.727701679603248,
                                      0.736546525283195, 0.482407830305287, 2.37417347335667, 0.828022142655586,
                                      0.753967973572438, 0.74478025473598, 0.752950465882046, 0.511526066766275,
                                      2.40618003363773, 0.851297794593357, "Question_01", 0.708956758905816,
                                      0.696200818345821, 0.706703449278422, 0.509789505745871, 2.74941005274127,
                                      0.923313042152352, 0.725416778506327, 0.714226029693282, 0.724004514548664,
                                      0.53784756469237, 2.78607545702061, 0.94854824592574, 0.741876798106837,
                                      0.732251241040743, 0.741305579818906, 0.564762630214891, 2.82274086129996,
                                      0.975211879276828, "Question_04", 0.72002148869523, 0.707408425935699,
                                      0.718536415916282, 0.459066914802536, 2.68499768018677, 0.939025770089386,
                                      0.735857963823041, 0.724655278339561, 0.734980705918487, 0.489045648668465,
                                      2.72228704784131, 0.964690420727732, 0.751694438950852, 0.741902130743423,
                                      0.751424995920693, 0.517911900604729, 2.75957641549584, 0.991807809628163,
                                      "Question_05", 0.720175760087626, 0.713554086936677, 0.719191916776892,
                                      0.442768393276605, 2.18377883264174, 1.09215249574698, 0.735979687155174,
                                      0.730026221541046, 0.735364197126949, 0.473324346566247, 2.22714896927266,
                                      1.12200227531637, 0.751783614222722, 0.746498356145415, 0.751536477477006,
                                      0.502782194317266, 2.27051910590357, 1.15354169085651, "Question_06",
                                      0.688053048555393, 0.685441617921256, 0.690274965780428, 0.546079746066424,
                                      2.88115419103322, 1.07303281846464, 0.705632108445892, 0.704421830534527,
                                      0.708784316995843, 0.572648610122145, 2.92376507195644, 1.10236003533831,
                                      0.72321116833639, 0.723402043147798, 0.727293668211258, 0.598066631097355,
                                      2.96637595287965, 1.13334730871044, "Question_07", 0.732655729480166,
                                      0.724887277530403, 0.733553557322497, 0.37469777923857, 2.20314427176157,
                                      0.849356548276965, 0.747791822403911, 0.741268949281306, 0.749193530651702,
                                      0.407453391591579, 2.23687281213536, 0.872570436301406, 0.762927915327655,
                                      0.757650621032208, 0.764833503980906, 0.439193128552637, 2.27060135250914,
                                      0.897098338057034, "Question_08", 0.737228558349435, 0.732686270853151,
                                      0.739809442164025, 0.335318331090338, 2.24691982021273, 0.853794175979056,
                                      0.752113195750575, 0.748493326739936, 0.75495617610098, 0.369190204284386,
                                      2.28082458187476, 0.877129349455149, 0.766997833151715, 0.76430038262672,
                                      0.770102910037935, 0.402108844508378, 2.31472934353679, 0.901785401981528,
                                      "Question_10"))
})

