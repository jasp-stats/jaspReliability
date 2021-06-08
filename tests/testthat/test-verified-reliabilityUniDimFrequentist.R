context("Unidimensional Reliability Frequentist -- Verification project")

options <- analysisOptions("reliabilityUniDimFrequentist")
options$omegaScale <- TRUE
options$alphaScale <- TRUE
options$lambda2Scale <- TRUE
options$lambda6Scale <- TRUE
options$meanScale <- TRUE
options$omegaItem <- TRUE
options$alphaItem <- TRUE
options$lambda2Item <- TRUE
options$lambda6Item <- TRUE
options$itemRestCor <- TRUE
options$meanItem <- TRUE
options$sdItem <- TRUE
options$omegaInterval <- "omegaAnalytic"
options$scoresMethod <- "meanScores"

options$variables <- c(paste("Question", c(1, 4:8), sep="_0"), 
                       paste("Question", 10, sep="_"))

set.seed(1)

result <- jaspTools::runAnalysis("reliabilityUniDimFrequentist", 
                                 "Reliability.csv", options)

test_that("Main (scale) table results match R, SPSS, SAS and MiniTab", { 
  tempResult <- result$results$stateContainer$collection$stateContainer_scaleTable
  resultTable <- tempResult$data
  
  jaspTools::expect_equal_tables(
    "test"=resultTable, 
    "ref"=list(0.762830469875067, 0.757381263545318, 0.764552908800385, 0.744068664329276,
               2.50730677335111, "Point estimate", 0.748803468173598, 0.742939526669925,
               0.750007365828022, 0.728936908337022, 2.48349947134054, "95% CI lower bound",
               0.776857471576537, 0.771189756266616, 0.778194560298748, 0.759434989460345,
               2.53111407536169, "95% CI upper bound"))
})

test_that("Main (item) table results match R, SPSS, SAS and MiniTab", { 
  tempResult <- result$results$stateContainer$collection$stateContainer_itemTable
  resultTable <- tempResult$data
  
  jaspTools::expect_equal_tables(
    "test"=resultTable, 
    "ref"=list(0.738318723638981, 0.72770167960325, 0.736546525283196, 0.702841729674047,
               0.482407830305298, 2.37417347335667, 0.828022142655599, "Question_01",
               0.725416778506326, 0.714226029693281, 0.724004514548663, 0.691632709918766,
               0.537847564692385, 2.78607545702061, 0.948548245925725, "Question_04",
               0.73585796382304, 0.72465527833956, 0.734980705918487, 0.704639145198906,
               0.489045648668466, 2.72228704784131, 0.964690420727732, "Question_05",
               0.735979687155172, 0.730026221541044, 0.735364197126947, 0.702040245237608,
               0.473324346566253, 2.22714896927266, 1.12200227531636, "Question_06",
               0.705632108445891, 0.704421830534527, 0.708784316995843, 0.681886282490119,
               0.572648610122144, 2.92376507195644, 1.10236003533831, "Question_07",
               0.747791822403911, 0.741268949281306, 0.749193530651701, 0.723243504328744,
               0.407453391591579, 2.23687281213536, 0.872570436301405, "Question_08",
               0.752113195750575, 0.748493326739936, 0.75495617610098, 0.731801981849852,
               0.369190204284374, 2.28082458187476, 0.877129349455171, "Question_10"
    ))
})