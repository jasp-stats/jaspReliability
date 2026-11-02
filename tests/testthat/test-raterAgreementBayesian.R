# ==== Bayesian Cohen's kappa and Fleiss' kappa, pairwise, with observed/chance agreement and posterior plots ====

options <- analysisOptions("raterAgreementBayesian")
options$variables <- c("V1", "V2", "V3")
options$cohensKappa <- TRUE
options$fleissKappa <- TRUE
options$observedAndChanceAgreement <- TRUE
options$posteriorPlot <- TRUE
options$setSeed <- TRUE
set.seed(1)
results <- runAnalysis("raterAgreementBayesian", testthat::test_path("binaryTestDt.csv"), options)

test_that("Bayesian Cohen's kappa table results match", {
  table <- results[["results"]][["cohensKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.497309168071882, 0.101456951312119, 0.2248625405901, 0.610314099076576,
                                      "V1 - V2", 0.361389077347582, 0.549270719208197, 0.00831861661275699,
                                      0.101477744164802, 0.594926522649455, "V1 - V3", 0.212873342164137,
                                      0.48120823704866, -0.0152519056512471, 0.065599743656698, 0.515139967514635,
                                      "V2 - V3", 0.153836357074565))
})

test_that("Bayesian Fleiss' kappa table results match", {
  table <- results[["results"]][["fleissKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.502247078591049, 0.0847726200041091, 0.217112915919341, 0.610314099076576,
                                      "V1 - V2", 0.350679562118738, 0.599705993446189, -0.143144301852982,
                                      -0.011904994892567, 0.594926522649455, "V1 - V3", 0.135149049441701,
                                      0.563938371006704, -0.23736733795919, -0.111789382618539, 0.515139967514635,
                                      "V2 - V3", 0.023679609250933))
})

test_that("Bayesian Cohen's kappa posterior plot matches", {
  plotName <- results[["results"]][["cohensKappaPosteriorPlots"]][["collection"]][["cohensKappaPosteriorPlots_pair1"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "cohens-kappa-posterior-v1-v2")
})

test_that("Bayesian Fleiss' kappa posterior plot matches", {
  plotName <- results[["results"]][["fleissKappaPosteriorPlots"]][["collection"]][["fleissKappaPosteriorPlots_pair1"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "fleiss-kappa-posterior-v1-v2")
})

# with a negligible prior the posterior mean is close to the classical estimates:
# Cohen's kappa (psych) and, for two raters, Fleiss' kappa (irr)
test_that("Posterior means are close to the classical kappas", {
  ratings   <- read.csv(testthat::test_path("binaryTestDt.csv"))[, c("V1", "V2")]
  cohen     <- psych::cohen.kappa(ratings)$kappa
  fleiss    <- irr::kappam.fleiss(ratings)$value
  cohenRow  <- results[["results"]][["cohensKappa"]][["data"]][[1]]
  fleissRow <- results[["results"]][["fleissKappa"]][["data"]][[1]]
  expect_equal(cohenRow[["mean"]], cohen, tolerance = 0.01)
  expect_equal(fleissRow[["mean"]], fleiss, tolerance = 0.01)
})


# ==== Weighted Cohen's kappa: declared level order is used, not alphabetical label order ====
test_that("Bayesian weighted Cohen's kappa respects ordered factor levels", {
  lv <- c("low", "medium", "high") # alphabetical order would be high < low < medium
  df <- data.frame(
    r1 = factor(c("low", "medium", "high", "low", "medium", "high", "low", "high"), levels = lv, ordered = TRUE),
    r2 = factor(c("low", "high", "medium", "low", "medium", "high", "medium", "high"), levels = lv, ordered = TRUE)
  )
  options <- analysisOptions("raterAgreementBayesian")
  options$variables       <- c("r1", "r2")
  options$variables.types <- c("ordinal", "ordinal")
  options$cohensKappa     <- TRUE
  options$cohensKappaType <- "weighted"
  options$weightType      <- "quadratic"
  options$setSeed         <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreementBayesian", df, options)
  table <- results[["results"]][["cohensKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.354036681660039, 0.683649818542856, "r1 - r2", 0.949148468430387))
})

test_that("Bayesian weighted Cohen's kappa requires ordinal variables", {
  options <- analysisOptions("raterAgreementBayesian")
  options$variables       <- c("V1", "V2")
  options$variables.types <- c("nominal", "nominal")
  options$cohensKappa     <- TRUE
  options$cohensKappaType <- "weighted"
  results <- runAnalysis("raterAgreementBayesian", testthat::test_path("binaryTestDt.csv"), options)
  expect_identical(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
                   "Weighted Cohen's kappa requires ordinal variables. Remove nominal variables or change their type.")
})


# ==== Missing values: every pair uses its own complete cases and reports its n ====
test_that("Bayesian Cohen's kappa uses pairwise complete cases", {
  df <- data.frame(r1 = c("a", "b", "a", "b", "a", NA, "b", "a"),
                   r2 = c("a", "b", "b", "b", NA, "a", "b", "a"),
                   r3 = c("a", "a", "a", "b", "a", "b", NA, "a"))
  options <- analysisOptions("raterAgreementBayesian")
  options$variables   <- c("r1", "r2", "r3")
  options$cohensKappa <- TRUE
  options$ci          <- FALSE
  options$setSeed     <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreementBayesian", df, options)
  table <- results[["results"]][["cohensKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.646968065111668, 6, "r1 - r2", 0.537569858350333, 6, "r1 - r3",
                                      0.0105306190761147, 6, "r2 - r3"))
})

test_that("Bayesian Cohen's kappa errors cleanly when rater pairs have no overlap", {
  df <- data.frame(r1 = c("a", "b", "a", NA, NA, NA), r2 = c(NA, NA, NA, "a", "b", "a"))
  options <- analysisOptions("raterAgreementBayesian")
  options$variables   <- c("r1", "r2")
  options$cohensKappa <- TRUE
  results <- runAnalysis("raterAgreementBayesian", df, options)
  expect_identical(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
                   "Cohen's kappa could not be computed for any rater pair: fewer than 3 jointly rated subjects/items.")
})


# ==== A pair whose ratings do not vary is skipped with a footnote; the seed is set per pair ====
test_that("Bayesian Fleiss' kappa skips pairs whose ratings do not vary", {
  df <- data.frame(r1 = c("a", "b", "a", "b", "a"),
                   r2 = c("a", "a", "a", "a", "a"),
                   r3 = c("a", "a", "a", "a", "a"))
  options <- analysisOptions("raterAgreementBayesian")
  options$variables   <- c("r1", "r2", "r3")
  options$fleissKappa <- TRUE
  options$setSeed     <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreementBayesian", df, options)
  table <- results[["results"]][["fleissKappa"]][["data"]]
  # r2 and r3 are identical, so with the seed set per pair r1 - r2 and r1 - r3 have identical posteriors
  jaspTools::expect_equal_tables(table,
                                 list(-0.599865167476499, -0.272496200226045, "r1 - r2", -0.0135617605383569,
                                      -0.599865167476499, -0.272496200226045, "r1 - r3", -0.0135617605383569,
                                      "", "", "r2 - r3", ""))
  footnotes <- vapply(results[["results"]][["fleissKappa"]][["footnotes"]], `[[`, character(1L), "text")
  expect_true("Some rater pairs could not be computed: r2 - r3 (the ratings do not vary)." %in% footnotes)
})


# ==== Raters in rows give the same posteriors as raters in columns ====
test_that("Raters in rows match raters in columns", {
  ratings <- read.csv(testthat::test_path("binaryTestDt.csv"))[, c("V1", "V2", "V3")]
  colMode <- as.data.frame(lapply(ratings, factor, levels = c(0, 1)))
  rowMode <- as.data.frame(t(ratings))
  rowMode[] <- lapply(rowMode, factor, levels = c(0, 1))

  options <- analysisOptions("raterAgreementBayesian")
  options$cohensKappa <- TRUE
  options$setSeed     <- TRUE

  options$variables <- colnames(colMode)
  resultsColumns    <- runAnalysis("raterAgreementBayesian", colMode, options)
  options$variables     <- colnames(rowMode)
  options$dataStructure <- "ratersInRows"
  resultsRows           <- runAnalysis("raterAgreementBayesian", rowMode, options)

  estimates <- function(results) lapply(results[["results"]][["cohensKappa"]][["data"]], function(row) unlist(row[c("mean", "lower", "upper")]))
  expect_equal(estimates(resultsRows), estimates(resultsColumns))
})
