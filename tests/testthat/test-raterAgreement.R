# ==== Ensure results are unchanged on JASP debug data ====

####Cohen's unweighted kappa and Fleiss' kappa and Krippendorff's alpha####

# Set options
options <- analysisOptions("raterAgreement")
options$variables <- paste0("V", 1:5)
options$dataStructure <- "ratersInColumns"
options$setSeed <- TRUE
options$bootstrapSamples <- 200
options$fleissKappa <- TRUE
options$krippendorffsAlpha <- TRUE
options$cohensKappa <- TRUE
set.seed(1)
results <- runAnalysis("raterAgreement", testthat::test_path("binaryTestDt.csv"), options, makeTests = F)


test_that("Cohen's kappa table results match", {
  # pairs are enumerated in combn() order (V1-V2, V1-V3, V1-V4, ...); the kappas are
  # unchanged, each pair is computed on its own pairwise-complete data
  table <- results[["results"]][["cohensKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("", "", "", 0.114984433765595, "Average kappa", 0.0924679737811612,
                                      0.357606541867125, 0.0676386326935962, 0.225037257824143, "V1 - V2",
                                      -0.00101486089827026, 0.203411665159256, 0.052150582273454,
                                      0.101198402130493, "V1 - V3", 0.0576342926736923, 0.309887992795108,
                                      0.0643516161804914, 0.1837611427344, "V1 - V4", -0.0240223485028134,
                                      0.221356704625165, 0.0625978474766621, 0.0986671780611757, "V1 - V5",
                                      -0.0200019787491091, 0.151023173354119, 0.0436296670378263,
                                      0.0655105973025049, "V2 - V3", -0.0882775571852617, 0.1399093155876,
                                      0.0582120065911339, 0.025815879201169, "V2 - V4", 0.0762065929175585,
                                      0.337621571800286, 0.0666887200338209, 0.206914082358922, "V2 - V5",
                                      -0.0325092696700619, 0.257559504097391, 0.0739984959048937,
                                      0.112525117213664, "V3 - V4", -0.0317077880237179, 0.0969460355453043,
                                      0.0328204560348627, 0.0326191237607932, "V3 - V5", 0.00907597751537416,
                                      0.186515136621999, 0.0452659233807975, 0.0977955570686868, "V4 - V5"))
})

test_that("Fleiss' kappa table results match", {
  table <- results[["results"]][["fleissKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0321990251552644, 0.119851279212922, 0.0223606797749979, 0.0760251521840934,
                                      "Overall", 0.0321990251552635, 0.119851279212922, 0.0223606797749979,
                                      0.0760251521840926, 0, 0.0321990251552635, 0.119851279212922, 0.0223606797749979,
                                      0.0760251521840925, 1))
})

test_that("Krippendorff's alpha table results match", {
  table <- results[["results"]][["krippendorffsAlpha"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0271069334114449, 0.124379082494376, 0.0257751838160253, 0.0762561458960466,
                                      "Nominal"))
})


#### Cohen's weighted kappa and Fleiss' kappa and Krippendorf's alpha with different CI range(99%)####

# Set options
options <- analysisOptions("raterAgreement")
options$variables <- c("facGender", "facExperim", "debBinMiss20")
options$ciLevel <- 0.99
options$cohensKappaType <- "weighted"
options$bootstrapSamples <- 200
options$dataStructure <- "ratersInColumns"
options$setSeed <- TRUE
options$fleissKappa <- TRUE
options$krippendorffsAlpha <- TRUE
options$cohensKappa <- TRUE
set.seed(1)
results <- runAnalysis("raterAgreement", "test.csv", options, makeTests = F)


test_that("Weighted Cohen's kappa refuses raters without a common ordinal scale", {
  # facGender {f,m}, facExperim {control,experimental} and debBinMiss20 {0,1} have
  # disjoint categories, so no common ordinal scale exists and the distances that
  # weighted kappa needs are undefined
  expect_match(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
})

test_that("Fleiss' kappa table results match", {
  table <- results[["results"]][["fleissKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(-0.276327852798572, -0.127177258822542, 0.0289519561273982, -0.201752555810557,
                                      "Overall", -0.360298917334901, -0.0277607841576364, 0.0645497224367903,
                                      -0.194029850746269, "f", -0.372299217342401, -0.0397610841651366,
                                      0.0645497224367903, -0.206030150753769, "m", -0.384543178263759,
                                      -0.0520050450864946, 0.0645497224367903, -0.218274111675127, "control",
                                      -0.348535076440849, -0.0159969432635844, 0.0645497224367903,
                                      -0.182266009852217, "experimental", -0.337000773905705, -0.00446264072844088,
                                      0.0645497224367903, -0.170731707317073, 0, -0.397038297357863,
                                      -0.0645001641805985, 0.0645497224367903, -0.230769230769231, 1))
})

test_that("Krippendorff's alpha table results match", {
  table <- results[["results"]][["krippendorffsAlpha"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(-0.209716461340596, -0.197507041090973, 0.00269939839460176, -0.199079048349962,
                                      "Nominal"))
})


# ==== Verify results of Fleiss' kappa with data set from Fleiss (1971) NOTE: Only verifying kappa values, not CIs====
test_that("Fleiss' kappa table results match", {
  options <- analysisOptions("raterAgreement")
  options$variables <- c("V1", "V2", "V3", "V4", "V5", "V6")
  options$dataStructure <- "ratersInColumns"
  options$fleissKappa <- TRUE
  options$ci <- FALSE
  set.seed(1)
  results <- runAnalysis("raterAgreement", testthat::test_path("Fleiss1971.csv"), options)
  table <- results[["results"]][["fleissKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.430244520060141, "Overall", 0.244755244755245, 1, 0.244755244755245,
                                      2, 0.52, 3, 0.471127272727273, 4, 0.566117806823969, 5))
})


#### Kendall's W ####

options <- analysisOptions("raterAgreement")
options$variables                    <- c("contNormal", "contGamma", "contcor1")
options$dataStructure                <- "ratersInColumns"
options$cohensKappa                  <- FALSE
options$fleissKappa                  <- FALSE
options$krippendorffsAlpha           <- FALSE
options$kendallW                     <- TRUE
options$correctForTies               <- FALSE # tie correction is on by default; keep testing the uncorrected path
options$bootstrapSamples             <- 200
options$setSeed                      <- TRUE
set.seed(1)
results <- runAnalysis("raterAgreement", "debug.csv", options, makeTests = F)

test_that("Kendall's W table results match", {
  # CI requested but tie correction off: no CI columns, explanatory footnote instead
  table <- results[["results"]][["kendallW"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.926742171157389, 0.316646331299797, 94.0439603960396, 99,
                                      98.3333333333333, 196.666666666667, 0.621972149059366, 0.660376274732495))
  footnotes <- sapply(results[["results"]][["kendallW"]][["footnotes"]], `[[`, "text")
  expect_true(any(grepl("Bootstrap CIs are only available with the tie correction", footnotes)))
})

test_that("Kendall's W with tie correction and no CI results match", {
  options2 <- analysisOptions("raterAgreement")
  options2$variables          <- c("contNormal", "contGamma", "contcor1", "contcor2")
  options2$dataStructure      <- "ratersInColumns"
  options2$cohensKappa        <- FALSE
  options2$fleissKappa        <- FALSE
  options2$krippendorffsAlpha <- FALSE
  options2$kendallW           <- TRUE
  options2$correctForTies     <- TRUE
  options2$ci                 <- FALSE
  set.seed(1)
  results2 <- runAnalysis("raterAgreement", "debug.csv", options2)
  table <- results2[["results"]][["kendallW"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list(1.37810466030173, 0.31477197719772, 124.649702970297, 99, 98.5, 295.5,
         0.0416297450235598, 0.0214991915973132))
})


# ==== Verify Kendall's W against published reference (DescTools anxiety dataset) ====
# Source: DescTools::KendallW documentation, 3 raters x 20 subjects, 1-6 scale (with ties)
# Expected: W = 0.5397, chi2 = 30.76, df = 19, p = 0.04288
test_that("Kendall's W matches DescTools reference (anxiety ratings, tie correction)", {
  options <- analysisOptions("raterAgreement")
  options$variables       <- c("rater1", "rater2", "rater3")
  options$dataStructure   <- "ratersInColumns"
  options$kendallW        <- TRUE
  options$correctForTies  <- TRUE
  options$ci              <- FALSE
  results <- runAnalysis("raterAgreement", testthat::test_path("anxietyRatings.csv"), options)
  table <- results[["results"]][["kendallW"]][["data"]]
  # F test cross-validated against vegan::kendall.global (Legendre, 2005)
  jaspTools::expect_equal_tables(table,
    list(2.34458536585366, 0.53965687595437, 30.7604419293991, 19, 18.3333333333333,
         36.6666666666667, 0.0428834731269479, 0.013806204775388))
})

test_that("Cohen's kappa table results match with linear weighting", {
  options <- analysisOptions("raterAgreement")
  options$variables <- c("V1", "V2")
  options$fleissKappa <- FALSE
  options$krippendorffsAlpha <- FALSE
  options$cohensKappa <- TRUE
  options$ci <- FALSE
  options$cohensKappaType <- "weighted"
  options$dataStructure <- "ratersInColumns"
  options$weightType <- "linear"
  set.seed(1)
  results <- runAnalysis("raterAgreement", testthat::test_path("Fleiss1971.csv"), options)
  table <- results[["results"]][["cohensKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.633093525179856, "Average kappa", 0.633093525179856, "V1 - V2"
                                 ))
})


# ==== Ordered factors: declared level order must be used, not alphabetical label order ====
test_that("Kendall's W and Krippendorff's alpha respect ordered factor levels", {
  lv <- c("low", "medium", "high") # alphabetical order would be high < low < medium
  df <- data.frame(
    r1 = factor(c("low", "medium", "high", "low",    "medium", "high"), levels = lv, ordered = TRUE),
    r2 = factor(c("low", "high",   "medium", "low",  "medium", "high"), levels = lv, ordered = TRUE),
    r3 = factor(c("medium", "medium", "high", "low", "low",    "high"), levels = lv, ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables                <- c("r1", "r2", "r3")
  options$variables.types          <- c("ordinal", "ordinal", "ordinal")
  options$dataStructure            <- "ratersInColumns"
  options$kendallW                 <- TRUE
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "ordinal"
  options$ci                       <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # reference: irr::kendall / irr::kripp.alpha on the level codes (1 = low, 2 = medium, 3 = high)
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(7, 0.777777777777778, 11.6666666666667, 5, 4.33333333333333, 8.66666666666667,
         0.0396519759960316, 0.00776230190368293))
  jaspTools::expect_equal_tables(results[["results"]][["krippendorffsAlpha"]][["data"]],
    list(0.675925925925926, "Ordinal"))
})

# ==== Ties without correction: result matches irr and the ties warning is shown ====
test_that("Uncorrected Kendall's W on tied data matches irr and warns about ties", {
  options <- analysisOptions("raterAgreement")
  options$variables      <- c("rater1", "rater2", "rater3")
  options$dataStructure  <- "ratersInColumns"
  options$kendallW       <- TRUE
  options$correctForTies <- FALSE
  options$ci             <- FALSE
  results <- runAnalysis("raterAgreement", testthat::test_path("anxietyRatings.csv"), options)
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(2.01543106340154, 0.501921470342523, 28.6095238095238, 19, 18.3333333333333,
         36.6666666666667, 0.0723803546937571, 0.0349217756891712))
  footnotes <- sapply(results[["results"]][["kendallW"]][["footnotes"]], `[[`, "text")
  expect_true(any(grepl("Ties are present", footnotes)))
})

# ==== Sparse incomplete data: bootstrap must not abort the analysis ====
test_that("Kendall's W bootstrap CI handles incomplete rows and failed replicates", {
  df <- data.frame(
    r1 = c(1, 2, 3, NA),
    r2 = c(1, 3, 2, 4),
    r3 = c(2, 1, 3, 4)
  )
  options <- analysisOptions("raterAgreement")
  options$variables        <- c("r1", "r2", "r3")
  options$dataStructure    <- "ratersInColumns"
  options$kendallW         <- TRUE
  options$ci               <- TRUE
  options$bootstrapSamples <- 100
  options$setSeed          <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(0.111111111111111, 1, 1.6, 0.350241068029438, 0.444444444444444,
         2.66666666666667, 2, 1.33333333333333, 2.66666666666667,
         0.263597138115727, 0.325271609713338))
  footnotes <- sapply(results[["results"]][["kendallW"]][["footnotes"]], `[[`, "text")
  expect_true(any(grepl("bootstrap samples could not be computed", footnotes)))
})

# ==== Raters in rows ====
test_that("Kendall's W with raters in rows matches raters-in-columns reference", {
  df <- data.frame( # 3 raters (rows) assessing 4 subjects (columns)
    s1 = c(1, 2, 1.5),
    s2 = c(2, 3, 2.5),
    s3 = c(3, 1, 3.5),
    s4 = c(4, 4, 1.0)
  )
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("s1", "s2", "s3", "s4")
  options$dataStructure <- "ratersInRows"
  options$kendallW      <- TRUE
  options$ci            <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  # reference: irr::kendall(t(df), correct = TRUE)
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(0.5, 0.2, 1.8, 3, 2.33333333333333, 4.66666666666667,
         0.614934935782537, 0.660869392065187))
})

# ==== Coefficient-specific variable type validation ====
test_that("Invalid variable types show table errors", {
  options <- analysisOptions("raterAgreement")
  options$variables       <- c("contNormal", "contGamma")
  options$variables.types <- c("scale", "scale")
  options$dataStructure   <- "ratersInColumns"
  options$cohensKappa     <- TRUE
  options$fleissKappa     <- TRUE
  results <- runAnalysis("raterAgreement", "debug.csv", options)
  expect_match(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
               "Cohen's kappa requires nominal or ordinal variables")
  expect_match(results[["results"]][["fleissKappa"]][["error"]][["errorMessage"]],
               "Fleiss' kappa requires nominal or ordinal variables")

  options <- analysisOptions("raterAgreement")
  options$variables       <- c("facGender", "facExperim")
  options$variables.types <- c("nominal", "nominal")
  options$dataStructure   <- "ratersInColumns"
  options$cohensKappa     <- TRUE
  options$cohensKappaType <- "weighted"
  options$kendallW        <- TRUE
  results <- runAnalysis("raterAgreement", "debug.csv", options)
  expect_match(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
               "Weighted Cohen's kappa requires ordinal variables")
  expect_match(results[["results"]][["kendallW"]][["error"]][["errorMessage"]],
               "Kendall's W requires ordinal or scale variables")
})

# ==== Fleiss' kappa: category present only in a listwise-deleted row ====
test_that("Fleiss' kappa drops categories that only occur in incomplete rows", {
  df <- data.frame(
    r1 = c("A", "B", "A", "C"),
    r2 = c("A", "B", "B", NA),
    r3 = c("A", "B", "A", "C")
  )
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("r1", "r2", "r3")
  options$dataStructure <- "ratersInColumns"
  options$fleissKappa   <- TRUE
  options$ci            <- TRUE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  jaspTools::expect_equal_tables(results[["results"]][["fleissKappa"]][["data"]],
    list(-0.103321328180018, 1.20332132818002, 0.333333333333333, 0.55, "Overall",
         -0.103321328180018, 1.20332132818002, 0.333333333333333, 0.55, "A",
         -0.103321328180018, 1.20332132818002, 0.333333333333333, 0.55, "B"))
})

# ==== Fleiss' kappa: zero kappa must yield the analytic SE, not NaN ====
test_that("Fleiss' kappa SEs are exact for zero kappa", {
  df <- data.frame(
    r1 = c("A", "A", "B", "B"),
    r2 = c("A", "B", "A", "B")
  )
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("r1", "r2")
  options$dataStructure <- "ratersInColumns"
  options$fleissKappa   <- TRUE
  options$ci            <- TRUE
  results <- runAnalysis("raterAgreement", df, options)
  # analytic SE = sqrt(2 / (ns * nr * (nr - 1))) = sqrt(2/8) = 0.5
  jaspTools::expect_equal_tables(results[["results"]][["fleissKappa"]][["data"]],
    list(-0.979981992270027, 0.979981992270027, 0.5, 0, "Overall",
         -0.979981992270027, 0.979981992270027, 0.5, 0, "A",
         -0.979981992270027, 0.979981992270027, 0.5, 0, "B"))
})

# ==== Weighted Cohen's kappa: declared ordinal order, not alphabetical label order ====
test_that("Weighted Cohen's kappa respects ordered factor levels", {
  lv <- c("low", "medium", "high")
  df <- data.frame(
    r1 = factor(c("low","low","medium","medium","high","high","low","medium","high","high"), levels = lv, ordered = TRUE),
    r2 = factor(c("low","medium","medium","high","high","high","low","low","medium","high"), levels = lv, ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables       <- c("r1", "r2")
  options$variables.types <- c("ordinal", "ordinal")
  options$dataStructure   <- "ratersInColumns"
  options$cohensKappa     <- TRUE
  options$cohensKappaType <- "weighted"
  options$ci              <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # reference: psych::cohen.kappa on the level codes, w.exp = 2 (alphabetical labels give 0.275)
  jaspTools::expect_equal_tables(results[["results"]][["cohensKappa"]][["data"]],
    list(0.710144927536232, "Average kappa", 0.710144927536232, "r1 - r2"))
})

# ==== Raters in rows with labelled ordered factors: level order must survive the transpose ====
test_that("Row-mode ordered factors keep their declared level order", {
  lv <- c("low", "medium", "high")
  # 3 raters (rows) x 6 subjects (columns); same ratings as the column-mode ordered-factor test
  rowMode <- data.frame(
    s1 = factor(c("low", "low", "medium"),    levels = lv, ordered = TRUE),
    s2 = factor(c("medium", "high", "medium"), levels = lv, ordered = TRUE),
    s3 = factor(c("high", "medium", "high"),  levels = lv, ordered = TRUE),
    s4 = factor(c("low", "low", "low"),       levels = lv, ordered = TRUE),
    s5 = factor(c("medium", "medium", "low"), levels = lv, ordered = TRUE),
    s6 = factor(c("high", "high", "high"),    levels = lv, ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables                <- names(rowMode)
  options$variables.types          <- rep("ordinal", ncol(rowMode))
  options$dataStructure            <- "ratersInRows"
  options$kendallW                 <- TRUE
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "ordinal"
  options$ci                       <- FALSE
  results <- runAnalysis("raterAgreement", rowMode, options)
  expect_identical(results[["status"]], "complete")
  # references identical to the column-mode ordered-factor test above
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(7, 0.777777777777778, 11.6666666666667, 5, 4.33333333333333, 8.66666666666667,
         0.0396519759960316, 0.00776230190368293))
  jaspTools::expect_equal_tables(results[["results"]][["krippendorffsAlpha"]][["data"]],
    list(0.675925925925926, "Ordinal"))
})

# ==== Degenerate inputs show table errors instead of crashing ====
test_that("Cohen's kappa errors cleanly when rater pairs have no overlap", {
  df <- data.frame(r1 = c("a", "b", "a", NA, NA, NA), r2 = c(NA, NA, NA, "a", "b", "a"))
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("r1", "r2")
  options$dataStructure <- "ratersInColumns"
  options$cohensKappa   <- TRUE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  expect_match(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
               "jointly rated")
})

test_that("Krippendorff's alpha errors cleanly on all-missing data", {
  df <- data.frame(r1 = c(NA, NA, NA), r2 = c(NA, NA, NA))
  options <- analysisOptions("raterAgreement")
  options$variables          <- c("r1", "r2")
  options$dataStructure      <- "ratersInColumns"
  options$krippendorffsAlpha <- TRUE
  options$ci                 <- TRUE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  expect_match(results[["results"]][["krippendorffsAlpha"]][["error"]][["errorMessage"]],
               "no pairable observations")
})

test_that("One rater in row mode errors cleanly for all coefficients", {
  df <- data.frame(s1 = 1, s2 = 2, s3 = 3, s4 = 4)
  options <- analysisOptions("raterAgreement")
  options$variables          <- c("s1", "s2", "s3", "s4")
  options$dataStructure      <- "ratersInRows"
  options$cohensKappa        <- TRUE
  options$krippendorffsAlpha <- TRUE
  options$kendallW           <- TRUE
  options$ci                 <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  for (tbl in c("cohensKappa", "krippendorffsAlpha", "kendallW"))
    expect_match(results[["results"]][[tbl]][["error"]][["errorMessage"]],
                 "at least 2 raters")
})

test_that("Kendall's W suppresses the F test outside its domain (n = 2, m = 2)", {
  df <- data.frame(r1 = c(1, 2), r2 = c(2, 1))
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("r1", "r2")
  options$dataStructure <- "ratersInColumns"
  options$kendallW      <- TRUE
  options$ci            <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list("", 0, 0, 1, "", "", 1, ""))
  footnotes <- sapply(results[["results"]][["kendallW"]][["footnotes"]], `[[`, "text")
  expect_true(any(grepl("F test is not available", footnotes)))
})

test_that("Fleiss' kappa errors cleanly when all ratings are one category", {
  df <- data.frame(r1 = c("A", "A", "A"), r2 = c("A", "A", "A"), r3 = c("A", "A", "A"))
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("r1", "r2", "r3")
  options$dataStructure <- "ratersInColumns"
  options$fleissKappa   <- TRUE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  expect_match(results[["results"]][["fleissKappa"]][["error"]][["errorMessage"]],
               "single category")
})

test_that("Kendall's W errors cleanly when rankings do not vary", {
  df <- data.frame(r1 = c(2, 2, 2, 2), r2 = c(3, 3, 3, 3), r3 = c(1, 1, 1, 1))
  options <- analysisOptions("raterAgreement")
  options$variables      <- c("r1", "r2", "r3")
  options$dataStructure  <- "ratersInColumns"
  options$kendallW       <- TRUE
  options$correctForTies <- TRUE
  options$ci             <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  expect_match(results[["results"]][["kendallW"]][["error"]][["errorMessage"]],
               "do not vary")
})

# ==== Bootstrap CIs only exist for the tie-corrected coefficient ====
test_that("Uncorrected Kendall's W with CI shows no CI and explains why", {
  df <- data.frame(r1 = 1:5, r2 = 1:5, r3 = 1:5) # identical tie-free rankings, W = 1
  options <- analysisOptions("raterAgreement")
  options$variables        <- c("r1", "r2", "r3")
  options$dataStructure    <- "ratersInColumns"
  options$kendallW         <- TRUE
  options$correctForTies   <- FALSE
  options$ci               <- TRUE
  options$bootstrapSamples <- 1000
  options$setSeed          <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreement", df, options)
  # no CI columns (resampling introduces ties, incompatible with the uncorrected W);
  # F test undefined at perfect concordance (W = 1)
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list("", 1, 12, 4, "", "", 0.0173512652366645, ""))
  footnotes <- sapply(results[["results"]][["kendallW"]][["footnotes"]], `[[`, "text")
  expect_true(any(grepl("Bootstrap CIs are only available with the tie correction", footnotes)))
  expect_true(any(grepl("undefined for perfect concordance", footnotes)))
})

test_that("Corrected Kendall's W bootstrap CI is not distorted by resampling ties", {
  df <- data.frame(r1 = 1:5, r2 = 1:5, r3 = 1:5) # identical tie-free rankings, W = 1
  options <- analysisOptions("raterAgreement")
  options$variables        <- c("r1", "r2", "r3")
  options$dataStructure    <- "ratersInColumns"
  options$kendallW         <- TRUE
  options$correctForTies   <- TRUE
  options$ci               <- TRUE
  options$bootstrapSamples <- 1000
  options$setSeed          <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreement", df, options)
  # every corrected replicate keeps W = 1: SE = 0, CI = [1, 1]
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(1, 1, "", 0, 1, 12, 4, "", "", 0.0173512652366645, ""))
})

# ==== Krippendorff interval method uses actual numeric label values, not level codes ====
test_that("Krippendorff's alpha interval method parses numeric labels as scores", {
  sc <- c("0", "10", "100")
  df <- data.frame(
    r1 = factor(c("0", "10", "100", "10", "0", "100"),   levels = sc),
    r2 = factor(c("0", "100", "100", "10", "10", "100"), levels = sc),
    r3 = factor(c("10", "10", "100", "0", "0", "100"),   levels = sc)
  )
  options <- analysisOptions("raterAgreement")
  options$variables                <- c("r1", "r2", "r3")
  options$dataStructure            <- "ratersInColumns"
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "interval"
  options$ci                       <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # reference: irr::kripp.alpha on the numeric scores (level codes would give 0.6698)
  jaspTools::expect_equal_tables(results[["results"]][["krippendorffsAlpha"]][["data"]],
    list(0.787939988459319, "Interval"))
})

test_that("Krippendorff's alpha interval method rejects non-numeric labels", {
  df <- data.frame(r1 = c("a", "b", "c"), r2 = c("a", "c", "b"))
  options <- analysisOptions("raterAgreement")
  options$variables                <- c("r1", "r2")
  options$dataStructure            <- "ratersInColumns"
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "interval"
  results <- runAnalysis("raterAgreement", df, options)
  expect_match(results[["results"]][["krippendorffsAlpha"]][["error"]][["errorMessage"]],
               "requires numeric ratings")
})

# ==== Raters with mismatched level sets: union codes must respect numeric label order ====
test_that("Union level codes are ordered by value when a rater misses a category", {
  df <- data.frame(
    r1 = factor(c("1", "3", "1", "3", "3", "1"), levels = c("1", "3")), # never used "2"
    r2 = factor(c("1", "2", "3", "2", "3", "1"), levels = c("1", "2", "3")),
    r3 = factor(c("1", "3", "2", "2", "3", "1"), levels = c("1", "2", "3"))
  )
  options <- analysisOptions("raterAgreement")
  options$variables                <- c("r1", "r2", "r3")
  options$variables.types          <- c("ordinal", "ordinal", "ordinal")
  options$dataStructure            <- "ratersInColumns"
  options$kendallW                 <- TRUE
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "ordinal"
  options$ci                       <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # references: irr on the numeric label values (appearance-order union "1","3","2"
  # would scramble the codes and give different results)
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(6.80645161290323, 0.772893772893773, 11.5934065934066, 5, 4.33333333333333,
         8.66666666666667, 0.0408043931372446, 0.00848102504588204))
  jaspTools::expect_equal_tables(results[["results"]][["krippendorffsAlpha"]][["data"]],
    list(0.652777777777778, "Ordinal"))
})

# ==== Fleiss' kappa with mixed factor and numeric columns ====
test_that("Fleiss' kappa handles mixed factor/numeric columns and multi-width labels", {
  df <- data.frame(
    r1 = factor(c("7", "10", "7", "10")),
    r2 = c(7, 10, 7, 10), # numeric column: unlist() would leak factor codes, as.matrix() would pad "7" to " 7"
    r3 = factor(c("7", "10", "10", "7"))
  )
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("r1", "r2", "r3")
  options$dataStructure <- "ratersInColumns"
  options$fleissKappa   <- TRUE
  options$ci            <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  jaspTools::expect_equal_tables(results[["results"]][["fleissKappa"]][["data"]],
    list(0.333333333333333, "Overall", 0.333333333333333, "7", 0.333333333333333, "10"))
})

# ==== Raters in rows with mixed column types is refused ====
test_that("Row mode with mixed categorical/continuous columns gives a validation error", {
  df <- data.frame(s1 = c(1.5, 2.5), s2 = factor(c("a", "b")), s3 = c(7.1, 3.2))
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("s1", "s2", "s3")
  options$dataStructure <- "ratersInRows"
  options$kendallW      <- TRUE
  options$ci            <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "validationError")
})

# ==== Constant ratings error regardless of the tie-correction setting ====
test_that("Kendall's W errors on constant rankings also without tie correction", {
  df <- data.frame(r1 = c(2, 2, 2, 2), r2 = c(3, 3, 3, 3), r3 = c(1, 1, 1, 1))
  options <- analysisOptions("raterAgreement")
  options$variables      <- c("r1", "r2", "r3")
  options$dataStructure  <- "ratersInColumns"
  options$kendallW       <- TRUE
  options$correctForTies <- FALSE # uncorrected W would be a meaningless finite 0
  options$ci             <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_match(results[["results"]][["kendallW"]][["error"]][["errorMessage"]],
               "do not vary")
})

# ==== Round 3: partial textual ordinal level sets ====
test_that("Partial textual level sets merge to one order, invariant to rater order", {
  # no rater declares the full scale: low<high, medium<high, low<medium
  mk <- function(cols) {
    df <- data.frame(
      r1 = factor(c("low", "high", "low", "high", "low", "high"),          levels = c("low", "high"),    ordered = TRUE),
      r2 = factor(c("medium", "high", "medium", "high", "high", "medium"), levels = c("medium", "high"), ordered = TRUE),
      r3 = factor(c("low", "medium", "low", "medium", "low", "high"),      levels = c("low", "medium"),  ordered = TRUE)
    )
    df[, cols, drop = FALSE]
  }
  runFor <- function(cols) {
    options <- analysisOptions("raterAgreement")
    options$variables                <- cols
    options$variables.types          <- rep("ordinal", 3)
    options$dataStructure            <- "ratersInColumns"
    options$cohensKappa              <- TRUE
    options$cohensKappaType          <- "weighted"
    options$krippendorffsAlpha       <- TRUE
    options$krippendorffsAlphaMethod <- "ordinal"
    options$ci                       <- FALSE
    runAnalysis("raterAgreement", mk(cols), options)
  }
  a <- runFor(c("r1", "r2", "r3"))
  b <- runFor(c("r2", "r3", "r1"))
  c <- runFor(c("r3", "r1", "r2"))

  # the constraints imply low < medium < high, whatever order the raters appear in
  avgKappa <- function(r) r[["results"]][["cohensKappa"]][["data"]][[1]][["cKappa"]]
  alpha    <- function(r) r[["results"]][["krippendorffsAlpha"]][["data"]][[1]][["kAlpha"]]
  expect_equal(avgKappa(a), 0.364923747276688)
  expect_equal(avgKappa(b), avgKappa(a))
  expect_equal(avgKappa(c), avgKappa(a))
  expect_equal(alpha(a), 0.294117647058823)
  expect_equal(alpha(b), alpha(a))
  expect_equal(alpha(c), alpha(a))
})

test_that("Ambiguous ordinal schemas are refused by order-sensitive coefficients", {
  # low<high and medium<high are declared, but low vs medium is undetermined
  df <- data.frame(
    a = factor(c("low", "high", "low", "high"),       levels = c("low", "high"),    ordered = TRUE),
    b = factor(c("medium", "high", "medium", "high"), levels = c("medium", "high"), ordered = TRUE),
    c = factor(c("low", "high", "medium", "high"),    levels = c("low", "high"),    ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables                <- c("a", "b", "c")
  options$variables.types          <- rep("ordinal", 3)
  options$dataStructure            <- "ratersInColumns"
  options$cohensKappa              <- TRUE
  options$cohensKappaType          <- "weighted"
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "ordinal"
  options$fleissKappa              <- TRUE
  options$ci                       <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_match(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
  expect_match(results[["results"]][["krippendorffsAlpha"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
  # Fleiss' kappa is order-free, so it still computes
  expect_gt(length(results[["results"]][["fleissKappa"]][["data"]]), 0)
})

# ==== Round 3: Kendall ranks within raters, so columns keep their own scale ====
test_that("Kendall's W with mixed ordinal and scale raters is subject-order invariant", {
  df <- data.frame(
    ord = factor(c("low", "medium", "high", "medium", "low", "high", "high", "low", "medium"),
                 levels = c("low", "medium", "high"), ordered = TRUE),
    sc1 = c(2.5, 7.1, 9.9, 4.0, 1.2, 8.8, 9.1, 0.5, 5.5),
    sc2 = c(1.0, 6.0, 9.0, 5.0, 2.0, 8.0, 9.5, 0.1, 4.4)
  )
  options <- analysisOptions("raterAgreement")
  options$variables       <- c("ord", "sc1", "sc2")
  options$variables.types <- c("ordinal", "scale", "scale")
  options$dataStructure   <- "ratersInColumns"
  options$kendallW        <- TRUE
  options$ci              <- FALSE
  a <- runAnalysis("raterAgreement", df, options)
  b <- runAnalysis("raterAgreement", df[c(5, 2, 9, 1, 7, 3, 8, 6, 4), ], options) # permute subjects

  # reference: irr::kendall on the per-column codes (numeric columns keep their values)
  expect_equal(a[["results"]][["kendallW"]][["data"]][[1]][["W"]], 0.96551724137931)
  expect_equal(b[["results"]][["kendallW"]][["data"]][[1]][["W"]],
               a[["results"]][["kendallW"]][["data"]][[1]][["W"]])
})

# ==== Round 3: Krippendorff validates the effective coincidence data ====
test_that("Krippendorff's alpha rejects data whose pairable ratings do not vary", {
  df <- data.frame(               # one pairable A/A item, the rest are singleton B ratings
    r1 = c("A", "B",  NA,  NA),
    r2 = c("A",  NA, "B",  NA),
    r3 = c(NA,   NA,  NA, "B")
  )
  options <- analysisOptions("raterAgreement")
  options$variables          <- c("r1", "r2", "r3")
  options$dataStructure      <- "ratersInColumns"
  options$krippendorffsAlpha <- TRUE
  options$ci                 <- TRUE
  options$bootstrapSamples   <- 200
  options$setSeed            <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreement", df, options)
  # previously: NaN estimate with a bootstrap CI of [1, 1]
  expect_match(results[["results"]][["krippendorffsAlpha"]][["error"]][["errorMessage"]],
               "pairable ratings do not vary")
  expect_length(results[["results"]][["krippendorffsAlpha"]][["data"]], 0)
})

# ==== Round 3: Cohen validates each analyzed rater pair ====
test_that("Cohen's kappa validates each rater pair on its own complete cases", {
  runFor <- function(df, vars = c("r1", "r2")) {
    options <- analysisOptions("raterAgreement")
    options$variables     <- vars
    options$dataStructure <- "ratersInColumns"
    options$cohensKappa   <- TRUE
    options$ci            <- FALSE
    runAnalysis("raterAgreement", df, options)
  }
  # both raters constant -> not estimable (was: blank kappa, no error)
  r <- runFor(data.frame(r1 = rep("a", 5), r2 = rep("a", 5)))
  expect_match(r[["results"]][["cohensKappa"]][["error"]][["errorMessage"]], "do not vary")

  # one jointly rated subject -> not enough (was: NaN)
  r <- runFor(data.frame(r1 = c("a", "b", NA, NA, NA), r2 = c("a", NA, "b", NA, NA)))
  expect_match(r[["results"]][["cohensKappa"]][["error"]][["errorMessage"]], "fewer than 3")

  # two jointly rated subjects -> still below the stated minimum (was: kappa = 1)
  r <- runFor(data.frame(r1 = c("a", "b", NA, NA, NA), r2 = c("a", "b", NA, NA, NA)))
  expect_match(r[["results"]][["cohensKappa"]][["error"]][["errorMessage"]], "fewer than 3")

  # a single unusable pair does not blank the usable ones
  df <- data.frame(
    r1 = c("a", "b", "a", "b", "a", "b"),
    r2 = c("a", "b", "b", "b", "a", "a"),
    r3 = c("a",  NA,  NA,  NA,  NA,  NA)
  )
  r <- runFor(df, vars = c("r1", "r2", "r3"))
  expect_identical(r[["status"]], "complete")
  rows <- r[["results"]][["cohensKappa"]][["data"]]
  expect_equal(rows[[1]][["cKappa"]], 0.333333333333333) # average over the valid pair only
  expect_equal(rows[[2]][["cKappa"]], 0.333333333333333) # r1 - r2
  expect_identical(rows[[3]][["cKappa"]], "")            # r1 - r3, not computable
  footnotes <- sapply(r[["results"]][["cohensKappa"]][["footnotes"]], `[[`, "text")
  expect_true(any(grepl("Some rater pairs could not be computed", footnotes)))
})

test_that("Cohen's kappa pair labels survive variable names containing spaces", {
  df <- data.frame(`Rater one` = c("a", "b", "a", "b"), `Rater two` = c("a", "b", "b", "a"),
                   check.names = FALSE)
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("Rater one", "Rater two")
  options$dataStructure <- "ratersInColumns"
  options$cohensKappa   <- TRUE
  options$ci            <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # was "Rater - one Rater two" (psych's pair names, first space replaced)
  expect_identical(results[["results"]][["cohensKappa"]][["data"]][[2]][["ratings"]],
                   "Rater one - Rater two")
})

# ==== Round 3: weighted Cohen with CI on a valid shared ordinal scale ====
test_that("Weighted Cohen's kappa with CI matches psych on a shared ordinal scale", {
  lv <- c("low", "medium", "high")
  df <- data.frame(
    r1 = factor(c("low","low","medium","medium","high","high","low","medium","high","high"), levels = lv, ordered = TRUE),
    r2 = factor(c("low","medium","medium","high","high","high","low","low","medium","high"), levels = lv, ordered = TRUE),
    r3 = factor(c("high","low","low","medium","medium","high","medium","low","high","low"),  levels = lv, ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables       <- c("r1", "r2", "r3")
  options$variables.types <- rep("ordinal", 3)
  options$dataStructure   <- "ratersInColumns"
  options$cohensKappa     <- TRUE
  options$cohensKappaType <- "weighted"
  options$ci              <- TRUE
  options$ciLevel         <- 0.95
  results <- runAnalysis("raterAgreement", df, options)
  # reference: psych::cohen.kappa(w.exp = 2) on the level codes, per pair
  jaspTools::expect_equal_tables(results[["results"]][["cohensKappa"]][["data"]],
    list("", "", "", 0.293053004014425, "Average kappa",
         0.43904800204193, 0.981241853030534, 0.138317299518093, 0.710144927536232, "r1 - r2",
         -0.455696679422334, 0.765555834351911, 0.311549733415341, 0.154929577464788, "r1 - r3",
         -0.547150588633365, 0.575319602717872, 0.286349698312096, 0.0140845070422532, "r2 - r3"))
})

# ==== Round 4: factor levels stay authoritative even when labels look numeric ====
test_that("Weighted Cohen's kappa uses declared factor order, not label value, for numeric-looking levels", {
  lv      <- c("low", "medium", "high")
  relabel <- c(low = "1", medium = "3", high = "2") # non-monotonic labels, same declared order
  r1chr <- c("low","low","medium","medium","high","high","low","medium","high","high")
  r2chr <- c("low","medium","medium","high","high","high","low","low","medium","high")
  df <- data.frame(
    r1 = factor(relabel[r1chr], levels = relabel[lv], ordered = TRUE),
    r2 = factor(relabel[r2chr], levels = relabel[lv], ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables       <- c("r1", "r2")
  options$variables.types <- c("ordinal", "ordinal")
  options$dataStructure   <- "ratersInColumns"
  options$cohensKappa     <- TRUE
  options$cohensKappaType <- "weighted"
  options$ci              <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # identical declared order to the low/medium/high test above -> identical average kappa
  # (0.710144927536232), even though the labels ("1","3","2") look numeric. Sorting by
  # label value instead would give 0.166666666666667.
  jaspTools::expect_equal_tables(results[["results"]][["cohensKappa"]][["data"]],
    list(0.710144927536232, "Average kappa", 0.710144927536232, "r1 - r2"))
})

test_that("Krippendorff's alpha and Kendall's W use declared factor order for numeric-looking levels", {
  lv      <- c("low", "medium", "high")
  relabel <- c(low = "1", medium = "3", high = "2")
  mkCol <- function(x) factor(relabel[x], levels = relabel[lv], ordered = TRUE)
  df <- data.frame(
    r1 = mkCol(c("low", "medium", "high", "low",    "medium", "high")),
    r2 = mkCol(c("low", "high",   "medium", "low",  "medium", "high")),
    r3 = mkCol(c("medium", "medium", "high", "low", "low",    "high"))
  )
  options <- analysisOptions("raterAgreement")
  options$variables                <- c("r1", "r2", "r3")
  options$variables.types          <- c("ordinal", "ordinal", "ordinal")
  options$dataStructure            <- "ratersInColumns"
  options$kendallW                 <- TRUE
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "ordinal"
  options$ci                       <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # same references as "Kendall's W and Krippendorff's alpha respect ordered factor levels"
  # above: Kendall ranks within raters (label-value-independent), Krippendorff uses the
  # declared factor order regardless of the numeric-looking labels
  jaspTools::expect_equal_tables(results[["results"]][["kendallW"]][["data"]],
    list(7, 0.777777777777778, 11.6666666666667, 5, 4.33333333333333, 8.66666666666667,
         0.0396519759960316, 0.00776230190368293))
  jaspTools::expect_equal_tables(results[["results"]][["krippendorffsAlpha"]][["data"]],
    list(0.675925925925926, "Ordinal"))
})

test_that("Contradictory numeric-looking factor orders are refused, not silently sorted", {
  # r1 declares 1 < 2 < 3; r2 declares 1 < 3 < 2 -- contradictory on 2 vs 3
  df <- data.frame(
    r1 = factor(c("1","2","3","1","2","3"), levels = c("1","2","3"), ordered = TRUE),
    r2 = factor(c("1","3","2","1","3","2"), levels = c("1","3","2"), ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables                <- c("r1", "r2")
  options$variables.types          <- c("ordinal", "ordinal")
  options$dataStructure             <- "ratersInColumns"
  options$cohensKappa              <- TRUE
  options$cohensKappaType          <- "weighted"
  options$krippendorffsAlpha       <- TRUE
  options$krippendorffsAlphaMethod <- "ordinal"
  options$ci                       <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  expect_match(results[["results"]][["cohensKappa"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
  expect_match(results[["results"]][["krippendorffsAlpha"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
})

# ==== Round 4: pairwise weighted Cohen must use the full common scale, not just the ====
# ==== categories a given pair happens to observe ====
test_that("Weighted Cohen's kappa keeps an unused interior category in the scale for each pair", {
  lv <- c("A", "B", "C", "D")
  # r1, r2 never rate "B", but the declared common scale still has 4 categories
  r1 <- c("A","A","C","C","D","D","A","C","D","A")
  r2 <- c("A","C","D","A","D","C","A","D","D","C")
  df <- data.frame(
    r1 = factor(r1, levels = lv, ordered = TRUE),
    r2 = factor(r2, levels = lv, ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables       <- c("r1", "r2")
  options$variables.types <- c("ordinal", "ordinal")
  options$dataStructure   <- "ratersInColumns"
  options$cohensKappa     <- TRUE
  options$cohensKappaType <- "weighted"
  options$ci              <- FALSE
  results <- runAnalysis("raterAgreement", df, options)
  # reference: psych::cohen.kappa(codes, w.exp = 2, levels = 1:4); without the levels
  # argument, psych collapses the confusion matrix to only {A,C,D} and gives 0.577464788732394
  jaspTools::expect_equal_tables(results[["results"]][["cohensKappa"]][["data"]],
    list(0.545454545454546, "Average kappa", 0.545454545454546, "r1 - r2"))
})

# ==== Round 4: row-mode Kendall must respect an ambiguous merged ordinal scale ====
test_that("Row-mode Kendall's W refuses an ambiguous common ordinal scale, invariant to column order", {
  # subject columns declare partial constraints A<C, B<C, C<D: A vs B is undetermined
  mk <- function(cols) {
    df <- data.frame(
      s1 = factor(c("A","C","A"), levels = c("A","C"), ordered = TRUE),
      s2 = factor(c("C","B","C"), levels = c("B","C"), ordered = TRUE),
      s3 = factor(c("D","C","D"), levels = c("C","D"), ordered = TRUE)
    )
    df[, cols, drop = FALSE]
  }
  runFor <- function(cols) {
    options <- analysisOptions("raterAgreement")
    options$variables       <- cols
    options$variables.types <- rep("ordinal", 3)
    options$dataStructure   <- "ratersInRows"
    options$kendallW        <- TRUE
    options$ci              <- FALSE
    runAnalysis("raterAgreement", mk(cols), options)
  }
  a <- runFor(c("s1", "s2", "s3"))
  b <- runFor(c("s2", "s1", "s3")) # merely reordering the subject columns must not change the outcome
  expect_match(a[["results"]][["kendallW"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
  expect_match(b[["results"]][["kendallW"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
})

test_that("Row-mode Kendall's W bootstrap is skipped when the common ordinal scale is ambiguous", {
  df <- data.frame(
    s1 = factor(c("A","C","A"), levels = c("A","C"), ordered = TRUE),
    s2 = factor(c("C","B","C"), levels = c("B","C"), ordered = TRUE),
    s3 = factor(c("D","C","D"), levels = c("C","D"), ordered = TRUE)
  )
  options <- analysisOptions("raterAgreement")
  options$variables        <- c("s1", "s2", "s3")
  options$variables.types  <- rep("ordinal", 3)
  options$dataStructure    <- "ratersInRows"
  options$kendallW         <- TRUE
  options$ci               <- TRUE
  options$bootstrapSamples <- 100
  options$setSeed          <- TRUE
  set.seed(1)
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  expect_match(results[["results"]][["kendallW"]][["error"]][["errorMessage"]],
               "requires a common ordinal scale")
})

# ==== Round 3: Kendall bootstrap must not run on invalid input ====
test_that("Kendall's W bootstrap is skipped for input the table rejects", {
  df <- data.frame(r1 = c(2, 2, 2, 2), r2 = c(3, 3, 3, 3), r3 = c(1, 1, 1, 1))
  options <- analysisOptions("raterAgreement")
  options$variables        <- c("r1", "r2", "r3")
  options$dataStructure    <- "ratersInColumns"
  options$kendallW         <- TRUE
  options$correctForTies   <- TRUE
  options$ci               <- TRUE
  options$bootstrapSamples <- 100000 # would take a long time if the guard were missing
  options$setSeed          <- TRUE
  elapsed <- system.time(results <- runAnalysis("raterAgreement", df, options))[["elapsed"]]
  expect_match(results[["results"]][["kendallW"]][["error"]][["errorMessage"]], "do not vary")
  expect_lt(elapsed, 20)
})
