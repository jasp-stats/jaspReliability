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
  table <- results[["results"]][["cohensKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("", "", "", 0.114984433765595, "Average kappa", 0.0924679737811612,
                                      0.357606541867125, 0.0676386326935962, 0.225037257824143, "V1 - V2",
                                      -0.00101486089827027, 0.203411665159256, 0.052150582273454,
                                      0.101198402130493, "V1 - V3", -0.0200019787491091, 0.151023173354119,
                                      0.0436296670378263, 0.0655105973025049, "V2 - V3", 0.0576342926736923,
                                      0.309887992795108, 0.0643516161804914, 0.1837611427344, "V1 - V4",
                                      -0.0882775571852617, 0.1399093155876, 0.0582120065911339, 0.025815879201169,
                                      "V2 - V4", -0.032509269670062, 0.257559504097391, 0.0739984959048937,
                                      0.112525117213664, "V3 - V4", -0.0240223485028135, 0.221356704625165,
                                      0.0625978474766621, 0.0986671780611757, "V1 - V5", 0.0762065929175585,
                                      0.337621571800286, 0.0666887200338209, 0.206914082358922, "V2 - V5",
                                      -0.0317077880237179, 0.0969460355453044, 0.0328204560348627,
                                      0.0326191237607932, "V3 - V5", 0.00907597751537415, 0.186515136621999,
                                      0.0452659233807975, 0.0977955570686868, "V4 - V5"))
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


test_that("Cohen's Weighted kappa table results match", {
  table <- results[["results"]][["cohensKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list("", "", "", -0.00184386638316336, "", "Average kappa", -0.0105967225936691,
                                      0.0461522781492244, 0.0110156757407617, 0.0177777777777777,
                                      100, "facGender - facExperim", -0.0393585273220352, 0.0269993851919001,
                                      0.0128808831436364, -0.00617957106506752, 80, "facGender - debBinMiss20",
                                      -0.0516244390801062, 0.0173648273557058, 0.013391661151762,
                                      -0.0171298058622002, 80, "facExperim - debBinMiss20"))
})

test_that("Fleiss' kappa table results match", {
  table <- results[["results"]][["fleissKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(-0.276327852798572, -0.127177258822542, 0.0289519561273983, -0.201752555810557,
                                      "Overall", -0.337000773905705, -0.00446264072844088, 0.0645497224367903,
                                      -0.170731707317073, 0, -0.397038297357863, -0.0645001641805985,
                                      0.0645497224367903, -0.230769230769231, 1, -0.384543178263759,
                                      -0.0520050450864946, 0.0645497224367903, -0.218274111675127, "control",
                                      -0.348535076440849, -0.0159969432635844, 0.0645497224367903,
                                      -0.182266009852217, "experimental", -0.360298917334901, -0.0277607841576364,
                                      0.0645497224367903, -0.194029850746269, "f", -0.372299217342401,
                                      -0.0397610841651366, 0.0645497224367903, -0.206030150753769, "m"))
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
