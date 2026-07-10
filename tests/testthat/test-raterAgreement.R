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
                                 list("", "", "", -0.00184386638316336, "Average kappa", -0.0105967225936691,
                                      0.0461522781492244, 0.0110156757407617, 0.0177777777777777,
                                      "facGender - facExperim", -0.0393585273220366, 0.026999385191902,
                                      0.012880883143637, -0.00617957106506728, "facGender - debBinMiss20",
                                      -0.0516244390801057, 0.0173648273557048, 0.0133916611517617,
                                      -0.0171298058622005, "facExperim - debBinMiss20"))
})

test_that("Fleiss' kappa table results match", {
  table <- results[["results"]][["fleissKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(-0.276327852798572, -0.127177258822542, 0.0289519561273982, -0.201752555810557,
                                      "Overall", -0.360298917334901, -0.0277607841576364, 0.0645497224367903,
                                      -0.194029850746269, "f", -0.372299217342401, -0.0397610841651366, 0.0645497224367903,
                                      -0.206030150753769, "m", -0.384543178263759, -0.0520050450864946, 0.0645497224367903,
                                      -0.218274111675127, "control", -0.348535076440849, -0.0159969432635844,
                                      0.0645497224367903, -0.182266009852217, "experimental", -0.337000773905705,
                                      -0.00446264072844088, 0.0645497224367903, -0.170731707317073, 0, -0.397038297357863,
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
  table <- results[["results"]][["kendallW"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.235537187052039, 0.405722388905557, 0.0428257441483345, 0.316646331299797,
                                      94.0439603960396, 99, 0.621972149059366))
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
    list(0.31477197719772, 124.649702970297, 99, 0.0416297450235598))
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
  jaspTools::expect_equal_tables(table,
    list(0.539656900212835, 30.7604444444444, 19, 0.0428834698290932))
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
    list(0.777777777777778, 11.6666666666667, 5, 0.0396519759960316))
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
    list(0.501921470342523, 28.6095238095238, 19, 0.0723803546937571))
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
    list(0.111111111111111, 1, 0.350241068029438, 0.444444444444444,
         2.66666666666667, 2, 0.263597138115727))
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
    list(0.2, 1.8, 3, 0.614934935782537))
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
    r1 = c("A", "B", "C"),
    r2 = c("A", "B", NA),
    r3 = c("A", "B", "C")
  )
  options <- analysisOptions("raterAgreement")
  options$variables     <- c("r1", "r2", "r3")
  options$dataStructure <- "ratersInColumns"
  options$fleissKappa   <- TRUE
  options$ci            <- TRUE
  results <- runAnalysis("raterAgreement", df, options)
  expect_identical(results[["status"]], "complete")
  jaspTools::expect_equal_tables(results[["results"]][["fleissKappa"]][["data"]],
    list(0.199848053940782, 1.80015194605922, 0.408248290463863, 1, "Overall",
         0.199848053940782, 1.80015194605922, 0.408248290463863, 1, "A",
         0.199848053940782, 1.80015194605922, 0.408248290463863, 1, "B"))
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
