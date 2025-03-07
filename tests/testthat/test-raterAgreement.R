# ==== Ensure results are unchanged on JASP debug data ====

####Cohen's unweighted kappa and Fleiss' kappa and Krippendorff's alpha####

# Set options
options <- analysisOptions("raterAgreement")
options$variables <- paste0("V", 1:5)
options$dataStructure <- "ratersInColumns"
options$setSeed <- TRUE
options$krippendorffsAlphaBootstrapSamplesForCI <- 200
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
                                 list(0.0321990251552643, 0.119851279212922, 0.0223606797749979, 0.0760251521840934,
                                      "Overall", 0.0321890403455753, 0.119810959654425, 0.0223529411764706,
                                      0.076, 0, 0.0321890403455753, 0.119810959654425, 0.0223529411764706,
                                      0.076, 1))
})

test_that("Krippendorff's alpha table results match", {
  table <- results[["results"]][["krippendorffsAlpha"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.0271069334114449, 0.124379082494376, 0.0257751838160253, 0.0762561458960466,
                                      "Nominal"))
})


####Cohen's weighted kappa and Fleiss' kappa and Krippendorf's alpha with different CI range(99%)####

# Set options
options <- analysisOptions("raterAgreement")
options$variables <- c("facGender", "facExperim", "debBinMiss20")
options$ciLevel <- 0.99
options$cohensKappaType <- "weighted"
options$krippendorffsAlphaBootstrapSamplesForCI <- 200
options$dataStructure <- "ratersInColumns"
options$setSeed <- TRUE
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
                                 list(-0.276327852798572, -0.127177258822542, 0.0289519561273983, -0.201752555810557,
                                      "Overall", -0.360528094860817, -0.0274719051391827, 0.0646502835538752,
                                      -0.194, "f", -0.37243820115239, -0.0395617988476095, 0.0646153846153846,
                                      -0.206, "m", -0.384084231935422, -0.0519157680645785, 0.0644779650990831,
                                      -0.218, "control", -0.348005996191891, -0.0159940038081091,
                                      0.0644475920679887, -0.182, "experimental", -0.33723781932418,
                                      -0.00476218067581949, 0.0645375914836993, -0.171, 0, -0.397234597910737,
                                      -0.0647654020892627, 0.0645363408521303, -0.231, 1))
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
  options$cohensKappa <- FALSE
  options$krippendorffsAlpha <- FALSE
  options$ci <- FALSE
  set.seed(1)
  results <- runAnalysis("raterAgreement", testthat::test_path("Fleiss1971.csv"), options)
  table <- results[["results"]][["fleissKappa"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(0.430244520060141, "Overall", 0.245, 1, 0.245, 2, 0.52, 3, 0.471,
                                      4, 0.566, 5))
})


test_that("Cohen's kappa table results match with linear weighting", {
  options <- analysisOptions("raterAgreement")
  options$variables <- c("V1", "V2")
  options$fleissKappa <- FALSE
  options$krippendorffsAlpha <- FALSE
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

