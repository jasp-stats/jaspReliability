

# test the sem functions for ordinal data

options <- list(
  .meta = list(
    variables = list(
      hasTypes = TRUE,
      shouldEncode = TRUE
    )
  ),
  anova = TRUE,
  combinedPointPlot = TRUE,
  feldt = TRUE,
  feldtNumberOfSplits = "5",
  hideTable = FALSE,
  histogramCounts = FALSE,
  irt = TRUE,
  keats = FALSE,
  lord = FALSE,
  lord2 = FALSE,
  lord2NumberOfSplits = "5",
  minimumGroupSize = 10,
  mollenkopfFeldt = TRUE,
  mollenkopfFeldtNumberOfSplits = "5",
  mollenkopfFeldtPolyDegree = 3,
  plotHeight = 320,
  plotWidth = 480,
  pointPlots = TRUE,
  reliabilityValue = 0.5,
  thorndike = TRUE,
  userReliability = FALSE,
  sumScoreCiTable = FALSE,
  sumScoreCiPlots = FALSE,
  variables.types = c("ordinal", "ordinal", "ordinal", "ordinal", "ordinal"),
  variables = c("asrm_1", "asrm_2", "asrm_3", "asrm_4", "asrm_5")
)

results <- runAnalysis("standardErrorOfMeasurement", dataset = testthat::test_path("asrm.csv"), options = options, makeTests = FALSE)

test_that("Standard Error of Measurement table results match", {
  table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(1.29559693908693, 3, 1.33767168953407, 0.229854412468377, 0.810972506211268,
                                  5, 0.925820099772551, 1.29559693908693, 0, 1.33767168953407,
                                  0.982935317522817, 1.07055114141815, 6, 0.925820099772551, 1.29559693908693,
                                  2, 1.33767168953407, 1.266173971192005, 1.23298708781765, 7,
                                  0.925820099772551, 1.29559693908693, 4, 1.33767168953407, 1.383637134534919,
                                  1.34041735804632, 8, 0.925820099772551, 1.29559693908693, 6,
                                  1.33767168953407, 1.406286664699453, 1.41137323284577, 9, 0.925820099772551,
                                  1.15961069747436, 7, 1.13338336703573, 1.428717889549439, 1.45689201232324,
                                  10, 0.866025403784439, 1.15961069747436, 5, 1.13338336703573,
                                  1.497031064514674, 1.48489344540488, 11, 0.866025403784439, 1.70538363230486,
                                  7, 1.70918401391582, 1.545108940646190, 1.50186786343546, 12,
                                  1.64189930669738, 1.70538363230486, 9, 1.70918401391582, 1.539869728173880,
                                  1.5136512119818, 13, 1.64189930669738, 1.62117399689396, 13,
                                  1.60299994073407, 1.518672162069087, 1.52578416955645, 14, 1.53589529557661,
                                  1.43811745632331, 11, 1.41939979368434, 1.524246191291583, 1.54361561010075,
                                  15, 1.57249078621379, 1.65968233543213, 6, 1.62928079660054,
                                  1.523552374802350, 1.57221038715823, 16, 1.84390889145858, 1.65968233543213,
                                  1, 1.62928079660054, 1.488853997223532, 1.61611151602076, 17, 1.84390889145858,
                                  1.65968233543213, 2, 1.62928079660054, 1.449888652774629, 1.67903900932886,
                                  18, 1.84390889145858, 1.65968233543213, 1, 1.62928079660054,
                                  1.462503108363163, 1.76363892004825, 19, 1.84390889145858, 1.65968233543213,
                                  0, 1.62928079660054, 1.472112424372880, 1.87138377130906, 20,
                                  1.84390889145858, 1.65968233543213, 0, 1.62928079660054, 1.400888085158652,
                                  2.00265997698042, 21, 1.84390889145858, 1.65968233543213, 1,
                                  1.62928079660054, 1.228420561649999, 2.1569982806834, 22, 1.84390889145858
                                 ))
})


# binary data
options <- list(.meta = list(variables = list(hasTypes = TRUE,
                                              shouldEncode = TRUE)), anova = FALSE, ciLevelPlots = 0.95,
                ciLevelTable = 0.95, combinedPointPlot = FALSE, feldt = FALSE,
                feldtNumberOfSplits = "2", hideTable = FALSE, histogramCounts = FALSE,
                irt = TRUE, keats = TRUE, lord = TRUE, lord2 = TRUE, lord2NumberOfSplits = "2",
                minimumGroupSize = 10, mollenkopfFeldt = FALSE, mollenkopfFeldtNumberOfSplits = "2",
                mollenkopfFeldtPolyDegree = 3, plotHeight = 320, plotWidth = 480,
                pointPlots = TRUE, reliabilityValue = 0.5, sumScoreCiPlots = FALSE,
                sumScoreCiPlotsCutoff = FALSE, sumScoreCiPlotsCutoffValue = 0,
                sumScoreCiTable = FALSE, thorndike = FALSE, userReliability = FALSE,
                variables = c("V1", "V2", "V3", "V4"), variables.types = c("nominal",
                                                                           "nominal", "nominal", "nominal"))

results <- runAnalysis("standardErrorOfMeasurement", dataset = testthat::test_path("binaryTestDt.csv"), options = options, makeTests = FALSE)

test_that("Standard error of measurement table results match", {
  table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
  # handle slightly different values on macOS compared to Windows and Linux
  options("jaspRoundToPrecision" = function(x) signif(round(x, digits = 3), digits = 3))
  jaspTools::expect_equal_tables(table,
                                 list(
                                     7.0, 0.49433301750863673, 0.0, 0.0, 0.84852813742385702, 0,
                                     18.0, 0.76104647303067507, 0.9151894370278959, 1.0, 0.84852813742385702, 1,
                                     56.0, 0.84058341505800005, 1.0567697356551156, 1.1547005383792515, 1.2955969390869324, 2,
                                     69.0, 0.81285834219580733, 0.9151894370278959, 1.0, 1.0, 3,
                                     50.0, 0.199948272704165, 0.0, 0.0, 0.0, 4))
  options("jaspRoundToPrecision" = NULL)
})


# and some extra plots and stuff

options <- list(.meta = list(variables = list(hasTypes = TRUE,
                                              shouldEncode = TRUE)), anova = TRUE, ciLevelPlots = 0.99,
                ciLevelTable = 0.1, combinedPointPlot = FALSE, feldt = TRUE,
                feldtNumberOfSplits = "5", hideTable = FALSE, histogramCounts = FALSE,
                irt = FALSE, keats = TRUE, lord = FALSE, lord2 = FALSE, lord2NumberOfSplits = "5",
                minimumGroupSize = 10, mollenkopfFeldt = FALSE, mollenkopfFeldtNumberOfSplits = "5",
                mollenkopfFeldtPolyDegree = 3, plotHeight = 320, plotWidth = 480,
                pointPlots = FALSE, reliabilityValue = 0.5, sumScoreCiPlots = TRUE,
                sumScoreCiPlotsCutoff = TRUE, sumScoreCiPlotsCutoffValue = 3,
                sumScoreCiTable = TRUE, thorndike = FALSE, userReliability = FALSE,
                variables = c("V1", "V2", "V3", "V4", "V5"), variables.types = c("nominal",
                                                                                 "nominal", "nominal", "nominal", "nominal"))

results <- runAnalysis("standardErrorOfMeasurement", data = testthat::test_path("binaryTestDt.csv"), options = options, makeTests = F)


test_that("Sum Score CI Table results match", {
  table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(-0.0793452462508615, -0.0854278968697352, "", 0, 0.0793452462508615,
                                  0.0854278968697352, "", 0.920654753749139, 0.914572103130265,
                                  0.888491541244813, 1, 1.07934524625086, 1.08542789686974, 1.11150845875519,
                                  1.87104202115792, 1.86828037011968, 1.8634305870228, 2, 2.12895797884208,
                                  2.13171962988032, 2.1365694129772, 2.87297101344401, 2.87024411426061,
                                  2.8634305870228, 3, 3.12702898655599, 3.12975588573939, 3.1365694129772,
                                  3.88445146704636, 3.88331172201698, 3.88849154124481, 4, 4.11554853295364,
                                  4.11668827798302, 4.11150845875519, "", 4.93833944832921, "",
                                  5, "", 5.06166055167079, ""))
})

test_that("Sum score CI plots matches for ", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciPlotsContainer"]][["collection"]][["semMainContainer_ciPlotsContainer_anovaPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "anova-ci-plot")
})

test_that("Sum score CI plots matches for ", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciPlotsContainer"]][["collection"]][["semMainContainer_ciPlotsContainer_feldtPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "feldt-ci-plot")
})

test_that("Sum score CI plots matches for ", {
  plotName <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_ciPlotsContainer"]][["collection"]][["semMainContainer_ciPlotsContainer_keatsPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "keats-ci-plot")
})



# regression test for the merged-group term of the ANOVA method:
# Emons (2023) notes that the ANOVA method applied to a single, unconditional sample
# must reproduce the SEM that coefficient alpha implies. Merging every score into one
# group is exactly that case, so the ANOVA column has to equal the unconditional SEM.
options <- list(.meta = list(variables = list(hasTypes = TRUE, shouldEncode = TRUE)),
                anova = TRUE, ciLevelPlots = 0.95, ciLevelTable = 0.95,
                combinedPointPlot = FALSE, feldt = FALSE, feldtNumberOfSplits = "2",
                hideTable = FALSE, histogramCounts = FALSE, irt = FALSE, keats = FALSE,
                lord = FALSE, lord2 = FALSE, lord2NumberOfSplits = "2",
                minimumGroupSize = 78, mollenkopfFeldt = FALSE, # 78 == nrow(asrm.csv)
                mollenkopfFeldtNumberOfSplits = "2", mollenkopfFeldtPolyDegree = 3,
                plotHeight = 320, plotWidth = 480, pointPlots = FALSE,
                reliabilityValue = 0.5, sumScoreCiPlots = FALSE,
                sumScoreCiPlotsCutoff = FALSE, sumScoreCiPlotsCutoffValue = 0,
                sumScoreCiTable = FALSE, thorndike = FALSE, userReliability = FALSE,
                variables = paste0("asrm_", 1:5), variables.types = rep("ordinal", 5))

results <- runAnalysis("standardErrorOfMeasurement", dataset = testthat::test_path("asrm.csv"),
                       options = options, makeTests = FALSE)

test_that("ANOVA method reduces to the unconditional SEM when all scores form one group", {
  dat       <- as.matrix(read.csv(testthat::test_path("asrm.csv")))
  nit       <- ncol(dat)
  sumScores <- rowSums(dat)
  alpha     <- nit / (nit - 1) * (1 - sum(apply(dat, 2, var)) / var(sumScores))
  semAlpha  <- sd(sumScores) * sqrt(1 - alpha)

  table    <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
  anovaSem <- unique(vapply(table, function(row) row[["anova"]], numeric(1)))

  expect_equal(length(anovaSem), 1)
  expect_equal(anovaSem, semAlpha, tolerance = 1e-6)
})


# the items must be numerically coded; items with text labels give an error that says how to fix it
asrmNumeric <- read.csv(testthat::test_path("asrm.csv"), fileEncoding = "UTF-8-BOM")
semOptions  <- modifyList(options, list(minimumGroupSize = 10, thorndike = TRUE))

coefficientsOf <- function(results) {
  results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]][["data"]]
}

test_that("Items with text labels give an error that says how to fix it", {
  textLevels <- c("never", "rarely", "sometimes", "often", "always")
  asrmText   <- as.data.frame(lapply(asrmNumeric, function(x) factor(textLevels[x], levels = textLevels)))
  results <- runAnalysis("standardErrorOfMeasurement", dataset = asrmText, options = semOptions, makeTests = FALSE)

  expect_identical(results[["status"]], "validationError")
  expect_match(results[["results"]][["errorMessage"]], "label editor")
  expect_match(results[["results"]][["errorMessage"]], "asrm_1")
})

test_that("Only the items with text labels are named in the error", {
  asrmMixed <- asrmNumeric
  asrmMixed[[2]] <- factor(c("low", "mid", "high", "higher", "highest")[asrmMixed[[2]]])
  results <- runAnalysis("standardErrorOfMeasurement", dataset = asrmMixed, options = semOptions, makeTests = FALSE)

  expect_identical(results[["status"]], "validationError")
  expect_match(results[["results"]][["errorMessage"]], "asrm_2")
  expect_false(grepl("asrm_1", results[["results"]][["errorMessage"]]))
})


# binomial methods: Lord, Keats and Lord generalized treat the sum score as the number of correct
# responses, so they need 0/1 items, and the IRT method puts its SEM on a scale that starts at 0.
# Other codes of dichotomous items (e.g. 1/2) must give the same SEM as 0/1 coding, on the scale
# of the original sum score.
binaryItems     <- read.csv(testthat::test_path("binaryTestDt.csv"), fileEncoding = "UTF-8-BOM")[, 1:4]
binomialOptions <- modifyList(options, list(anova = FALSE, thorndike = FALSE, lord = TRUE, keats = TRUE, lord2 = TRUE, irt = TRUE,
                                            lord2NumberOfSplits = "2", minimumGroupSize = 10,
                                            variables = colnames(binaryItems)))

binomialSems <- function(data, types = rep("nominal", 4), onlyObserved = FALSE) {
  results <- runAnalysis("standardErrorOfMeasurement", dataset = data,
                         options = modifyList(binomialOptions, list(variables.types = types)), makeTests = FALSE)
  rows <- coefficientsOf(results)
  if (onlyObserved) # codes with gaps between the sum scores add empty score rows
    rows <- Filter(function(row) row[["counts"]] > 0, rows)
  do.call(rbind, lapply(rows, function(row) unlist(row[c("lord", "keats", "lord2", "irt")])))
}

binomialReference <- binomialSems(binaryItems)

test_that("Binomial methods give the same SEM for items coded 1/2 as for items coded 0/1", {
  expect_false(anyNA(binomialReference))
  expect_equal(binomialSems(binaryItems + 1), binomialReference)
})

test_that("Binomial methods scale the SEM with the distance between the two item values", {
  observed <- c("lord", "keats", "irt")
  expect_equal(binomialSems(binaryItems * 2, onlyObserved = TRUE)[, observed],
               2 * binomialReference[, observed])
})

test_that("IRT method gives the same SEM for polytomous items whatever the lowest code is", {
  irtOptions  <- modifyList(semOptions, list(anova = FALSE, thorndike = FALSE, irt = TRUE))
  irtSems <- function(data) {
    results <- runAnalysis("standardErrorOfMeasurement", dataset = data, options = irtOptions, makeTests = FALSE)
    vapply(coefficientsOf(results), function(row) row[["irt"]], numeric(1))
  }
  # same answers, coded 1-5 and 0-4: the SEM must belong to the same respondents
  expect_equal(irtSems(asrmNumeric), irtSems(asrmNumeric - 1))
})

test_that("IRT method does not depend on the direction in which dichotomous items are coded", {
  irtSems <- function(data) {
    results <- runAnalysis("standardErrorOfMeasurement", dataset = data,
                           options = modifyList(binomialOptions, list(lord = FALSE, keats = FALSE, lord2 = FALSE)), makeTests = FALSE)
    vapply(coefficientsOf(results), function(row) row[["irt"]], numeric(1))
  }
  # every item flipped: the same respondents, only the sum score is mirrored, so the SEM at the
  # lowest and the highest sum score must swap places as well
  # (mirt fits the two data sets separately, so they agree to its convergence tolerance only)
  expect_equal(irtSems(binaryItems), rev(irtSems(1 - binaryItems)), tolerance = 1e-4)
})

test_that("Binary items that are coded with different numbers give an error", {
  mixedCodes <- binaryItems
  mixedCodes[, 3:4] <- mixedCodes[, 3:4] + 1 # items 1 and 2 use 0/1, items 3 and 4 use 1/2
  results <- runAnalysis("standardErrorOfMeasurement", dataset = mixedCodes, options = binomialOptions, makeTests = FALSE)

  expect_identical(results[["status"]], "validationError")
  expect_match(results[["results"]][["errorMessage"]], "same two numbers")
  expect_match(results[["results"]][["errorMessage"]], "V1, V2 \\(values 0, 1\\)")
  expect_match(results[["results"]][["errorMessage"]], "V3, V4 \\(values 1, 2\\)")
})

test_that("The table says how the items entered the sum score", {
  footnotesOf <- function(results) {
    table <- results[["results"]][["semMainContainer"]][["collection"]][["semMainContainer_coefficientTable"]]
    paste(vapply(table[["footnotes"]], function(footnote) footnote[["text"]], character(1)), collapse = " ")
  }
  polytomous <- runAnalysis("standardErrorOfMeasurement", dataset = asrmNumeric, options = semOptions, makeTests = FALSE)
  dichotomous <- runAnalysis("standardErrorOfMeasurement", dataset = binaryItems, options = binomialOptions, makeTests = FALSE)

  expect_match(footnotesOf(polytomous), "observed values from 1 to 5")
  expect_match(footnotesOf(dichotomous), "dichotomous with the values 0 and 1")
})

test_that("IRT method puts its SEM on the scale of the item values", {
  irtOptions <- modifyList(semOptions, list(anova = FALSE, thorndike = FALSE, irt = TRUE))
  irtSems <- function(data) {
    results <- runAnalysis("standardErrorOfMeasurement", dataset = data, options = irtOptions, makeTests = FALSE)
    rows <- Filter(function(row) row[["counts"]] > 0, coefficientsOf(results)) # codes with gaps add empty score rows
    vapply(rows, function(row) row[["irt"]], numeric(1))
  }
  # same answers with the labels 2 * x + 3: the sum scores spread out by a factor 2, and so must the SEM
  expect_equal(irtSems(2 * asrmNumeric + 3), 2 * irtSems(asrmNumeric), tolerance = 1e-6)
})
