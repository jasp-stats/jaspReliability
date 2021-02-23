options <- analysisOptions("RaterAgreement")

# Dataset from Shrout & Fleiss (1979)
sf_data <- matrix(c(
  9,    2,   5,    8,
  6,    1,   3,    2,
  8,    4,   6,    8,
  7,    1,   2,    6,
  10,   5,   6,    9,
  6,   2,   4,    7
), ncol=4, byrow=TRUE)
dataset <- as.data.frame(sf_data)

# Set options
options$confidenceIntervalValue <- 0.95
options$variables <- paste0("V", 1:4)


test_that("Intraclass Correlation 1 table results are unchanged", {
  options$iccType <- "icc1"
  options$iccRatingAverage <- FALSE

  results <- runAnalysis("RaterAgreement", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.165742286184985, -0.0967218809492438, "ICC1,1", 0.64339874385296)
  )
})

test_that("Intraclass Correlation 2 table results are unchanged", {
  options$iccType <- "icc2"
  options$iccRatingAverage <- FALSE

  results <- runAnalysis("RaterAgreement", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.28976419371026, 0.0429013074043754, "ICC2,1", 0.691071005905065)
  )
})

test_that("Intraclass Correlation 3 table results are unchanged", {
  options$iccType <- "icc3"
  options$iccRatingAverage <- FALSE

  results <- runAnalysis("RaterAgreement", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.714841511511933, 0.411835299105756, "ICC3,1", 0.925833056596294)
  )
})

test_that("Intraclass Correlation 1k table results are unchanged", {
  options$iccType <- "icc1"
  options$iccRatingAverage <- TRUE

  results <- runAnalysis("RaterAgreement", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.442798057590153, -0.54503916, "ICC1,k", 0.878301237198259)
  )
})

test_that("Intraclass Correlation 2k table results are unchanged", {
  options$iccType <- "icc2"
  options$iccRatingAverage <- TRUE

  results <- runAnalysis("RaterAgreement", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.620051021729246, 0.15203742, "ICC2,k", 0.899476869219885)
  )
})

test_that("Intraclass Correlation 3k table results are unchanged", {
  options$iccType <- "icc3"
  options$iccRatingAverage <- TRUE

  results <- runAnalysis("RaterAgreement", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.909315864654189, 0.7368986135979, "ICC3,k", 0.980366125822908)
  )
})

test_that("ICC coefficients match the ones published in Shrout & Fleiss (1979)", {
  # ICC coefficients from the paper
  sf_coefs <- c(
    icc1 = .17,
    icc2 = .29,
    icc3 = .71,
    icc1k = .44,
    icc2k = .62,
    icc3k = .91
  )

  jasp_coefs <- c()
  for (icc_type in names(sf_coefs)) {
    # Remove the k if it is part of icc_type
    options$iccType <- substring(icc_type, 0, 4)
    options$iccRatingAverage <- substring(icc_type, 5) == "k"

    results <- runAnalysis("RaterAgreement", dataset, options)
    jasp_coefs[icc_type] <- results[["results"]][["table"]][["data"]][[1]]$ICC
  }

  # Tolerance is set to just above the rounding difference
  # because R would round 0.715 to 0.71 instead of 0.72
  expect_equal(sf_coefs, jasp_coefs, tolerance = .0051)
})
