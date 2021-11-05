options <- analysisOptions("intraclassCorrelation")
# ==== Ensure results are unchanged on JASP debug data ====

# Set options
options$confidenceIntervalValue <- 0.95
options$variables <- c("contcor1", "contcor2", "contNormal")
dataset <- "test.csv"

test_that("Intraclass Correlation 1 table results are unchanged", {
  options$iccType <- "icc1"
  options$iccRatingAverage <- FALSE

  results <- runAnalysis("intraclassCorrelation", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.27073667055252, 0.166492336930429, "ICC1,1", 0.380050661447667)
  )
})

test_that("Intraclass Correlation 2 table results are unchanged", {
  options$iccType <- "icc2"
  options$iccRatingAverage <- FALSE

  results <- runAnalysis("intraclassCorrelation", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.273772717560708, 0.170415127396369, "ICC2,1", 0.382274954775931)
  )
})

test_that("Intraclass Correlation 3 table results are unchanged", {
  options$iccType <- "icc3"
  options$iccRatingAverage <- FALSE

  results <- runAnalysis("intraclassCorrelation", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.277235249748641, 0.17271511380534, "ICC3,1", 0.386496184673786)
  )
})

test_that("Intraclass Correlation 1k table results are unchanged", {
  options$iccType <- "icc1"
  options$iccRatingAverage <- TRUE

  results <- runAnalysis("intraclassCorrelation", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.526905000559601, 0.37470574162312, "ICC1,k", 0.647776335096136)
  )
})

test_that("Intraclass Correlation 2k table results are unchanged", {
  options$iccType <- "icc2"
  options$iccRatingAverage <- TRUE

  results <- runAnalysis("intraclassCorrelation", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.530723127116256, 0.38129015985557, "ICC2,k", 0.649924866460166)
  )
})

test_that("Intraclass Correlation 3k table results are unchanged", {
  options$iccType <- "icc3"
  options$iccRatingAverage <- TRUE

  results <- runAnalysis("intraclassCorrelation", dataset, options)
  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(
    table,
    list(0.535041192171159, 0.385114984621821, "ICC3,k", 0.653972670197125)
  )
})

# ==== Replicate Shrout & Fleiss ====

# Dataset from Shrout & Fleiss (1979)
sf_data <- matrix(c(
  9,    2,   5,    8,
  6,    1,   3,    2,
  8,    4,   6,    8,
  7,    1,   2,    6,
  10,   5,   6,    9,
  6,   2,   4,    7
), ncol = 4, byrow = TRUE)
dataset <- as.data.frame(sf_data)

# Set options
options$confidenceIntervalValue <- 0.95
options$variables <- paste0("V", 1:4)

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

    results <- runAnalysis("intraclassCorrelation", dataset, options)
    jasp_coefs[icc_type] <- results[["results"]][["table"]][["data"]][[1]]$ICC
  }

  # Tolerance is set to just above the rounding difference
  # because R would round 0.715 to 0.71 instead of 0.72
  expect_equal(sf_coefs, jasp_coefs, tolerance = .0051)
})
