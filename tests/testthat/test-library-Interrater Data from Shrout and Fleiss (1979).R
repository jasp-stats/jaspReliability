context("Library: Interrater Data from Shrout and Fleiss (1979)")

# This test file was auto-generated from a JASP example file.
# The JASP file is stored in tests/testthat/jaspfiles/library/.

test_that("intraclassCorrelation results match", {

  # Load from JASP example file
  jaspFile <- testthat::test_path("jaspfiles", "library", "Interrater Data from Shrout and Fleiss (1979).jasp")
  opts <- jaspTools::analysisOptions(jaspFile)
  dataset <- jaspTools::extractDatasetFromJASPFile(jaspFile)

  # Encode and run analysis
  encoded <- jaspTools:::encodeOptionsAndDataset(opts, dataset)
  set.seed(1)
  results <- jaspTools::runAnalysis("intraclassCorrelation", encoded$dataset, encoded$options, encodedDataset = TRUE)

  table <- results[["results"]][["table"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list(0.714841511510497, 0.342465949421523, "ICC3,1", 0.945858444552497
    ))

})

