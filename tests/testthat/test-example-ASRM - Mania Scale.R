context("Example: ASRM - Mania Scale")

# This test file was auto-generated from a JASP example file.
# The JASP file is stored in the module's examples/ folder.

test_that("unidimensionalReliabilityBayesian results match", {

  # Load from JASP example file
  jaspFile <- testthat::test_path("..", "..", "examples", "ASRM - Mania Scale.jasp")
  opts <- jaspTools::analysisOptions(jaspFile)
  dataset <- jaspTools::extractDatasetFromJASPFile(jaspFile)

  # Encode and run analysis
  encoded <- jaspTools:::encodeOptionsAndDataset(opts, dataset)
  set.seed(1)
  results <- jaspTools::runAnalysis("unidimensionalReliabilityBayesian", encoded$dataset, encoded$options, encodedDataset = TRUE)

  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_itemTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list(0.58987759834391, 0.614342019249597, 0.703045796133977, 0.715955981194311,
     0.796335009943826, 0.808618599153813, "jaspColumn1", 0.606907150313666,
     0.630389241283377, 0.709072993213719, 0.725391162288912, 0.810580370121484,
     0.819269024306521, "jaspColumn2", 0.689128334512277, 0.708826535712507,
     0.777015008669417, 0.787492448326537, 0.846543601821102, 0.858728964076031,
     "jaspColumn3", 0.65428181416818, 0.675553141405965, 0.747343941503956,
     0.759267052305755, 0.834464603163858, 0.847777120079131, "jaspColumn4",
     0.612152559989754, 0.626104718610868, 0.7120001343938, 0.726470848396046,
     0.814362060526522, 0.821369974053776, "jaspColumn5"))

  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_probabilityTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list(0.945614035087719, 0.137022782278187, "McDonald's <unicode>",
     0.977894736842105, 0.234019851705111, "Cronbach's <unicode>"
    ))

  table <- results[["results"]][["stateContainer"]][["collection"]][["stateContainer_scaleTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
    list("McDonald's <unicode>", 0.771766917760584, 0.689713315172189,
     0.999795841484768, 0.850929327827215, "Cronbach's <unicode>",
     0.784672940230205, 0.705832968348325, 1.00291308581205, 0.852091839974414
    ))

})

