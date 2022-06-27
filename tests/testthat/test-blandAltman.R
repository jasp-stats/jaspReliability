
options <- analysisOptions("blandAltman")
options$variables <- c("contcor1", "contcor2")
options$pairs <- list(c("contcor1", "contcor2"))
options$blandAltmanTable <- TRUE
options$ciDisplay <- TRUE
options$ciShading <- TRUE
results <- jaspTools::runAnalysis("blandAltman", "test.csv", options)

test_that("Bland-Altman plot matches", {
  plotName <- results[["results"]][["plotsBlandAltman"]][["collection"]][["plotsBlandAltman_contcor1 - contcor2"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "baplot")
})
# plot was verified with 'blandr::blandr.draw(debug$contcor1, debug$contcor2)'

test_that("Bland-Altman table results match", {
  table <- results[["results"]][["tabBlandAltman"]][["collection"]][["tabBlandAltman_contcor1 - contcor2"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(1.61924361312177, 1.33231181919255, "Mean difference + 1.96 SD",
                                      1.90617540705099, -0.01713939797, -0.182799546434097, "Mean difference",
                                      0.148520750494097, -1.65352240906177, -1.94045420299099, "Mean difference - 1.96 SD",
                                      -1.36659061513255))
})
