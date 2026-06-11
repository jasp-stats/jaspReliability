library(jaspTools)
library(testthat)

pkgload::load_all(quiet = TRUE)

jaspTools::runTestsTravis(module = getwd())
