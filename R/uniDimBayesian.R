

reliabilityUniDimBayesian <- function(jaspResults, dataset, options) {

  dataset <- .readData(dataset, options)

  if (length(options[["reverseScaledItems"]]) > 0L) {
    dataset <- .reverseScoreItems(dataset, options)
  }


  .checkErrors(dataset, options, Bayes = TRUE)



  model <- .BayesianPreCalc(jaspResults, dataset, options)
  options <- .scaleItemBoxAlign(options)

  model[["derivedOptions"]] <- .BayesianDerivedOptions(options)
  model[["omegaScale"]] <- .BayesianOmegaScale(jaspResults, dataset, options, model)
  model[["omegaItem"]] <- .BayesianOmegaItem(jaspResults, dataset, options, model)
  model[["alphaScale"]] <- .BayesianAlphaScale(jaspResults, dataset, options, model)
  model[["alphaItem"]] <- .BayesianAlphaItem(jaspResults, dataset, options, model)
  model[["lambda2Scale"]] <- .BayesianLambda2Scale(jaspResults, dataset, options, model)
  model[["lambda2Item"]] <- .BayesianLambda2Item(jaspResults, dataset, options, model)
  model[["lambda6Scale"]] <- .BayesianLambda6Scale(jaspResults, dataset, options, model)
  model[["lambda6Item"]] <- .BayesianLambda6Item(jaspResults, dataset, options, model)
  model[["glbScale"]] <- .BayesianGlbScale(jaspResults, dataset, options, model)
  model[["glbItem"]] <- .BayesianGlbItem(jaspResults, dataset, options, model)

  model[["averageInterItemCor"]] <- .BayesianAverageCor(jaspResults, dataset, options, model)
  model[["meanScale"]] <- .BayesianMean(jaspResults, dataset, options, model)
  model[["sdScale"]] <- .BayesianStdDev(jaspResults, dataset, options, model)
  model[["itemRestCor"]] <- .BayesianItemRestCor(jaspResults, dataset, options, model)
  model[["meanItem"]] <- .BayesianMeanItem(jaspResults, dataset, options, model)
  model[["sdItem"]] <- .BayesianSdItem(jaspResults, dataset, options, model)

  .BayesianScaleTable(jaspResults, model, options)
  .BayesianItemTable(jaspResults, model, options)
  .BayesianProbTable(jaspResults, model, options)
  .BayesianPosteriorPlot(jaspResults, model, options)
  .BayesianIfItemPlot(jaspResults, model, options)
  .omegaPosteriorPredictive(jaspResults, model, options)
  .BayesianTracePlot(jaspResults, model, options)
  return()

}

.BayesianDerivedOptions <- function(options) {

  # order of appearance in Bayesrel
  derivedOptions <- list(
    selectedEstimators  = unlist(options[c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale", "glbScale",
                                           "averageInterItemCor", "meanScale", "sdScale")]),
    selectedEstimatorsPlots  = unlist(options[c("omegaScale", "alphaScale", "lambda2Scale", "lambda6Scale",
                                                "glbScale")]),
    itemDroppedSelected = unlist(options[c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item", "glbItem",
                                           "itemRestCor", "meanItem", "sdItem")]),
    itemDroppedSelectedItem = unlist(options[c("omegaItem", "alphaItem", "lambda2Item", "lambda6Item",
                                               "glbItem")]),

    namesEstimators     = list(
      tables = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                 "Greatest Lower Bound", "Average interitem correlation", "mean", "sd"),
      tables_item = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                      gettext("Greatest Lower Bound"), gettext("Item-rest correlation"), gettext("mean"), gettext("sd")),
      coefficients = c("McDonald's \u03C9", "Cronbach's \u03B1", "Guttman's \u03BB2", "Guttman's \u03BB6",
                       gettext("Greatest Lower Bound"), gettext("Item-rest correlation")),
      plots = list(expression("McDonald's"~omega), expression("Cronbach\'s"~alpha), expression("Guttman's"~lambda[2]),
                   expression("Guttman's"~lambda[6]), gettext("Greatest Lower Bound")),
      plotsNoGreek = c("omega", "alpha", "lambda2", "lambda6", "glb")
    )

  )
  return(derivedOptions)
}


.getStateContainerB <- function(jaspResults) {
  if (!is.null(jaspResults[["stateContainer"]]))
    return(jaspResults[["stateContainer"]])

  jaspResults[["stateContainer"]] <- createJaspContainer(dependencies = c("variables", "reverseScaledItems",
                                                                          "noSamples", "noBurnin", "noThin",
                                                                          "noChains", "missingValues", "setSeed",
                                                                          "seed", "disableSampleSave")
  )

  return(jaspResults[["stateContainer"]])
}

.summarizePosteriorStats <- function(samples, ciValue) {
  return(list(
    mean(samples),
    coda::HPDinterval(coda::mcmc(as.vector(samples)), prob = ciValue)
  ))
}

.summarizePosteriorItems <- function(samples, ciValue) {
  return(list(
    colMeans(samples),
    coda::HPDinterval(coda::mcmc(samples), prob = ciValue)
  ))
}


.samplePrior <- function(k, estimate, callback = function(){}) {

  n_samp <- 2e3
  v0 <- k
  k0 <- 1e-10
  t <- diag(k)
  T0 <- solve(t / k0)
  m <- array(0, c(n_samp, k, k))

  for (i in seq_len(n_samp)) {
    m[i, , ] <- LaplacesDemon::rinvwishart(v0, T0)
    callback()
  }

  if (estimate == "alphaScale") {
    prioralpha <- apply(m, MARGIN = 1, Bayesrel:::applyalpha, callback)
    out <- density(prioralpha, from = 0, to = 1, n = 512)
  }
  if (estimate == "lambda2Scale") {
    priorlambda2 <- apply(m, MARGIN = 1, Bayesrel:::applylambda2, callback)
    out <- density(priorlambda2, from = 0, to = 1, n = 512)
  }
  if (estimate == "lambda6Scale") {
    priorlambda6 <- apply(m, MARGIN = 1, Bayesrel:::applylambda6, callback)
    out <- density(priorlambda6, from = 0, to = 1, n = 512)
  }
  if (estimate == "glbScale") {
    priorglb <- Bayesrel:::glbOnArrayCustom(m, callback)
    out <- density(priorglb, from = 0, to = 1, n = 512)
  }

  if (estimate == "omegaScale") {
    H0 <- 1 # prior multiplier matrix for lambdas variance
    l0k <- rep(0, k) # prior lambdas
    a0k <- 1 # prior gamma function for psis
    b0k <- 2 # prior gamma for psi
    prioromega <- numeric(n_samp)
    for (i in seq_len(n_samp)) {
      invpsi <- rgamma(k, a0k, b0k)
      psi <- 1 / invpsi
      lambda <- rnorm(k, l0k, sqrt(psi * H0))
      prioromega[i] <- Bayesrel:::omegaBasic(lambda, psi)
    }
    out <- density(prioromega, from = 0, to = 1, n = 512)
  }

  return(out)
}

.BayesItemDroppedStats <- function(cov_samp, f1 = function(){}, callback = function(){}) {

  dd <- dim(cov_samp)
  out <- matrix(0, dd[1] * dd[2], dd[3])
  cov_samp <- array(cov_samp, c(dd[1] * dd[2], dd[3], dd[3]))
  for (i in seq_len(dd[3])) {
    out[, i] <- apply(cov_samp[, -i, -i], c(1), f1, callback)
  }

  return(out)
}
