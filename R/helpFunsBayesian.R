


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


.summarizePosteriorItems <- function(samples, ciValue) {
  return(list(
    colMeans(samples),
    coda::HPDinterval(coda::mcmc(samples), prob = ciValue)
  ))
}


.samplePrior <- function(k, estimate, callback = function(){}, k0, df0, a0, b0, m0) {

  n_samp <- 2e3

  if (estimate == "omegaScale") {
    H0 <- 1 # prior multiplier matrix for lambdas variance
    l0k <- rep(m0, k) # prior lambdas
    a0k <- a0 # prior gamma function for psis
    b0k <- b0 # prior gamma for psi
    prioromega <- numeric(n_samp)
    for (i in seq_len(n_samp)) {
      invpsi <- rgamma(k, a0k, b0k)
      psi <- 1 / invpsi
      lambda <- rnorm(k, l0k, sqrt(psi * H0))
      prioromega[i] <- Bayesrel:::omegaBasic(lambda, psi)
      callback()
    }
    out <- density(prioromega, from = 0, to = 1, n = 512)
    return(out)
  }

  v0 <- df0
  k0 <- k0
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
    return(out)

  }
  if (estimate == "lambda2Scale") {
    priorlambda2 <- apply(m, MARGIN = 1, Bayesrel:::applylambda2, callback)
    out <- density(priorlambda2, from = 0, to = 1, n = 512)
    return(out)

  }
  if (estimate == "lambda6Scale") {
    priorlambda6 <- apply(m, MARGIN = 1, Bayesrel:::applylambda6, callback)
    out <- density(priorlambda6, from = 0, to = 1, n = 512)
    return(out)

  }
  if (estimate == "glbScale") {
    priorglb <- Bayesrel:::glbOnArrayCustom(m, callback)
    out <- density(priorglb, from = 0, to = 1, n = 512)
    return(out)

  }

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



.itemRestCor <- function(dataset, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0) {

  ircor_samp <- matrix(0, n.chains * length(seq(1, n.iter - n.burnin, thin)), ncol(dataset))
  for (i in seq(ncol(dataset))) {
    help_dat <- cbind(as.matrix(dataset[, i]), rowMeans(as.matrix(dataset[, -i]), na.rm = TRUE))
    ircor_samp[, i] <- .WishartCorTransform(help_dat, n.iter = n.iter, n.burnin = n.burnin, thin = thin,
                                            n.chains = n.chains, pairwise = pairwise, callback = callback, k0)
  }

  return(ircor_samp)
}

.WishartCorTransform <- function(x, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0) {
  tmp_cov <- Bayesrel:::covSamp(x, n.iter, n.burnin, thin, n.chains, pairwise, callback, k0 = k0, df0 = NULL)$cov_mat
  dd <- dim(tmp_cov)
  tmp_cov <- array(tmp_cov, c(dd[1] * dd[2], dd[3], dd[4]))
  tmp_cor <- apply(tmp_cov, c(1), cov2cor)
  out <- tmp_cor[2, ]
  callback()
  return(out)
}


.colMedians <- function(x) {
  return(apply(x, 2, median))
}


.stdFactorLoads <- function(ll, ee, callback = function(){}) {
  ds <- dim(ll)
  out <- ll
  for (i in seq_len(ds[1])) {
    for (j in seq_len(ds[2])) {
      implV <- diag(ll[i, j, ] %*% t(ll[i, j, ]) + diag(ee[i, j, ]))
      out[i, j, ] <- ll[i, j, ] / sqrt(implV)
      callback()
    }
  }
  return(out)
}

.omegaOnArray <- function(ll, ee, callback = function(){}) {
  ds <- dim(ll)
  out <- matrix(0, ds[1], ds[2])
  for (i in seq_len(ds[1])) {
    for (j in seq_len(ds[2])) {
      out[i, j] <- Bayesrel:::omegaBasic(ll[i, j, ], ee[i, j, ])
      callback()
    }
  }
  return(out)
}

.implCovs <- function(ll, ee, callback = function(){}) {
  ds <- dim(ll)
  out <- array(0, c(ds[1], ds[2], ds[3], ds[3]))
  for (i in seq_len(ds[1])) {
    for (j in seq_len(ds[2])) {
      out[i, j, , ] <- ll[i, j, ] %*% t(ll[i, j, ]) + diag(ee[i, j, ])
      callback()
    }
  }
  return(out)
}


.SRMR <- function(cdat, impl) {
  nvar <- ncol(cdat)
  e <- nvar * (nvar + 1) / 2
  sqrt.d <- 1 / sqrt(diag(cdat))
  D <- diag(sqrt.d, ncol = length(sqrt.d))
  R <- D %*% (cdat - impl) %*% D
  srmr <- sqrt(sum(R[lower.tri(R, diag = TRUE)]^2) / e)
  return(srmr)
}

.LR <- function(cdat, impl, n) {
  k <- ncol(cdat)

  # general function
  LR <- (n-1) * (log(det(impl)) + sum(diag(solve(impl) %*% cdat)) - log(det(cdat)) - k)

  return(max(0, LR))
}

.lavLR <- function(cdat, impl, n) {
  k <- ncol(cdat)

  # lavaan likelihoods
  LOG.2PI <- log(2 * pi)
  logdet_cdat <- log(det(cdat))
  logdet_impl <- log(det(impl))

  loglikh1 <- -n/2 * (k * LOG.2PI + logdet_cdat + k)
  loglikh0 <- -n/2 * (k * LOG.2PI + logdet_impl + k)

  LR <- 2 * (loglikh1 - loglikh0)

  return(max(0, LR))
}


.modelDeviance <- function(cdat, implieds, n) {

  implieds_v <- apply(implieds, c(3, 4), as.vector)
  k <- ncol(cdat)

  srmr_obs <- numeric(nrow(implieds_v))
  LR_obs <- numeric(nrow(implieds_v))

  # srmr_rep <- numeric(nrow(implieds_v))
  # LR_rep <- numeric(nrow(implieds_v))

  for (i in 1:nrow(implieds_v)) {

    srmr_obs[i] <- .SRMR(cdat, implieds_v[i, , ])
    LR_obs[i] <- .LR(cdat, implieds_v[i, , ], n)

    # tmp_dat <- MASS::mvrnorm(n, integer(k), implieds_v[i, , ])

    # # compare the post model predicted cov matrix with original cov matrix:
    # srmr_rep[i] <- .SRMR(.rescale(cov(tmp_dat), n), cdat)
    # LR_rep[i] <- .LR(.rescale(cov(tmp_dat), n), cdat, n)

    # # compare the post model predicted cov matrix with the post model implied cov matrix
    # srmr_rep[i] <- .SRMR(.rescale(cov(tmp_dat), n), implieds_v[i, , ])
    # LR_rep[i] <- .LR(.rescale(cov(tmp_dat), n), implieds_v[i, , ], n)

  }
  return(list(srmr_obs = srmr_obs, LR_obs = LR_obs
              # ,srmr_rep = srmr_rep, LR_rep = LR_rep
              ))
}

.RMSEA <- function(chisq, df, n) {
  rmsea <- sqrt(max(0, (chisq - df) / (df * n)))
  return(rmsea)
}

.BRMSEA <- function(chisq, p, pD, n) {
  dChisq <- (chisq - p)
  dChisq[dChisq < 0] <- 0
  rmsea <- sqrt(dChisq / ((p - pD) * n))
  return(rmsea)
}

.rescale <- function(cc, n) {
  return(cc * ((n - 1) / n))
}
