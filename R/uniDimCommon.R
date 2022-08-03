# This is a temporary fix
# TODO: remove it when R will solve this problem!
gettextf <- function(fmt, ..., domain = NULL)  {
  return(sprintf(gettext(fmt, domain = domain), ...))
}

.readData <- function(dataset, options) {

  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(
      columns.as.numeric  = variables,
      exclude.na.listwise = if (options[["missingValues"]] == "excludeCasesListwise") variables else NULL
    )
  }
  return(dataset)
}


.checkErrors <- function(dataset, options, Bayes = FALSE) {

  if (Bayes) {
    # check for sensible MCMC values
    .checkMCMCBounds <- function() {
      # are the remaining samples after burnin and thinning larger than 2?
      if (ceiling((options[["noSamples"]] - options[["noBurnin"]]) / options[["noThin"]]) < 2) {
        return(gettext("Too few MCMC samples will remain after the removal of the burnin samples and thinning,
                     please increase No.samples."))
      }
      return(NULL)
    }


  } else {
    .checkMCMCBounds <- NULL
  }

  jaspBase::.hasErrors(dataset = dataset,
             type = c("infinity", "variance", "observations"),
             observations.amount = " < 3",
             infinity.target = options$variables,
             variance.target = options$variables,
             observations.target = options$variables,
             custom = .checkMCMCBounds,
             exitAnalysisIfErrors = TRUE)
}


.checkLoadings <- function(dataset, variables) {
  if (ncol(dataset >= 2)) {
    # check for negative loadings:
    prin <- psych::principal(dataset)
    idx <- prin[["loadings"]] < 0
    sidx <- sum(idx)
    if (sidx == 0) {
      footnote <- ""
    } else {
      footnote <- sprintf(ngettext(sidx,
                                   "The following item correlated negatively with the scale: %s. ",
                                   "The following items correlated negatively with the scale: %s. "),
                          paste(variables[idx], collapse = ", "))
    }

    # check for perfect correlations:
    cr <- cor(dataset, use = "pairwise.complete.obs")
    cr[lower.tri(cr, diag = TRUE)] <- 0
    pos <- which(round(cr, 3) == 1, arr.ind = TRUE)
    if (length(pos) == 0) {
      footnote <- gettextf("%s", footnote)
    } else {
      for (i in seq_len(nrow(pos))) {
        footnote <- gettextf("%s Variables %s and %s correlated perfectly. ",
                             footnote, variables[pos[i, 1]], variables[pos[i, 2]])
      }
    }

    return(footnote)

  } else {
    return(.atLeast2Variables())
  }
}

.reverseScoreItems <- function(dataset, options) {
  dataset_rev <- as.matrix(dataset) # fails for string factors!
  cols <- match(unlist(options[["reverseScaledItems"]]), colnames(dataset))
  total <- apply(as.matrix(dataset[, cols]), 2, min, na.rm = TRUE) +
    apply(as.matrix(dataset[, cols]), 2, max, na.rm = TRUE)
  dataset_rev[, cols] <- matrix(rep(total, nrow(dataset)), nrow(dataset), length(cols), byrow = TRUE) -
    as.matrix(dataset[, cols])
  return(as.data.frame(dataset_rev))
}


.cov2cor.callback <- function(C, callback) {
  callback()
  return(cov2cor(C))
}

# calculate the kullback leibler distance between two samples
.KLD.statistic <- function(x, y) {
  # transform the samples to PDFs:
  xdf <- .get_approx_density(x)
  ydf <- .get_approx_density(y)

  xx <- seq(0, 1, length.out = 1e4)
  t <- LaplacesDemon::KLD(xdf(xx), ydf(xx))
  t$sum.KLD.py.px
}

# calculate the kolomogorov smirnov distances between some samples and the original sample
.ks.test.statistic <- function(x, y) {
  t <- stats::ks.test(x, y)
  t$statistic
}

# konvert empirical samples to cumulative density functions
.get_approx_density <- function(x) {
  d <- density(x, n = 2^12)
  f <- approxfun(d$x, d$y, yleft = 0, yright = 0)
  c <- integrate(f, 0, 1)$value
  return(
    function(x) {
      return(f(x) / c)
    }
  )
}



# change options when scale box is unchecked
.scaleItemBoxAlign <- function(options) {
  opts <- options
  if (!options[["omegaScale"]])
    opts[["omegaItem"]] <- FALSE
  if (!options[["alphaScale"]])
    opts[["alphaItem"]] <- FALSE
  if (!options[["lambda2Scale"]])
    opts[["lambda2Item"]] <- FALSE
  if (!options[["lambda6Scale"]])
    opts[["lambda6Item"]] <- FALSE
  if (!options[["glbScale"]])
    opts[["glbItem"]] <- FALSE

  return(opts)

}

.addFootnoteReverseScaledItems <- function(options) {
  out <- sprintf(ngettext(length(options[["reverseScaledItems"]]),
                          "The following item was reverse scaled: %s. ",
                          "The following items were reverse scaled: %s. "),
                 paste(options[["reverseScaledItems"]], collapse = ", "))
  return(out)
}

.atLeast2Variables <- function() {
  return(gettext("Please enter at least 2 variables to do an analysis"))
}

.is.empty <- function(model) {
  !is.null(model[["empty"]])
}


.sampleListHelper <- function(ll, llName) {
  samps <- lapply(ll, `[[`, llName)
  samps[lengths(samps) == 0] <- NULL
  return(samps)
}
