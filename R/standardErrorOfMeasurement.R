#
# Copyright (C) 2021 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' @export
standardErrorOfMeasurement <- function(jaspResults, dataset, options) {


  sink(file = "~/Downloads/log.txt")
  on.exit(sink(NULL))

  ready <- length(options[["variables"]]) > 1

  dataset <- .semReadData(dataset, options)

  .semErrorCheck(dataset, options, ready)

  .semCreateMainContainer(jaspResults, options)

  .semComputeCoefficients(jaspResults, dataset, options, ready)

  .semCoefficientsTable(jaspResults, dataset, options, ready)

  .semHistPlot(jaspResults, dataset, options, ready)

  .semPointPlots(jaspResults, dataset, options, ready)

  return()
}

# Read in the dataset
.semReadData <- function(dataset, options) {

  variables <- unlist(options[["variables"]])

  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.numeric = variables, exclude.na.listwise = variables)
  }

  return(dataset)
}


# check errors
.semErrorCheck <- function(dataset, options, ready) {

  if (!ready) return()
  .semCheckCategories <- function() {

    nc <- length(unique(c(as.matrix(dataset))))
    ncSupplied <- options[["responseCategories"]]

    if (nc > ncSupplied) {
      return(gettext("It seems the number of response categories in the data is higher than the ones you specified."))
    } else {
      return(NULL)
    }

  }

  .hasErrors(dataset = dataset,
             type = c('missingValues', "variance", "infinity", "variance"),
             missingValues.target = options$variables,
             exitAnalysisIfErrors = TRUE,
             custom = .semCheckCategories)
}



# create Main container
.semCreateMainContainer <- function(jaspResults, options) {

  if (!is.null(jaspResults[["semMainContainer"]])) return()

  semMainContainer <- createJaspContainer(dependencies = c("variables", "responseCategories"))
  jaspResults[["semMainContainer"]] <- semMainContainer


  return()
}


##### Computation #####
.semComputeCoefficients <- function(jaspResults, dataset, options, ready) {

  if (!ready) return()

  nc <- options[["responseCategories"]]
  counts <- .counts(dataset, nc)
  countsState <- createJaspState(counts, dependencies = NULL)
  jaspResults[["semMainContainer"]][["countsState"]] <- countsState

  # first the average sem
  if (is.null(jaspResults[["semMainContainer"]][["averageState"]])) {
    if (options[["userReliability"]]) {
      rel <- options[["reliabilityValue"]]
    } else {
      rel <- Bayesrel:::applyalpha(cov(dataset))
    }
    average <- as.numeric(sd(rowSums(dataset)) * sqrt(1 - rel))
    averageState <- createJaspState(average, dependencies = NULL)
    jaspResults[["semMainContainer"]][["averageState"]] <- averageState
  }

  if (options[["thorndike"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["thorndikeState"]])) {
      out <- .semThorn(dataset, nc = nc, caseMin = options[["minimumGroupSize"]], splits = NULL, counts)
      thorndikeState <- createJaspState(out, dependencies = c("thorndike", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["thorndikeState"]] <- thorndikeState
    }
  }

  if (options[["feldt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["feldtState"]])) {
      out <- .semFeldt(dataset, K = options[["feldtNumberOfSplits"]], nc = nc, caseMin = options[["minimumGroupSize"]],
                       splits = NULL, counts)
      feldtState <- createJaspState(out, dependencies = c("feldt", "feldtNumberOfSplits", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["feldtState"]] <- feldtState
    }
  }

  if (options[["mollenkopfFeldt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["mfState"]])) {
      out <- .semMF(dataset, K = options[["mollenkopfFeldtNumberOfSplits"]], nc = nc,
                    n_poly = options[["mollenkopfFeldtPolyDegree"]], splits = NULL, counts)
      mfState <- createJaspState(out, dependencies = c("mollenkopfFeldt", "mollenkopfFeldtNumberOfSplits",
                                                          "mollenkopfFeldtPolyDegree"))
      jaspResults[["semMainContainer"]][["mfState"]] <- mfState
    }
  }

  if (options[["anova"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["anovaState"]])) {
      out <- .semAnova(dataset, nc = nc, caseMin = options[["minimumGroupSize"]], counts)
      anovaState <- createJaspState(out, dependencies = c("anova", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["anovaState"]] <- anovaState
    }
  }

  # works for both data types
  if (options[["irt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["irtState"]])) {
      out <- .semIrt(dataset, nc = nc)
      irtState <- createJaspState(out, dependencies = "irt")
      jaspResults[["semMainContainer"]][["irtState"]] <- irtState
    }
  }

  # only for binary data
  if (options[["lord"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["lordState"]])) {
      out <- .semLord(ncol(dataset))
      lordState <- createJaspState(out, dependencies = "lord")
      jaspResults[["semMainContainer"]][["lordState"]] <- lordState
    }
  }

  if (options[["keats"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["keatsState"]])) {
      out <- .semKeats(dataset, options)
      keatsState <- createJaspState(out, dependencies = "keats")
      jaspResults[["semMainContainer"]][["keatsState"]] <- keatsState
    }
  }

  if (options[["lord2"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["lord2State"]])) {
      out <- .semLord2(dataset, options[["lord2NumberOfSplits"]], counts)
      lord2State <- createJaspState(out, dependencies = c("lord2", "lord2NumberOfSplits"))
      jaspResults[["semMainContainer"]][["lord2State"]] <- lord2State
    }
  }

 return()
}

##### Output tables #####
.semCoefficientsTable <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["coefficientTable"]])) return()

  coefficientTable <- createJaspTable(gettext("Standard error of measurement"),
                                      dependencies = c("thorndike", "feldt", "mollenkopfFeldt",
                                                       "anova", "irt", "lord", "keats", "lord2"))
  jaspResults[["semMainContainer"]][["coefficientTable"]] <- coefficientTable

  coefficientTable$addColumnInfo(name = "score", title = gettext("Score"), type = "string")
  coefficientTable$addColumnInfo(name = "average", title = gettext("Average"), type = "number")

  if (!ready || jaspResults[["semMainContainer"]]$getError()) return()

  dtFill <- data.frame(score = "all")
  average <- jaspResults[["semMainContainer"]][["averageState"]]$object
  dtFill$average <- average

  if (any(c(options[["thorndike"]], options[["feldt"]], options[["mollenkopfFeldt"]], options[["anova"]],
            options[["irt"]], options[["lord"]], options[["keats"]], options[["lord2"]]))) {

    # the repetition of this seems annoying...
    coefficientTable <- createJaspTable(gettext("Standard error of measurement"),
                                        dependencies = c("thorndike", "feldt", "mollenkopfFeldt",
                                                         "anova", "irt", "feldtNumberOfSplits",
                                                         "mollenkopfFeldtNumberOfSplits", "mollenkopfFeldtPolyDegree",
                                                         "minimumGroupSize",
                                                         "lord", "keats", "lord2", "lord2NumberOfSplits"))
    jaspResults[["semMainContainer"]][["coefficientTable"]] <- coefficientTable

    coefficientTable$addColumnInfo(name = "score", title = gettext("Score"), type = "string")
    coefficientTable$addColumnInfo(name = "counts", title = gettext("Counts"), type = "string")
    coefficientTable$addColumnInfo(name = "average", title = gettext("Average"), type = "number")

    dtFill <- data.frame(score = 0:(ncol(dataset) * (options[["responseCategories"]] - 1)))
    counts <- jaspResults[["semMainContainer"]][["countsState"]]$object
    dtFill$counts <- counts
    dtFill$average <- average

    if (options[["thorndike"]]) {
      out <- jaspResults[["semMainContainer"]][["thorndikeState"]]$object
      coefficientTable$addColumnInfo(name = "thorndike", title = gettext("Thorndike"), type = "number")
      dtFill$thorndike <- out[, 2]
    }
    if (options[["feldt"]]) {
      out <- jaspResults[["semMainContainer"]][["feldtState"]]$object
      coefficientTable$addColumnInfo(name = "feldt", title = gettext("Feldt"), type = "number")
      dtFill$feldt <- out[, 2]
    }
    if (options[["mollenkopfFeldt"]]) {
      out <- jaspResults[["semMainContainer"]][["mfState"]]$object
      coefficientTable$addColumnInfo(name = "moll", title = gettext("Mollenkopf-Feldt"), type = "number")
      dtFill$moll <- out[, 2]
    }
    if (options[["anova"]]) {
      out <- jaspResults[["semMainContainer"]][["anovaState"]]$object
      coefficientTable$addColumnInfo(name = "anova", title = gettext("ANOVA"), type = "number")
      dtFill$anova <- out[, 2]
    }

    if (options[["irt"]]) {
      out <- jaspResults[["semMainContainer"]][["irtState"]]$object
      # because the irt method does not produce score groups but a continous score along the scale, we need to group it
      outDt <- data.frame(scores = out$scores, sems = out$sems)
      outOrdered <- outDt[order(outDt$scores, decreasing = FALSE), ]
      uniqueOutOrdered <- unique(outOrdered)
      print(str(out))
      irtBinned <- approx(uniqueOutOrdered$scores, uniqueOutOrdered$sems, xout = out$discScores)$y

      coefficientTable$addColumnInfo(name = "irt", title = gettext("IRT"), type = "number")
      dtFill$irt <- irtBinned
    }

    if (options[["lord"]]) {
      out <- jaspResults[["semMainContainer"]][["lordState"]]$object
      coefficientTable$addColumnInfo(name = "lord", title = gettext("Lord"), type = "number")
      dtFill$lord <- out[, 2]
    }

    if (options[["keats"]]) {
      out <- jaspResults[["semMainContainer"]][["keatsState"]]$object
      coefficientTable$addColumnInfo(name = "keats", title = gettext("Keats"), type = "number")
      dtFill$keats <- out[, 2]
    }

    if (options[["lord2"]]) {
      out <- jaspResults[["semMainContainer"]][["lord2State"]]$object
      coefficientTable$addColumnInfo(name = "lord2", title = gettext("Lord-2"), type = "number")
      dtFill$lord2 <- out[, 2]
    }
  }

  coefficientTable$setData(dtFill)

  return()

}


#### Output plots ####

.semHistPlot <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["histPlot"]]) || !options[["histogramCounts"]]) return()

  if (!ready || jaspResults[["semMainContainer"]]$getError()) return()

  ss <- data.frame(scores = rowSums(dataset))
  p <- ggplot2::ggplot(ss) +
    ggplot2::geom_histogram(ggplot2::aes(scores), bins = nrow(unique(ss)), binwidth = .5) +
    ggplot2::ylab(gettext("Counts")) +
    ggplot2::xlab(gettext("Scores"))
  p <- jaspGraphs::themeJasp(p)

  histPlot <- createJaspPlot(plot = p, title = gettext("Histogram of counts per score group"))
  jaspResults[["semMainContainer"]][["histPlot"]] <- histPlot

  return()

}

.semPointPlots <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["pointPlots"]]) || !options[["pointPlots"]]) return()

  if (!ready || jaspResults[["semMainContainer"]]$getError()) return()

  pointPlotsContainer <- createJaspContainer(title = gettext("Point plots"))
  pointPlotsContainer$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]])
  jaspResults[["semMainContainer"]][["pointPlotsContainer"]] <- pointPlotsContainer

  if (options[["thorndike"]]) {
    pl <- .makeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["thorndikeState"]]$object,
                               title = gettext("Thorndike method point plot"), dep = "thorndike")

    pointPlotsContainer[["thorndikePlot"]] <- pl
  }


  return()

}

.makeSinglePointPlot <- function(resultsObject, title, dep) {

  dat <- as.data.frame(resultsObject)
  colnames(dat) <- c("score", "sem", "counts", "merged")
  pl <- ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(x = score, y = sem), size = 2.5) +
    ggplot2::labs(x = "Test Score", y = "sem")

  outPlot <- createJaspPlot(jaspGraphs::themeJasp(pl), title = title, dependencies = dep)

  return(outPlot)

}


##### sem compute functions #####

.counts <- function(X, nc) {
  nit <- ncol(X)
  S <- rowSums(X)
  scores <- 0:(nit * (nc - 1))
  # create a matrix that counts the scores and also checks the caseMins
  counts <- numeric(length(scores))
  # first count all the different sum score occurences
  for (i in seq_len(length(counts))) {
    counts[i] <- sum(S == scores[i])
  }
  return(counts)
}


.semIrt <- function(X, nc) {

  if (nc == 2) {
    ity <- "2PL"
  } else {
    ity <- "graded"
  }

  res <- mirt::mirt(X, 1, itemtype=ity)
  # x <- mirt::fscores(res)
  x <- seq(-5, 5, by = 0.1)

  cofs <- mirt::coef(res, IRTpars=TRUE)
  cofs$GroupPars <- NULL
  cof_mat <- t(sapply(cofs, function(x) x))
  cof_mat <- cof_mat[, c(2:nc, 1)]
  colnames(cof_mat) <- c("a", paste0("b", 1:(nc-1)))

  a <- cof_mat[, "a"]
  nit <- length(a)
  outIRT <- list(scores = numeric(length(x)), sems = numeric(length(x)))

  if (nc == 2) {
    b <- cof_mat[, "b1"]
    for (i in 1:length(x)) {
      outIRT$scores[i] <- sum(plogis(a * (x[i] - b)))
      outIRT$sems[i] <- sqrt(sum(plogis(a * (x[i] - b)) * (1 - plogis(a * (x[i] - b)))))
    }

  } else {
    b <- cbind(-Inf, cof_mat[, 1:(nc - 1)], Inf)
    probs <- matrix(NA, nit, nc)
    for (i in 1:length(x)) {
      for (j in 1:(nc)) {
        probs[, j] <- plogis(a * (x[i] - b[, j])) - plogis(a * (x[i] - b[, j + 1]))
      }
      ev <- rowSums(t(t(probs) * (0:(nc - 1))))
      dev <- (matrix(0:(nc - 1), nit, nc, TRUE) - ev)^2

      outIRT$scores[i] <- sum(ev)
      outIRT$sems[i] <- sqrt(sum(dev * probs))
    }
  }


  outIRT$discScores <- 0:(nit * (nc - 1))

  return(outIRT)
}



.semThorn <- function(X, nc, caseMin = 3, splits = NULL, counts) {

  nit <- ncol(X)
  N <- nrow(X)

  S <- rowSums(X)

  if (is.null(splits)) {
    # when no specific splits are defined, we just evenly distribute the items among the number of splits
    k <- split(seq_len(nit), 1:2)
  }

  partSUMS1 <- rowSums(X[, k[[1]], drop = FALSE])
  partSUMS2 <- rowSums(X[, k[[2]], drop = FALSE])

  out <- matrix(NA, (nit * (nc - 1)) + 1, 4)
  scores <- 0:(nit * (nc - 1))
  out[, 1] <- scores
  out[, 3] <- counts
  out[, 4] <- FALSE
  ii <- 1

  while (ii <= nrow(out)) {
    if (counts[ii] >= caseMin) {
      ind <- which(S == scores[ii])
      out[ii, 2] <- sd(partSUMS1[ind] - partSUMS2[ind])
      ii <- ii + 1
    } else {
      csum <- cumsum(counts[ii:length(counts)])
      # index of the first score that surpasses the min size in the submatrix with the cumulative sums
      firstInd <- which(csum >= caseMin)
      # at the higher end of the score groups, we might encounter that merging groups will not reach the
      # minsize anymore, so we need to check
      if (length(firstInd) > 0) {
        firstInd <- firstInd[1]
        # actual index in the scores and out matrix
        nextInd <- ii + (firstInd - 1)
        ind <- which(S %in% scores[ii:nextInd])
        out[ii:nextInd, 2] <- sd(partSUMS1[ind] - partSUMS2[ind])
        out[ii:nextInd, 4] <- TRUE
        ii <- nextInd + 1
      } else {
        break
      }
    }
  }

  # now lets also do kind of the same from the back to capture the cases at the end that dont reach the min size
  backInd <- nrow(out)
  if (counts[backInd] < caseMin) {
    csum <- cumsum(counts[backInd:1])
    # index of the first score that surpasses the min size in the submatrix with the cumulative sums
    firstInd <- which(csum >= caseMin)
    firstInd <- firstInd[1]
    # actual index in the scores and out matrix
    nextInd <- backInd - (firstInd - 1)
    ind <- which(S %in% scores[backInd:nextInd])
    out[backInd:nextInd, 2] <- sd(partSUMS1[ind] - partSUMS2[ind])
    out[backInd:nextInd, 4] <- TRUE
  }

  return(out)
}


.semFeldt <- function(X, K, nc, caseMin, splits = NULL, counts) {

  nit <- ncol(X)
  N <- nrow(X)

  S <- rowSums(X)
  partSUMS <- matrix(NA, N, K)

  if (K < nit) {
    if (is.null(splits)) {
      # when no specific splits are defined, we just evenly distribute the items among the number of splits
      k <- split(seq_len(nit), 1:K)
      for (i in 1:K) {
        partSUMS[, i] <- rowSums(X[, k[[i]], drop = FALSE])
      }
    } else {
      # TODO: implement the specified splits
    }

    d <- nit / mean(sapply(k, length))
  }

  if (K == nit) {
    for (i in 1:K) {
      partSUMS[, i] <- X[, i]
    }
    d <- nit
  }

  out <- matrix(NA, (nit * (nc - 1)) + 1, 4)
  scores <- 0:(nit * (nc - 1))
  out[, 1] <- scores
  out[, 3] <- counts
  out[, 4] <- FALSE
  ii <- 1

  while (ii <= nrow(out)) {
    if (counts[ii] >= caseMin) {
      ind <- which(S == scores[ii])
      mean_diff <- partSUMS[ind, ] - rowMeans(partSUMS[ind, ]) - matrix(colMeans(partSUMS[ind, ]), length(ind), K, TRUE) + mean(partSUMS[ind, ])
      out[ii, 2] <- sqrt(d * sum(rowSums(mean_diff^2) / (K - 1)) / length(ind))
      ii <- ii + 1
    } else {
      csum <- cumsum(counts[ii:length(counts)])
      # index of the first score that surpasses the min size in the submatrix with the cumulative sums
      firstInd <- which(csum >= caseMin)
      # at the higher end of the score groups, we might encounter that merging groups will not reach the
      # minsize anymore, so we need to check
      if (length(firstInd) > 0) {
        firstInd <- firstInd[1]
        # actual index in the scores and out matrix
        nextInd <- ii + (firstInd - 1)
        ind <- which(S %in% scores[ii:nextInd])
        mean_diff <- partSUMS[ind, ] - rowMeans(partSUMS[ind, ]) - matrix(colMeans(partSUMS[ind, ]), length(ind), K, TRUE) + mean(partSUMS[ind, ])
        out[ii:nextInd, 2] <- sqrt(d * sum(rowSums(mean_diff^2) / (K - 1)) / length(ind))
        out[ii:nextInd, 4] <- TRUE
        ii <- nextInd + 1
      } else {
        break
      }
    }
  }

  # now lets also do kind of the same from the back to capture the cases at the end that dont reach the min size
  backInd <- nrow(out)
  if (counts[backInd] < caseMin) {
    csum <- cumsum(counts[backInd:1])
    # index of the first score that surpasses the min size in the submatrix with the cumulative sums
    firstInd <- which(csum >= caseMin)
    firstInd <- firstInd[1]
    # actual index in the scores and out matrix
    nextInd <- backInd - (firstInd - 1)
    ind <- which(S %in% scores[backInd:nextInd])
    mean_diff <- partSUMS[ind, ] - rowMeans(partSUMS[ind, ]) - matrix(colMeans(partSUMS[ind, ]), length(ind), K, TRUE) + mean(partSUMS[ind, ])
    out[backInd:nextInd, 2] <- sqrt(d * sum(rowSums(mean_diff^2) / (K - 1)) / length(ind))
    out[backInd:nextInd, 4] <- TRUE
  }

  return(out)
}


.semMF <- function(X, K, nc, n_poly, splits = NULL, counts) {


  nit <- ncol(X)
  N <- nrow(X)

  S <- rowSums(X)
  partSUMS <- matrix(NA, N, K)

  if (K < nit) {
    if (is.null(splits)) {
      # when no specific splits are defined, we just evenly distribute the items among the number of splits
      k <- split(seq_len(nit), 1:K)
      for (i in 1:K) {
        partSUMS[, i] <- rowSums(X[, k[[i]], drop = FALSE])
      }
    } else {
      # TODO: implement the specified splits
    }

    d <- nit / mean(sapply(k, length))
  }

  if (K == nit) {
    for (i in 1:K) {
      partSUMS[, i] <- X[, i]
    }
    d <- nit
  }

  out <- matrix(NA, (nit * (nc - 1)) + 1, 4)
  scores <- 0:(nit * (nc - 1))
  out[, 1] <- scores
  out[, 3] <- counts

  rawDiffK <- d *
    rowSums((partSUMS - matrix(colMeans(partSUMS), N, K, TRUE) - rowMeans(partSUMS) + mean(partSUMS))^2) /
    (K - 1)
  betaK <- coef(lm(rawDiffK ~ poly(S, n_poly, raw = TRUE)))
  scrs <- sqrt(betaK[1] + rowSums(matrix(betaK[-1], length(scores), n_poly, TRUE) * poly(scores, n_poly, raw = TRUE)))
  out[, 2] <- scrs

  return(out)
}


.semAnova <- function(X, nc, caseMin, counts) {

  nit <- ncol(X)
  N <- nrow(X)

  S <- rowSums(X)

  out <- matrix(NA, (nit * (nc - 1)) + 1, 4)
  scores <- 0:(nit * (nc - 1))
  out[, 1] <- scores
  out[, 3] <- counts
  out[, 4] <- FALSE
  ii <- 1

  while (ii <= nrow(out)) {
    if (counts[ii] >= caseMin) {
      ind <- which(S == scores[ii])
      out[ii, 2] <- sqrt(nit / (nit - 1) * sum(diag(cov(X[ind, ]))))
      ii <- ii + 1
    } else {
      csum <- cumsum(counts[ii:length(counts)])
      # index of the first score that surpasses the min size in the submatrix with the cumulative sums
      firstInd <- which(csum >= caseMin)
      # at the higher end of the score groups, we might encounter that merging groups will not reach the
      # minsize anymore, so we need to check
      if (length(firstInd) > 0) {
        firstInd <- firstInd[1]
        # actual index in the scores and out matrix
        nextInd <- ii + (firstInd - 1)
        ind <- which(S %in% scores[ii:nextInd])
        out[ii:nextInd, 2] <- sqrt(nit / (nit - 1) * sum(diag(cov(X[ind, ]))))
        out[ii:nextInd, 4] <- TRUE
        ii <- nextInd + 1
      } else {
        break
      }
    }
  }

  # now lets also do kind of the same from the back to capture the cases at the end that dont reach the min size
  backInd <- nrow(out)
  if (counts[backInd] < caseMin) {
    csum <- cumsum(counts[backInd:1])
    # index of the first score that surpasses the min size in the submatrix with the cumulative sums
    firstInd <- which(csum >= caseMin)
    firstInd <- firstInd[1]
    # actual index in the scores and out matrix
    nextInd <- backInd - (firstInd - 1)
    ind <- which(S %in% scores[backInd:nextInd])
    out[backInd:nextInd, 2] <- sqrt(nit / (nit - 1) * sum(diag(cov(X[ind, ]))))
    out[backInd:nextInd, 4] <- TRUE
  }

  return(out)
}


.semLord <- function(nit) {
  x <- 0:nit
  out_binom <- sqrt((x) * (nit - (x)) / (nit - 1))
  return(cbind(0:nit, out_binom))
}

.semKeats <- function(X, options) {

  nit <- ncol(X)
  x <- 0:nit
  S <- rowSums(X)
  KR21 <- nit / (nit - 1) * (1 - (mean(S) * (nit - mean(S))) / (nit * var(S)))
  if (options[["userReliability"]]) {
    rel <- options[["reliabilityValue"]]
  } else {
    rel <- Bayesrel:::applyalpha(cov(X))
  }
  correction <- (1 - rel) / (1 - KR21)
  ou_binom2 <- sqrt((x) * (nit - (x)) / (nit - 1) * correction)

  return(cbind(0:nit, ou_binom2))
}


# X = dataset, K = number of splits, counts = counts per score group
.semLord2 <- function(X, K, counts, splits = NULL) {

  nit <- ncol(X)
  N <- nrow(X)
  S <- rowSums(X)

  partSUMS <- matrix(NA, N, K)

  if (K < nit) {
    if (is.null(splits)) {
      # when no specific splits are defined, we just evenly distribute the items among the number of splits
      k <- split(seq_len(nit), 1:K)
      for (i in 1:K) {
        partSUMS[, i] <- rowSums(X[, k[[i]], drop = FALSE])
      }
      cc <- sapply(k, length)
    } else {
      # TODO: implement the specified splits
    }
  }

  if (K == nit) {
    for (j in 1:K) {
      partSUMS[, j] <- X[, j]
    }
    cc <- rep(1, K)
  }

  scores <- 0:nit
  out2 <- numeric(length(scores))
  ii <- 1
  while (ii <= length(scores)) {
    ind <- which(S == scores[ii])
    ccmat <- matrix(cc, length(ind), length(cc), byrow = TRUE)
    out2[ii] <- sqrt(mean(rowSums(partSUMS[ind, ] * (ccmat - partSUMS[ind, ]) / (ccmat - 1))))
    ii <- ii + 1
  }

  return(cbind(0:nit, out2))

}



