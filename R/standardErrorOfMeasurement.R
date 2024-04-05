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
      thorndikeState <- createJaspState(out, dependencies = "thorndike")
      jaspResults[["semMainContainer"]][["thorndikeState"]] <- thorndikeState
    }
  }

  if (options[["feldt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["feldtState"]])) {
      out <- .semFeldt(dataset, nc = nc, caseMin = options[["minimumGroupSize"]], splits = NULL, counts)
      feldtState <- createJaspState(out, dependencies = "thorndike")
      jaspResults[["semMainContainer"]][["feldtState"]] <- feldtStatte
    }
  }

 return()
}

##### Output tables #####
.semCoefficientsTable <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["coefficientTable"]])) return()

  coefficientTable <- createJaspTable(gettext("Standard error of measurement"),
                                      dependencies = c("thorndike", "feldt", "mollenkopf", "mollenkopfFeldt",
                                                       "anova", "irt"))
  jaspResults[["semMainContainer"]][["coefficientTable"]] <- coefficientTable

  coefficientTable$addColumnInfo(name = "score", title = gettext("Score"), type = "string")
  coefficientTable$addColumnInfo(name = "average", title = gettext("Average"), type = "number")

  if (!ready || jaspResults[["semMainContainer"]]$getError()) return()

  dtFill <- data.frame(score = "all")
  average <- jaspResults[["semMainContainer"]][["averageState"]]$object
  dtFill$average <- average

  if (any(c(options[["thorndike"]], options[["feldt"]], options[["mollenkopf"]], options[["mollenkopfFeldt"]],
          options[["anova"]], options[["irt"]]))) {

    # the repetition of this seems annoying...
    coefficientTable <- createJaspTable(gettext("Standard error of measurement"),
                                        dependencies = c("thorndike", "feldt", "mollenkopf", "mollenkopfFeldt",
                                                         "anova", "irt"))
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


.semGRM <- function(X) {

  # res <- ltm::grm(X)
  res <- mirt::mirt(X, itemtype = "graded")
  x <- mirt::fscores(res)
  # x <- seq(-5, 5, length.out = nrow(X))
  cofs <- mirt::coef(res, IRTpars=TRUE)
  cofs$GroupPars <- NULL
  cnames <- colnames(cofs[[1]])
  cof_mat <- t(sapply(cofs, function(x) x))
  nc <- dim(cof_mat)[2]  # only works for equal categories over items
  colnames(cof_mat) <- c("a", paste0("b", 1:(nc-1)))
  cof_mat <- cof_mat[, c(2:nc, 1)]

  a <- cof_mat[, "a"]
  nit <- length(a)
  b <- cbind(-Inf, cof_mat[, 1:(nc - 1)], Inf)
  ouIRT <- matrix(NA, length(x), 2)
  probs <- matrix(NA, nit, nc)
  for (i in 1:length(x)) {
    for (j in 1:(nc)) {
      probs[, j] <- plogis(a * (x[i] - b[, j])) - plogis(a * (x[i] - b[, j + 1]))
    }
    ev <- rowSums(t(t(probs) * (0:(nc - 1))))
    dev <- (matrix(0:(nc - 1), nit, nc, TRUE) - ev)^2

    ouIRT[i, 1] <- sum(ev)
    ouIRT[i, 2] <- sqrt(sum(dev * probs))
  }

  return(ouIRT)
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

  # create a matrix that counts the scores and also checks the caseMins
  counts <- numeric(length(scores))
  # first count all the different sum score occurences
  for (i in seq_len(nrow(out))) {
    counts[i] <- sum(S == scores[i])
  }
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


.semFeldt <- function(X, K, nc, caseMin = 3, splits = NULL) {

  nit <- ncol(X)
  N <- nrow(X)

  S <- rowSums(X)
  partSUMS <- matrix(NA, N, K)

  if (K < nit) {
    if (is.null(splits)) {
      # when no specific splits are defined, we just evenly distribute the items among the number of splits
      k <- split(seq_len(nit), 1:K)
    }
    for (i in 1:K) {
      partSUMS[, i] <- rowSums(X[, k[[i]]])
    }
    d <- nit / mean(sapply(k, length))
  }

  if (K == nit) {
    for (i in 1:K) {
      partSUMS[, i] <- X[, i]
    }
    d <- nit
  }

  out <- matrix(NA, (nit * (nc - 1)) + 1, 3)
  out[, 1] <- 0:(nit * (nc - 1))

  ss <- 0
  while (ss[length(ss)] < (nit * (nc - 1))) {

    ind <- which(S == ss)

    if (length(ind) >= caseMin) {
      mean_diff <- partSUMS[ind, ] - rowMeans(partSUMS[ind, ]) - matrix(colMeans(partSUMS[ind, ]), length(ind), K, TRUE) + mean(partSUMS[ind, ])
      out[ss + 1, 2] <- sqrt(d * sum(rowSums(mean_diff^2) / (K - 1)) / length(ind))
      out[ss + 1, 3] <- length(ind)
      ss <- ss + 1
    } else { # ind < caseMin, too few observations for the single score, build group
      out[ss + 1, 3] <- length(ind)
      while (length(ind) < caseMin) {
        ss <- c(ss, ss[length(ss)] + 1)
        ind <- which(S %in% ss)
        # so what to do if we are at the end of the score range, and we just never reach the caseMin anymore?
        # we just use the sem from the previous score (-group).
        if (ss[length(ss)] + 1 > nrow(out)) {
          ss <- ss[-length(ss)]
          out[ss + 1, 2] <- out[ss[1], 2]
          breakout <- TRUE
          break
        }
        breakout <- FALSE
        out[ss[length(ss)] + 1, 3] <- sum(S == ss[length(ss)])
      }
      if (!breakout) {
        mean_diff <- partSUMS[ind, ] - matrix(colMeans(partSUMS[ind, ]), length(ind), K, TRUE) - rowMeans(partSUMS[ind, ]) + mean(partSUMS[ind, ])
        out[ss + 1, 2] <- sqrt(d * sum(rowSums(mean_diff^2) / (K - 1)) / length(ind))
        ss <- ss[length(ss)] + 1
      }
    }
  }

  return(out)
}


.semMF <- function(X, nc, splits = NULL, n_poly = 3, K = 2) {

  N <- nrow(X)
  nit <- ncol(X)
  partSUMS <- matrix(NA, N, K)
  S <- rowSums(X)

  if (K < nit) {
    if (is.null(splits)) {
      # when no specific splits are defined, we just evenly distribute the items among the number of splits
      k <- split(seq_len(nit), 1:K)
    }
    for (i in 1:K) {
      partSUMS[, i] <- rowSums(X[, k[[i]]])
    }
    d <- nit / mean(sapply(k, length))
  }

  if (K == nit) {
    for (i in 1:K) {
      partSUMS[, i] <- X[, i]
    }
    d <- nit
  }

  x <- 0:(nit * (nc - 1))
  out <- matrix(NA, length(x), 2)
  out[, 1] <- x

  # scored <- c()
  # for (j in 0:(length(x)-1)) {
  #   ind <- which(S == j)
  #   scored[j + 1] <- ifelse(length(ind) < caseMin, FALSE, TRUE)
  # }

  rawDiffK <- d *
    rowSums((partSUMS - matrix(colMeans(partSUMS), N, K, TRUE) - rowMeans(partSUMS) + mean(partSUMS))^2) /
    (K - 1)
  betaK <- coef(lm(rawDiffK ~ poly(S, n_poly, raw = TRUE)))
  scrs <- sqrt(betaK[1] + rowSums(matrix(betaK[-1], length(x), n_poly, TRUE) * poly(x, n_poly, raw = TRUE)))
  out[, 2] <- scrs

  return(out)
}


.semAnova <- function(X, nc, caseMin = 5) {

  out <- c()
  S <- rowSums(X)
  nit <- ncol(X)

  out <- matrix(NA, (nit * (nc - 1)) + 1, 3)
  out[, 1] <- 0:(nit * (nc - 1))

  ss <- 0
  while (ss[length(ss)] < (nit * (nc - 1))) {

    ind <- which(S == ss)

    if (length(ind) >= caseMin) {
      out[ss + 1, 2] <- sqrt(nit / (nit - 1) * sum(diag(cov(X[ind, ]))))
      out[ss + 1, 3] <- length(ind)
      ss <- ss + 1
    } else { # ind < caseMin, too few observations for the single score, build group
      out[ss + 1, 3] <- length(ind)
      while (length(ind) < caseMin) {
        ss <- c(ss, ss[length(ss)] + 1)
        ind <- which(S %in% ss)
        # so what to do if we are at the end of the score range, and we just never reach the caseMin anymore?
        # we just use the sem from the previous score (group).
        if (ss[length(ss)] + 1 > nrow(out)) {
          ss <- ss[-length(ss)]
          out[ss + 1, 2] <- out[ss[1], 2]
          breakout <- TRUE
          break
        }
        breakout <- FALSE
        out[ss[length(ss)] + 1, 3] <- sum(S == ss[length(ss)])
      }
      if (!breakout) {
        out[ss + 1, 2] <- sqrt(nit / (nit - 1) * sum(diag(cov(X[ind, ]))))
        ss <- ss[length(ss)] + 1
      }
    }
  }

  return(out)
}


