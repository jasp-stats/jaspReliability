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

  .semCombinedPointPlot(jaspResults, dataset, options, ready)

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

  .hasErrors(dataset = dataset,
             type = c('missingValues', "variance", "infinity", "variance"),
             missingValues.target = options$variables,
             exitAnalysisIfErrors = TRUE)
}



# create Main container
.semCreateMainContainer <- function(jaspResults, options) {

  if (!is.null(jaspResults[["semMainContainer"]])) return()

  semMainContainer <- createJaspContainer(dependencies = "variables")
  jaspResults[["semMainContainer"]] <- semMainContainer

  return()
}


##### Computation #####
.semComputeCoefficients <- function(jaspResults, dataset, options, ready) {

  if (!ready) return()

  nc <- length(unique(c(as.matrix(dataset))))
  scrs <- .semCounts(dataset, nc)
  counts <- scrs$counts
  scores <- scrs$scores
  countsState <- createJaspState(scrs, dependencies = NULL)
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
      out <- .semThorn(dataset, K = 2, nc = nc, caseMin = options[["minimumGroupSize"]], splits = NULL, scrs)
      thorndikeState <- createJaspState(out, dependencies = c("thorndike", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["thorndikeState"]] <- thorndikeState
    }
  }

  if (options[["feldt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["feldtState"]])) {
      out <- .semFeldt(dataset, K = options[["feldtNumberOfSplits"]], nc = nc, caseMin = options[["minimumGroupSize"]],
                       splits = NULL, scrs)
      feldtState <- createJaspState(out, dependencies = c("feldt", "feldtNumberOfSplits", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["feldtState"]] <- feldtState
    }
  }

  if (options[["mollenkopfFeldt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["mfState"]])) {
      out <- .semMF(dataset, K = options[["mollenkopfFeldtNumberOfSplits"]], nc = nc,
                    n_poly = options[["mollenkopfFeldtPolyDegree"]], splits = NULL, scrs)
      mfState <- createJaspState(out, dependencies = c("mollenkopfFeldt", "mollenkopfFeldtNumberOfSplits",
                                                          "mollenkopfFeldtPolyDegree"))
      jaspResults[["semMainContainer"]][["mfState"]] <- mfState
    }
  }

  if (options[["anova"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["anovaState"]])) {
      out <- .semAnova(dataset, nc = nc, caseMin = options[["minimumGroupSize"]], scrs)
      anovaState <- createJaspState(out, dependencies = c("anova", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["anovaState"]] <- anovaState
    }
  }

  # works for both data types
  if (options[["irt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["irtState"]])) {
      out <- .semIRT(dataset, nc = nc, scores)
      irtState <- createJaspState(out, dependencies = "irt")
      jaspResults[["semMainContainer"]][["irtState"]] <- irtState
    }
  }

  if (any(options[["lord"]], options[["keats"]], options[["lord2"]]) && nc > 2) {
    .quitAnalysis(gettext("The Lord, Keats, and Lord's compound methods are only available for binary data."))
  }

  # only for binary data
  if (options[["lord"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["lordState"]])) {
      out <- .semLord(ncol(dataset), scrs)
      lordState <- createJaspState(out, dependencies = "lord")
      jaspResults[["semMainContainer"]][["lordState"]] <- lordState
    }
  }

  if (options[["keats"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["keatsState"]])) {
      out <- .semKeats(dataset, options, scrs)
      keatsState <- createJaspState(out, dependencies = "keats")
      jaspResults[["semMainContainer"]][["keatsState"]] <- keatsState
    }
  }

  if (options[["lord2"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["lord2State"]])) {
      out <- .semLord2(dataset, options[["lord2NumberOfSplits"]], scrs, options[["minimumGroupSize"]])
      lord2State <- createJaspState(out, dependencies = c("lord2", "lord2NumberOfSplits", "minimumGroupSize"))
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
                                                       "anova", "irt", "lord", "keats", "lord2", "hideTable",
                                                       "feldtNumberOfSplits", "mollenkopfFeldtNumberOfSplits",
                                                       "mollenkopfFeldtPolyDegree", "minimumGroupSize",
                                                       "lord2NumberOfSplits"))
  jaspResults[["semMainContainer"]][["coefficientTable"]] <- coefficientTable

  coefficientTable$addColumnInfo(name = "score", title = gettext("Sum score"), type = "string")
  coefficientTable$addColumnInfo(name = "average", title = gettext("Traditional"), type = "number")

  if (!ready || jaspResults[["semMainContainer"]]$getError()) return()

  dtFill <- data.frame(score = "all")
  average <- jaspResults[["semMainContainer"]][["averageState"]]$object
  dtFill$average <- average

  if (!options[["hideTable"]]) {
    if (any(c(options[["thorndike"]], options[["feldt"]], options[["mollenkopfFeldt"]], options[["anova"]],
              options[["irt"]], options[["lord"]], options[["keats"]], options[["lord2"]]))) {

      # the repetition of this seems annoying...
      coefficientTable <- createJaspTable(gettext("Standard error of measurement"))
      coefficientTable$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]])
      jaspResults[["semMainContainer"]][["coefficientTable"]] <- coefficientTable

      coefficientTable$addColumnInfo(name = "score", title = gettext("Sum score"), type = "string")
      coefficientTable$addColumnInfo(name = "counts", title = gettext("Counts"), type = "string")
      coefficientTable$addFootnote(message = gettextf("The traditional sem value equals %1$1.3f.", average))

      scrs <- jaspResults[["semMainContainer"]][["countsState"]]$object
      counts <- scrs$counts
      dtFill <- data.frame(score = scrs$scores)
      dtFill$counts <- counts

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
        out <- jaspResults[["semMainContainer"]][["irtState"]]$object$binned
        coefficientTable$addColumnInfo(name = "irt", title = gettext("IRT"), type = "number")
        dtFill$irt <- out[, 2]
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
        coefficientTable$addColumnInfo(name = "lord2", title = gettext("Lord's compound"), type = "number")
        dtFill$lord2 <- out[, 2]
      }
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
    ggplot2::xlab(gettext("Sum scores"))
  p <- jaspGraphs::themeJasp(p)

  histPlot <- createJaspPlot(plot = p, title = gettext("Histogram of counts per sum score group"))
  jaspResults[["semMainContainer"]][["histPlot"]] <- histPlot

  return()

}

.semPointPlots <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["pointPlots"]]) || !options[["pointPlots"]]
      || !ready || jaspResults[["semMainContainer"]]$getError()) {return()}


  pointPlotsContainer <- createJaspContainer(title = gettext("Plots"))
  pointPlotsContainer$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]], options = "pointPlots")
  jaspResults[["semMainContainer"]][["pointPlotsContainer"]] <- pointPlotsContainer

  if (options[["thorndike"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["thorndikeState"]],
                                  title = "Thorndike")

    pointPlotsContainer[["thorndikePlot"]] <- pl
  }

  if (options[["feldt"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["feldtState"]],
                               title = "Feldt")

    pointPlotsContainer[["feldtPlot"]] <- pl
  }

  if (options[["mollenkopfFeldt"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["mfState"]],
                               title = "Mollenkopf-Feldt")

    pointPlotsContainer[["mfPlot"]] <- pl
  }

  if (options[["anova"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["anovaState"]],
                               title = gettext("ANOVA"))

    pointPlotsContainer[["anovaPlot"]] <- pl
  }

  if (options[["irt"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["irtState"]],
                               title = gettext("IRT"), irt = TRUE)

    pointPlotsContainer[["irtPlot"]] <- pl
  }

  if (options[["lord"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["lordState"]],
                                  title = "Lord")

    pointPlotsContainer[["lordPlot"]] <- pl
  }

  if (options[["keats"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["keatsState"]],
                                  title = "Keats")

    pointPlotsContainer[["keatsPlot"]] <- pl
  }

  if (options[["lord2"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["lord2State"]],
                                  title = gettext("Lord's compound"))

    pointPlotsContainer[["lordPlot"]] <- pl
  }



  return()

}

.semMakeSinglePointPlot <- function(resultsObject, title, irt = FALSE) {

  if (!irt) {
    dat <- as.data.frame(resultsObject$object)
  } else {
    dat <- as.data.frame(resultsObject$object$unbinned)
  }
  colnames(dat) <- c("score", "sem")
  pl <- ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(x = score, y = sem), size = 2.5) +
    ggplot2::labs(x = "Sum Score", y = "sem")

  outPlot <- createJaspPlot(jaspGraphs::themeJasp(pl), title = title, width = 400)
  outPlot$dependOn(optionsFromObject = resultsObject)

  return(outPlot)

}

.semCombinedPointPlot <- function(jaspResults, dataset, options, ready) {

  if (!options[["combinedPointPlot"]]
      || !is.null(jaspResults[["semMainContainer"]][["combinedPlot"]])
      || !ready
      || jaspResults[["semMainContainer"]]$getError()) {
    return()
  }



  if (any(c(options[["thorndike"]], options[["feldt"]], options[["mollenkopfFeldt"]], options[["anova"]],
            options[["irt"]], options[["lord"]], options[["keats"]], options[["lord2"]]))) {

    nit <- ncol(dataset)
    nc <- length(unique(c(as.matrix(dataset))))
    scores <- jaspResults[["semMainContainer"]][["countsState"]]$object$scores
    dt <- data.frame(score = integer(), sem = double(), Type = character())
    dtIRT <- NULL

    if (options[["thorndike"]]) {
      out <- jaspResults[["semMainContainer"]][["thorndikeState"]]$object
      dtBind <- data.frame(score = scores, sem = out[, 2], Type = "Thorndike")
      dt <- rbind(dt, dtBind)
    }

    if (options[["feldt"]]) {
      out <- jaspResults[["semMainContainer"]][["feldtState"]]$object
      dtBind <- data.frame(score = scores, sem = out[, 2], Type = "Feldt")
      dt <- rbind(dt, dtBind)
    }
    if (options[["mollenkopfFeldt"]]) {
      out <- jaspResults[["semMainContainer"]][["mfState"]]$object
      dtBind <- data.frame(score = scores, sem = out[, 2], Type = "MF")
      dt <- rbind(dt, dtBind)
    }
    if (options[["anova"]]) {
      out <- jaspResults[["semMainContainer"]][["anovaState"]]$object
      dtBind <- data.frame(score = scores, sem = out[, 2], Type = gettext("ANOVA"))
      dt <- rbind(dt, dtBind)
    }

    if (options[["irt"]]) {
      out <- jaspResults[["semMainContainer"]][["irtState"]]$object$unbinned
      dtIRT <- data.frame(score = out[, 1], sem = out[, 2], Type = gettext("IRT"))
    }

    if (options[["lord"]]) {
      out <- jaspResults[["semMainContainer"]][["lordState"]]$object
      dtBind <- data.frame(score = scores, sem = out[, 2], Type = gettext("Lord"))
      dt <- rbind(dt, dtBind)
    }

    if (options[["keats"]]) {
      out <- jaspResults[["semMainContainer"]][["keatsState"]]$object
      dtBind <- data.frame(score = scores, sem = out[, 2], Type = gettext("Keats"))
      dt <- rbind(dt, dtBind)
    }

    if (options[["lord2"]]) {
      out <- jaspResults[["semMainContainer"]][["lord2State"]]$object
      dtBind <- data.frame(score = scores, sem = out[, 2], Type = gettext("Lord's compound"))
      dt <- rbind(dt, dtBind)
    }


    # so that the legend does not become alphabetically ordered:
    levs <- unique(dt$Type)
    dt$Type <- factor(dt$Type, levels = levs, ordered = FALSE, labels = levs)

    pl <- ggplot2::ggplot() +
      ggplot2::geom_point(data = dt, ggplot2::aes(x = score, y = sem, shape = Type), size = 2.5) +
      ggplot2::labs(x = "Sum Score", y = "sem") +
      ggplot2::theme(plot.margin = ggplot2::margin(t=1,r=5,b=1.5,l=1, "cm"))

    if (!is.null(dtIRT)) {
      pl <- pl + ggplot2::geom_line(data = dtIRT, ggplot2::aes(x = score, y = sem, linetype = ""), color = "black", linewidth = 1) +
        ggplot2::scale_linetype_discrete(name = "", labels = dtIRT[1, 3])
      # move the legend of IRT underneath the other:
      pl <- pl + ggplot2::guides(linetype = ggplot2::guide_legend(order = 0),
             shape = ggplot2::guide_legend(order = 1))
    }

    plot <- createJaspPlot(jaspGraphs::themeJasp(pl, legend.position = "right"), title = gettext("Combined plot"),
                           width = 600)
    plot$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]])
    jaspResults[["semMainContainer"]][["combinedPlot"]] <- plot
  }

  return()
}



##### Sem compute functions #####


.semIRT <- function(X, nc, scores) {

  if (nc == 2) {
    ity <- "2PL"
  } else {
    ity <- "graded"
  }

  res <- mirt::mirt(X, model = 1, itemtype = ity)
  # x <- mirt::fscores(res)
  x <- seq(-5, 5, by = 0.1)

  cofs <- mirt::coef(res, IRTpars=TRUE)
  cofs$GroupPars <- NULL
  cof_mat <- t(sapply(cofs, function(x) x))
  # make data frame just so that all columns are named and we dont have to be explicit
  cof_mat <- as.data.frame(cof_mat)
  colnames(cof_mat) <- c("a", paste0("b", 1:(nc-1)))

  a <- cof_mat[, "a"]
  nit <- length(a)
  outIRT <- matrix(NA, length(x), 2)

  if (nc == 2) {
    b <- cof_mat[, "b1"]
    for (i in 1:length(x)) {
      outIRT[i, 1] <- sum(plogis(a * (x[i] - b)))
      outIRT[i, 2] <- sqrt(sum(plogis(a * (x[i] - b)) * (1 - plogis(a * (x[i] - b)))))
    }

  } else {
    b <- cbind(-Inf, cof_mat[, paste0("b", 1:(nc-1))], Inf)
    probs <- matrix(NA, nit, nc)
    for (i in 1:length(x)) {
      for (j in 1:nc) {
        probs[, j] <- plogis(a * (x[i] - b[, j])) - plogis(a * (x[i] - b[, j + 1]))
      }
      ev <- rowSums(t(t(probs) * (0:(nc - 1))))
      dev <- (matrix(0:(nc - 1), nit, nc, TRUE) - ev)^2

      outIRT[i, 1] <- sum(ev)
      outIRT[i, 2] <- sqrt(sum(dev * probs))
    }
  }

  discScores <- scores

  # because the irt method does not produce score groups but a continuous score along the scale, we need to group it
  outDt <- as.data.frame(outIRT)
  outOrdered <- outDt[order(outDt[, 1], decreasing = FALSE), ]
  uniqueOutOrdered <- unique(outOrdered)
  irtBinned <- approx(uniqueOutOrdered[, 1], uniqueOutOrdered[, 2], xout = discScores)$y
  outBinned <- data.frame(scores = discScores, sem = irtBinned)

  outIRT <- list(unbinned = outIRT, binned = outBinned)
  return(outIRT)

}


.semThorn <- function(X, K, nc, caseMin, splits = NULL, scrs) {

  prep <- .semPrepareSumScores(X, K, splits = NULL)
  partSUMS <- prep$partSUMS
  S <- prep$S
  out <- .semPrepareOutMatrix(ncol(X), nc, scrs)
  fun <- function(partSUMS, ind, cc) {
    return(sd(partSUMS[ind, 1] - partSUMS[ind, 2]))
  }

  out <- .semComputeWithCaseMin(out, S, caseMin, partSUMS, fun)

  return(out)
}


.semFeldt <- function(X, K, nc, caseMin, splits = NULL, scrs) {

  prep <- .semPrepareSumScores(X, K, splits = NULL)
  partSUMS <- prep$partSUMS
  S <- prep$S
  d <- prep$d
  out <- .semPrepareOutMatrix(ncol(X), nc, scrs)
  fun <- function(partSUMS, ind, cc) {
    K <- ncol(partSUMS)
    mean_diff <- partSUMS[ind, ] - rowMeans(partSUMS[ind, ]) - matrix(colMeans(partSUMS[ind, ]), length(ind), K, TRUE) + mean(partSUMS[ind, ])
    ret <- sqrt(d * sum(rowSums(mean_diff^2) / (K - 1)) / length(ind))
    return(ret)
  }
  out <- .semComputeWithCaseMin(out, S, caseMin, partSUMS, fun)

  return(out)
}


.semMF <- function(X, K, nc, n_poly, splits = NULL, scrs) {

  prep <- .semPrepareSumScores(X, K, splits = NULL)
  partSUMS <- prep$partSUMS
  S <- prep$S
  d <- prep$d
  N <- nrow(X)
  out <- .semPrepareOutMatrix(ncol(X), nc, scrs)
  scores <- out[, 1]

  rawDiffK <- d *
    rowSums((partSUMS - matrix(colMeans(partSUMS), N, K, TRUE) - rowMeans(partSUMS) + mean(partSUMS))^2) /
    (K - 1)
  betaK <- coef(lm(rawDiffK ~ poly(S, n_poly, raw = TRUE)))
  scrs <- sqrt(betaK[1] + rowSums(matrix(betaK[-1], length(scores), n_poly, TRUE) * poly(scores, n_poly, raw = TRUE)))
  out[, 2] <- scrs

  return(out)
}


.semAnova <- function(X, nc, caseMin, scrs) {

  prep <- .semPrepareSumScores(X, K = 1, splits = NULL)
  S <- prep$S

  out <- .semPrepareOutMatrix(ncol(X), nc, scrs)
  fun <- function(X, ind, cc) {
    nit <- ncol(X)
    return(sqrt(nit / (nit - 1) * sum(diag(cov(X[ind, ])))))
  }
  out <- .semComputeWithCaseMin(out, S, caseMin, X, fun)

  return(out)
}


.semLord <- function(nit, scrs) {
  x <- scrs$scores
  out_binom <- sqrt((x) * (nit - (x)) / (nit - 1))
  return(cbind(x, out_binom))
}

.semKeats <- function(X, options, scrs) {

  nit <- ncol(X)
  x <- scrs$scores
  S <- rowSums(X)
  KR21 <- nit / (nit - 1) * (1 - (mean(S) * (nit - mean(S))) / (nit * var(S)))
  if (options[["userReliability"]]) {
    rel <- options[["reliabilityValue"]]
  } else {
    rel <- Bayesrel:::applyalpha(cov(X))
  }
  correction <- (1 - rel) / (1 - KR21)
  ou_binom2 <- sqrt((x) * (nit - (x)) / (nit - 1) * correction)

  return(cbind(x, ou_binom2))
}


# X = dataset, K = number of splits, counts = counts per score group
.semLord2 <- function(X, K, scrs, caseMin, splits = NULL) {

  nc <- 2
  prep <- .semPrepareSumScores(X, K, splits = NULL)
  partSUMS <- prep$partSUMS
  S <- prep$S
  cc <- prep$cc

  out <- .semPrepareOutMatrix(ncol(X), nc, scrs)

  fun <- function(partSUMS, ind, cc) {
    ccmat <- matrix(cc, length(ind), length(cc), byrow = TRUE)
    ret <- sqrt(mean(rowSums(partSUMS[ind, ] * (ccmat - partSUMS[ind, ]) / (ccmat - 1))))
    return(ret)
  }
  out <- .semComputeWithCaseMin(out, S, caseMin, partSUMS, fun, cc)

  return(out)

}



#### Helper functions ####
.semCounts <- function(X, nc) {
  nit <- ncol(X)
  S <- rowSums(X)
  scoreRange <- range(S)
  scores <- scoreRange[1]:scoreRange[2]
  # create a matrix that counts the scores and also checks the caseMins
  counts <- numeric(length(scores))
  # first count all the different sum score occurences
  for (i in seq_len(length(counts))) {
    counts[i] <- sum(S == scores[i])
  }
  return(list(counts = counts, scores = scores))
}

.semPrepareSumScores <- function(X, K, splits) {

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
    cc <- sapply(k, length)
    d <- nit / mean(cc)

  }

  if (K == nit) {
    for (i in 1:K) {
      partSUMS[, i] <- X[, i]
    }
    cc <- rep(1, K)
    d <- nit
  }

  return(list(S = S, partSUMS = partSUMS, d = d, cc = cc))
}


.semPrepareOutMatrix <- function(nit, nc, scrs) {
  scores <- scrs$scores
  counts <- scrs$counts
  out <- matrix(NA, length(scores), 4)
  out[, 1] <- scores
  out[, 3] <- counts
  out[, 4] <- FALSE

  return(out)
}


.semComputeWithCaseMin <- function(out, S, caseMin, partSUMS, fun, cc = NULL) {

  scores <- out[, 1]
  counts <- out[, 3]
  ii <- 1
  while (ii <= nrow(out)) {
    if (counts[ii] >= caseMin) {
      ind <- which(S == scores[ii])
      out[ii, 2] <- fun(partSUMS, ind, cc)
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
        out[ii:nextInd, 2] <- fun(partSUMS, ind)
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
    out[backInd:nextInd, 2] <- fun(partSUMS, ind)
    out[backInd:nextInd, 4] <- TRUE
  }

  return(out)
}
