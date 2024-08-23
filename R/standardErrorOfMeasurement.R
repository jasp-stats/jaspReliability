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



  ready <- length(options[["variables"]]) > 1

  dataset <- .semReadData(dataset, options)

  .semErrorCheck(dataset, options, ready)

  .semCreateMainContainer(jaspResults, options)

  options <- .semOptionsHelper(options, dataset)

  .semComputeCoefficients(jaspResults, dataset, options, ready)

  .semComputeSumScoresCi(jaspResults, dataset, options, ready)

  .semCoefficientsTable(jaspResults, dataset, options, ready)

  .semSumScoreCiTable(jaspResults, dataset, options, ready)

  .semHistPlot(jaspResults, dataset, options, ready)

  .semPointPlots(jaspResults, dataset, options, ready)

  .semCombinedPointPlot(jaspResults, dataset, options, ready)

  .semSumScoreCiPlots(jaspResults, dataset, options, ready)

  return()
}


#### the common functions ####
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

  if (options[["lord2"]] && options[["lord2NumberOfSplits"]] == "") {
    .quitAnalysis(gettext("For the Lord's compound method, the test could not be split in equally sized parts with more than 1 item per part. Consider adding or removing items."))
  }
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

  scrs <- .semCounts(dataset, nc)
  counts <- scrs$counts
  scores <- scrs$scores
  countsState <- createJaspState(scrs, dependencies = NULL)
  nc <- length(unique(c(as.matrix(dataset)))) # may be needed for IRT

  jaspResults[["semMainContainer"]][["countsState"]] <- countsState

  # first the average sem
  if (is.null(jaspResults[["semMainContainer"]][["averageState"]])) {
    if (options[["userReliability"]]) {
      rel <- options[["reliabilityValue"]]
    } else {
      rel <- Bayesrel:::applyalpha(cov(dataset))
    }

    average <- list(est = as.numeric(sd(rowSums(dataset)) * sqrt(1 - rel)))

    averageState <- createJaspState(average, dependencies = c("userReliability", "reliabilityValue"))
    jaspResults[["semMainContainer"]][["averageState"]] <- averageState
  }

  selected <- options[["selected"]]

  if (length(selected) > 0) {
    method <- names(selected)
    # at least one method is selected
    for (i in 1:length(selected)) {
      if (is.na(selected[[i]][["name"]])) {
        .quitAnalysis(gettext("The Lord, Keats, and Lord's compound method are only available for binary data."))
      }
      if (is.null(jaspResults[["semMainContainer"]][[paste0(method[i], "State")]])) {
        out <- eval(parse(text = selected[[i]][["funString"]]))
        state <- createJaspState(out, dependencies = selected[[i]][["dependencies"]])
        jaspResults[["semMainContainer"]][[paste0(method[i], "State")]] <- state
      }
    }
  }

 return()
}

.semComputeSumScoresCi <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["sumScoreState"]]) ||
      !ready) return()

  if (!options[["sumScoreCiTable"]] && !options[["sumScoreCiPlots"]]) return()

  scrs <- jaspResults[["semMainContainer"]][["countsState"]]$object

  dtFill <- list(table = data.frame(score = scrs$scores), plots = data.frame(score = scrs$scores))

  selected <- options[["selected"]]
  if (length(selected) > 0) {
    method <- names(selected)
    # at least one method is selected
    for (i in 1:length(selected)) {
      out <- jaspResults[["semMainContainer"]][[paste0(method[i], "State")]]$object

      if (!is.null(names(out))) { # then we have IRT and actually need to decide which of the scores to take
        out <- out$binned
      }

      cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
      dtFill[["table"]][[paste0("lower", method[i])]] <- cis[, 1]
      dtFill[["table"]][[paste0("upper", method[i])]] <- cis[, 2]
      dtFill[["plots"]][[paste0("lower", method[i])]] <- cis[, 3]
      dtFill[["plots"]][[paste0("upper", method[i])]] <- cis[, 4]
    }
  }

  ciDataState <- createJaspState(dtFill,
                                 dependencies = c("thorndike", "feldt", "mollenkopfFeldt",
                                                  "anova", "irt", "lord", "keats", "lord2", "hideTable",
                                                  "feldtNumberOfSplits", "mollenkopfFeldtNumberOfSplits",
                                                  "mollenkopfFeldtPolyDegree", "minimumGroupSize",
                                                  "lord2NumberOfSplits", "userReliability", "reliabilityValue",
                                                  "sumScoreCiTable", "sumScoreCiPlots",
                                                  "ciLevelTable", "ciLevelPlots"))
  jaspResults[["semMainContainer"]][["ciDataState"]] <- ciDataState

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
                                                       "lord2NumberOfSplits", "userReliability", "reliabilityValue"))
  coefficientTable$position <- 1
  jaspResults[["semMainContainer"]][["coefficientTable"]] <- coefficientTable

  coefficientTable$addColumnInfo(name = "score", title = gettext("Sum score"), type = "string")
  coefficientTable$addColumnInfo(name = "average", title = gettext("Traditional"), type = "number")

  if (!ready || jaspResults[["semMainContainer"]]$getError()) return()

  dtFill <- data.frame(score = "all")
  average <- jaspResults[["semMainContainer"]][["averageState"]]$object
  dtFill$average <- average$est

  if (!options[["hideTable"]]) {

    selected <- options[["selected"]]
    if (length(selected) > 0) {  # at least one method is selected

      coefficientTable <- createJaspTable(gettext("Standard error of measurement"))
      coefficientTable$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]])
      coefficientTable$position <- 1
      jaspResults[["semMainContainer"]][["coefficientTable"]] <- coefficientTable

      coefficientTable$addColumnInfo(name = "score", title = gettext("Sum score"), type = "string")
      coefficientTable$addColumnInfo(name = "counts", title = gettext("Counts"), type = "string")
      coefficientTable$addFootnote(message = gettextf("The traditional sem value equals %1$1.3f.", average$est))

      scrs <- jaspResults[["semMainContainer"]][["countsState"]]$object
      counts <- scrs$counts
      dtFill <- data.frame(score = scrs$scores)
      dtFill$counts <- counts
      ci <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)

      method <- names(selected)

      for (i in 1:length(selected)) {
        out <- jaspResults[["semMainContainer"]][[paste0(method[i], "State")]]$object

        if (!is.null(names(out))) { # then we have IRT and actually need to decide which of the scores to take
          out <- out$binned
        }
        coefficientTable$addColumnInfo(name = method[i], title = selected[[i]][["name"]], type = "number")
        dtFill[[method[i]]] <- out[, 2]
      }
    }
  }

  coefficientTable$setData(dtFill)

  return()

}

.semSumScoreCiTable <- function(jaspResults, dataset, options, ready) {
  if (!is.null(jaspResults[["semMainContainer"]][["ciTable"]]) ||
      !ready) return()

  if (!options[["sumScoreCiTable"]] || is.null(jaspResults[["semMainContainer"]][["ciDataState"]])) return()

  ciTable <- createJaspTable(gettext("Sum score CI table"))
  ciTable$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]],
                   options = c("ciLevelTable", "sumScoreCiTable"))
  ciTable$position <- 2
  jaspResults[["semMainContainer"]][["ciTable"]] <- ciTable

  ciTable$addColumnInfo(name = "score", title = gettext("Sum score"), type = "string")
  ci <- format(100 * options[["ciLevelTable"]], digits = 3, drop0trailing = TRUE)

  ciData <- jaspResults[["semMainContainer"]][["ciDataState"]]$object$table

  selected <- options[["selected"]]

  if (length(selected) > 0) { # at least one method is selected
    method <- names(selected)
    for (i in 1:length(selected)) {
      ciTable$addColumnInfo(name = paste0("lower", method[i]), title = gettextf("Lower %s%% CI", ci), type = "number",
                            overtitle = selected[[i]][["name"]])
      ciTable$addColumnInfo(name = paste0("upper", method[i]), title = gettextf("Upper %s%% CI", ci), type = "number",
                            overtitle = selected[[i]][["name"]])
    }
  }

  ciTable$setData(ciData)

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
  p <- p + jaspGraphs::themeJaspRaw() + jaspGraphs::geom_rangeframe()

  histPlot <- createJaspPlot(plot = p, title = gettext("Histogram of counts per sum score group"), width = 500)
  jaspResults[["semMainContainer"]][["histPlot"]] <- histPlot

  return()

}

.semPointPlots <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["pointPlots"]]) || !options[["pointPlots"]]
      || !ready || jaspResults[["semMainContainer"]]$getError()) {return()}

  pointPlotsContainer <- createJaspContainer(title = gettext("Plots"))
  pointPlotsContainer$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]], options = "pointPlots")
  jaspResults[["semMainContainer"]][["pointPlotsContainer"]] <- pointPlotsContainer

  selected <- options[["selected"]]
  if (length(selected) > 0) {
    method <- names(selected)
    # at least one method is selected
    for (i in 1:length(selected)) {
      out <- jaspResults[["semMainContainer"]][[paste0(method[i], "State")]]
      irt <- FALSE
      if (!is.null(names(out$object))) { # then we have IRT and actually need to decide which of the scores to take
        irt <- TRUE
      }

      pl <- .semMakeSinglePointPlot(resultsObject = out,
                                    title = selected[[i]][["name"]], irt = irt)
      pointPlotsContainer[[paste0(method[i], "Plot")]] <- pl
    }
  }

  return()

}

.semMakeSinglePointPlot <- function(resultsObject, title, irt = FALSE) {

  if (!irt) {
    dat <- as.data.frame(resultsObject$object)
    colnames(dat)[1:2] <- c("score", "sem")
    pl <- ggplot2::ggplot(dat) +
      ggplot2::geom_point(ggplot2::aes(x = score, y = sem), size = 2.5) +
      ggplot2::labs(x = gettext("Sum Score"), y = gettext("sem"))
  } else {
    dat <- as.data.frame(resultsObject$object$unbinned)
    colnames(dat)[1:2] <- c("score", "sem")
    pl <- ggplot2::ggplot(dat) +
      ggplot2::geom_line(ggplot2::aes(x = score, y = sem), size = 2.5) +
      ggplot2::labs(x = gettext("Sum Score"), y = gettext("sem"))
  }

  pl <- pl + jaspGraphs::themeJaspRaw() + jaspGraphs::geom_rangeframe()
  outPlot <- createJaspPlot(pl, title = title, width = 400)
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

  selected <- options[["selected"]]
  if (length(selected) > 0) { # at least one method is selected

    nit <- ncol(dataset)
    nc <- length(unique(c(as.matrix(dataset))))
    scores <- jaspResults[["semMainContainer"]][["countsState"]]$object$scores
    dt <- data.frame(score = integer(), sem = double(), Type = character())
    dtIRT <- NULL

    method <- names(selected)
    for (i in 1:length(selected)) {
      out <- jaspResults[["semMainContainer"]][[paste0(method[i], "State")]]$object
      if (!is.null(names(out))) { # then we have IRT and actually need to decide which of the scores to take
        out <- out$unbinned
        dtIRT <- data.frame(score = out[, 1], sem = out[, 2], Type = selected[[i]][["name"]])
      } else {
        dtBind <- data.frame(score = scores, sem = out[, 2], Type = selected[[i]][["name"]])
        dt <- rbind(dt, dtBind)
      }
    }

    # so that the legend does not become alphabetically ordered:
    levs <- unique(dt$Type)
    dt$Type <- factor(dt$Type, levels = levs, ordered = FALSE, labels = levs)

    pl <- ggplot2::ggplot() +
      ggplot2::geom_point(data = dt, ggplot2::aes(x = score, y = sem, shape = Type), size = 2.5) +
      ggplot2::labs(x = "Sum Score", y = "sem") +
      ggplot2::theme(plot.margin = ggplot2::margin(t=1,r=5,b=1.5,l=1, "cm")) +
      jaspGraphs::themeJaspRaw(legend.position = "right") + jaspGraphs::geom_rangeframe()

    if (!is.null(dtIRT)) {
      pl <- pl + ggplot2::geom_line(data = dtIRT, ggplot2::aes(x = score, y = sem, linetype = ""), color = "black", linewidth = 1) +
        ggplot2::scale_linetype_discrete(name = "", labels = dtIRT[1, 3])
      # move the legend of IRT underneath the other:
      pl <- pl + ggplot2::guides(linetype = ggplot2::guide_legend(order = 0),
                                 shape = ggplot2::guide_legend(order = 1))
    }

    plot <- createJaspPlot(pl, title = gettext("Combined plot"),
                           width = 600)
    plot$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]],
                  options = "combinedPointPlot")
    jaspResults[["semMainContainer"]][["combinedPlot"]] <- plot
  }

  return()
}

.semSumScoreCiPlots <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["ciPlots"]]) || !options[["sumScoreCiPlots"]]
      || !ready || jaspResults[["semMainContainer"]]$getError()) {return()}

  nc <- length(unique(c(as.matrix(dataset)))) # may be needed for IRT

  ciPlotsContainer <- createJaspContainer(title = gettext("Sum Score CI Plots"))
  ciPlotsContainer$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]],
                            options = c("sumScoreCiPlots", "sumScoreCiPlotsCutoff", "sumScoreCiPlotsCutoffValue", "ciLevelPlots"))
  jaspResults[["semMainContainer"]][["ciPlotsContainer"]] <- ciPlotsContainer

  ciData <- jaspResults[["semMainContainer"]][["ciDataState"]]$object$plots
  scores <- jaspResults[["semMainContainer"]][["countsState"]]$object$scores

  selected <- options[["selected"]]
  if (length(selected) > 0) { # at least one method is selected

    nit <- ncol(dataset)
    nc <- length(unique(c(as.matrix(dataset))))
    scores <- jaspResults[["semMainContainer"]][["countsState"]]$object$scores
    dt <- data.frame(score = integer(), sem = double(), Type = character())
    dtIRT <- NULL

    method <- names(selected)
    for (i in 1:length(selected)) {

      if (method[i] == "irt") { # then we have IRT and actually need to decide which of the scores to take
        outUnbinned <- jaspResults[["semMainContainer"]][["irtState"]]$object$unbinned
        cis <- .semComputeCis(outUnbinned, outUnbinned[, 1], ciLevelTable = options[["ciLevelTable"]],
                              ciLevelPlots = options[["ciLevelPlots"]])
        lowerIrtUnbinned <- cis[, 3]
        upperIrtUnbinned <- cis[, 4]
        ciDataUse <- cbind(outUnbinned[, 1], lowerIrtUnbinned, upperIrtUnbinned)
        irt <- TRUE
      } else {
        ciDataUse <- ciData[, c("score", paste0("lower", method[i]), paste0("upper", method[i]))]
        irt <- FALSE
      }
      pl <- .semMakeSingleCiPlot(title = selected[[i]][["name"]], ciData = ciDataUse, irt = irt,
                                 cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
      ciPlotsContainer[[paste0(method[i], "Plot")]] <- pl

    }
  }

  return()

}

.semMakeSingleCiPlot <- function(ciData, title, irt = FALSE, cutoff = NA) {

  if (!irt) {
    dat <- as.data.frame(ciData)
    colnames(dat)[2:3] <- c("lower", "upper")
    dat$tscore <- dat$score
    pl <- ggplot2::ggplot(dat) +
      ggplot2::geom_point(ggplot2::aes(x = score, y = tscore), size = 2.5) +
      ggplot2::geom_errorbar(ggplot2::aes(x = score, ymin = lower, ymax = upper), width = 0.5) +
      ggplot2::labs(x = gettext("Sum Score"), y = gettext("True Score"))
  } else {
    dat <- as.data.frame(ciData)
    colnames(dat) <- c("score", "lower", "upper")
    dat$tscore <- dat$score
    pl <- ggplot2::ggplot(dat) +
      ggplot2::geom_ribbon(ggplot2::aes(x = score, ymin = lower, ymax = upper), fill = "grey80") +
      ggplot2::geom_line(ggplot2::aes(x = score, y = tscore)) +
      ggplot2::labs(x = gettext("Sum Score"), y = gettext("True Score"))
  }

  if (!is.na(cutoff)) {
    pl <- pl + ggplot2::geom_hline(yintercept = cutoff, linetype = "dashed", color = "red")
  }
  pl <- pl + jaspGraphs::themeJaspRaw() + jaspGraphs::geom_rangeframe()
  outPlot <- createJaspPlot(pl, title = title, width = 400, height = 400)

  return(outPlot)

}


##### the functions to compute the sems #####

#' if not stated otherwise:
#' X is the data matrix
#' K is the number of splits
#' nc is the number of response categories
#' caseMin is the minimum number of cases in a group
#' scrs is a vector of sum scores

.semThorn <- function(X, nc, caseMin, scoresObj) {

  prep <- .semPrepareSumScores(X, K = 2)
  partSUMS <- prep$partSUMS
  S <- prep$S
  out <- .semPrepareOutMatrix(ncol(X), nc, scoresObj)
  fun <- function(partSUMS, ind, cc) {
    return(sd(partSUMS[ind, 1] - partSUMS[ind, 2]))
  }

  out <- .semComputeWithCaseMin(out, S, caseMin, partSUMS, fun)
  return(out)
}

.semFeldt <- function(X, K, nc, caseMin, scoresObj) {

  prep <- .semPrepareSumScores(X, K)
  partSUMS <- prep$partSUMS
  S <- prep$S
  d <- prep$d
  out <- .semPrepareOutMatrix(ncol(X), nc, scoresObj)
  fun <- function(partSUMS, ind, cc) {
    K <- ncol(partSUMS)
    mean_diff <- partSUMS[ind, ] - rowMeans(partSUMS[ind, ]) - matrix(colMeans(partSUMS[ind, ]), length(ind), K, TRUE) + mean(partSUMS[ind, ])
    ret <- sqrt(d * sum(rowSums(mean_diff^2) / (K - 1)) / length(ind))
    return(ret)
  }
  out <- .semComputeWithCaseMin(out, S, caseMin, partSUMS, fun)
  return(out)
}

.semMollenkopfFeldt <- function(X, K, nc, n_poly, scoresObj) {

  prep <- .semPrepareSumScores(X, K)
  partSUMS <- prep$partSUMS
  S <- prep$S
  d <- prep$d
  N <- nrow(X)
  out <- .semPrepareOutMatrix(ncol(X), nc, scoresObj)
  scores <- out[, 1]

  rawDiffK <- d *
    rowSums((partSUMS - matrix(colMeans(partSUMS), N, K, TRUE) - rowMeans(partSUMS) + mean(partSUMS))^2) /
    (K - 1)
  betaK <- coef(lm(rawDiffK ~ poly(S, n_poly, raw = TRUE)))
  scrs <- sqrt(betaK[1] + rowSums(matrix(betaK[-1], length(scores), n_poly, TRUE) * poly(scores, n_poly, raw = TRUE)))
  out[, 2] <- scrs

  return(out)
}

.semAnova <- function(X, nc, caseMin, scoresObj) {

  prep <- .semPrepareSumScores(X, K = 1)
  S <- prep$S

  out <- .semPrepareOutMatrix(ncol(X), nc, scoresObj)
  fun <- function(X, ind, cc) {
    nit <- ncol(X)
    return(sqrt(nit / (nit - 1) * sum(diag(cov(X[ind, ])))))
  }
  out <- .semComputeWithCaseMin(out, S, caseMin, X, fun)

  return(out)
}

.semIrt <- function(X, nc, scores) {

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

  # we need to go with lists here because with ordinal items, we may encounter that
  # some items have a reduced number of answering categories because not all were chosen
  as <- sapply(cofs, function(x) x[, "a"])
  bs <- lapply(cofs, function(x) x[, grep("b", colnames(x))])
  bMulti <- lapply(bs, function(x) 0:length(x))

  nit <- length(as)
  outIRT <- matrix(NA, length(x), 2)

  if (nc == 2) {
    # we dont assume that for binary data items will end up with different number of answering categories
    # that would mean an item with no variance
    bs <- unlist(bs)
    for (i in 1:length(x)) {
      outIRT[i, 1] <- sum(plogis(as * (x[i] - bs)))
      outIRT[i, 2] <- sqrt(sum(plogis(as * (x[i] - bs)) * (1 - plogis(as * (x[i] - bs)))))
    }

  } else {

    bsInf <- lapply(bs, function(x) c(-Inf, x, Inf))
    probs <- vector("list", length = nit)

    for (i in 1:length(x)) {
      for (k in 1:nit) {
        blen <- length(bsInf[[k]]) - 1
        probs[[k]] <- numeric(blen)
        for (j in 1:blen) {
          probs[[k]][j] <- plogis(as[k] * (x[i] - bsInf[[k]][j])) - plogis(as[k] * (x[i] - bsInf[[k]][j + 1]))
        }
      }
      pprobs <- Map("*", probs, bMulti)
      ev <- sapply(pprobs, sum)
      dev <- Map("-", bMulti, ev)
      dev <- lapply(dev, function(x) x^2)
      devPr <- Map("*", dev, probs)

      outIRT[i, 1] <- sum(ev)
      outIRT[i, 2] <- sqrt(sum(unlist(devPr)))
    }
  }
  colnames(outIRT) <- c("scores", "sem")
  discScores <- scores

  # because the irt method does not produce score groups but a continuous score along the scale, we need to group it
  outDt <- as.data.frame(outIRT)
  outOrdered <- outDt[order(outDt[, 1], decreasing = FALSE), ]
  uniqueOutOrdered <- unique(outOrdered)
  irtBinned <- approx(uniqueOutOrdered[, 1], uniqueOutOrdered[, 2], xout = discScores)$y
  # take the value from the continous sems for the smallest latent trait value
  irtBinned[is.na(irtBinned)] <- uniqueOutOrdered[1, 2]
  outBinned <- data.frame(scores = discScores, sem = irtBinned)

  outIRT <- list(unbinned = outIRT, binned = outBinned)

  return(outIRT)

}

.semLord <- function(nit, scoresObj) {
  x <- scoresObj$scores
  out_binom <- sqrt((x) * (nit - (x)) / (nit - 1))
  return(cbind(x, out_binom))
}

.semKeats <- function(X, options, scoresObj) {

  nit <- ncol(X)
  x <- scoresObj$scores
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
.semLord2 <- function(X, K, scoresObj, caseMin) {

  nc <- 2
  prep <- .semPrepareSumScores(X, K)
  partSUMS <- prep$partSUMS
  S <- prep$S
  cc <- prep$cc

  out <- .semPrepareOutMatrix(ncol(X), nc, scoresObj)

  fun <- function(partSUMS, ind, cc) {
    ccmat <- matrix(cc, length(ind), length(cc), byrow = TRUE)
    ret <- sqrt(mean(rowSums(partSUMS[ind, ] * (ccmat - partSUMS[ind, ]) / (ccmat - 1))))
    return(ret)
  }
  out <- .semComputeWithCaseMin(out, S, caseMin, partSUMS, fun, cc)

  return(out)

}



#### Helper functions ####
#' X is the data matrix
#' nc is the number of response categories
#' K is the number of splits

#' this function counts the sum scores and returns the scores and their counts
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

#' this function prepares the sum scores for the SEM computations:
#' for some functions it implements the split and saves the sum scores of the test parts
.semPrepareSumScores <- function(X, K) {

  nit <- ncol(X)
  N <- nrow(X)

  S <- rowSums(X)
  partSUMS <- matrix(NA, N, K)

  if (K < nit) {

    k <- split(seq_len(nit), 1:K)

    for (i in 1:K) {
      partSUMS[, i] <- rowSums(X[, k[[i]], drop = FALSE])
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

#' this function prepares the output matrix for the SEM computations
.semPrepareOutMatrix <- function(nit, nc, scrs) {
  scores <- scrs$scores
  counts <- scrs$counts
  out <- matrix(NA, length(scores), 4)
  out[, 1] <- scores
  out[, 3] <- counts
  out[, 4] <- FALSE

  return(out)
}

#' this function computes the SEMs for the different score groups
#' there is a special treatment to collapse some scores into groups if they do not reach the min size
#' out is the output matrix from the previous function
#' S is a vector with the sum scores
#' partSUMS is a matrix with the sum scores of the test parts
#' fun is the function that computes the SEM, this one differs for each method
#' cc is a vector with the number of items in each test part, this is only needed for the Lord's compound method
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
        # actual index in the scores and out-matrix
        firstInd <- firstInd[1]
        nextInd <- ii + (firstInd - 1)
        ind <- which(S %in% scores[ii:nextInd])
        out[ii:nextInd, 2] <- fun(partSUMS, ind, cc)
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
    out[backInd:nextInd, 2] <- fun(partSUMS, ind, cc)
    out[backInd:nextInd, 4] <- TRUE
  }

  colnames(out) <- c("score", "sem", "counts", "collapsed")

  return(out)
}

#' this function computes the confidence intervals for the sum scores
#' computes two CIs, one for the table and one for the plot in JASP
#' out is the output matrix from the sem functions, that is, it contains the sem values
#' scores is a vector with the sum scores
#' ciLevelTable is the confidence level for the table
#' ciLevelPlots is the confidence level for the plots
.semComputeCis <- function(out, scores, ciLevelTable, ciLevelPlots) {

  cis <- matrix(NA, nrow(out), 4)
  zValueTable <- qnorm(1 - (1 - ciLevelTable) / 2)
  zValuePlots <- qnorm(1 - (1 - ciLevelPlots) / 2)
  for (i in 1:nrow(out)) {
    cis[i, ] <- c(scores[i] - zValueTable * out[i, 2], scores[i] + zValueTable * out[i, 2],
                   scores[i] - zValuePlots * out[i, 2], scores[i] + zValuePlots * out[i, 2])
  }

  return(cis)
}

#' this might be the most important function :-
.semOptionsHelper <- function(options, dataset) {

  nc <- length(unique(c(as.matrix(dataset))))

  # Define the options and their corresponding pretty names
  optionsFull <- list(
    thorndike = list(name = "Thorndike",
                     funString = ".semThorn(dataset, nc = nc, caseMin = options$minimumGroupSize, scrs)",
                     dependencies = c("thorndike", "minimumGroupSize")),
    feldt = list(name = "Feldt",
                 funString = ".semFeldt(dataset, K = as.numeric(options$feldtNumberOfSplits), nc = nc, caseMin = options$minimumGroupSize, scrs)",
                 dependencies = c("feldt", "feldtNumberOfSplits", "minimumGroupSize")),
    mollenkopfFeldt = list(name = "Mollenkopf-Feldt",
                           funString = ".semMollenkopfFeldt(dataset, K = as.numeric(options$mollenkopfFeldtNumberOfSplits), nc = nc, n_poly = options$mollenkopfFeldtPolyDegree, scrs)",
                           dependencies = c("mollenkopfFeldt", "mollenkopfFeldtNumberOfSplits", "mollenkopfFeldtPolyDegree")),
    anova = list(name = "ANOVA",
                 funString = ".semAnova(dataset, nc = nc, caseMin = options$minimumGroupSize, scrs)",
                 dependencies = c("anova", "minimumGroupSize")),
    irt = list(name = ifelse(nc == 2, "IRT-2PLM", "IRT-GRM"),
               funString = ".semIrt(dataset, nc = nc, scores)",
               dependencies = "irt"),
    lord = switch(as.character(nc),
                  "2" = list(name = "Lord",
                             funString = ".semLord(ncol(dataset), scrs)",
                             dependencies = "lord"),
                  list(name = NA,
                       funString = NA,
                       dependencies = NA)),
    keats = switch(as.character(nc),
                   "2" = list(name = "Keats",
                              funString = ".semKeats(dataset, options, scrs)",
                              dependencies = "keats"),
                   list(name = NA,
                        funString = NA,
                        dependencies = NA)),
    lord2 = switch(as.character(nc),
                   "2" = list(name = "Lord's compound",
                              funString = ".semLord2(dataset, as.numeric(options$lord2NumberOfSplits), scrs, options$minimumGroupSize)",
                              dependencies = c("lord2", "lord2NumberOfSplits", "minimumGroupSize")),
                   list(name = NA,
                        funString = NA,
                        dependencies = NA))

  )

  # Initialize the output list
  outList <- list()
  # Loop through the options and append to the list if the option is selected
  for (option in names(optionsFull)) {
    if (options[[option]]) {
      outList[[option]] <- optionsFull[[option]]
    }
  }

  options[["selected"]] <- outList
  return(options)
}

