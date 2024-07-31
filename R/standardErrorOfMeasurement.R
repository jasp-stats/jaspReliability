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


# sink(file = "~/Downloads/log.txt")
# on.exit(sink(NULL))

  ready <- length(options[["variables"]]) > 1

  dataset <- .semReadData(dataset, options)

  .semErrorCheck(dataset, options, ready)

  .semCreateMainContainer(jaspResults, options)

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

    average <- list(est = as.numeric(sd(rowSums(dataset)) * sqrt(1 - rel)))

    averageState <- createJaspState(average, dependencies = c("userReliability", "reliabilityValue"))
    jaspResults[["semMainContainer"]][["averageState"]] <- averageState
  }

  if (options[["thorndike"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["thorndikeState"]])) {
      out <- .semThorn(dataset, K = 2, nc = nc, caseMin = options[["minimumGroupSize"]], scrs)
      thorndikeState <- createJaspState(out, dependencies = c("thorndike", "minimumGroupSize", "ciLevel"))
      jaspResults[["semMainContainer"]][["thorndikeState"]] <- thorndikeState
    }
  }

  if (options[["feldt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["feldtState"]])) {
      out <- .semFeldt(dataset, K = as.numeric(options[["feldtNumberOfSplits"]]), nc = nc, caseMin = options[["minimumGroupSize"]],
                       scrs)
      feldtState <- createJaspState(out, dependencies = c("feldt", "feldtNumberOfSplits", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["feldtState"]] <- feldtState
    }
  }

  if (options[["mollenkopfFeldt"]]) {
    if (is.null(jaspResults[["semMainContainer"]][["mfState"]])) {
      out <- .semMF(dataset, K = as.numeric(options[["mollenkopfFeldtNumberOfSplits"]]), nc = nc,
                    n_poly = options[["mollenkopfFeldtPolyDegree"]], scrs)
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
      out <- .semIrt(dataset, nc = nc, scores)
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
      out <- .semLord2(dataset, as.numeric(options[["lord2NumberOfSplits"]]), scrs, options[["minimumGroupSize"]])
      lord2State <- createJaspState(out, dependencies = c("lord2", "lord2NumberOfSplits", "minimumGroupSize"))
      jaspResults[["semMainContainer"]][["lord2State"]] <- lord2State
    }
  }

 return()
}


.semComputeSumScoresCi <- function(jaspResults, dataset, options, ready) {

  if (!is.null(jaspResults[["semMainContainer"]][["sumScoreState"]]) ||
      !ready) return()

  if (!options[["sumScoreCiTable"]] && !options[["sumScoreCiPlots"]]) return()

  if (!any(c(options[["thorndike"]], options[["feldt"]], options[["mollenkopfFeldt"]], options[["anova"]],
            options[["irt"]], options[["lord"]], options[["keats"]], options[["lord2"]]))) return()

  scrs <- jaspResults[["semMainContainer"]][["countsState"]]$object
  dtFill <- list(table = data.frame(score = scrs$scores), plots = data.frame(score = scrs$scores))

  # all of this should be so fast we do not need dependency checks...
  if (options[["thorndike"]]) {
    out <- jaspResults[["semMainContainer"]][["thorndikeState"]]$object
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerThorn <- cis[, 1]
    dtFill$table$upperThorn <- cis[, 2]
    dtFill$plots$lowerThorn <- cis[, 3]
    dtFill$plots$upperThorn <- cis[, 4]
  }

  if (options[["feldt"]]) {
    out <- jaspResults[["semMainContainer"]][["feldtState"]]$object
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerFeldt <- cis[, 1]
    dtFill$table$upperFeldt <- cis[, 2]
    dtFill$plots$lowerFeldt <- cis[, 3]
    dtFill$plots$upperFeldt <- cis[, 4]
  }

  if (options[["mollenkopfFeldt"]]) {
    out <- jaspResults[["semMainContainer"]][["mfState"]]$object
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerMoll <- cis[, 1]
    dtFill$table$upperMoll <- cis[, 2]
    dtFill$plots$lowerMoll <- cis[, 3]
    dtFill$plots$upperMoll <- cis[, 4]
  }

  if (options[["anova"]]) {
    out <- jaspResults[["semMainContainer"]][["anovaState"]]$object
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerAnova <- cis[, 1]
    dtFill$table$upperAnova <- cis[, 2]
    dtFill$plots$lowerAnova <- cis[, 3]
    dtFill$plots$upperAnova <- cis[, 4]
  }

  if (options[["irt"]]) {
    nc <- length(unique(c(as.matrix(dataset)))) # needed for IRT
    out <- jaspResults[["semMainContainer"]][["irtState"]]$object$binned
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerIrt <- cis[, 1]
    dtFill$table$upperIrt <- cis[, 2]
    dtFill$plots$lowerIrt <- cis[, 3]
    dtFill$plots$upperIrt <- cis[, 4]
  }

  if (options[["lord"]]) {
    out <- jaspResults[["semMainContainer"]][["lordState"]]$object
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerLord <- cis[, 1]
    dtFill$table$upperLord <- cis[, 2]
    dtFill$plots$lowerLord <- cis[, 3]
    dtFill$plots$upperLord <- cis[, 4]
  }

  if (options[["keats"]]) {
    out <- jaspResults[["semMainContainer"]][["keatsState"]]$object
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerKeats <- cis[, 1]
    dtFill$table$upperKeats <- cis[, 2]
    dtFill$plots$lowerKeats <- cis[, 3]
    dtFill$plots$upperKeats <- cis[, 4]
  }

  if (options[["lord2"]]) {
    out <- jaspResults[["semMainContainer"]][["lord2State"]]$object
    cis <- .semComputeCis(out, scrs$scores, options[["ciLevelTable"]], options[["ciLevelPlots"]])
    dtFill$table$lowerLord2 <- cis[, 1]
    dtFill$table$upperLord2 <- cis[, 2]
    dtFill$plots$lowerLord2 <- cis[, 3]
    dtFill$plots$upperLord2 <- cis[, 4]
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

  nc <- length(unique(c(as.matrix(dataset)))) # needed for IRT

  coefficientTable$addColumnInfo(name = "score", title = gettext("Sum score"), type = "string")
  coefficientTable$addColumnInfo(name = "average", title = gettext("Traditional"), type = "number")

  if (!ready || jaspResults[["semMainContainer"]]$getError()) return()

  dtFill <- data.frame(score = "all")
  average <- jaspResults[["semMainContainer"]][["averageState"]]$object
  dtFill$average <- average$est

  if (!options[["hideTable"]]) {
    if (any(c(options[["thorndike"]], options[["feldt"]], options[["mollenkopfFeldt"]], options[["anova"]],
              options[["irt"]], options[["lord"]], options[["keats"]], options[["lord2"]]))) {

      # the repetition of this seems annoying...
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
        irtTitle <- ifelse(nc == 2, gettext("IRT-2PL"), gettext("IRT-GRM"))
        out <- jaspResults[["semMainContainer"]][["irtState"]]$object$binned
        coefficientTable$addColumnInfo(name = "irt", title = irtTitle, type = "number")
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
  if (options[["thorndike"]]) {
    ciTable$addColumnInfo(name = "lowerThorn", title = gettextf("Lower %s%% CI", ci), type = "number",
                                   overtitle = gettext("Thorndike"))
    ciTable$addColumnInfo(name = "upperThorn", title = gettextf("Upper %s%% CI", ci), type = "number",
                                   overtitle = gettext("Thorndike"))
  }

  if (options[["feldt"]]) {
    ciTable$addColumnInfo(name = "lowerFeldt", title = gettextf("Lower %s%% CI", ci), type = "number",
                          overtitle = gettext("Feldt"))
    ciTable$addColumnInfo(name = "upperFeldt", title = gettextf("Upper %s%% CI", ci), type = "number",
                          overtitle = gettext("Feldt"))
  }

  if (options[["mollenkopfFeldt"]]) {
    ciTable$addColumnInfo(name = "lowerMoll", title = gettextf("Lower %s%% CI", ci), type = "number",
                          overtitle = gettext("Mollenkopf-Feldt"))
    ciTable$addColumnInfo(name = "upperMoll", title = gettextf("Upper %s%% CI", ci), type = "number",
                          overtitle = gettext("Mollenkopf-Feldt"))
  }

  if (options[["anova"]]) {
    ciTable$addColumnInfo(name = "lowerAnova", title = gettextf("Lower %s%% CI", ci), type = "number",
                          overtitle = gettext("Anova"))
    ciTable$addColumnInfo(name = "upperAnova", title = gettextf("Upper %s%% CI", ci), type = "number",
                          overtitle = gettext("Anova"))
  }

  if (options[["irt"]]) {
    nc <- length(unique(c(as.matrix(dataset)))) # needed for IRT
    irtTitle <- ifelse(nc == 2, gettext("IRT-2PL"), gettext("IRT-GRM"))
    ciTable$addColumnInfo(name = "lowerIrt", title = gettextf("Lower %s%% CI", ci), type = "number",
                          overtitle = irtTitle)
    ciTable$addColumnInfo(name = "upperIrt", title = gettextf("Upper %s%% CI", ci), type = "number",
                          overtitle = irtTitle)
  }

  if (options[["lord"]]) {
    ciTable$addColumnInfo(name = "lowerLord", title = gettextf("Lower %s%% CI", ci), type = "number",
                          overtitle = gettext("Lord"))
    ciTable$addColumnInfo(name = "upperLord", title = gettextf("Upper %s%% CI", ci), type = "number",
                          overtitle = gettext("Lord"))
  }

  if (options[["keats"]]) {
    ciTable$addColumnInfo(name = "lowerKeats", title = gettextf("Lower %s%% CI", ci), type = "number",
                          overtitle = gettext("Keats"))
    ciTable$addColumnInfo(name = "upperKeats", title = gettextf("Upper %s%% CI", ci), type = "number",
                          overtitle = gettext("Keats"))
  }

  if (options[["lord2"]]) {
    ciTable$addColumnInfo(name = "lowerLord2", title = gettextf("Lower %s%% CI", ci), type = "number",
                          overtitle = gettext("Lord's compound"))
    ciTable$addColumnInfo(name = "upperLord2", title = gettextf("Upper %s%% CI", ci), type = "number",
                          overtitle = gettext("Lord's compound"))
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

  nc <- length(unique(c(as.matrix(dataset)))) # may be needed for IRT

  pointPlotsContainer <- createJaspContainer(title = gettext("Plots"))
  pointPlotsContainer$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]], options = "pointPlots")
  jaspResults[["semMainContainer"]][["pointPlotsContainer"]] <- pointPlotsContainer
  if (options[["thorndike"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["thorndikeState"]],
                                  title = gettext("Thorndike"))

    pointPlotsContainer[["thorndikePlot"]] <- pl
  }

  if (options[["feldt"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["feldtState"]],
                               title = gettext("Feldt"))

    pointPlotsContainer[["feldtPlot"]] <- pl
  }

  if (options[["mollenkopfFeldt"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["mfState"]],
                               title = gettext("Mollenkopf-Feldt"))

    pointPlotsContainer[["mfPlot"]] <- pl
  }

  if (options[["anova"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["anovaState"]],
                               title = gettext("ANOVA"))

    pointPlotsContainer[["anovaPlot"]] <- pl
  }

  if (options[["irt"]]) {
    irtTitle <- ifelse(nc == 2, gettext("IRT-2PL"), gettext("IRT-GRM"))
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["irtState"]],
                               title = irtTitle, irt = TRUE)

    pointPlotsContainer[["irtPlot"]] <- pl
  }

  if (options[["lord"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["lordState"]],
                                  title = gettext("Lord"))
    pointPlotsContainer[["lordPlot"]] <- pl
  }

  if (options[["keats"]]) {
    pl <- .semMakeSinglePointPlot(resultsObject = jaspResults[["semMainContainer"]][["keatsState"]],
                                  title = gettext("Keats"))

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
    colnames(dat) <- c("score", "sem")
    pl <- ggplot2::ggplot(dat) +
      ggplot2::geom_point(ggplot2::aes(x = score, y = sem), size = 2.5) +
      ggplot2::labs(x = gettext("Sum Score"), y = gettext("sem"))
  } else {
    dat <- as.data.frame(resultsObject$object$unbinned)
    colnames(dat) <- c("score", "sem")
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
      irtTitle <- ifelse(nc == 2, gettext("IRT-2PL"), gettext("IRT-GRM"))
      out <- jaspResults[["semMainContainer"]][["irtState"]]$object$unbinned
      dtIRT <- data.frame(score = out[, 1], sem = out[, 2], Type = irtTitle)
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

    plot <- createJaspPlot(jaspGraphs::themeJaspRaw(pl, legend.position = "right"), title = gettext("Combined plot"),
                           width = 600)
    plot$dependOn(optionsFromObject = jaspResults[["semMainContainer"]][["coefficientTable"]])
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

  if (options[["thorndike"]]) {
    ciDataUse <- ciData[, c("score", "lowerThorn", "upperThorn")]
    pl <- .semMakeSingleCiPlot(title = gettext("Thorndike"), ciData = ciDataUse,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["thorndikePlot"]] <- pl
  }

  if (options[["feldt"]]) {
    ciDataUse <- ciData[, c("score", "lowerFeldt", "upperFeldt")]
    pl <- .semMakeSingleCiPlot(title = gettext("Feldt"), ciData = ciDataUse,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["feldtPlot"]] <- pl
  }

  if (options[["mollenkopfFeldt"]]) {
    ciDataUse <- ciData[, c("score", "lowerMoll", "upperMoll")]
    pl <- .semMakeSingleCiPlot(title = gettext("Mollenkopf-Feldt"), ciData = ciDataUse,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["mfPlot"]] <- pl
  }

  if (options[["anova"]]) {
    ciDataUse <- ciData[, c("score", "lowerAnova", "upperAnova")]
    pl <- .semMakeSingleCiPlot(title = gettext("Anova"), ciData = ciDataUse,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["anovaPlot"]] <- pl
  }

  if (options[["irt"]]) {
    irtTitle <- ifelse(nc == 2, gettext("IRT-2PL"), gettext("IRT-GRM"))
    outUnbinned <- jaspResults[["semMainContainer"]][["irtState"]]$object$unbinned
    cis <- .semComputeCis(outUnbinned, outUnbinned[, 1], options[["ciLevel"]])
    lowerIrtUnbinned <- cis[, 1]
    upperIrtUnbinned <- cis[, 2]
    ciDataUse <- cbind(outUnbinned[, 1], lowerIrtUnbinned, upperIrtUnbinned)
    pl <- .semMakeSingleCiPlot(title = irtTitle, ciData = ciDataUse, irt = TRUE,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["irtPlot"]] <- pl
  }

  if (options[["lord"]]) {
    ciDataUse <- ciData[, c("score", "lowerLord", "upperLord")]
    pl <- .semMakeSingleCiPlot(title = gettext("Lord"), ciData = ciDataUse,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["lordPlot"]] <- pl
  }

  if (options[["keats"]]) {
    ciDataUse <- ciData[, c("score", "lowerKeats", "upperKeats")]
    pl <- .semMakeSingleCiPlot(title = gettext("Keats"), ciData = ciDataUse,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["keatsPlot"]] <- pl
  }

  if (options[["lord2"]]) {
    ciDataUse <- ciData[, c("score", "lowerLord2", "upperLord2")]
    pl <- .semMakeSingleCiPlot(title = gettext("Lord2"), ciData = ciDataUse,
                               cutoff = ifelse(options[["sumScoreCiPlotsCutoff"]], options[["sumScoreCiPlotsCutoffValue"]], NA))
    ciPlotsContainer[["lord2Plot"]] <- pl
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


##### Sem compute functions #####


.semThorn <- function(X, K, nc, caseMin, scrs) {

  prep <- .semPrepareSumScores(X, K)
  partSUMS <- prep$partSUMS
  S <- prep$S
  out <- .semPrepareOutMatrix(ncol(X), nc, scrs)
  fun <- function(partSUMS, ind, cc) {
    return(sd(partSUMS[ind, 1] - partSUMS[ind, 2]))
  }

  out <- .semComputeWithCaseMin(out, S, caseMin, partSUMS, fun)

  return(out)
}


.semFeldt <- function(X, K, nc, caseMin, scrs) {

  prep <- .semPrepareSumScores(X, K)
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


.semMF <- function(X, K, nc, n_poly, scrs) {

  prep <- .semPrepareSumScores(X, K)
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

  prep <- .semPrepareSumScores(X, K = 1)
  S <- prep$S

  out <- .semPrepareOutMatrix(ncol(X), nc, scrs)
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
.semLord2 <- function(X, K, scrs, caseMin) {

  nc <- 2
  prep <- .semPrepareSumScores(X, K)
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

.semPrepareSumScores <- function(X, K) {

  nit <- ncol(X)
  N <- nrow(X)

  S <- rowSums(X)
  partSUMS <- matrix(NA, N, K)

  if (K < nit) {

    # firstSplit <- ceiling(nit / K)

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

  return(out)
}


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
