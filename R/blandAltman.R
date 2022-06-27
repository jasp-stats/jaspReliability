
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
#' @export
blandAltman <- function(jaspResults, dataset, options) {

  dataset <- .blandAltmanReadData(dataset, options)

  .blandAltmanPlots(jaspResults, dataset, options)
  .blandAltmanTable(jaspResults, dataset, options)

  return()
}

# Read in the dataset (copied from .reliabilityReadData)
.blandAltmanReadData <- function(dataset, options) {
  variables <- unlist(options[["pairs"]])
  variables <- variables[variables != ""]
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.numeric = variables)
  }
  return(dataset)
}

.blandAltmanPlots <- function(jaspResults, dataset, options) {

  if (length(options[["pairs"]]) < 1) return()

  if (is.null(jaspResults[["plotsBlandAltman"]])) {
    jaspResults[["plotsBlandAltman"]] <- createJaspContainer(gettext("Bland-Altman Plots"))
    jaspResults[["plotsBlandAltman"]]$dependOn(c("ciDisplay", "ciShading", "ciValue"))
    subcontainer <- jaspResults[["plotsBlandAltman"]]
  } else {
    subcontainer <- jaspResults[["plotsBlandAltman"]]
  }

  for (pair in options[["pairs"]]) {
    title <- paste(pair, collapse = " - ")
    if (!is.null(subcontainer[[title]]))
      next
    blandAltmanPlots <- createJaspPlot(title = title, width = 600, height = 420)
    blandAltmanPlots$dependOn(optionContainsValue = list(pairs = pair))
    subcontainer[[title]] <- blandAltmanPlots

    if (pair[[1]] == "" || pair[[2]] == "") {
      blandAltmanPlots$setError(gettext("Please provide another variable"))
      next
    }

    subData <- data.frame(dataset[, pair])
    errorMessage <- .baCheckErrors(subData, c(pair[[1]], pair[[2]]), obsAmount = "< 2")
    if (!is.null(errorMessage)) {
      blandAltmanPlots$setError(gettextf("Plotting not possible: %s", errorMessage))
      next
    }

    subData <- na.omit(subData)
    blandAltmanPlots$plotObject <- .createBlandAltmanPlot(subData, options)
  }

  return()

}

.createBlandAltmanPlot <- function(dataset, options) {

  ba <- .blandAltmanStats(dataset = dataset, options = options)
  values <- data.frame(m = ba[["means"]], d = ba[["diffs"]])

  # Bland-Altman Plot
  yBreaks <- if (options[["ciDisplay"]]) {
    jaspGraphs::getPrettyAxisBreaks(c(ba[["diffs"]], ba[["CiLines"]]))
  } else {
    jaspGraphs::getPrettyAxisBreaks(c(ba[["diffs"]], ba[["lines"]]))
  }
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(ba[["means"]])

  p <- ggplot2::ggplot(values, ggplot2::aes(x = m, y = d)) +
    ggplot2::geom_hline(yintercept = ba[["lines"]], linetype = 2, size = 1) +
    ggplot2::xlab("Mean of Measurements") +
    ggplot2::ylab("Difference of Measurements")

  if (options[["ciDisplay"]]) {
    p <- p + ggplot2::geom_hline(yintercept = ba[["CiLines"]], linetype = 2, size = 0.5)

    if (options[["ciShading"]]) {
      p <- p + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                                 ymin = ba[["CiLines"]][3],
                                 ymax = ba[["CiLines"]][4],
                                 fill = "blue", alpha = 0.3) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = ba[["CiLines"]][5],
                          ymax = ba[["CiLines"]][6],
                          fill = "green", alpha = 0.3) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = ba[["CiLines"]][1],
                          ymax = ba[["CiLines"]][2],
                          fill = "red", alpha = 0.3)
    }
  }

  p <- p + ggplot2::geom_point(shape = 21, size = 3, colour = "black", fill = "grey") +
    ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
    ggplot2::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks))+
    jaspGraphs::geom_rangeframe(sides = "rbl") +
    jaspGraphs::themeJaspRaw() +
    ggplot2::theme(plot.margin = ggplot2::margin(5))

  return(p)
}

.blandAltmanTable <- function(jaspResults, dataset, options) {

  if (!options[["blandAltmanTable"]])
    return()

  ready <- length(options[["pairs"]]) > 0

  if (is.null(jaspResults[["tabBlandAltman"]])) {
    jaspResults[["tabBlandAltman"]] <- createJaspContainer(gettext("Bland-Altman Tables"))
    jaspResults[["tabBlandAltman"]]$dependOn(c("blandAltmanTable", "ciDisplay", "ciValue"))
    subcontainer <- jaspResults[["tabBlandAltman"]]
  } else {
    subcontainer <- jaspResults[["tabBlandAltman"]]
  }

  if (ready) {
    for (pair in options[["pairs"]]) {
      title <- paste(pair, collapse = " - ")
      if (!is.null(subcontainer[[title]]))
        next
      tablesBlandAltman <- createJaspTable(title = title)
      tablesBlandAltman$dependOn(optionContainsValue = list(pairs = pair))
      tablesBlandAltman$addColumnInfo(name = "names", title = gettext("Bias & Limits"), type = "string")
      tablesBlandAltman$addColumnInfo(name = "agree", title = gettext("Point Value"), type = "number")
      subcontainer[[title]] <- tablesBlandAltman

      intervalLow <- gettextf("Mean difference - 1.96 SD")
      intervalUp <- gettextf("Mean difference + 1.96 SD")
      allData <- data.frame(names = c(intervalUp, gettext("Mean difference"), intervalLow))
      formattedCIPercent <- format(
        100 * options[["ciValue"]],
        digits = 3,
        drop0trailing = TRUE
      )

      if (options[["ciDisplay"]]) {
        tablesBlandAltman$addColumnInfo(
          name = "lowerBound",
          title = gettextf("Lower %s%% CI", formattedCIPercent),
          type = "number"
        )
        tablesBlandAltman$addColumnInfo(
          name = "upperBound",
          title = gettextf("Upper %s%% CI", formattedCIPercent),
          type = "number"
        )
      }

      if (pair[[1]] == "" || pair[[2]] == "") {
        tablesBlandAltman$addFootnote(gettext("Please provide another variable"))
        next
      }

      subData <- data.frame(dataset[, unlist(pair)])
      errorMessage <- .baCheckErrors(subData, c(pair[[1]], pair[[2]]), obsAmount = "< 2")
      if (!is.null(errorMessage)) {
        tablesBlandAltman$addFootnote(gettextf("Calculations not possible: %s", errorMessage))
        next
      }

      subData <- na.omit(subData)
      ba <- .blandAltmanStats(dataset = subData, options = options)

      linesData <- data.frame(agree = c(ba[["lines"]][3], ba[["lines"]][2], ba[["lines"]][1]))
      allData <- cbind(allData, linesData)

      if (options[["ciDisplay"]]) {
        lowCiData <- data.frame(lowerBound = c(ba[["CiLines"]][5], ba[["CiLines"]][3], ba[["CiLines"]][1]))
        topCiData <- data.frame(upperBound = c(ba[["CiLines"]][6], ba[["CiLines"]][4], ba[["CiLines"]][2]))
        allData <- cbind(allData, lowCiData, topCiData)
      }

      tablesBlandAltman$setData(allData)
      tablesBlandAltman$addFootnote(gettextf("%s subjects and 2 measurements.", nrow(dataset)))
    }
  }
  return()

}

.blandAltmanStats <- function(dataset, options) {

  # Calculating main components
  ci <- options[["ciValue"]]
  diffs <- dataset[[1]] - dataset[[2]]
  means <- (dataset[[1]] + dataset[[2]])/2

  obs <- length(dataset[[1]])
  meanDiffs <- mean(diffs)
  criticalDiff <- 1.96 * sd(diffs) #1.96 used instead of suggested 2 by Bland & Altman
  lowerLimit <- meanDiffs - criticalDiff
  upperLimit <- meanDiffs + criticalDiff
  lines <- c(lowerLimit, meanDiffs, upperLimit)

  t1 <- qt((1 - ci)/2, df = obs - 1)
  t2 <- qt((ci + 1)/2, df = obs - 1)

  # CI calculations based on Bland & Altman, 1986
  lowerLimitCiLower <- lowerLimit + t1 * sqrt(sd(diffs)^2 * 3/obs)
  lowerLimitCiUpper <- lowerLimit + t2 * sqrt(sd(diffs)^2 * 3/obs)
  meanDiffCiLower <- meanDiffs + t1 * sd(diffs)/sqrt(obs)
  meanDiffCiUpper <- meanDiffs + t2 * sd(diffs)/sqrt(obs)
  upperLimitCiLower <- upperLimit + t1 * sqrt(sd(diffs)^2 * 3/obs)
  upperLimitCiUpper <- upperLimit + t2 * sqrt(sd(diffs)^2 * 3/obs)
  CiLines <- c(lowerLimitCiLower, lowerLimitCiUpper, meanDiffCiLower, meanDiffCiUpper, upperLimitCiLower, upperLimitCiUpper)

  return(list(means = means, diffs = diffs, lines = lines, CiLines = CiLines))
}

.baCheckErrors <- function(dataset, vars, obsAmount) {
  errors <- .hasErrors(dataset, all.target = vars, message = "short", type = c("infinity", "observations"), observations.amount = obsAmount)
  if (!isFALSE(errors))
    return(errors$message)

  return(NULL)
}

