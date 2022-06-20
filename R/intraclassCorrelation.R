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
intraclassCorrelation <- function(jaspResults, dataset, options) {

  dataset <- .intraclassCorrelationReadData(dataset, options)

  jaspResults[["table"]] <- .handleIntraclassCorrelation(dataset, options)

  .descriptivesBlandAltman(     jaspResults, dataset, options)
  .descriptivesBlandAltmanTable(jaspResults, dataset, options)

  return()
}

# Read in the dataset (copied from .reliabilityReadData)
.intraclassCorrelationReadData <- function(dataset, options) {
  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.numeric = variables)
  }
  return(dataset)
}

.handleIntraclassCorrelation <- function(dataset, options) {

  # Check for errors using JASPs internal convenience function
  jaspBase::.hasErrors(
    dataset = dataset,
    type = c("infinity", "observations"),
    observations.amount = c("< 2"),
    exitAnalysisIfErrors = TRUE
  )

  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettext("Intraclass Correlation"))
  jaspTable$addColumnInfo(name = "type", title = gettext("Type"), type = "string")
  jaspTable$addColumnInfo(name = "ICC", title = gettext("Point Estimate"), type = "number")
  formattedCIPercent <- format(
    100 * options[["confidenceIntervalValue"]],
    digits = 3,
    drop0trailing = TRUE
  )

  if (options[["intervalOn"]]) {
    jaspTable$addColumnInfo(
      name = "lower bound",
      title = gettextf("Lower %s%% CI", formattedCIPercent),
      type = "number"
    )
    jaspTable$addColumnInfo(
      name = "upper bound",
      title = gettextf("Upper %s%% CI", formattedCIPercent),
      type = "number"
    )
  }
  jaspTable$dependOn(
    options = c(
      "variables",
      "intervalOn",
      "confidenceIntervalValue",
      "iccType",
      "iccRatingAverage"
    )
  )

  # Only run ICC analysis if conditions are met,
  # but still return table elsewise for better feedback
  if (ncol(dataset) >= 2) {
    # Get the ICC type e.g. "ICC1k"
    type <- toupper(options["iccType"])
    if (options[["iccRatingAverage"]]) {
      type <- paste0(type, ",k")
    } else {
      type <- paste0(type, ",1")
    }

    # Compute the ICC
    full_results <- psych::ICC(
      dataset,
      alpha = 1 - options[["confidenceIntervalValue"]]
    )
    icc_results <- full_results$results
    icc_results$type <- c("ICC1,1", "ICC2,1", "ICC3,1", "ICC1,k", "ICC2,k", "ICC3,k")
    # Select correct ICC
    icc <- icc_results[icc_results$type == type, ]
    rownames(icc) <- NULL

    # Only report relevant columns
    cols <- c("type", "ICC")
    # Report CI
    if (options[["intervalOn"]]) {
      cols <- c(cols, "lower bound", "upper bound")
    }
    icc <- icc[, cols]

    # Add results to table
    jaspTable$setData(icc)
    jaspTable$addFootnote(
      gettextf(
        "%s subjects and %s raters/measurements. ICC type as referenced by Shrout & Fleiss (1979).",
        full_results$n.obs,
        full_results$n.judge
      )
    )
  } else {
    jaspTable$addFootnote(
      gettext("Please select at least 2 columns")
    )
  }

  return(jaspTable)
}

.descriptivesBlandAltman <- function(jaspResults, dataset, options) {

  if (!options[["descriptivesBlandAltman"]])
    return()

  ready <- length(options[["pairs"]]) > 0

  if (is.null(jaspResults[["plotsBlandAltman"]])) {
    jaspResults[["plotsBlandAltman"]] <- createJaspContainer(gettext("Bland-Altman Plots"))
    jaspResults[["plotsBlandAltman"]]$dependOn(c("descriptivesBlandAltman", "ciDisplay", "ciShading", "ciValue"))
    subcontainer <- jaspResults[["plotsBlandAltman"]]
  } else {
    subcontainer <- jaspResults[["plotsBlandAltman"]]
  }

  if (ready) {
    for (pair in options[["pairs"]]) {
      title <- paste(pair, collapse = " - ")
      if (!is.null(subcontainer[[title]]))
        next
      descriptivesBlandAltman <- createJaspPlot(title = title, width = 600, height = 420)
      descriptivesBlandAltman$dependOn(optionContainsValue = list(pairs = pair))
      subcontainer[[title]] <- descriptivesBlandAltman

      if (pair[[1]] == "" || pair[[2]] == "") {
        descriptivesBlandAltman$setError(gettext("Please provide another variable"))
        next
      }

      subData <- data.frame(dataset[, unlist(pair)])
      errorMessage <- .baCheckErrors(subData, c(pair[[1]], pair[[2]]), obsAmount = "< 2")
      if (!is.null(errorMessage)) {
        descriptivesBlandAltman$setError(gettextf("Plotting not possible: %s", errorMessage))
        next
      }

      subData <- na.omit(subData)
      descriptivesBlandAltman$plotObject <- .descriptivesBlandAltmanPlot(subData, options)
    }
  }
  return()

}

.descriptivesBlandAltmanPlot <- function(dataset, options) {

  ba <- .descriptivesBlandAltmanStats(dataset = dataset, options = options)
  values <- data.frame(m = ba$means, d = ba$diffs)

  # Bland-Altman Plot
  if (options[["ciDisplay"]]) {
    sum <- c(ba$diffs, ba$CiLines)
    yBreaks <- jaspGraphs::getPrettyAxisBreaks(sum)
  } else {
    sum <- c(ba$diffs, ba$lines)
    yBreaks <- jaspGraphs::getPrettyAxisBreaks(sum)
  }
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(ba$means)

  p <- ggplot2::ggplot(values, ggplot2::aes(x = m, y = d)) +
    ggplot2::geom_hline(yintercept = ba$lines, linetype = 2, size = 1) +
    ggplot2::xlab("Mean of Measurements") +
    ggplot2::ylab("Difference of Measurements")

  if (options[["ciDisplay"]]) {
    p <- p + ggplot2::geom_hline(yintercept = ba$CiLines, linetype = 2, size = 0.5)

    if (options[["ciShading"]]) {
      p <- p + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                                 ymin = ba$CiLines[3],
                                 ymax = ba$CiLines[4],
                                 fill = "blue", alpha = 0.3) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = ba$CiLines[5],
                          ymax = ba$CiLines[6],
                          fill = "green", alpha = 0.3) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = ba$CiLines[1],
                          ymax = ba$CiLines[2],
                          fill = "red", alpha = 0.3)
    }
  }

  p <- p + ggplot2::geom_point(shape = 21, size = 3, colour = "black", fill = "grey") +
    ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
    ggplot2::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks))+
    jaspGraphs::geom_rangeframe(sides = "rbl") +
    jaspGraphs::themeJaspRaw() +
    ggplot2::theme(plot.margin = ggplot2::margin(5))
}

.descriptivesBlandAltmanTable <- function(jaspResults, dataset, options) {

  if (!options[["descriptivesBlandAltmanTable"]])
    return()

  ready <- length(options[["pairs"]]) > 0

  if (is.null(jaspResults[["tabBlandAltman"]])) {
    jaspResults[["tabBlandAltman"]] <- createJaspContainer(gettext("Bland-Altman Tables"))
    jaspResults[["tabBlandAltman"]]$dependOn(c("descriptivesBlandAltmanTable", "ciDisplay", "ciValue"))
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
      ba <- .descriptivesBlandAltmanStats(dataset = subData, options = options)

      linesData <- data.frame(agree = c(ba$lines[3], ba$lines[2], ba$lines[1]))
      allData <- cbind(allData, linesData)

      if (options[["ciDisplay"]]) {
        lowCiData <- data.frame(lowerBound = c(ba$CiLines[5], ba$CiLines[3], ba$CiLines[1]))
        topCiData <- data.frame(upperBound = c(ba$CiLines[6], ba$CiLines[4], ba$CiLines[2]))
        allData <- cbind(allData, lowCiData, topCiData)
      }

      tablesBlandAltman$setData(allData)
      tablesBlandAltman$addFootnote(gettextf("%s subjects and 2 measurements.", nrow(dataset)))
    }
  }
  return()

}

.descriptivesBlandAltmanStats <- function(dataset, options) {

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





