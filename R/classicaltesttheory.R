#
# Copyright (C) 2013-2023 University of Amsterdam
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

classicalTestTheory <- function(jaspResults, dataset, options, ...) {
  dataset <- .cttReadData(dataset, options)
  ready <- .cttReady(options)

  # Create the summary table
  .cttSummaryTable(dataset, options, jaspResults, ready, position = 1)

  # Create the descriptives table
  .cttDescriptivesTable(dataset, options, jaspResults, ready, position = 3)

  # Create the Cronbachs alpha table
  .cttAlphaTable(dataset, options, jaspResults, ready, position = 5)

  # Create the item statistics table
  .cttItemStatisticsTable(dataset, options, jaspResults, ready, position = 7)

  # Create the histogram of test scores
  .cttHistogram(dataset, options, jaspResults, ready, position = 9)

  # Create the difficulty and discrimination plot
  .cttDiscrimination(dataset, options, jaspResults, ready, position = 11)

  # Create the correlation heatmap
  .cttHeatmap(dataset, options, jaspResults, ready, position = 13)
}

.cttReadData <- function(dataset, options) {
  variables <- unlist(options[["items"]])
  variables <- variables[variables != ""]
  dataset <- jaspBase::excludeNaListwise(dataset, variables)
  .hasErrors(dataset,
    type = c("infinity", "observations"),
    all.target = c(options[["items"]], options[["covariates"]]),
    observations.amount = paste0("< ", nrow(dataset))
  )
  return(dataset)
}

.cttReady <- function(options) {
  ready <- length(options[["items"]]) > 1 && all(unlist(options[["items"]]) != "")
  return(ready)
}

.cttDeps <- function(type = "none") {
  deps <- c("items", "customMaxScore", "tableCronbachsAlphaCI")
  return(deps)
}

.cttState <- function(dataset, options, jaspResults) {
  if (!is.null(jaspResults[["state"]])) {
    return(jaspResults[["state"]]$object)
  }
  p <- try({
    items <- dataset[, unlist(options[["items"]]), drop = FALSE]
    if (!options[["customMaxScore"]]) {
      maxScores <- apply(items, 2, max)
    } else {
      maxScores <- items[1, ]
      items <- items[-1, ]
    }
    sumScores <- rowSums(items)
    maxScore <- sum(maxScores)
    descriptives <- data.frame(
      items = length(options[["items"]]),
      n = nrow(items),
      min = min(sumScores, na.rm = TRUE),
      max = max(sumScores, na.rm = TRUE),
      mean = mean(sumScores, na.rm = TRUE),
      median = median(sumScores, na.rm = TRUE),
      sd = sd(sumScores, na.rm = TRUE),
      skew = moments::skewness(sumScores, na.rm = TRUE),
      kurt = moments::kurtosis(sumScores, na.rm = TRUE)
    )
    rir <- alpha <- numeric(length(options[["items"]]))
    for (i in seq_along(options[["items"]])) {
      rir[i] <- cor(items[, i], rowSums(items[, -i, drop = FALSE]))
      if (ncol(items) > 2) {
        alpha[i] <- ltm::cronbach.alpha(items[, -i, drop = FALSE])[["alpha"]]
      } else {
        alpha[i] <- NA
      }
    }
    itemInfo <- data.frame(
      item = unlist(options[["items"]]),
      score = as.numeric(colMeans(items)),
      sd = as.numeric(apply(items, 2, sd)),
      rit = as.numeric(apply(items, 2, function(x) cor(x, rowSums(items)))),
      rir = rir,
      p = as.numeric(colMeans(items) / maxScores),
      alpha = alpha
    )
    result <- list()
    result[["items"]] <- items
    result[["avgP"]] <- mean(itemInfo[["p"]], na.rm = TRUE)
    result[["avgRit"]] <- mean(itemInfo[["rit"]], na.rm = TRUE)
    result[["avgRir"]] <- mean(itemInfo[["rir"]], na.rm = TRUE)
    result[["maxScore"]] <- maxScore
    result[["maxScores"]] <- maxScores
    result[["sumScores"]] <- sumScores
    result[["descriptives"]] <- descriptives
    result[["itemInfo"]] <- itemInfo
    result[["cronbach"]] <- ltm::cronbach.alpha(items, CI = TRUE, probs = c((1 - options[["tableCronbachsAlphaCI"]]) / 2, options[["tableCronbachsAlphaCI"]] + (1 - options[["tableCronbachsAlphaCI"]]) / 2))
  })
  if (jaspBase:::isTryError(p)) {
    jaspBase:::.quitAnalysis(gettextf("An error occurred in the analysis: %1$s", jaspBase:::.extractErrorMessage(p)))
  }
  jaspResults[["state"]] <- createJaspState(result)
  jaspResults[["state"]]$dependOn(options = .cttDeps())
  return(jaspResults[["state"]]$object)
}

.cttSummaryTable <- function(dataset, options, jaspResults, ready, position) {
  if (options[["explanatoryText"]]) {
    text <- createJaspHtml(gettext("<h3>Explanatory Text: Classical Test Theory</h3>Classical Test Theory (CTT) is a framework for assessing and understanding the reliability and validity of psychological and educational tests. It suggests that an observed test score is composed of two components: the true score (the actual ability being measured) and noise (random variability in test scores). CTT focuses on quantifying and separating these components to evaluate the quality and consistency of a test."))
    text$position <- position
    text$dependOn(options = "explanatoryText")
    jaspResults[["tableSummaryText"]] <- text
  }
  if (!is.null(jaspResults[["tableSummary"]])) {
    return()
  }
  # Table
  tb <- createJaspTable(title = gettext("Test Summary"))
  tb$position <- position + 1
  tb$addColumnInfo(name = "items", title = gettext("Items"), type = "integer")
  tb$addColumnInfo(name = "n", title = gettext("Respondents"), type = "integer")
  tb$addColumnInfo(name = "avgp", title = gettext("Avg. difficulty (P)"), type = "number")
  tb$addColumnInfo(name = "avgrit", title = gettext("Avg. r<sub>item-total</sub>"), type = "number")
  tb$addColumnInfo(name = "avgrir", title = gettext("Avg. r<sub>item-rest</sub>"), type = "number")
  tb$dependOn(options = .cttDeps())
  jaspResults[["tableSummary"]] <- tb
  if (!ready) {
    tb$addFootnote(gettext("The maximum possible test score is ...."))
    return()
  }
  state <- .cttState(dataset, options, jaspResults)
  tb[["items"]] <- state[["descriptives"]][["items"]]
  tb[["n"]] <- state[["descriptives"]][["n"]]
  tb[["avgp"]] <- state[["avgP"]]
  tb[["avgrit"]] <- state[["avgRit"]]
  tb[["avgrir"]] <- state[["avgRir"]]
  tb$addFootnote(gettextf("The maximum possible test score is %1$s.", state[["maxScore"]]))
}

.cttDescriptivesTable <- function(dataset, options, jaspResults, ready, position) {
  if (options[["explanatoryText"]] && options[["tableDescriptives"]]) {
    text <- createJaspHtml(gettext("<h3>Explanatory Text: Descriptive Statistics</h3>The table below provides descriptive statistics the test scores (i.e., sum of item scores). It includes information about the lowest and highest recorded test score, the average (mean), the middle point (median), a measure of how spread out the scores are (standard deviation), and two measures that describe the shape of the score distribution (skewness and kurtosis)."))
    text$position <- position
    text$dependOn(options = c("explanatoryText", "tableDescriptives"))
    jaspResults[["tableDescriptivesText"]] <- text
  }
  if (!is.null(jaspResults[["tableDescriptives"]]) || !options[["tableDescriptives"]]) {
    return()
  }
  tb <- createJaspTable(title = gettext("Descriptive Statistics of Sum Scores"))
  tb$position <- position + 1
  tb$addColumnInfo(name = "min", title = gettext("Min."), type = "number")
  tb$addColumnInfo(name = "max", title = gettext("Max."), type = "number")
  tb$addColumnInfo(name = "mean", title = gettext("Mean"), type = "number")
  tb$addColumnInfo(name = "median", title = gettext("Median"), type = "number")
  tb$addColumnInfo(name = "sd", title = gettext("Std. dev"), type = "number")
  tb$addColumnInfo(name = "skew", title = gettext("Skewness"), type = "number")
  tb$addColumnInfo(name = "kurt", title = gettext("Kurtosis"), type = "number")
  tb$dependOn(options = c(.cttDeps(), "tableDescriptives"))
  jaspResults[["tableDescriptives"]] <- tb
  if (!ready) {
    return()
  }
  state <- .cttState(dataset, options, jaspResults)
  tb[["min"]] <- state[["descriptives"]][["min"]]
  tb[["max"]] <- state[["descriptives"]][["max"]]
  tb[["mean"]] <- state[["descriptives"]][["mean"]]
  tb[["median"]] <- state[["descriptives"]][["median"]]
  tb[["sd"]] <- state[["descriptives"]][["sd"]]
  tb[["skew"]] <- state[["descriptives"]][["skew"]]
  tb[["kurt"]] <- state[["descriptives"]][["kurt"]]
}

.cttAlphaTable <- function(dataset, options, jaspResults, ready, position) {
  if (options[["explanatoryText"]] && options[["tableCronbachsAlpha"]]) {
    text <- createJaspHtml(gettextf("<h3>Explanatory Text: Test Reliability</h3>Cronbach's %1$s is a way to quantify how reliable (i.e., consistent) a test is. It tells us if all the questions in the test are measuring the same construct. This helps to understand if the test is consistent and if it is a good measure of what it is supposed to measure. When interpreting Cronbach's %1$s for a test, a high %1$s value (i.e., close to 1) indicates good internal consistency, suggesting that the test items are reliably measuring the same construct. On the other hand, a low %1$s value (i.e., close to 0) may suggest that some items in the assessment are not contributing to the overall measurement consistency and might need revision or removal. A rule of thumb is that a 'good' test has at least an %1$s value of 0.8.", "\u03B1"))
    text$position <- position
    text$dependOn(options = c("explanatoryText", "tableCronbachsAlpha"))
    jaspResults[["tableCronbachsAlphaText"]] <- text
  }
  if (!is.null(jaspResults[["tableCronbachsAlpha"]]) || !options[["tableCronbachsAlpha"]]) {
    return()
  }
  tb <- createJaspTable(title = gettext("Test Reliability"))
  tb$position <- position + 1
  tb$addColumnInfo(name = "alpha", title = gettextf("Cronbach's %1$s", "\u03B1"), type = "number")
  overtitle <- gettextf("%1$s%% Confidence interval", round(options[["tableCronbachsAlphaCI"]] * 100, 3))
  tb$addColumnInfo(name = "alpha_lower", title = gettext("Lower"), type = "number", overtitle = overtitle)
  tb$addColumnInfo(name = "alpha_upper", title = gettext("Upper"), type = "number", overtitle = overtitle)
  tb$dependOn(options = c(.cttDeps(), "tableCronbachsAlpha"))
  jaspResults[["tableCronbachsAlpha"]] <- tb
  if (!ready) {
    return()
  }
  if (length(options[["items"]]) < 2) {
    tb$addFootnote(gettext("There must be at least two items to compute Cronbach's alpha."))
    return()
  }
  state <- .cttState(dataset, options, jaspResults)
  tb[["alpha"]] <- state[["cronbach"]][["alpha"]]
  tb[["alpha_lower"]] <- as.numeric(state[["cronbach"]][["ci"]][1])
  tb[["alpha_upper"]] <- as.numeric(state[["cronbach"]][["ci"]][2])
}

.cttItemStatisticsTable <- function(dataset, options, jaspResults, ready, position) {
  if (options[["explanatoryText"]] && options[["tableItemStatistics"]]) {
    text <- createJaspHtml(gettextf("<h3>Explanatory Text: Item Information</h3>In educational testing and assessment, several statistics can evaluate the performance and quality of test items. These statistics help us understand how well each question or item is working and contributing to the overall quality of the test. Let's break down these important concepts in a straightforward manner:\n\n<b>1. Avg. Score:</b> This is the typical score that people get on a specific test item. It gives us an idea of how well test-takers are performing on that particular question.\n<b>2. Std. dev:</b> Standard deviation is a measure of how much individual scores on a test item vary from the average score. It helps us see how consistent or spread out the scores are for that item.\n<b>3. Difficulty (P):</b> Difficulty refers to how hard or easy a test item is for the test-takers. We estimate this by looking at the average score that people achieve on the item, and then we divide it by the range of possible scores. Essentially, it helps us understand if the item is generally difficult or not.\n<b>4. r<sub>item-total</sub>, r<sub>item-rest</sub>:</b> These correlations indicate how well a specific test item relates to the overall test score and how it relates to the rest of the test items. In simple terms, they show how closely linked the item is to the overall performance and to the other questions in the test.\n<b>5. %1$s if removed:</b> Cronbach's %1$s is a measure of the internal consistency or reliability of a test. The column represents what happens to the reliability of the test when a particular item is removed. It helps us assess how important or unimportant a specific question is for the overall reliability of the test.", "\u03B1"))
    text$position <- position
    text$dependOn(options = c("explanatoryText", "tableItemStatistics"))
    jaspResults[["tableItemStatisticsText"]] <- text
  }
  if (!is.null(jaspResults[["tableItemStatistics"]]) || !options[["tableItemStatistics"]]) {
    return()
  }
  tb <- createJaspTable(title = gettext("Item Information"))
  tb$position <- position + 1
  tb$addColumnInfo(name = "item", title = gettext("Item"), type = "string")
  tb$addColumnInfo(name = "score", title = gettext("Avg. score"), type = "number")
  tb$addColumnInfo(name = "sd", title = gettext("Std. dev"), type = "number")
  tb$addColumnInfo(name = "p", title = gettext("Difficulty (P)"), type = "number")
  tb$addColumnInfo(name = "rit", title = gettext("r<sub>item-total</sub>"), type = "number")
  tb$addColumnInfo(name = "rir", title = gettext("r<sub>item-rest</sub>"), type = "number")
  tb$addColumnInfo(name = "alpha", title = gettextf("%1$s if removed", "\u03B1"), type = "number")
  tb$dependOn(options = c(.cttDeps(), "tableItemStatistics"))
  jaspResults[["tableItemStatistics"]] <- tb
  if (!ready) {
    return()
  }
  state <- .cttState(dataset, options, jaspResults)
  tb[["item"]] <- state[["itemInfo"]][["item"]]
  tb[["score"]] <- state[["itemInfo"]][["score"]]
  tb[["sd"]] <- state[["itemInfo"]][["sd"]]
  tb[["p"]] <- state[["itemInfo"]][["p"]]
  tb[["rit"]] <- state[["itemInfo"]][["rit"]]
  tb[["rir"]] <- state[["itemInfo"]][["rir"]]
  tb[["alpha"]] <- state[["itemInfo"]][["alpha"]]
}

.cttHistogram <- function(dataset, options, jaspResults, ready, position) {
  if (options[["explanatoryText"]] && options[["plotHistogram"]]) {
    text <- createJaspHtml(gettext("<h3>Explanatory Text: Histogram of Sum Scores</h3>The figure below displays a histogram of the test (i.e., sum) scores. It shows the frequencies of different total scores or sums of responses, which can help you understand the central tendency and variability of the data. Additionally, you can assess patterns or trends in the distribution, such as skewness or bimodality,"))
    text$position <- position
    text$dependOn(options = c("explanatoryText", "plotHistogram"))
    jaspResults[["plotHistogramText"]] <- text
  }
  if (!is.null(jaspResults[["plotHistogram"]]) || !options[["plotHistogram"]]) {
    return()
  }
  fg <- createJaspPlot(title = gettext("Histogram of Sum Scores"), height = 320, width = 480)
  fg$position <- position + 1
  fg$dependOn(options = c(.cttDeps(), "plotHistogram"))
  jaspResults[["plotHistogram"]] <- fg
  if (!ready) {
    return()
  }
  state <- .cttState(dataset, options, jaspResults)
  plotdata <- data.frame(x = seq_len(state[["maxScore"]]), y = rep(NA, state[["maxScore"]]))
  for (i in seq_len(state[["maxScore"]])) {
    plotdata$y[i] <- length(which(state[["sumScores"]] < i & state[["sumScores"]] >= (i - 1)))
  }
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(plotdata$x, min.n = 4)
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(plotdata$y, min.n = 4)
  p <- ggplot2::ggplot(data = plotdata, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_col(fill = "lightgray", col = "black") +
    ggplot2::scale_x_continuous(name = gettext("Sum Score"), breaks = xBreaks, limits = range(xBreaks)) +
    ggplot2::scale_y_continuous(name = gettext("Counts"), breaks = yBreaks, limits = range(yBreaks)) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw()
  fg$plotObject <- p
}

.cttDiscrimination <- function(dataset, options, jaspResults, ready, position) {
  if (options[["explanatoryText"]] && options[["plotItems"]]) {
    text <- createJaspHtml(gettext("<h3>Explanatory Text: Difficulty and Discrimination Plot</h3>The figure below displays a bar chart of the item difficulty (measured by P) and discrimination (measured by r<sub>item-total</sub> or Rit). This figure shows which items were the easiest to answer (left) and which items were the most difficult to answer (right)."))
    text$position <- position
    text$dependOn(options = c("explanatoryText", "plotItems"))
    jaspResults[["plotItemsText"]] <- text
  }
  if (!is.null(jaspResults[["plotItems"]]) || !options[["plotItems"]]) {
    return()
  }
  fg <- createJaspPlot(title = gettext("Item Difficulty and Discrimination Plot"), height = 320, width = 750)
  fg$position <- position
  fg$dependOn(options = c(.cttDeps(), "plotItems"))
  jaspResults[["plotItems"]] <- fg
  if (!ready) {
    return()
  }
  state <- .cttState(dataset, options, jaspResults)
  plotdata <- data.frame(
    x = factor(rep(options[["items"]], 2)),
    y = c(state[["itemInfo"]][["p"]], state[["itemInfo"]][["rit"]]),
    type = c(rep(gettext("Difficulty (P)"), length(options[["items"]])), rep(gettext("Discrimination (Rit)"), length(options[["items"]])))
  )
  plotdata$x <- factor(plotdata$x, levels = plotdata$x[order(subset(plotdata$y, plotdata$type == gettext("Difficulty (P)")))])
  xxName <- gettext("difficulty")
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, plotdata$y), min.n = 4)
  xName <- gettextf("Item (ordered by %1$s)", xxName)
  p <- ggplot2::ggplot(data = plotdata, mapping = ggplot2::aes(x = x, y = y, fill = type)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.75), col = "black", width = 0.75) +
    ggplot2::scale_x_discrete(name = xName) +
    ggplot2::scale_y_continuous(name = NULL, breaks = yBreaks, limits = range(yBreaks)) +
    ggplot2::scale_fill_manual(name = NULL, values = c("firebrick", "dodgerblue")) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "top") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 12))
  fg$plotObject <- p
}

.cttHeatmap <- function(dataset, options, jaspResults, ready, position) {
  if (options[["explanatoryText"]] && options[["plotCorrelationHeatmap"]]) {
    text <- createJaspHtml(gettext("<h3>Explanatory Text: Correlation Heatmap</h3>A correlation heat map displays Pearson correlations between item scores. The Pearson correlation coefficient describes linear correlation between two items. Positive correlations are shown in blue, while negative correlations are shown in red."))
    text$position <- position
    text$dependOn(options = c("explanatoryText", "plotCorrelationHeatmap"))
    jaspResults[["plotCorrelationHeatmapText"]] <- text
  }
  if (!is.null(jaspResults[["plotCorrelationHeatmap"]]) || !options[["plotCorrelationHeatmap"]]) {
    return()
  }
  fg <- createJaspPlot(title = gettext("Item Correlation Heatmap"), height = 40 * length(options[["items"]]), width = 40 * length(options[["items"]]) + 100)
  fg$position <- position + 1
  fg$dependOn(options = c(.cttDeps(), "plotCorrelationHeatmap", "plotCorrelationHeatmapShowValues"))
  jaspResults[["plotCorrelationHeatmap"]] <- fg
  if (!ready) {
    return()
  }
  state <- .cttState(dataset, options, jaspResults)
  plotdata <- as.data.frame(cor(state[["items"]]))
  plotdata$row <- rownames(plotdata)
  plotdata <- reshape2::melt(plotdata, id.var = "row")
  colnames(plotdata) <- c("row", "col", "value")
  plotdata$row <- factor(plotdata$row, levels = options[["items"]])
  plotdata$col <- factor(plotdata$col, levels = options[["items"]])
  p <- ggplot2::ggplot(data = plotdata, mapping = ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_x_discrete(name = NULL) +
    ggplot2::scale_y_discrete(name = NULL) +
    ggplot2::scale_fill_gradient2(name = NULL, low = "#721A1D", mid = "white", high = "#322C89", limits = c(-1, 1)) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right")
  if (options[["plotCorrelationHeatmapShowValues"]]) {
    p <- p + ggplot2::geom_text(label = round(plotdata$value, 2))
  }
  fg$plotObject <- p
}
