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

RaterAgreement <- function(jaspResults, dataset, options) {

  dataset <- .raterAgreementReadData(dataset, options)

  jaspResults[["table"]] <- .handleIntraclassCorrelation(dataset, options)

  return()
}

# Read in the dataset (copied from .reliabilityReadData)
.raterAgreementReadData <- function(dataset, options) {
  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.numeric = variables, columns.as.factor = NULL, exclude.na.listwise = NULL)
  }
  return(dataset)
}

.handleIntraclassCorrelation <- function (dataset, options) {

  # Check for errors using JASPs internal convenience function
  .hasErrors(
    dataset = dataset,
    type = c("infinity", "negativeValues", "observations"),
    observations.amount = c("< 3"),
    exitAnalysisIfErrors = TRUE
  )

  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettext("Intraclass Correlation"))
  jaspTable$addColumnInfo(name = "type", title = gettext("Type"), type = "string")
  jaspTable$addColumnInfo(name = "ICC", title = gettext("Estimate"), type = "number")
  formattedCIPercent <- format(
    100 * options[["confidenceIntervalValue"]],
    digits = 3,
    drop0trailing = TRUE
  )
  jaspTable$addColumnInfo(
    name = "lower bound",
    title = gettextf("lower %s%% CI", formattedCIPercent),
    type = "number"
  )
  jaspTable$addColumnInfo(
    name = "upper bound",
    title = gettextf("upper %s%% CI", formattedCIPercent),
    type = "number"
  )
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
      type <- paste0(type, "k")
    }

    # Compute the ICC
    full_results <- psych::ICC(
      dataset,
      alpha = 1 - options[["confidenceIntervalValue"]]
    )
    icc_results <- full_results$results

    # Select correct ICC
    icc <- icc_results[icc_results$type == type, ]
    rownames(icc) <- NULL

    # Round all numeric columns
    numeric_columns <- unlist(lapply(icc, is.numeric))
    icc[ , numeric_columns] <- round(icc[ , numeric_columns], digits = 3)

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
        "%s subjects and %s judges. ICC type as referenced by Shrout & Fleiss (1979).",
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
