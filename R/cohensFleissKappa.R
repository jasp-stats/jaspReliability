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
cohensFleissKappa <- function(jaspResults, dataset, options) {
  
  ready <- length(options[["variables"]]) > 1
  
  dataset <- .readDataCohensFleissKappa(dataset, options)
    
  if (options[["cohensKappa"]])
    jaspResults[["cohensKappa"]] <- .computeCohensKappaTable(dataset, options, ready)
  
  if (options[["fleissKappa"]])
    jaspResults[["fleissKappa"]] <- .computeFleissKappaTable(dataset, options, ready)
  
  return()
}

# Read in the dataset (copied from .reliabilityReadData)
.readDataCohensFleissKappa <- function(dataset, options) {
  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.factor = variables)
  }
  return(dataset)
}

.computeCohensKappaTable <- function(dataset, options, ready) {
  
  weighted <- options[["cohensWeightedOrNot"]] == "cohensWeighted"
  
  weightedString <- ifelse(weighted, "Weighted", "Unweighted")
  
  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Cohen's %s Kappa", casefold(weightedString)))
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "cKappa", title = gettextf("%s Kappa", weightedString), type = "number")
  jaspTable$position <- 1
  
  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "cohensKappa",
      "cohensWeightedOrNot",
      "kappaIntervalOn",
      "kappaConfidenceIntervalValue"
    )
  )
  
  
  formattedCIPercent <- format(
    100 * options[["kappaConfidenceIntervalValue"]],
    digits = 3,
    drop0trailing = TRUE
  )
  
  
  if (ready) {
    #calculate Cohen's Kappas
    nPairs <- ncol(dataset) * (ncol(dataset) - 1) / 2
    out_kappa <- psych::cohen.kappa(dataset, alpha = 1 - options[["kappaConfidenceIntervalValue"]])
    if (nPairs == 1) {
      allKappaData <- list(out_kappa)
      allPairStrings <- paste(options[["variables"]], collapse = " - ")
    } else {
      allKappaData <- out_kappa[2:(nPairs + 1)]
      allPairStrings <- sub(" ", " - ", names(out_kappa[2:(nPairs + 1)]))
    }

    #Extract Kappas and CIs
    allKappas <- c()
    allLowerBounds <- c()
    allUpperBounds <- c()
    
    for (i in allKappaData) {
      kappaData <- i$confid
      k <- ifelse(weighted, 2, 1)
      allKappas <- c(allKappas, kappaData[k, 2])
      allLowerBounds <- c(allLowerBounds, kappaData[k, 1])
      allUpperBounds <- c(allUpperBounds, kappaData[k, 3])
      averageKappa <- ifelse(weighted, out_kappa$av.wt, out_kappa$av.kappa)
    }
    
    tableData <- list("ratings" = c("Average Kappa", allPairStrings),
                      "cKappa" = c(averageKappa, allKappas))
    
    if (options[["kappaIntervalOn"]]) {
      jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
      jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
      tableData[["CIL"]] <- c(NA, allLowerBounds)
      tableData[["CIU"]] <- c(NA, allUpperBounds)
    }
    
    
    jaspTable$setData(tableData)
    jaspTable$addFootnote(gettextf('%i subjects/items and %i raters/measurements.', nrow(dataset), ncol(dataset)))
    
    #if weighted kappa option is on but data only has 2 levels
    if (weighted && length(levels(unlist(dataset))) < 3)
      jaspTable$addFootnote(gettext('If there are only 2 levels, weighted kappa is equal to unweighted kappa.'))
  }
  return(jaspTable)
}

.computeFleissKappaTable <- function(dataset, options, ready) {
  
  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Fleiss' Kappa"))
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "fKappa", title = gettext("Fleiss' Kappa"), type = "number")
  jaspTable$position <- 2
  
  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "fleissKappa",
      "kappaIntervalOn",
      "kappaConfidenceIntervalValue"
    )
  )

  
  formattedCIPercent <- format(
    100 * options[["kappaConfidenceIntervalValue"]],
    digits = 3,
    drop0trailing = TRUE
  )
  
  if (ready) {
    #calculate Fleiss' Kappa
    allKappaData <- .fleissKappaMod(dataset, detail = TRUE)
    overallKappa <- allKappaData$value
    categoryKappas <- allKappaData$detail[, 1]
    overallSE <- allKappaData$se
    categorySE <- allKappaData$se_cat
    alpha <- 1 - options[["kappaConfidenceIntervalValue"]]
    
    # for nominal text data we want the rating text to be displayed:
    # if the transformation to numeric goes wrong:
    if (anyNA(as.numeric(as.character(unique(unlist(dataset)))))) {
      categories <- sort(unique(unlist(dataset)))
    } else { # if data is ordinal (can be transformed to numeric)
      categories <- sort(as.numeric(as.character(unique(unlist(dataset)))))
    }
    ratings <- c("Overall", as.character(categories))
    # when the data was read with columns as factors the kappa function
    # would not order the categories numerically, this is a fix
    categoryKappas <- categoryKappas[as.character(categories)]
    
    tableData <- list("ratings" = ratings,
                      "fKappa" = c(overallKappa, categoryKappas))
    if (options[["kappaIntervalOn"]]) {
      overallCI <- overallKappa + c(-1,1)*qnorm(1 - alpha/2) * overallSE 
      categoryCIL <- categoryKappas - qnorm(1 - alpha/2) * categorySE
      categoryCIU <- categoryKappas + qnorm(1 - alpha/2) * categorySE
      jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
      jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
      tableData[["CIL"]] <- c(overallCI[1], categoryCIL)
      tableData[["CIU"]] <- c(overallCI[2], categoryCIU)
      
      jaspTable$addFootnote(gettext('Confidence intervals are asymptotic.'))
    }
    
    jaspTable$setData(tableData)
    
    jaspTable$addFootnote(gettextf('%i subjects/items and %i raters/measurements.', nrow(dataset), ncol(dataset)))
  }
  return(jaspTable)
}
