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
cohenFleissKappa <- function(jaspResults, dataset, options) {
  
  ready <- length(options[["variables"]]) > 1
  
  dataset <- .cohensFleissKappaReadData(dataset, options)
    
  if(options[["cohensKappa"]])
    jaspResults[["cohensKappa"]] <- .cohensKappa(dataset, options, ready)
  
  if(options[["fleissKappa"]])
    jaspResults[["fleissKappa"]] <- .fleissKappa(dataset, options, ready)
  
  return()
}

# Read in the dataset (copied from .reliabilityReadData)
.cohensFleissKappaReadData <- function(dataset, options) {
  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.factor = variables)
  }
  return(dataset)
}

.cohensKappa <- function(dataset, options, ready){
  
  weighted <- options[["cohensWeightedOrNot"]] == "cohensWeighted"
  
  if(weighted){
    weightedString <- "Weighted"
  }else{
    weightedString <- "Unweighted"
  }
  
  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Cohen's %s Kappa", casefold(weightedString, upper = F)))
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
  
  
  if(ready){
    #Create all pairs
    variables <- colnames(dataset)
    allPairs <- .getPairs(variables)
    allDataframes <- .getDataframes(dataset, allPairs)
    allPairStrings <- paste(allPairs$x, allPairs$y, sep = " - ")
    
    #calculate Cohen's Kappas
    allKappaData <- lapply(allDataframes, psych::cohen.kappa, alpha = 1 - options[["kappaConfidenceIntervalValue"]])
    
    #Extract Kappas and CIs
    allKappas <- c()
    allLowerBounds <- c()
    allUpperBounds <- c()
    
    for(i in allKappaData){
      kappaData <- i$confid
      if(weighted){
        k <- 2
      }else{
        k <-1
      }
      allKappas <- c(allKappas, kappaData[k,2])
      allLowerBounds <- c(allLowerBounds, kappaData[k,1])
      allUpperBounds <- c(allUpperBounds, kappaData[k,3])
    }
    
    tableData <- list("ratings" = allPairStrings,
                      "cKappa" = allKappas)
    
    if (options[["kappaIntervalOn"]]) {
      jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
      jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
      tableData[["CIL"]] <- allLowerBounds
      tableData[["CIU"]] <- allUpperBounds
    }
    
    
    jaspTable$setData(tableData)
    jaspTable$addFootnote(gettextf('%i subjects/items and %i judges/measurements.', nrow(dataset), ncol(dataset)))
    
    #if weighted kappa option is on but data only has 2 levels
    if(weighted && length(levels(unlist(dataset))) < 3)
      jaspTable$addFootnote(gettext('Weighted Kappa is only meaningful for more than 2 levels. Calculated weighted Kappa is the same as unweighted Kappa.'))
    
  }
  
  return(jaspTable)
}

.fleissKappa <- function(dataset, options, ready){
  
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
  
  
  # if (options[["kappaIntervalOn"]]) {
  #   jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
  #   jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI for Kappa", formattedCIPercent))
  # }
  
  if(ready){
    #calculate Fleiss' Kappa
    allKappaData <- irr::kappam.fleiss(dataset, detail = T)
    overallKappa <- allKappaData$value
    categoryKappas <- allKappaData$detail[,1]
    
    cateogries <- sort(unique(unlist(dataset)))
    ratings <- c("Overall", as.character(cateogries))
    
    jaspTable$setData(list("ratings" = ratings,
                           "fKappa" = c(overallKappa, categoryKappas)))
    
    # "CIL"    = c(overallKappa, categoryKappas),
    # "CIU"    = c(overallKappa, categoryKappas))
    
    jaspTable$addFootnote(gettextf('%i subjects/items and %i judges/measurements.', nrow(dataset), ncol(dataset)))
  }

  return(jaspTable)
  
  
}

.getPairs <- function(variable){
  df <- data.frame()
  for(i in 1:length(variable)){
    for(j in variable[i:length(variable)]){
      firstVariable <- variable[i]
      secondVariable <- j
      if(firstVariable != secondVariable){
        row <- list(x = firstVariable, y = secondVariable)
        df <- rbind(df, row)
      }
    }
  }
  return(df)
}


.getDataframes <- function(dataset, pairs){
  l <- list()
  for(i in 1:nrow(pairs)){
    l[[i]] <- data.frame(x = dataset[pairs$x[i]], y = dataset[pairs$y[i]])
  }
  return(l)
}
