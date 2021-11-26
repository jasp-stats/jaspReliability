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

  dataset <- .intraclassCorrelationReadData(dataset, options)
  
  if(options[["cohensKappa"]])
    jaspResults[["cohensKappa"]] <- .cohensKappa(dataset, options)
  
  if(options[["fleissKappa"]])
    jaspResults[["fleissKappa"]] <- .fleissKappa(dataset, options)
  
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

.handleIntraclassCorrelation <- function(dataset, options) {

  # Check for errors using JASPs internal convenience function
  .hasErrors(
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

  return(jaspTable)
}


.cohensKappa <- function(dataset, options){
  
  if(options[["cohensWeightedOrNot"]] == "cohensWeighted"){
    weightedString <- "Weighted"
  }else{
    weightedString <- "Unweighted"
  }
  
  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Cohen's %s Kappa", casefold(weightedString, upper = F)))
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "cKappa", title = gettextf("%s Kappa", weightedString), type = "number")
  
  
  formattedCIPercent <- format(
    100 * options[["kappaConfidenceIntervalValue"]],
    digits = 3,
    drop0trailing = TRUE
  )
  
  
  if (options[["kappaIntervalOn"]]) {
    jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "integer", overtitle = gettext("Confidence interval of 95%"))
    jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "integer", overtitle = gettext("Confidence interval of 95%"))
  }
  
  
  
  
  jaspTable$setData(list("ratings" = c("contNormal - contGamma", "contNormal - V1", "contGamma - V1"),
                         "cKappa" = c(1,1,1),
                         "CIL"    = c(0.9, 0.9,0.9),
                         "CIU"    = c(1.1,1.1,1.1)))
  
  
  return(jaspTable)
  
  
}



.fleissKappa <- function(dataset, options){
  
  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Fleiss' Kappa"))
  
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "fKappa", title = gettext("Fleiss' Kappa"), type = "number")

  
  formattedCIPercent <- format(
    100 * options[["kappaConfidenceIntervalValue"]],
    digits = 3,
    drop0trailing = TRUE
  )
  
  
  if (options[["kappaIntervalOn"]]) {
    jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "integer", overtitle = gettext("Confidence interval of 95%"))
    jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "integer", overtitle = gettext("Confidence interval of 95%"))
  }
  
  
  jaspTable$setData(list("ratings" = c("Overall", "Item 1", "Item 2"),
                         "fKappa" = c(1,1,1),
                         "CIL"    = c(0.9, 0.9,0.9),
                         "CIU"    = c(1.1,1.1,1.1)))

  
  
  return(jaspTable)
  
  
}
