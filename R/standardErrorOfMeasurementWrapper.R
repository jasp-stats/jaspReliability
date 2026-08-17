#
# Copyright (C) 2013-2025 University of Amsterdam
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

# This is a generated file. Don't change it!

#' Standard Error of Measurement
#'
#' @param anova, Estimates conditional SEM from a repeated-measures ANOVA using the ICC(3,k) approach (Emons, 2023).
#'    Defaults to \code{FALSE}.
#' @param combinedPointPlot, Combined plot displaying all methods' SEM values in a single panel for easy comparison.
#'    Defaults to \code{FALSE}.
#' @param feldt, Splits items into multiple equal parts to estimate error variance per score group.
#'    Defaults to \code{FALSE}.
#' @param feldtNumberOfSplits, How many equal parts to split the test into. Must be a divisor of the number of items.
#' @param hideTable, Hide the main SEM table showing conditional error estimates per sum score and method.
#'    Defaults to \code{FALSE}.
#' @param histogramCounts, Histogram of the number of respondents per sum score group.
#'    Defaults to \code{FALSE}.
#' @param irt, IRT-based SEM using the 2PLM for dichotomous items or the graded response model for polytomous items. Assumes a single latent variable.
#'    Defaults to \code{FALSE}.
#' @param keats, Keats' correction of Lord's method; uses a reliability coefficient to reduce bias.
#'    Defaults to \code{FALSE}.
#' @param lord, Lord's binomial method based on the number of correct and incorrect responses.
#'    Defaults to \code{FALSE}.
#' @param lord2, Lord's method generalised to multiple test parts.
#'    Defaults to \code{FALSE}.
#' @param lord2NumberOfSplits, How many parts to split the test into. Must be a divisor of the number of items (excluding the last value).
#' @param minimumGroupSize, Score groups with fewer observations are merged with adjacent groups before estimating the SEM (default: 20). Applies to the Thorndike, Feldt, ANOVA, and Lord generalized methods.
#' @param mollenkopfFeldt, Splits items into parts and models score differences with polynomial regression to estimate the conditional SEM.
#'    Defaults to \code{FALSE}.
#' @param mollenkopfFeldtNumberOfSplits, How many equal parts to split the test into. Must be a divisor of the number of items.
#' @param mollenkopfFeldtPolyDegree, Degree of the polynomial used to predict score differences, e.g., 3 fits Y = X + X² + X³.
#' @param pointPlots, Separate point plot per method showing the conditional SEM values across sum scores.
#'    Defaults to \code{FALSE}.
#' @param sumScoreCiPlots, Plot sum scores with confidence interval bands per method. Intervals are normal-theory CIs using the estimated SEM.
#'    Defaults to \code{FALSE}.
#' @param sumScoreCiPlotsCutoff, Add a horizontal reference line at the specified cut score.
#'    Defaults to \code{FALSE}.
#' @param sumScoreCiTable, Table of sum scores with normal-theory confidence intervals computed from the estimated SEM.
#'    Defaults to \code{FALSE}.
#' @param thorndike, Splits items into odd- and even-numbered halves. Reorder variables in the list to change the split.
#'    Defaults to \code{FALSE}.
#' @param userReliability, Override the reliability estimate used for the unconditional SEM and the Keats method.
#'    Defaults to \code{FALSE}.
#' @param variables, Items to include in the SEM analysis. Must be ordinally or nominally scaled (dichotomous or polytomous).
#' @export
standardErrorOfMeasurement <- function(
          data = NULL,
          version = "0.97.1",
          anova = FALSE,
          ciLevelPlots = 0.95,
          ciLevelTable = 0.95,
          combinedPointPlot = FALSE,
          feldt = FALSE,
          feldtNumberOfSplits = "2",
          hideTable = FALSE,
          histogramCounts = FALSE,
          irt = FALSE,
          keats = FALSE,
          lord = FALSE,
          lord2 = FALSE,
          lord2NumberOfSplits = "2",
          minimumGroupSize = 20,
          mollenkopfFeldt = FALSE,
          mollenkopfFeldtNumberOfSplits = "2",
          mollenkopfFeldtPolyDegree = 2,
          plotHeight = 320,
          plotWidth = 480,
          pointPlots = FALSE,
          reliabilityValue = 0.5,
          sumScoreCiPlots = FALSE,
          sumScoreCiPlotsCutoff = FALSE,
          sumScoreCiPlotsCutoffValue = 0,
          sumScoreCiTable = FALSE,
          thorndike = FALSE,
          userReliability = FALSE,
          variables = list(types = list(), value = list())) {

   defaultArgCalls <- formals(jaspReliability::standardErrorOfMeasurement)
   defaultArgs <- lapply(defaultArgCalls, eval)
   options <- as.list(match.call())[-1L]
   options <- lapply(options, eval)
   defaults <- setdiff(names(defaultArgs), names(options))
   options[defaults] <- defaultArgs[defaults]
   options[["data"]] <- NULL
   options[["version"]] <- NULL


   if (!jaspBase::jaspResultsCalledFromJasp() && !is.null(data)) {
      jaspBase::storeDataSet(data)
   }

   optionsWithFormula <- c("feldtNumberOfSplits", "lord2NumberOfSplits", "mollenkopfFeldtNumberOfSplits", "variables")
   for (name in optionsWithFormula) {
      if ((name %in% optionsWithFormula) && inherits(options[[name]], "formula")) options[[name]] = jaspBase::jaspFormula(options[[name]], data)   }

   return(jaspBase::runWrappedAnalysis("jaspReliability", "standardErrorOfMeasurement", "StandardErrorOfMeasurement.qml", options, version, TRUE))
}
