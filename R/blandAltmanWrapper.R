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

#' Bland-Altman Plots
#'
#' @param blandAltmanTable, Table showing the mean difference, the upper and lower limits of agreement, and their confidence intervals.
#'    Defaults to \code{FALSE}.
#' @param ci, Display confidence intervals around the mean difference and its upper and lower limits of agreement in the plot.
#'    Defaults to \code{FALSE}.
#' @param ciShading, Shade the confidence interval regions around the mean difference and limits of agreement.
#'    Defaults to \code{FALSE}.
#' @param ciShadingWithColour, Use colour (rather than grey) to fill the shaded confidence regions.
#'    Defaults to \code{FALSE}.
#' @param pairs, Pairs of scale variables to compare. Each pair produces a Bland-Altman plot of the mean versus the difference between measurements.
#' @export
blandAltman <- function(
          data = NULL,
          version = "0.97.1",
          blandAltmanTable = FALSE,
          ci = FALSE,
          ciLevel = 0.95,
          ciShading = FALSE,
          ciShadingWithColour = FALSE,
          pairs = list(),
          plotHeight = 320,
          plotWidth = 480) {

   defaultArgCalls <- formals(jaspReliability::blandAltman)
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

   optionsWithFormula <- c("pairs")
   for (name in optionsWithFormula) {
      if ((name %in% optionsWithFormula) && inherits(options[[name]], "formula")) options[[name]] = jaspBase::jaspFormula(options[[name]], data)   }

   return(jaspBase::runWrappedAnalysis("jaspReliability", "blandAltman", "BlandAltman.qml", options, version, TRUE))
}
