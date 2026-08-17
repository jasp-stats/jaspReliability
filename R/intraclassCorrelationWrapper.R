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

#' Intraclass Correlation
#'
#' @param averagedRating, When ratings are averaged across raters, the k-rater version of the selected ICC is computed (e.g. ICC(2,k) instead of ICC(2,1)) per Shrout & Fleiss (1979). Averaging generally raises the reliability estimate.
#'    Defaults to \code{FALSE}.
#' @param ci, Report a confidence interval for the ICC estimate based on the F-distribution.
#'    Defaults to \code{TRUE}.
#' @param iccType, Determines the ICC model. ICC(1): each subject rated by a different random rater. ICC(2): all subjects rated by the same random sample of raters. ICC(3): all subjects rated by the same fixed raters. See Shrout & Fleiss (1979).
#' \itemize{
#'   \item \code{"icc3"}
#'   \item \code{"icc1"}
#'   \item \code{"icc2"}
#' }
#' @param variables, Rating variables to include. Each variable is one rater; each row is a subject being rated.
#' @export
intraclassCorrelation <- function(
          data = NULL,
          version = "0.97.1",
          averagedRating = FALSE,
          ci = TRUE,
          ciLevel = 0.95,
          iccType = "icc1",
          plotHeight = 320,
          plotWidth = 480,
          variables = list(types = list(), value = list())) {

   defaultArgCalls <- formals(jaspReliability::intraclassCorrelation)
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

   optionsWithFormula <- c("variables")
   for (name in optionsWithFormula) {
      if ((name %in% optionsWithFormula) && inherits(options[[name]], "formula")) options[[name]] = jaspBase::jaspFormula(options[[name]], data)   }

   return(jaspBase::runWrappedAnalysis("jaspReliability", "intraclassCorrelation", "IntraclassCorrelation.qml", options, version, TRUE))
}
