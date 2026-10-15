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

#' Multidimensional Reliability
#'
#' @param bootstrapSamples, Number of bootstrap replications. The factor model is refit for every replication, so higher values take proportionally longer.
#' @param ciLevel, Coverage of the confidence interval for the reliability coefficients.
#' @param factors, Assign the items to their factors. The analysis requires at least two factors with at least two items each. An item may be assigned to more than one factor (cross-loading); the bi-factor model does not support cross-loadings.
#' @param fitMeasures, Fit measures of the factor model: chi-square, CFI, TLI, RMSEA, SRMR, AIC and BIC.
#'    Defaults to \code{FALSE}.
#' @param intervalMethod, Analytic intervals are the Wald-type intervals of the factor model. Because omega is a ratio, these intervals are symmetric and can extend beyond zero or one; the bootstrapped interval resamples the data and refits the model, which is slower but does not make that assumption.
#' \itemize{
#'   \item \code{"analytic"}
#'   \item \code{"bootstrapped"}
#' }
#' @param itemDeletedOmegaH, McDonald's omega_h for the remaining items when this item is dropped. Not available for the correlated-factors model.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedOmegaT, McDonald's omega_t for the remaining items when this item is dropped. The model is refit once per item.
#'    Defaults to \code{FALSE}.
#' @param itemRestCorrelation, Pearson correlation between each item and the sum of the remaining items.
#'    Defaults to \code{FALSE}.
#' @param meanSdScoresMethod, Whether the mean and standard deviation in the scale table are based on participants' sum scores or mean scores across items.
#' \itemize{
#'   \item \code{"meanScores"}
#'   \item \code{"sumScores"}
#' }
#' @param modelType, The factor model used to estimate the reliability coefficients. McDonald's omega_h (general/group-common reliability) is only available for the second-order and bi-factor models.
#' @param naAction, Full information maximum likelihood uses all available observations; listwise deletion removes any row with a missing value.
#' \itemize{
#'   \item \code{"fiml"}
#'   \item \code{"listwise"}
#' }
#' @param reverseScaledItems, Items assigned here are recoded (reverse-scored) before the analysis.
#' @param samplesSavingDisabled, When checked, bootstrap samples are not stored in the output file. Reduces file size but changing the confidence level requires resampling.
#'    Defaults to \code{FALSE}.
#' @param setSeed, Fix the random number generator seed to make the bootstrap results reproducible.
#'    Defaults to \code{FALSE}.
#' @param standardizedLoadings, Table of standardized loadings of the items on their group factors. For the bi-factor model the loadings on the general factor are added as a column; for the second-order model the loadings of the group factors on the general factor are shown in a separate table.
#'    Defaults to \code{FALSE}.
#' @export
reliabilityMultidimensionalFrequentist <- function(
          data = NULL,
          version = "0.98.1",
          bootstrapSamples = 1000,
          ciLevel = 0.95,
          factors = list(list(indicators = list(), name = "Factor1", title = "Factor 1"), list(indicators = list(), name = "Factor2", title = "Factor 2")),
          fitMeasures = FALSE,
          intervalMethod = "analytic",
          itemDeletedOmegaH = FALSE,
          itemDeletedOmegaT = FALSE,
          itemRestCorrelation = FALSE,
          meanSdScoresMethod = "sumScores",
          modelType = "secondOrder",
          naAction = "fiml",
          plotHeight = 320,
          plotWidth = 480,
          reverseScaledItems = list(types = list(), value = list()),
          samplesSavingDisabled = FALSE,
          seed = 1234,
          setSeed = FALSE,
          standardizedLoadings = FALSE) {

   defaultArgCalls <- formals(jaspReliability::reliabilityMultidimensionalFrequentist)
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

   optionsWithFormula <- c("factors", "modelType", "reverseScaledItems")
   for (name in optionsWithFormula) {
      if ((name %in% optionsWithFormula) && inherits(options[[name]], "formula")) options[[name]] = jaspBase::jaspFormula(options[[name]], data)   }

   return(jaspBase::runWrappedAnalysis("jaspReliability", "reliabilityMultidimensionalFrequentist", "ReliabilityMultidimensionalFrequentist.qml", options, version, TRUE))
}
