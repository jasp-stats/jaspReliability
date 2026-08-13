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

#' Reliability
#'
#' @param averageInterItemCorrelation, Mean of all pairwise Pearson correlations between items.
#'    Defaults to \code{FALSE}.
#' @param bootstrapSamples, Number of bootstrap replications. Higher values yield more stable interval estimates.
#' @param bootstrapType, Non-parametric bootstrap resamples the data; parametric bootstrap samples from a multivariate normal with the estimated parameters.
#' \itemize{
#'   \item \code{"parametric"}
#'   \item \code{"nonParametric"}
#' }
#' @param ciLevel, Coverage of the confidence intervals for scale reliability statistics.
#' @param intervalMethod, Analytic intervals use normal-theory standard errors (van der Ark, 2024). Bootstrapped intervals use percentile resampling.
#' \itemize{
#'   \item \code{"bootstrapped"}
#'   \item \code{"analytic"}
#' }
#' @param intervalMethodVar, Chi-square-based intervals assume normality; non-parametric intervals use a two-step bootstrap procedure (van der Ark, 2024).
#' \itemize{
#'   \item \code{"twostep"}
#'   \item \code{"chisq"}
#' }
#' @param itemCiLevel, Coverage of the confidence intervals for item-level statistics.
#' @param itemDeletedAlpha, Alpha of the remaining items when this item is removed from the scale.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedLambda2, Lambda 2 of the remaining items when this item is removed from the scale.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedOmega, Omega of the remaining items when this item is removed from the scale.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedSplithalf, Split-half coefficient of the remaining items when this item is removed from the scale.
#'    Defaults to \code{FALSE}.
#' @param itemRestCorrelation, Correlation of each item with the sum of the remaining items.
#'    Defaults to \code{FALSE}.
#' @param meanSdScoresMethod, Whether the mean, variance, and SD are based on sum scores or mean scores across items.
#' \itemize{
#'   \item \code{"meanScores"}
#'   \item \code{"sumScores"}
#' }
#' @param naAction, Pairwise uses all available observations per covariance pair; listwise deletes any row with a missing value.
#' \itemize{
#'   \item \code{"pairwise"}
#'   \item \code{"listwise"}
#' }
#' @param omegaEstimationMethod, CFA fits the single-factor model via confirmatory factor analysis; PFA uses principal factor analysis.
#' \itemize{
#'   \item \code{"cfa"}
#'   \item \code{"pfa"}
#' }
#' @param omegaFitMeasures, Chi-square test, RMSEA, SRMR and other fit indices for the single-factor model underlying omega.
#'    Defaults to \code{FALSE}.
#' @param samplesSavingDisabled, When checked, bootstrap samples are not stored in the output file. This reduces file size but may slow down re-running the analysis, because samples are precomputed and cached for speed.
#'    Defaults to \code{FALSE}.
#' @param scaleAlpha, Cronbach's alpha. For binary items this equals KR-20.
#'    Defaults to \code{FALSE}.
#' @param scaleLambda2, Guttman's lambda 2, a lower bound for reliability.
#'    Defaults to \code{FALSE}.
#' @param scaleOmega, McDonald's omega for unidimensional data based on the single-factor model. The denominator uses model-implied total variance.
#'    Defaults to \code{TRUE}.
#' @param scaleSplithalf, Splits items into two halves (odd/even by default). Unstandardized: Flanagan-Rulon coefficient; Standardized: Spearman-Brown coefficient.
#'    Defaults to \code{FALSE}.
#' @param standardizedLoadings, Table of standardized loadings from the single-factor model.
#'    Defaults to \code{FALSE}.
#' @param variables, Items/variables to include in the reliability analysis. Must be scale variables.
#' @export
reliabilityUnidimensionalFrequentist <- function(
          data = NULL,
          version = "0.97.1",
          averageInterItemCorrelation = FALSE,
          bootstrapSamples = 1000,
          bootstrapType = "nonParametric",
          ciLevel = 0.95,
          coefficientType = "unstandardized",
          hiddenScaleThreshold = 0,
          intervalMethod = "analytic",
          intervalMethodVar = "chisq",
          itemCiLevel = 0.95,
          itemDeletedAlpha = FALSE,
          itemDeletedLambda2 = FALSE,
          itemDeletedOmega = FALSE,
          itemDeletedSplithalf = FALSE,
          itemMean = FALSE,
          itemRestCorrelation = FALSE,
          itemSd = FALSE,
          itemVar = FALSE,
          meanSdScoresMethod = "sumScores",
          naAction = "pairwise",
          omegaEstimationMethod = "cfa",
          omegaFitMeasures = FALSE,
          plotHeight = 320,
          plotWidth = 480,
          reverseScaledItems = list(types = list(), value = list()),
          samplesSavingDisabled = FALSE,
          scaleAlpha = FALSE,
          scaleLambda2 = FALSE,
          scaleMean = FALSE,
          scaleOmega = TRUE,
          scaleSd = FALSE,
          scaleSplithalf = FALSE,
          scaleVar = FALSE,
          seed = 1234,
          setSeed = FALSE,
          standardizedLoadings = FALSE,
          variables = list(types = list(), value = list())) {

   defaultArgCalls <- formals(jaspReliability::reliabilityUnidimensionalFrequentist)
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

   optionsWithFormula <- c("reverseScaledItems", "variables")
   for (name in optionsWithFormula) {
      if ((name %in% optionsWithFormula) && inherits(options[[name]], "formula")) options[[name]] = jaspBase::jaspFormula(options[[name]], data)   }

   return(jaspBase::runWrappedAnalysis("jaspReliability", "reliabilityUnidimensionalFrequentist", "ReliabilityUnidimensionalFrequentist.qml", options, version, TRUE))
}
