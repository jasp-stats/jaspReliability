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

#' Bayesian Reliability
#'
#' @param averageInterItemCorrelation, Mean of all pairwise Pearson correlations between items.
#'    Defaults to \code{FALSE}.
#' @param burnin, Initial samples discarded while the chain converges to the posterior.
#' @param chains, Number of independent MCMC chains. Multiple chains enable R-hat convergence diagnostics.
#' @param coefficientType, Unstandardized uses the covariance matrix; standardized uses the correlation matrix. Standardized alpha is contested in the literature.
#' \itemize{
#'   \item \code{"unstandardized"}
#'   \item \code{"standardized"}
#' }
#' @param effectiveSampleSize, Effective sample size: number of independent posterior samples after accounting for autocorrelation.
#'    Defaults to \code{FALSE}.
#' @param inverseGammaPriorScale, Scale parameter of the inverse gamma prior on residual variances.
#' @param inverseGammaPriorShape, Shape parameter of the inverse gamma prior on residual variances.
#' @param inverseWishartPriorDf, Degrees of freedom of the inverse Wishart prior. Minimum equals the number of items.
#' @param inverseWishartPriorScale, Precision values on the diagonal of the inverse Wishart scaling matrix.
#' @param itemCiLevel, Width of the credible interval for item-level statistics.
#' @param itemDeletedAlpha, Posterior distribution of alpha for the remaining items when this item is removed.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedLambda2, Posterior distribution of lambda 2 for the remaining items when this item is removed.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedOmega, Posterior distribution of omega for the remaining items when this item is removed.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedPlot, Displays posterior densities of the reliability of the remaining items when each item is dropped.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedPlotOrdered, Sort the densities by how much removing an item changes the posterior (by mean, KL-divergence, or KS-distance).
#'    Defaults to \code{FALSE}.
#' @param itemDeletedSplithalf, Posterior distribution of the split-half coefficient for the remaining items when this item is removed.
#'    Defaults to \code{FALSE}.
#' @param itemRestCorrelation, Correlation of each item with the sum of the remaining items.
#'    Defaults to \code{FALSE}.
#' @param meanSdScoresMethod, Whether the mean, variance, and SD are based on sum scores or mean scores across items.
#' \itemize{
#'   \item \code{"meanScores"}
#'   \item \code{"sumScores"}
#' }
#' @param naAction, Bayesian imputation treats missing values as unknown parameters sampled from the posterior; listwise deletion removes any row with a missing value.
#' \itemize{
#'   \item \code{"listwise"}
#'   \item \code{"imputation"}
#' }
#' @param normalPriorMean, Mean of the normal prior on factor loadings.
#' @param omegaFitMeasures, Bayesian fit indices (B-LR, B-RMSEA, B-CFI, B-TLI) for the single-factor model with probability statements relative to cutoffs.
#'    Defaults to \code{FALSE}.
#' @param omegaPosteriorPredictiveCheck, Graphical check for single-factor model fit: eigenvalues of the observed covariance matrix are compared against the model-implied posterior distribution.
#'    Defaults to \code{FALSE}.
#' @param pointEstimate, Whether to report the posterior mean or median as the point estimate in the tables.
#' \itemize{
#'   \item \code{"mean"}
#'   \item \code{"median"}
#' }
#' @param posteriorPlot, Display posterior density plots for the selected reliability coefficients.
#'    Defaults to \code{FALSE}.
#' @param posteriorPlotFixedRange, Fix the x-axis of the posterior plots to [0, 1] for easier comparison between coefficients.
#'    Defaults to \code{FALSE}.
#' @param posteriorPlotPriorDisplayed, Add the prior distribution to the posterior density plot.
#'    Defaults to \code{FALSE}.
#' @param posteriorPlotShaded, Shade the posterior region corresponding to the probability interval in the density plot.
#'    Defaults to \code{FALSE}.
#' @param probabilityTable, Report the prior and posterior probability that a reliability coefficient falls within the specified interval.
#'    Defaults to \code{FALSE}.
#' @param rHat, Potential scale reduction factor. Values close to 1 (< 1.1) indicate convergence across chains.
#'    Defaults to \code{FALSE}.
#' @param samples, Total number of MCMC iterations per chain, including burn-in. The number of retained posterior samples per chain is ceiling((samples - burn-in) / thinning).
#' @param samplesSavingDisabled, When checked, MCMC samples are not stored in the output file. Reduces file size but may slow down re-running the analysis, because samples are precomputed and cached for speed.
#'    Defaults to \code{FALSE}.
#' @param scaleAlpha, Cronbach's alpha. For binary items this equals KR-20.
#'    Defaults to \code{FALSE}.
#' @param scaleCiLevel, Width of the credible interval for scale reliability statistics.
#' @param scaleLambda2, Guttman's lambda 2, a lower bound for reliability.
#'    Defaults to \code{FALSE}.
#' @param scaleOmega, McDonald's omega for unidimensional data based on the single-factor model.
#'    Defaults to \code{TRUE}.
#' @param scaleSplithalf, Splits items into two halves (odd/even by default). Unstandardized: Flanagan-Rulon; Standardized: Spearman-Brown.
#'    Defaults to \code{FALSE}.
#' @param standardizedLoadings, Table of standardized loadings from the single-factor model (posterior mean or median).
#'    Defaults to \code{FALSE}.
#' @param thinning, Keep every k-th sample to reduce autocorrelation. A value of 1 keeps all samples.
#' @param tracePlot, Plot of sampled values per chain over iterations. Well-mixed chains indicate convergence.
#'    Defaults to \code{FALSE}.
#' @param variables, Items/variables to include in the reliability analysis. Must be scale variables.
#' @export
reliabilityUnidimensionalBayesian <- function(
          data = NULL,
          version = "0.97.1",
          averageInterItemCorrelation = FALSE,
          burnin = 50,
          chains = 3,
          coefficientType = "unstandardized",
          effectiveSampleSize = FALSE,
          inverseGammaPriorScale = 1,
          inverseGammaPriorShape = 2,
          inverseWishartPriorDf = 0,
          inverseWishartPriorScale = 1e-10,
          itemCiLevel = 0.95,
          itemDeletedAlpha = FALSE,
          itemDeletedLambda2 = FALSE,
          itemDeletedOmega = FALSE,
          itemDeletedPlot = FALSE,
          itemDeletedPlotOrdered = FALSE,
          itemDeletedPlotOrderedType = "mean",
          itemDeletedSplithalf = FALSE,
          itemMean = FALSE,
          itemRestCorrelation = FALSE,
          itemSd = FALSE,
          itemVar = FALSE,
          meanSdScoresMethod = "sumScores",
          naAction = "imputation",
          normalPriorMean = 0,
          omegaFitMeasures = FALSE,
          omegaFitMeasuresCiLevel = 0.9,
          omegaFitMeasuresCutoffCfiTli = 0.9,
          omegaFitMeasuresCutoffRmsea = 0.08,
          omegaPosteriorPredictiveCheck = FALSE,
          plotHeight = 320,
          plotWidth = 480,
          pointEstimate = "mean",
          posteriorPlot = FALSE,
          posteriorPlotFixedRange = FALSE,
          posteriorPlotPriorDisplayed = FALSE,
          posteriorPlotShaded = FALSE,
          probabilityTable = FALSE,
          probabilityTableLowerBound = 0.7,
          probabilityTableUpperBound = 0.9,
          rHat = FALSE,
          reverseScaledItems = list(types = list(), value = list()),
          samples = 1000,
          samplesSavingDisabled = FALSE,
          scaleAlpha = FALSE,
          scaleCiLevel = 0.95,
          scaleLambda2 = FALSE,
          scaleMean = FALSE,
          scaleOmega = TRUE,
          scaleSd = FALSE,
          scaleSplithalf = FALSE,
          scaleVar = FALSE,
          seed = 1234,
          setSeed = FALSE,
          standardizedLoadings = FALSE,
          thinning = 1,
          tracePlot = FALSE,
          variables = list(types = list(), value = list())) {

   defaultArgCalls <- formals(jaspReliability::reliabilityUnidimensionalBayesian)
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

   return(jaspBase::runWrappedAnalysis("jaspReliability", "reliabilityUnidimensionalBayesian", "ReliabilityUnidimensionalBayesian.qml", options, version, TRUE))
}
