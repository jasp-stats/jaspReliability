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

#' Bayesian Multidimensional Reliability
#'
#' @param burnin, Initial samples discarded while the chain converges to the posterior.
#' @param chains, Number of independent MCMC chains. Multiple chains enable R-hat convergence diagnostics.
#' @param ciLevel, Width of the credible interval for the reliability coefficients.
#' @param factors, Assign the items to their factors. The analysis requires at least two factors with at least two items each. An item may be assigned to more than one factor (cross-loading); the bi-factor model does not support cross-loadings.
#' @param fitMeasures, Bayesian fit measures for the factor model: B-LR, B-SRMR, and B-RMSEA with a posterior probability relative to the cutoff.
#'    Defaults to \code{FALSE}.
#' @param fitMeasuresCiLevel, Width of the credible interval for the fit measures.
#' @param fitMeasuresCutoffRmsea, Cutoff for the posterior probability that the B-RMSEA is below this value.
#' @param igScaleGFactor, Scale parameter of the inverse gamma prior on the general factor variance.
#' @param igScaleLatent, Scale parameter of the inverse gamma prior on the latent residual variances.
#' @param igScaleManifest, Scale parameter of the inverse gamma prior on the item residual variances.
#' @param igShapeGFactor, Shape parameter of the inverse gamma prior on the general factor variance.
#' @param igShapeLatent, Shape parameter of the inverse gamma prior on the latent residual variances.
#' @param igShapeManifest, Shape parameter of the inverse gamma prior on the item residual variances.
#' @param itemDeletedOmegaH, Posterior point estimate of ω_h for the remaining items when this item is dropped. Not available for the correlated-factors model.
#'    Defaults to \code{FALSE}.
#' @param itemDeletedOmegaT, Posterior point estimate of ω_t for the remaining items when this item is dropped. The model is refit once per item.
#'    Defaults to \code{FALSE}.
#' @param itemRestCorrelation, Pearson correlation between each item and the sum of the remaining items.
#'    Defaults to \code{FALSE}.
#' @param latentCorDf, Degrees of freedom of the inverse Wishart prior on the latent correlation matrix. Values below the number of factors are raised to the number of factors to keep the prior proper.
#' @param loadMeanLatent, Mean of the normal prior on the structural loadings of the group factors on the general factor.
#' @param loadMeanManifest, Mean of the normal prior on the item factor loadings.
#' @param loadScaleLatent, Scales the variance of the normal prior on the structural loadings of the group factors on the general factor. Larger values mean a less informative prior.
#' @param loadScaleManifest, Scales the variance of the normal prior on the item factor loadings. The prior is conditional on the item residual variance, so a loading has prior standard deviation sqrt(residual variance x variance scaling); larger values mean a less informative prior.
#' @param meanSdScoresMethod, Whether the mean and standard deviation in the scale table are based on participants' sum scores or mean scores across items.
#' \itemize{
#'   \item \code{"meanScores"}
#'   \item \code{"sumScores"}
#' }
#' @param modelType, The factor model used to estimate the reliability coefficients. McDonald's ω_h (general/group-common reliability) is only available for the second-order and bi-factor models.
#' @param naAction, Bayesian imputation treats missing values as unknown parameters sampled from the posterior; listwise deletion removes any row with a missing value.
#' \itemize{
#'   \item \code{"listwise"}
#'   \item \code{"imputation"}
#' }
#' @param pointEstimate, Whether to report the posterior mean or median as the point estimate in the tables.
#' \itemize{
#'   \item \code{"median"}
#'   \item \code{"mean"}
#' }
#' @param posteriorPlot, Display posterior density plots for the reliability coefficients.
#'    Defaults to \code{FALSE}.
#' @param posteriorPlotFixedRange, Fix the x-axis of the posterior plots to [0, 1] for easier comparison between coefficients.
#'    Defaults to \code{FALSE}.
#' @param posteriorPlotPriorDisplayed, Add the prior distribution to the posterior density plot.
#'    Defaults to \code{FALSE}.
#' @param posteriorPlotShaded, Shade the posterior region corresponding to the probability interval in the density plot.
#'    Defaults to \code{FALSE}.
#' @param posteriorPredictiveCheck, Graphical check of factor model fit: eigenvalues of the observed covariance matrix are compared against the model-implied posterior distribution.
#'    Defaults to \code{FALSE}.
#' @param probabilityTable, Report the posterior probability that a reliability coefficient falls within the specified interval.
#'    Defaults to \code{FALSE}.
#' @param probabilityTableLowerBound, Lower bound of the reliability interval.
#' @param probabilityTableUpperBound, Upper bound of the reliability interval.
#' @param rHat, Potential scale reduction factor. Values close to 1 (< 1.1) indicate convergence across chains.
#'    Defaults to \code{FALSE}.
#' @param reverseScaledItems, Items assigned here are recoded (reverse-scored) before the analysis.
#' @param samples, Total number of MCMC samples per chain, including burn-in.
#' @param samplesSavingDisabled, When checked, MCMC samples are not stored in the output file. Reduces file size but re-running the analysis or changing options requires resampling.
#'    Defaults to \code{FALSE}.
#' @param setSeed, Fix the random number generator seed to make the MCMC results reproducible.
#'    Defaults to \code{FALSE}.
#' @param thinning, Keep every k-th sample to reduce autocorrelation. A value of 1 keeps all samples.
#' @param tracePlot, Plot of sampled values per chain over iterations. Well-mixed chains indicate convergence.
#'    Defaults to \code{FALSE}.
#' @export
reliabilityMultidimensionalBayesian <- function(
          data = NULL,
          version = "0.97.1",
          burnin = 200,
          chains = 3,
          ciLevel = 0.95,
          factors = list(list(indicators = list(), name = "Factor1", title = "Factor 1"), list(indicators = list(), name = "Factor2", title = "Factor 2")),
          fitMeasures = FALSE,
          fitMeasuresCiLevel = 0.9,
          fitMeasuresCutoffRmsea = 0.08,
          igScaleGFactor = 1,
          igScaleLatent = 1,
          igScaleManifest = 1,
          igShapeGFactor = 2,
          igShapeLatent = 2,
          igShapeManifest = 2,
          itemDeletedOmegaH = FALSE,
          itemDeletedOmegaT = FALSE,
          itemRestCorrelation = FALSE,
          latentCorDf = 2,
          loadMeanLatent = 0,
          loadMeanManifest = 0,
          loadScaleLatent = 2.5,
          loadScaleManifest = 1,
          meanSdScoresMethod = "sumScores",
          modelType = "secondOrder",
          naAction = "imputation",
          plotHeight = 320,
          plotWidth = 480,
          pointEstimate = "mean",
          posteriorPlot = FALSE,
          posteriorPlotFixedRange = FALSE,
          posteriorPlotPriorDisplayed = FALSE,
          posteriorPlotShaded = FALSE,
          posteriorPredictiveCheck = FALSE,
          probabilityTable = FALSE,
          probabilityTableLowerBound = 0.7,
          probabilityTableUpperBound = 0.9,
          rHat = FALSE,
          reverseScaledItems = list(types = list(), value = list()),
          samples = 2000,
          samplesSavingDisabled = FALSE,
          seed = 1234,
          setSeed = FALSE,
          thinning = 1,
          tracePlot = FALSE) {

   defaultArgCalls <- formals(jaspReliability::reliabilityMultidimensionalBayesian)
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

   return(jaspBase::runWrappedAnalysis("jaspReliability", "reliabilityMultidimensionalBayesian", "ReliabilityMultidimensionalBayesian.qml", options, version, TRUE))
}
