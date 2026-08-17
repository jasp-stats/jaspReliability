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

#' Rater Agreement
#'
#' @param bootstrapSamples, Number of bootstrap replications used to compute confidence intervals. Higher values give more stable estimates. Only applies to Krippendorff's alpha and Kendall's W, whose CIs are bootstrap-based.
#' @param ci, Report a confidence interval for each agreement coefficient. Cohen's kappa and Fleiss' kappa use asymptotic CIs; Krippendorff's alpha and Kendall's W use bootstrap CIs.
#'    Defaults to \code{TRUE}.
#' @param ciLevel, Width of the confidence interval.
#' @param cohensKappa, Measures agreement between exactly two raters. When more than two raters are entered, all pairwise combinations are computed.
#'    Defaults to \code{FALSE}.
#' @param cohensKappaType, Unweighted kappa treats all disagreements equally. Weighted kappa accounts for the degree of disagreement and requires ordinal ratings.
#' \itemize{
#'   \item \code{"unweighted"}: All disagreements are treated as equal regardless of their magnitude.
#'   \item \code{"weighted"}: Disagreements are penalised according to their magnitude. Requires ordinal ratings.
#' }
#' @param correctForTies, Apply a correction to Kendall's W when tied ranks are present in the data. Recommended whenever ties can occur.
#'    Defaults to \code{TRUE}.
#' @param dataStructure, Specify whether raters are arranged in columns (default) or rows in the dataset.
#' \itemize{
#'   \item \code{"ratersInColumns"}: Each column is one rater and each row is one subject or item being rated.
#'   \item \code{"ratersInRows"}: Each row is one rater and each column is one subject or item being rated.
#' }
#' @param fleissKappa, Measures agreement among two or more raters on nominal categories. Generalises Cohen's kappa to multiple raters.
#'    Defaults to \code{FALSE}.
#' @param kendallW, Measures concordance among multiple raters on ordinal rankings. Ranges from 0 (no agreement) to 1 (perfect agreement).
#'    Defaults to \code{FALSE}.
#' @param krippendorffsAlpha, Measures agreement among two or more raters. Applicable to nominal, ordinal, interval, or ratio data.
#'    Defaults to \code{FALSE}.
#' @param krippendorffsAlphaMethod, Level of measurement determines how disagreements are quantified in the alpha calculation.
#' @param variables, Rating variables to include. Whether a variable represents a rater or a subject/item depends on the data structure setting.
#' @param weightType, Weighting scheme applied to disagreements between ordinal categories.
#' \itemize{
#'   \item \code{"quadratic"}: Penalises larger disagreements quadratically; sensitive to large discrepancies.
#'   \item \code{"linear"}: Penalises disagreements proportionally to their size.
#' }
#' @export
raterAgreement <- function(
          data = NULL,
          version = "0.97.1",
          bootstrapSamples = 1000,
          ci = TRUE,
          ciLevel = 0.95,
          cohensKappa = FALSE,
          cohensKappaType = "unweighted",
          correctForTies = TRUE,
          dataStructure = "ratersInColumns",
          fleissKappa = FALSE,
          kendallW = FALSE,
          krippendorffsAlpha = FALSE,
          krippendorffsAlphaMethod = "nominal",
          plotHeight = 320,
          plotWidth = 480,
          seed = 1,
          setSeed = FALSE,
          variables = list(types = list(), value = list()),
          weightType = "quadratic") {

   defaultArgCalls <- formals(jaspReliability::raterAgreement)
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

   optionsWithFormula <- c("krippendorffsAlphaMethod", "variables")
   for (name in optionsWithFormula) {
      if ((name %in% optionsWithFormula) && inherits(options[[name]], "formula")) options[[name]] = jaspBase::jaspFormula(options[[name]], data)   }

   return(jaspBase::runWrappedAnalysis("jaspReliability", "raterAgreement", "RaterAgreement.qml", options, version, TRUE))
}
