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

#' Bayesian Rater Agreement
#'
#' Bayesian Rater Agreement estimates the chance-corrected agreement between raters who assign subjects to nominal or ordinal categories. For every rater pair, the cell probabilities of the pair's agreement table get a Dirichlet posterior, and Cohen's kappa or Fleiss' kappa is computed from every posterior draw (Calle-Alonso & Pérez Sánchez, 2015; Pfadt et al., 2026).
#'
#' @param ci, Report the highest posterior density interval for each agreement coefficient.
#'    Defaults to \code{TRUE}.
#' @param ciLevel, Width of the credible interval.
#' @param cohensKappa, Chance-corrected agreement between two raters, with each rater's own category proportions. Use it when the same raters rated every subject. With more than two raters, every rater pair gets its own posterior.
#'    Defaults to \code{FALSE}.
#' @param cohensKappaType, Unweighted kappa treats all disagreements equally. Weighted kappa accounts for the degree of disagreement and requires ordinal ratings.
#' \itemize{
#'   \item \code{"unweighted"}: All disagreements are treated as equal regardless of their magnitude.
#'   \item \code{"weighted"}: Disagreements are penalised according to their magnitude. Requires ordinal ratings.
#' }
#' @param dataStructure, Specify whether raters are arranged in columns (default) or rows in the dataset.
#' \itemize{
#'   \item \code{"ratersInColumns"}: Each column is one rater and each row is one subject or item being rated.
#'   \item \code{"ratersInRows"}: Each row is one rater and each column is one subject or item being rated.
#' }
#' @param dirichletPriorConcentration, Concentration of the Dirichlet prior on each cell of a rater pair's agreement table, that is, the number of prior pseudo-observations per cell. Small values let the data dominate. Larger values pull kappa towards 0, especially with many categories, because the pseudo-observations also fall into the disagreement cells.
#' @param fleissKappa, Chance-corrected agreement with the category proportions pooled over both raters of a pair, so raters are treated as interchangeable. Use it when different raters rated different subjects. With two raters it equals Fleiss' kappa and Scott's pi; with more than two raters, every rater pair gets its own posterior.
#'    Defaults to \code{FALSE}.
#' @param observedAndChanceAgreement, Report the posterior means of the observed agreement and of the agreement expected by chance, from which kappa is computed.
#'    Defaults to \code{FALSE}.
#' @param posteriorPlot, Plot the posterior density of each coefficient for every rater pair, with the credible interval shaded.
#'    Defaults to \code{FALSE}.
#' @param samples, Number of draws from the posterior distribution of the cell probabilities. Higher values give more stable estimates.
#' @param variables, Rating variables to include. Whether a variable represents a rater or a subject/item depends on the data structure setting.
#' @param weightType, Weighting scheme applied to disagreements between ordinal categories.
#' \itemize{
#'   \item \code{"quadratic"}: Penalises larger disagreements quadratically; sensitive to large discrepancies.
#'   \item \code{"linear"}: Penalises disagreements proportionally to their size.
#' }
#' @export
raterAgreementBayesian <- function(
          data = NULL,
          version = "0.97.1",
          ci = TRUE,
          ciLevel = 0.95,
          cohensKappa = FALSE,
          cohensKappaType = "unweighted",
          dataStructure = "ratersInColumns",
          dirichletPriorConcentration = 0.001,
          fleissKappa = FALSE,
          observedAndChanceAgreement = FALSE,
          plotHeight = 320,
          plotWidth = 480,
          posteriorPlot = FALSE,
          samples = 5000,
          seed = 1,
          setSeed = FALSE,
          variables = list(types = list(), value = list()),
          weightType = "quadratic") {

   defaultArgCalls <- formals(jaspReliability::raterAgreementBayesian)
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

   return(jaspBase::runWrappedAnalysis("jaspReliability", "raterAgreementBayesian", "RaterAgreementBayesian.qml", options, version, TRUE))
}
