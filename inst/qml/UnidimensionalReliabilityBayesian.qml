//
// Copyright (C) 2013-2020 University of Amsterdam
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//

import QtQuick
import JASP.Controls
import QtQuick.Layouts

Form
{
	VariablesForm
	{
		height: 300
		AvailableVariablesList { name: "allVariablesList" }
		AssignedVariablesList
		{
			name:			"variables"
			title:			qsTr("Variables")
			allowedColumns:	["scale"]
			id:				vars
			info:			qsTr("Items/variables to include in the reliability analysis. Must be scale variables.")
		}
	}

	Section
	{
		title: qsTr("Analysis")

		Group
		{
			title: qsTr("Scale Statistics")
			info:  qsTr("Reliability coefficients and descriptive statistics computed for the full scale.")

			CIField
			{
				name:			"scaleCiLevel";
				label:			qsTr("Credible interval");
				defaultValue:	95
				info:			qsTr("Width of the credible interval for scale reliability statistics.")
			}

			CheckBox
			{
				id:			omega
				name:		"scaleOmega"
				label:		qsTr("McDonald's ω")
				checked:	true
				info:		qsTr("McDonald's omega for unidimensional data based on the single-factor model.")
			}

			CheckBox
			{
				id:		alpha
				name:	"scaleAlpha";
				label:	qsTr("Cronbach's α");
				info:	qsTr("Cronbach's alpha. For binary items this equals KR-20.")
			}

			CheckBox
			{
				id:		lambda2
				name:	"scaleLambda2";
				label:	qsTr("Guttman's λ2");
				info:	qsTr("Guttman's lambda 2, a lower bound for reliability.")
			}

			CheckBox
			{
				id:		splithalf
				name:	"scaleSplithalf";
				label:	qsTr("Split-half coefficient");
				info:	qsTr("Splits items into two halves (odd/even by default). Unstandardized: Flanagan-Rulon; Standardized: Spearman-Brown.")
			}

			CheckBox { name: "averageInterItemCorrelation";	label: qsTr("Average interitem correlation"); info: qsTr("Mean of all pairwise Pearson correlations between items.") }

			RowLayout {
				CheckBox { name: "scaleMean";	label: qsTr("Mean");	id: mean}
				CheckBox { name: "scaleVar";	label: qsTr("Variance"); id: variance}
				CheckBox { name: "scaleSd";		label: qsTr("SD");		id: sd}
			}
			RadioButtonGroup
			{
				indent:		true
				enabled:	mean.checked || sd.checked || variance.checked
				name:		"meanSdScoresMethod"
				info:		qsTr("Whether the mean, variance, and SD are based on sum scores or mean scores across items.")

				RadioButton { value: "sumScores";	label: qsTr("of participants' sum scores"); checked: true}
				RadioButton { value: "meanScores";	label: qsTr("of participants' mean scores")}
			}

		}

		Group
		{
			title: qsTr("Individual Item Statistics")
			info:  qsTr("Reliability coefficients and descriptive statistics per item, including posterior distributions when that item is removed.")

			CIField
			{
				name: 			"itemCiLevel";
				label: 			qsTr("Credible interval");
				defaultValue: 	95
				info:			qsTr("Width of the credible interval for item-level statistics.")
			}

			CheckBox
			{
				id: 		omegaItem
				name: 		"itemDeletedOmega";
				label: 		qsTr("McDonald's ω (if item dropped)");
				enabled: 	omega.checked
				info:		qsTr("Posterior distribution of omega for the remaining items when this item is removed.")
			}

			CheckBox
			{
				id: 		alphaItem
				name: 		"itemDeletedAlpha";
				label: 		qsTr("Cronbach's α (if item dropped)");
				enabled: 	alpha.checked
				info:		qsTr("Posterior distribution of alpha for the remaining items when this item is removed.")
			}

			CheckBox
			{
				id: 		lambda2Item
				name: 		"itemDeletedLambda2";
				label: 		qsTr("Guttman's λ2 (if item dropped)");
				enabled: 	lambda2.checked
				info:		qsTr("Posterior distribution of lambda 2 for the remaining items when this item is removed.")
			}

			CheckBox
			{
				id: 		splithalfItem
				name: 		"itemDeletedSplithalf";
				label: 		qsTr("Split-half coefficient (if item dropped)");
				enabled: 	splithalf.checked
				info:		qsTr("Posterior distribution of the split-half coefficient for the remaining items when this item is removed.")
			}

			CheckBox
			{
				id: 		itemPlot
				name: 		"itemDeletedPlot";
				label: 		qsTr("If item dropped plot");
				enabled: 	omegaItem.checked || alphaItem.checked || lambda2Item.checked || splithalfItem.checked
				info:		qsTr("Displays posterior densities of the reliability of the remaining items when each item is dropped.")

				CheckBox
				{
					name: 		"itemDeletedPlotOrdered";
					label: 		qsTr("Order items");
					enabled: 	itemPlot.checked
					info:		qsTr("Sort the densities by how much removing an item changes the posterior (by mean, KL-divergence, or KS-distance).")

					RadioButtonGroup
					{
						title: 	""
						name: 	"itemDeletedPlotOrderedType"

						RadioButton { value: "mean";				label: qsTr("Order items by mean");			checked: true	}
						RadioButton { value: "kullbackLeibler"; 	label: qsTr("Order items by KL-divergence")					}
						RadioButton { value: "kolmogorovSmirnov"; 	label: qsTr("Order items by KS-distance")					}
					}
				}
			}

			CheckBox { name: "itemRestCorrelation";				label: qsTr("Item-rest correlation"); info: qsTr("Correlation of each item with the sum of the remaining items.") }
			RowLayout {
				CheckBox { name: "itemMean";						label: qsTr("Mean")								}
				CheckBox { name: "itemVar";						label: qsTr("Variance")								}
				CheckBox { name: "itemSd";							label: qsTr("SD")				}
			}

		}

		Group
		{
			CheckBox
			{
				name: 	"posteriorPlot";
				label: 	qsTr("Plot Posteriors");
				id:		postPlot
				info:	qsTr("Display posterior density plots for the selected reliability coefficients.")

				CheckBox
				{
					name: 	"posteriorPlotFixedRange";
					label: 	qsTr("Fix range to 0-1")
					info:	qsTr("Fix the x-axis of the posterior plots to [0, 1] for easier comparison between coefficients.")
				}

				CheckBox
				{
					name: 	"posteriorPlotPriorDisplayed";
					label: 	qsTr("Display Priors")
					info:	qsTr("Add the prior distribution to the posterior density plot.")
				}

			}
		}

		Group
		{
			columns: 1
			CheckBox
			{
				id:					probTable
				name:				"probabilityTable"
				label:				qsTr("Probability for:")
				childrenOnSameRow:	true
				info:				qsTr("Report the prior and posterior probability that a reliability coefficient falls within the specified interval.")

				RowLayout
				{
					DoubleField
					{
						name:			"probabilityTableLowerBound"
						label:			""
						defaultValue:	0.70
						min:			0
						max:			1
						decimals:		2
						fieldWidth: 	40
					}

					Label
					{	text: qsTr("< Reliability <")}

					DoubleField
					{
						name:			"probabilityTableUpperBound"
						label:			""
						defaultValue:	.90
						min:			0
						max:			1
						decimals:		2
						fieldWidth: 	40
					}
				}
			}

			CheckBox
			{
				id:			shadePlots
				name:		"posteriorPlotShaded";
				indent:		true
				label:		qsTr("Shade posterior region in plot");
				enabled:	probTable.checked & postPlot.checked
				info:		qsTr("Shade the posterior region corresponding to the probability interval in the density plot.")
				Layout.leftMargin:	childControlsArea.anchors.leftMargin
			}
		}
	}

	Section
	{
		title: qsTr("Convergence")

		Group
		{
			title: qsTr("MCMC parameters");
			info:  qsTr("Parameters controlling the MCMC sampler used to estimate the posterior distributions.")

			IntegerField
			{
				name: 			"samples"
				label: 			qsTr("No. samples")
				defaultValue: 	1000
				fieldWidth: 	100
				min: 			100
				max: 			1e7
				info:			qsTr("Total number of MCMC iterations per chain, including burn-in. The number of retained posterior samples per chain is ceiling((samples - burn-in) / thinning).")
			}

			IntegerField
			{
				name: 			"burnin"
				label: 			qsTr("No. burnin samples")
				defaultValue: 	50
				fieldWidth: 	100
				min: 			1
				max:			1e6
				info:			qsTr("Initial samples discarded while the chain converges to the posterior.")
			}

			IntegerField
			{
				name: 			"thinning"
				label: 			qsTr("Thinning")
				defaultValue: 	1
				fieldWidth: 	100
				min: 			1
				max: 			1e5
				info:			qsTr("Keep every k-th sample to reduce autocorrelation. A value of 1 keeps all samples.")
			}

			IntegerField
			{
				name: 			"chains"
				label: 			qsTr("No. chains")
				defaultValue: 	3
				fieldWidth: 	100
				min: 			2
				max: 			100
				info:			qsTr("Number of independent MCMC chains. Multiple chains enable R-hat convergence diagnostics.")
			}
		}

		Group
		{
			title: qsTr("Diagnostics")
			info:  qsTr("MCMC convergence diagnostics to assess whether the sampler has converged to the posterior.")

			CheckBox {	name: "rHat";							label: qsTr("R-hat");		info: qsTr("Potential scale reduction factor. Values close to 1 (< 1.1) indicate convergence across chains.")}
			CheckBox {	name: "tracePlot";					label: qsTr("Traceplots");	info: qsTr("Plot of sampled values per chain over iterations. Well-mixed chains indicate convergence.")}
			CheckBox {	name: "effectiveSampleSize"; 	label: qsTr("ESS");			info: qsTr("Effective sample size: number of independent posterior samples after accounting for autocorrelation.")}
		}

		Group
		{
			title: qsTr("Repeatability")
			info:  qsTr("Set a random seed to reproduce identical posterior samples across runs.")

			CheckBox
			{
				name: 				"setSeed"
				label: 				qsTr("Set seed")
				childrenOnSameRow: 	true

				IntegerField
				{
					name: 			"seed"
					label: 			""
					defaultValue: 	1234
					fieldWidth: 	100
					min: 			1
					max: 			1e9
				}
			}
		}

		Group
		{
			title: qsTr("Samples")

			CheckBox
			{
				name:			"samplesSavingDisabled"
				label:			qsTr("Disable saving samples")
				checked:		false
				info:			qsTr("When checked, MCMC samples are not stored in the output file. Reduces file size but may slow down re-running the analysis, because samples are precomputed and cached for speed.")
			}

		}
	}




	Section
	{
		title: qsTr("Reverse-Scaled Items")

		VariablesForm
		{
			height: 150
			AvailableVariablesList { name: "normalScaledItems";	 title: qsTr("Normal-Scaled Items"); source: "variables" }
			AssignedVariablesList {  
				name: "reverseScaledItems"; 
				title: qsTr("Reverse-Scaled Items") 
				allowedColumns: ["scale"]
			}
		}
	}

	Section
	{
		title: qsTr("Priors")
		Group
		{
			title: qsTr("CTT-Coefficients (α, λ2)")
			info:  qsTr("The prior on alpha, lambda 2, and the split-half is induced by an inverse Wishart prior on the covariance matrix.")

			FormulaField
			{
				name:			"inverseWishartPriorScale"
				label:			qsTr("Inverse Wishart scale")
				defaultValue:	"1e-10"
				min:			0
				max:			100
				fieldWidth: 	40
				info:			qsTr("Precision values on the diagonal of the inverse Wishart scaling matrix.")
			}
			DoubleField
			{
				name:			"inverseWishartPriorDf"
				label:			qsTr("Inverse Wishart df")
				defaultValue:	vars.count
				min:			vars.count
				max:			vars.count + 100
				fieldWidth: 	40
				info:			qsTr("Degrees of freedom of the inverse Wishart prior. Minimum equals the number of items.")
			}

		}
		Group
		{
			title: qsTr("McDonald's ω")
			info:  qsTr("The prior on omega is induced by an inverse gamma on residual variances and a normal prior on factor loadings.")

			RowLayout
			{
				Label
				{	text: qsTr("Inverse gamma:")}

				IntegerField
				{
					name:			"inverseGammaPriorShape"
					label:			qsTr("shape")
					defaultValue:	2
					min:			0
					max:			100
					fieldWidth: 	40
					info:			qsTr("Shape parameter of the inverse gamma prior on residual variances.")
				}


				IntegerField
				{
					name:			"inverseGammaPriorScale"
					label:			qsTr("scale")
					defaultValue:	1
					min:			0
					max:			100
					fieldWidth: 	40
					info:			qsTr("Scale parameter of the inverse gamma prior on residual variances.")
				}
			}
			RowLayout
			{
				Label
				{	text: qsTr("Normal:")}

				DoubleField
				{
					name:			"normalPriorMean"
					label:			qsTr("mean")
					defaultValue:	0
					min:			-10
					max:			10
					fieldWidth: 	40
					info:			qsTr("Mean of the normal prior on factor loadings.")
				}
			}


		}
	}
	Section
	{
		title: qsTr("Advanced Options")

		Group
		{
			title: qsTr("Missing Values")
			RadioButtonGroup
			{
				name: 	"naAction"
				info:	qsTr("Bayesian imputation treats missing values as unknown parameters sampled from the posterior; listwise deletion removes any row with a missing value.")

				RadioButton { value: "imputation";	label: qsTr("Bayesian imputation"); checked: true	}
				RadioButton { value: "listwise";	label: qsTr("Exclude cases listwise")				}
			}
		}

		Group
		{
			title: qsTr("McDonald's ω Estimation")
			enabled: omega.checked
			info:  qsTr("Options for fitting the single-factor model underlying McDonald's omega.")
			CheckBox
			{
				name:		"omegaPosteriorPredictiveCheck"
				label:		qsTr("Posterior predictive check");
				info:		qsTr("Graphical check for single-factor model fit: eigenvalues of the observed covariance matrix are compared against the model-implied posterior distribution.")
			}
			CheckBox
			{
				name:		"omegaFitMeasures"
				label:		qsTr("Fit measures");
				info:		qsTr("Bayesian fit indices (B-LR, B-RMSEA, B-CFI, B-TLI) for the single-factor model with probability statements relative to cutoffs.")

				CIField
				{
					name:			"omegaFitMeasuresCiLevel";
					label:			qsTr("Credible interval");
					defaultValue:	90
				}

				DoubleField
				{
					name:			"omegaFitMeasuresCutoffRmsea"
					label:			qsTr("p(RMSEA <")
					defaultValue:	.08
					min:			0
					max:			1
					fieldWidth: 	40
					afterLabel:		qsTr(")")
				}
				DoubleField
				{
					name:			"omegaFitMeasuresCutoffCfiTli"
					label:			qsTr("p(CFI/TLI >")
					defaultValue:	.9
					min:			0
					max:			1
					fieldWidth: 	40
					afterLabel:		qsTr(")")
				}
			}
			CheckBox
			{
				name:		"standardizedLoadings"
				label:		qsTr("Standardized factor loadings");
				info:		qsTr("Table of standardized loadings from the single-factor model (posterior mean or median).")
			}
		}

		Group
		{
			RadioButtonGroup
			{
				title: qsTr("Coefficients")
				name: "coefficientType"
				info:  qsTr("Unstandardized uses the covariance matrix; standardized uses the correlation matrix. Standardized alpha is contested in the literature.")

				RadioButton{ value: "unstandardized"; label: qsTr("Unstandardized"); checked: true }
				RadioButton{ value: "standardized";	label: qsTr("Standardized");
				}

			}
		}
		Group
		{
			RadioButtonGroup
			{
				name: "pointEstimate"
				title: qsTr("Posterior Point Estimate")
				info:  qsTr("Whether to report the posterior mean or median as the point estimate in the tables.")
				RadioButton{ value: "mean";		label: qsTr("Mean"); checked: true	}
				RadioButton{ value: "median";	label: qsTr("Median")				}
			}
		}


	}
}
