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
import QtQuick.Layouts
import QtQuick.Controls
import JASP.Controls 
import JASP.Theme	
import JASP.Widgets 

Form
{
	FactorsForm
	{
		id: factors
		name: "factors"
		initNumberFactors: 2
		allowedColumns: ["scale"]
		info: qsTr("Assign the items to their factors. The analysis requires at least two factors with at least two items each.")
	}

	DropDown
	{
		id:				modelType
		name:			"modelType"
		label:			qsTr("Model")
		values:
		[
			{ label: qsTr("Second-order"),		value: "secondOrder"	},
			{ label: qsTr("Bi-factor"),			value: "biFactor"		},
			{ label: qsTr("Correlated factors"),	value: "correlated"		}
		]
		info: qsTr("The factor model used to estimate the reliability coefficients. McDonald's ω_h (general/group-common reliability) is only available for the second-order and bi-factor models.")
	}

	Section
	{
		title: qsTr("Analysis")
		Group 
		{
			CIField
			{
				name:			"credibleIntervalValue";
				label:			qsTr("Credible interval");
				defaultValue:	95;
				info:			qsTr("Width of the credible interval for the reliability coefficients.")
			}
			RadioButtonGroup
			{
				name:		"scoresMethod"
				title:		qsTr("Mean and standard deviation of")
				info:		qsTr("Whether the mean and standard deviation in the scale table are based on participants' sum scores or mean scores across items.")
				RadioButton { value: "sumScores";	label: qsTr("participants' sum scores"); checked: true}
				RadioButton { value: "meanScores";	label: qsTr("participants' mean scores")}
			}
		}

		Group
		{
			title: qsTr("Item statistics")
			CheckBox { name: "itemDeletedOmegaT";	label: qsTr("McDonald's ω_t (if item dropped)"); info: qsTr("Posterior point estimate of ω_t for the remaining items when this item is dropped. The model is refit once per item.") }
			CheckBox { name: "itemDeletedOmegaH";	label: qsTr("McDonald's ω_h (if item dropped)"); enabled: modelType.currentValue !== "correlated"; info: qsTr("Posterior point estimate of ω_h for the remaining items when this item is dropped. Not available for the correlated-factors model.") }
			CheckBox { name: "itemRestCor";	label: qsTr("Item-rest correlation"); info: qsTr("Pearson correlation between each item and the sum of the remaining items.") }
		}

		Group
		{
			CheckBox
			{
				name: 	"plotPosterior";
				label: 	qsTr("Plot Posteriors");
				id:		postPlot
				info:	qsTr("Display posterior density plots for the reliability coefficients.")

				CheckBox
				{
					name: 	"fixXRange";
					label: 	qsTr("Fix range to 0-1")
					info:	qsTr("Fix the x-axis of the posterior plots to [0, 1] for easier comparison between coefficients.")
				}

				CheckBox
				{
					name: 	"dispPrior";
					label: 	qsTr("Display Priors")
					info:	qsTr("Add the prior distribution to the posterior density plot.")
				}

			}
		}

		Group
		{
			CheckBox
			{
				id:					probTable
				name:				"probTable"
				label:				qsTr("Probability for:")
				childrenOnSameRow:	true
				info:				qsTr("Report the posterior probability that a reliability coefficient falls within the specified interval.")

				RowLayout
				{
					DoubleField
					{
						id:				probTableValueLow
						name:			"probTableValueLow"
						label:			""
						info:			qsTr("Lower bound of the reliability interval.")
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
						name:			"probTableValueHigh"
						label:			""
						info:			qsTr("Upper bound of the reliability interval.")
						defaultValue:	.90
						min:			0
						max:			1
						decimals:		2
						fieldWidth: 	40
					}
				}
			}


			Item
			{
				width:	shadePlots.width + Theme.subOptionOffset
				height: shadePlots.height

				CheckBox
				{
					id:			shadePlots
					name:		"shadePlots";
					indent:		true
					label:		qsTr("Shade posterior region in plot");
					info:		qsTr("Shade the posterior region corresponding to the probability interval in the density plot.")
					enabled:	probTable.checked & postPlot.checked
					x:			Theme.subOptionOffset
				}
			}
		}
		
		Group
		{
			title: qsTr("Model fit")
			CheckBox
			{
				name:		"dispPPC"
				label:		qsTr("Posterior predictive check");
				info:		qsTr("Graphical check of factor model fit: eigenvalues of the observed covariance matrix are compared against the model-implied posterior distribution.")
			}
			CheckBox
			{
				name:		"fitMeasures"
				label:		qsTr("Fit measures");
				info:		qsTr("Bayesian fit measures for the factor model: B-LR, B-SRMR, and B-RMSEA with a posterior probability relative to the cutoff.")

				CIField
				{
					name:			"credibleIntervalValueFitMeasures";
					label:			qsTr("Credible interval");
					defaultValue:	90
					info:			qsTr("Width of the credible interval for the fit measures.")
				}

				DoubleField
				{
					name:			"fitCutoffSat"
					label:			qsTr("p(RMSEA <")
					defaultValue:	.08
					info:			qsTr("Cutoff for the posterior probability that the B-RMSEA is below this value.")
					min:			0
					max:			1
					fieldWidth: 	40
					afterLabel:		qsTr(")")
				}
			}
		}

	}

	

	Section
	{
		title: qsTr("MCMC options")

		Group
		{
			title: qsTr("MCMC parameters");

			IntegerField
			{
				id:				noSamples
				name: 			"noSamples"
				label: 			qsTr("No. samples")
				defaultValue: 	2000
				info: 			qsTr("Total number of MCMC samples per chain, including burn-in.")
				fieldWidth: 	60
				min: 			100
				max: 			1e7
			}

			IntegerField
			{
				id:				noBurnin
				name: 			"noBurnin"
				label: 			qsTr("No. burnin samples")
				defaultValue: 	200
				info: 			qsTr("Initial samples discarded while the chain converges to the posterior.")
				fieldWidth: 	60
				min: 			1
				max:			1e6
			}

			IntegerField
			{
				id:				noThin
				name: 			"noThin"
				label: 			qsTr("Thinning")
				defaultValue: 	1
				info: 			qsTr("Keep every k-th sample to reduce autocorrelation. A value of 1 keeps all samples.")
				fieldWidth: 	40
				min: 			1
				max: 			1e5
			}

			IntegerField
			{
				name: 			"noChains"
				label: 			qsTr("No. chains")
				defaultValue: 	3
				info: 			qsTr("Number of independent MCMC chains. Multiple chains enable R-hat convergence diagnostics.")
				fieldWidth: 	40
				min: 			2
				max: 			100
			}
		}

		Group
		{
			title: qsTr("Diagnostics")

			CheckBox {	name: "rHat";		label: qsTr("R-hat");		info: qsTr("Potential scale reduction factor. Values close to 1 (< 1.1) indicate convergence across chains.")	}
			CheckBox {	name: "tracePlot";	label: qsTr("Traceplots");	info: qsTr("Plot of sampled values per chain over iterations. Well-mixed chains indicate convergence.")	}
		}

		Group
		{
			title: qsTr("Repeatability")

			CheckBox
			{
				name: 				"setSeed"
				label: 				qsTr("Set seed")
				childrenOnSameRow: 	true
				info: 				qsTr("Fix the random number generator seed to make the MCMC results reproducible.")

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
				name:				"disableSampleSave"
				label:				qsTr("Disable saving samples")
				checked:			false
				info:				qsTr("When checked, MCMC samples are not stored in the output file. Reduces file size but re-running the analysis or changing options requires resampling.")
			}

		}
	}

	Section
	{
		title: qsTr("Reverse-Scaled Items")

		VariablesForm
		{
			height: 150
			AvailableVariablesList { name: "normalScaledItems";	 title: qsTr("Normal-Scaled Items"); source: factors.name }
			AssignedVariablesList {  name: "reverseScaledItems"; title: qsTr("Reverse-Scaled Items"); info: qsTr("Items assigned here are recoded (reverse-scored) before the analysis.") }
		}
	}

	Section
	{
		title: qsTr("Priors")

		Group
		{
			title: qsTr("Measurement level")
			info:  qsTr("Priors on the measurement-level parameters: item residual variances and item factor loadings.")

			Group
			{
				title: qsTr("Residual variances: Inverse gamma")
				IntegerField
				{
					name:			"igShapeManifest"
					label:			qsTr("shape")
					info:			qsTr("Shape parameter of the inverse gamma prior on the item residual variances.")
					defaultValue:	2
					min:			0
					max:			100
					fieldWidth: 	40
				}
				IntegerField
				{
					name:			"igScaleManifest"
					label:			qsTr("scale")
					info:			qsTr("Scale parameter of the inverse gamma prior on the item residual variances.")
					defaultValue:	1
					min:			0
					max:			100
					fieldWidth: 	40
				}
			}

			Group
			{
				title: qsTr("Loadings: Normal")
				DoubleField
				{
					name:			"loadMeanManifest"
					label:			qsTr("mean")
					info:			qsTr("Mean of the normal prior on the item factor loadings.")
					defaultValue:	0
					min:			-10
					max:			10
					fieldWidth: 	40
				}
			}
		}

		Group
		{
			title: qsTr("Structural level")
			visible: modelType.currentValue !== "correlated"
			info:  qsTr("Priors on the structural parameters of the second-order and bi-factor models: latent residual variances, structural loadings, and the general factor variance.")

			Group
			{
				title: qsTr("Residual variances: Inverse gamma")
				IntegerField
				{
					name:			"igShapeLatent"
					label:			qsTr("shape")
					info:			qsTr("Shape parameter of the inverse gamma prior on the latent residual variances.")
					defaultValue:	2
					min:			0
					max:			100
					fieldWidth: 	40
				}
				IntegerField
				{
					name:			"igScaleLatent"
					label:			qsTr("scale")
					info:			qsTr("Scale parameter of the inverse gamma prior on the latent residual variances.")
					defaultValue:	1
					min:			0
					max:			100
					fieldWidth: 	40
				}
			}

			Group
			{
				title: qsTr("Loadings: Normal")
				DoubleField
				{
					name:			"loadMeanLatent"
					label:			qsTr("mean")
					info:			qsTr("Mean of the normal prior on the structural loadings of the group factors on the general factor.")
					defaultValue:	0
					min:			-10
					max:			10
					fieldWidth: 	40
				}
			}
			Group
			{
				title: qsTr("Factor variance: Inverse gamma")

				IntegerField
				{
					name:			"igShapeGFactor"
					label:			qsTr("shape")
					info:			qsTr("Shape parameter of the inverse gamma prior on the general factor variance.")
					defaultValue:	2
					min:			1
					max:			100
					fieldWidth: 	40
				}
				IntegerField
				{
					name:			"igScaleGFactor"
					label:			qsTr("scale")
					info:			qsTr("Scale parameter of the inverse gamma prior on the general factor variance.")
					defaultValue:	1
					min:			1
					max:			100
					fieldWidth: 	40
				}
			}
		}



		Group
		{
			title: qsTr("Latent correlations")
			visible: modelType.currentValue === "correlated"

			Group
			{
				title: qsTr("Correlation matrix: Inverse-Wishart")
				IntegerField
				{
					name:			"latentCorDf"
					label:			qsTr("degrees of freedom")
					info:			qsTr("Degrees of freedom of the inverse Wishart prior on the latent correlation matrix. Values below the number of factors are raised to the number of factors to keep the prior proper.")
					defaultValue:	2
					min:			1
					max:			1000
					fieldWidth: 	40
				}
			}
		}
	}

	Section
	{
		title: qsTr("Advanced Options")
		Group
		{
			title: qsTr("")

			RadioButtonGroup
			{
				name: "pointEst"
				title: qsTr("Posterior Point Estimate")
				info: qsTr("Whether to report the posterior mean or median as the point estimate in the tables.")
				RadioButton{ value: "mean"; label: qsTr("Mean"); checked: true }
				RadioButton{ value: "median"; label: qsTr("Median") }
			}
		}

		Group
		{
			title: qsTr("Missing Values")
			RadioButtonGroup
			{
				title: 	qsTr("")
				name: 	"missingValues"
				info: 	qsTr("Bayesian imputation treats missing values as unknown parameters sampled from the posterior; listwise deletion removes any row with a missing value.")

				RadioButton { value: "excludeCasesPairwise"; label: qsTr("Bayesian imputation"); checked: true}
				RadioButton { value: "excludeCasesListwise"; label: qsTr("Exclude cases listwise")}
			}
		}

	}

}
