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

import QtQuick 			2.8
import JASP.Controls 	1.0
import JASP.Theme 		1.0
import QtQuick.Layouts	1.3

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
			allowedColumns:	["scale", "ordinal"]
			id:				vars
		}
	}

	Section
	{
		title: qsTr("Analysis")

		Group
		{
			title: qsTr("Scale Statistics")

			CIField
			{
				name:			"scaleCiLevel";
				label:			qsTr("Credible interval");
				defaultValue:	95
			}

			CheckBox
			{
				id:			omega
				name:		"scaleOmega"
				label:		qsTr("McDonald's ω")
				checked:	true
			}

			CheckBox
			{
				id:		alpha
				name:	"scaleAlpha";
				label:	qsTr("Cronbach's α");
			}

			CheckBox
			{
				id:		lambda2
				name:	"scaleLambda2";
				label:	qsTr("Guttman's λ2");
			}

			CheckBox
			{
				id:		lambda6
				name:	"scaleLambda6";
				label:	qsTr("Guttman's λ6");
			}

			CheckBox
			{
				id:		glb
				name:	"scaleGreatestLowerBound";
				label:	qsTr("Greatest lower bound");
			}

			CheckBox { name: "averageInterItemCorrelation";	label: qsTr("Average interitem correlation")	}

			RowLayout {
				CheckBox { name: "scaleMean";	label: qsTr("Mean");	id: mean}
				CheckBox { name: "scaleSd";		label: qsTr("SD");		id: sd}

			}
			RadioButtonGroup
			{
				indent:		true
				enabled:	mean.checked || sd.checked
				name:		"meanSdScoresMethod"

				RadioButton { value: "sumScores";	label: qsTr("of participants' sum scores"); checked: true}
				RadioButton { value: "meanScores";	label: qsTr("of participants' mean scores")}
			}

		}

		Group
		{
			title: qsTr("Individual Item Statistics")

			CIField
			{
				name: 			"itemCiLevel";
				label: 			qsTr("Credible interval");
				defaultValue: 	95
			}

			CheckBox
			{
				id: 		omegaItem
				name: 		"itemDeletedOmega";
				label: 		qsTr("McDonald's ω  (if item dropped)");
				enabled: 	omega.checked
			}

			CheckBox
			{
				id: 		alphaItem
				name: 		"itemDeletedAlpha";
				label: 		qsTr("Cronbach's α (if item dropped)");
				enabled: 	alpha.checked
			}

			CheckBox
			{
				id: 		lambda2Item
				name: 		"itemDeletedLambda2";
				label: 		qsTr("Guttman's λ2 (if item dropped)");
				enabled: 	lambda2.checked
			}

			CheckBox
			{
				id: 		lambda6Item
				name: 		"itemDeletedLambda6"
				label: 		qsTr("Guttman's λ6 (if item dropped)");
				enabled: 	lambda6.checked
			}

			CheckBox
			{
				id: 		glbItem
				name: 		"itemDeletedGreatestLowerBound";
				label: 		qsTr("Greatest lower bound (if item dropped)");
				enabled: 	glb.checked
			}

			CheckBox
			{
				id: 		itemPlot
				name: 		"itemDeletedPlot";
				label: 		qsTr("If item dropped plot");
				enabled: 	omegaItem.checked || alphaItem.checked || lambda2Item.checked || lambda6Item.checked || glbItem.checked;

				CheckBox
				{
					name: 		"itemDeletedPlotOrdered";
					label: 		qsTr("Order items");
					enabled: 	itemPlot.checked

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

			CheckBox { name: "itemRestCorrelation";				label: qsTr("Item-rest correlation")			}
			CheckBox { name: "itemMean";						label: qsTr("Mean")								}
			CheckBox { name: "itemSd";							label: qsTr("Standard deviation")				}
		}

		Group
		{
			CheckBox
			{
				name: 	"posteriorPlot";
				label: 	qsTr("Plot Posteriors");
				id:		postPlot

				CheckBox
				{
					name: 	"posteriorPlotFixedRange";
					label: 	qsTr("Fix range to 0-1")
				}

				CheckBox
				{
					name: 	"posteriorPlotPriorDisplayed";
					label: 	qsTr("Display Priors")
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

			IntegerField
			{
				name: 			"samples"
				label: 			qsTr("No. samples")
				defaultValue: 	1000
				fieldWidth: 	100
				min: 			100
				max: 			1e7
			}

			IntegerField
			{
				name: 			"burnin"
				label: 			qsTr("No. burnin samples")
				defaultValue: 	50
				fieldWidth: 	100
				min: 			1
				max:			1e6
			}

			IntegerField
			{
				name: 			"thinning"
				label: 			qsTr("Thinning")
				defaultValue: 	1
				fieldWidth: 	100
				min: 			1
				max: 			1e5
			}

			IntegerField
			{
				name: 			"chains"
				label: 			qsTr("No. chains")
				defaultValue: 	3
				fieldWidth: 	100
				min: 			2
				max: 			100
			}
		}

		Group
		{
			title: qsTr("Diagnostics")

			CheckBox {	name: "rHat";		label: qsTr("R-hat");		}
			CheckBox {	name: "tracePlot";	label: qsTr("Traceplots");	}
		}

		Group
		{
			title: qsTr("Repeatability")

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
				name:				"samplesSavingDisabled"
				label:				qsTr("Disable saving samples")
				checked:			false
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
			AssignedVariablesList {  name: "reverseScaledItems"; title: qsTr("Reverse-Scaled Items") }
		}
	}

	Section
	{
		title: qsTr("Priors")
		Group
		{
			title: qsTr("CTT-Coefficients (α, λ2, λ6, glb)")

			FormulaField
			{
				name:			"inverseWishartPriorScale"
				label:			qsTr("Inverse Wishart scale")
				defaultValue:	"1e-10"
				min:			0
				max:			100
				fieldWidth: 	40
			}
			DoubleField
			{
				name:			"inverseWishartPriorDf"
				label:			qsTr("Inverse Wishart df")
				defaultValue:	vars.count
				min:			vars.count
				max:			vars.count + 100
				fieldWidth: 	40
			}

		}
		Group
		{
			title: qsTr("McDonald's ω")

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
				}


				IntegerField
				{
					name:			"inverseGammaPriorScale"
					label:			qsTr("scale")
					defaultValue:	1
					min:			0
					max:			100
					fieldWidth: 	40
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

				RadioButton { value: "imputation";	label: qsTr("Bayesian imputation"); checked: true	}
				RadioButton { value: "listwise";	label: qsTr("Exclude cases listwise")				}
			}
		}

		Group
		{
			title: qsTr("McDonald's ω Estimation")
			enabled: omega.checked
			CheckBox
			{
				name:		"omegaPosteriorPredictiveCheck"
				label:		qsTr("Posterior predictive check");
			}
			CheckBox
			{
				name:		"omegaFitMeasures"
				label:		qsTr("Fit measures");

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
			}
		}

		Group
		{
			RadioButtonGroup
			{
				title: qsTr("Coefficients")
				name: "coefficientType"

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
				RadioButton{ value: "mean";		label: qsTr("Mean"); checked: true	}
				RadioButton{ value: "median";	label: qsTr("Median")				}
			}
		}


	}
}
