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

	RadioButtonGroup
	{
		name:		"scoresMethod"
		title:		qsTr("Mean and standard deviation of")

		RadioButton { value: "sumScores";	label: qsTr("participants' sum scores"); checked: true}
		RadioButton { value: "meanScores";	label: qsTr("participants' mean scores")}
	}

	Group
	{
		title: qsTr("Item statistics")
		CheckBox { name: "itemDeletedOmegaT";	label: qsTr("McDonald's ω_t (if item dropped)") }
		CheckBox { name: "itemDeletedOmegaH";	label: qsTr("McDonald's ω_h (if item dropped)"); enabled: modelType.currentValue !== "correlated" }
		CheckBox { name: "itemRestCor";	label: qsTr("Item-rest correlation") }
	}

	Section
	{
		title: 	qsTr("Output options")
		Group
		{
			title: qsTr("Scale Statistics")
			CIField
			{
				name:			"credibleIntervalValue";
				label:			qsTr("Credible interval");
				defaultValue:	95;
			}
		}

		Group
		{
			CheckBox
			{
				name: 	"plotPosterior";
				label: 	qsTr("Plot Posteriors");
				id:		postPlot

				CheckBox
				{
					name: 	"fixXRange";
					label: 	qsTr("Fix range to 0-1")
				}

				CheckBox
				{
					name: 	"dispPrior";
					label: 	qsTr("Display Priors")
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

				RowLayout
				{
					DoubleField
					{
						id:				probTableValueLow
						name:			"probTableValueLow"
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
						name:			"probTableValueHigh"
						label:			""
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
					enabled:	probTable.checked & postPlot.checked
					x:			Theme.subOptionOffset
				}
			}
		}

		Group
		{
			title: qsTr("")

			RadioButtonGroup
			{
				name: "pointEst"
				title: qsTr("Posterior Point Estimate")
				RadioButton{ value: "mean"; label: qsTr("Mean"); checked: true }
				RadioButton{ value: "median"; label: qsTr("Median") }
			}
		}

		Group
		{
			title: qsTr("Model fit")
			CheckBox
			{
				name:		"dispPPC"
				label:		qsTr("Posterior predictive check");
			}
			CheckBox
			{
				name:		"fitMeasures"
				label:		qsTr("Fit measures");

				CIField
				{
					name:			"credibleIntervalValueFitMeasures";
					label:			qsTr("Credible interval");
					defaultValue:	90
				}

				DoubleField
				{
					name:			"fitCutoffSat"
					label:			qsTr("p(RMSEA <")
					defaultValue:	.08
					min:			0
					max:			1
					fieldWidth: 	40
					afterLabel:		qsTr(")")
				}
				DoubleField
				{
					name:			"fitCutoffNull"
					label:			qsTr("p(CFI/TLI >")
					defaultValue:	.9
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
				fieldWidth: 	40
				min: 			1
				max: 			1e5
			}

			IntegerField
			{
				name: 			"noChains"
				label: 			qsTr("No. chains")
				defaultValue: 	3
				fieldWidth: 	40
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
				name:				"disableSampleSave"
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
			AvailableVariablesList { name: "normalScaledItems";	 title: qsTr("Normal-Scaled Items"); source: factors.name }
			AssignedVariablesList {  name: "reverseScaledItems"; title: qsTr("Reverse-Scaled Items") }
		}
	}

	Section
	{
		title: qsTr("Priors")

		Group
		{
			title: qsTr("Measurement level")

			Group
			{
				title: qsTr("Residual variances: Inverse gamma")
				IntegerField
				{
					name:			"igShapeManifest"
					label:			qsTr("shape")
					defaultValue:	2
					min:			0
					max:			100
					fieldWidth: 	40
				}
				IntegerField
				{
					name:			"igScaleManifest"
					label:			qsTr("scale")
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

			Group
			{
				title: qsTr("Residual variances: Inverse gamma")
				IntegerField
				{
					name:			"igShapeLatent"
					label:			qsTr("shape")
					defaultValue:	2
					min:			0
					max:			100
					fieldWidth: 	40
				}
				IntegerField
				{
					name:			"igScaleLatent"
					label:			qsTr("scale")
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
					defaultValue:	2
					min:			1
					max:			100
					fieldWidth: 	40
				}
				IntegerField
				{
					name:			"igScaleGFactor"
					label:			qsTr("scale")
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
			title: qsTr("Missing Values")
			RadioButtonGroup
			{
				title: 	qsTr("")
				name: 	"missingValues"

				RadioButton { value: "excludeCasesPairwise"; label: qsTr("Bayesian imputation"); checked: true}
				RadioButton { value: "excludeCasesListwise"; label: qsTr("Exclude cases listwise")}
			}
		}

	}

}
