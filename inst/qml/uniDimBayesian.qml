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
import JASP.Widgets		1.0

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
				name:			"credibleIntervalValueScale";
				label:			qsTr("Credible interval");
				defaultValue:	95
			}

			CheckBox
			{
				id:			omega
				name:		"omegaScale"
				label:		qsTr("McDonald's ω")
				checked:	true

				CheckBox
				{
					name:		"dispPPC"
					label:		qsTr("Posterior predictive check");
					enabled:	omega.checked
				}

			}

			CheckBox
			{
				id:		alpha
				name:	"alphaScale";
				label:	qsTr("Cronbach's α");
			}

			CheckBox
			{
				id:		lambda2
				name:	"lambda2Scale";
				label:	qsTr("Guttman's λ2");
			}

			CheckBox
			{
				id:		lambda6
				name:	"lambda6Scale";
				label:	qsTr("Guttman's λ6");
			}

			CheckBox
			{
				id:		glb
				name:	"glbScale";
				label:	qsTr("Greatest lower bound");
			}

			CheckBox { name: "averageInterItemCor";	label: qsTr("Average interitem correlation")	}
			CheckBox
			{
				name: "meanScale"
				label: qsTr("Mean")

				RadioButtonGroup
				{
						title: 	qsTr("")
						name: 	"meanMethod"

						RadioButton { value: "sumScores"; label: qsTr("of sum scores"); checked: true}
						RadioButton { value: "meanScores"; label: qsTr("of mean scores")}
				}
			}
			CheckBox
			{
				name:	"sdScale"
				label:	qsTr("Standard deviation")

				RadioButtonGroup
				{
						title: 	qsTr("")
						name: 	"sdMethod"

						RadioButton { value: "sumScores"; label: qsTr("of sum scores"); checked: true}
						RadioButton { value: "meanScores"; label: qsTr("of mean scores")}
				}
			}
		}

		Group
		{
			title: qsTr("Individual Item Statistics")

			CIField
			{
				name: 			"credibleIntervalValueItem";
				label: 			qsTr("Credible interval");
				defaultValue: 	95
			}

			CheckBox
			{
				id: 		omegaItem
				name: 		"omegaItem";
				label: 		qsTr("McDonald's ω  (if item dropped)");
				enabled: 	omega.checked
			}

			CheckBox
			{
				id: 		alphaItem
				name: 		"alphaItem";
				label: 		qsTr("Cronbach's α (if item dropped)");
				enabled: 	alpha.checked
			}

			CheckBox
			{
				id: 		lambda2Item
				name: 		"lambda2Item";
				label: 		qsTr("Guttman's λ2 (if item dropped)");
				enabled: 	lambda2.checked
			}

			CheckBox
			{
				id: 		lambda6item
				name: 		"lambda6Item"
				label: 		qsTr("Guttman's λ6 (if item dropped)");
				enabled: 	lambda6.checked
			}

			CheckBox
			{
				id: 		glbItem
				name: 		"glbItem";
				label: 		qsTr("Greatest lower bound (if item dropped)");
				enabled: 	glb.checked
			}

			CheckBox
			{
				id: 		plotItem
				name: 		"plotItem";
				label: 		qsTr("If item dropped plot");
				enabled: 	omegaItem.checked || alphaItem.checked || lambda2Item.checked || lambda6Item.checked || glbItem.checked;

				CheckBox
				{
					name: 		"orderItem";
					label: 		qsTr("Order items");
					enabled: 	plotItem.checked

					RadioButtonGroup
					{
						title: 	""
						name: 	"orderType"

						RadioButton { value: "orderItemMean"; 	label: qsTr("Order items by mean");			checked: true	}
						RadioButton { value: "orderItemKL"; 	label: qsTr("Order items by KL-divergence")					}
						RadioButton { value: "orderItemKS"; 	label: qsTr("Order items by KS-distance")					}
					}
				}
			}

			CheckBox { name: "itemRestCor";						label: qsTr("Item-rest correlation")			}
			CheckBox { name: "itemMean";						label: qsTr("Mean")								}
			CheckBox { name: "itemSd";							label: qsTr("Standard deviation")				}
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
	}

	Section
	{
		title: qsTr("Convergence")

		Group
		{
			title: qsTr("MCMC parameters");

			IntegerField
			{
				id:				noSamples
				name: 			"noSamples"
				label: 			qsTr("No. samples")
				defaultValue: 	1000
				fieldWidth: 	100
				min: 			parseFloat(noBurnin.value) + parseFloat(noThin.value) + 1
				max: 			1e7
			}

			IntegerField
			{
				id:				noBurnin
				name: 			"noBurnin"
				label: 			qsTr("No. burnin samples")
				defaultValue: 	50
				fieldWidth: 	100
				min: 			1
				max:			parseFloat(noSamples.value) - parseFloat(noThin.value) -1
			}

			IntegerField
			{
				id:				noThin
				name: 			"noThin"
				label: 			qsTr("Thinning")
				defaultValue: 	1
				fieldWidth: 	100
				min: 			1
				max: 			parseFloat(noSamples.value) - parseFloat(noBurnin.value) -1
			}

			IntegerField
			{
				name: 			"noChains"
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

			RowLayout
			{
				CheckBox
				{
					name:				"disableSampleSave"
					label:				qsTr("Disable saving samples")
					checked:			false
				}
				HelpButton
				{
					toolTip: 						qsTr("Click to learn more about saving the samples.")
					helpPage:						"toolTip/sampleSavingBayes"
				}
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
		title: qsTr("Missing Data Handling")

		RadioButtonGroup
		{
			title: 	qsTr("Missing Values")
			name: 	"missingValues"

			RadioButton { value: "excludeCasesPairwise"; label: qsTr("Exclude cases pairwise"); checked: true}
			RadioButton { value: "excludeCasesListwise"; label: qsTr("Exclude cases listwise")}
		}
	}
}
