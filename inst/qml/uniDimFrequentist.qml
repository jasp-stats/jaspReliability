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
import JASP.Theme		1.0
import JASP.Widgets 	1.0

Form
{
	VariablesForm
	{
		height: 300

		AvailableVariablesList { name: "allVariablesList" }

		AssignedVariablesList
		{
			name: 			"variables"
			title: 			qsTr("Variables")
			allowedColumns: ["scale", "ordinal"]
		}
	}

	Section
	{
		title: 	qsTr("Analysis")

		Group
		{
			title: qsTr("Scale Statistics")

			CheckBox
			{
				name: 				"intervalOn"
				label:				qsTr("Confidence Interval")
				checked: 			true
				childrenOnSameRow: 	true
				id:					interval

				CIField
				{
					name:			"confidenceIntervalValue";
					label:			"";
					defaultValue:	95;
				}
			}

			CheckBox
			{
				id:			omega
				name:		"omegaScale"
				label:		qsTr("McDonald's ω")
				checked:	true
			}

			CheckBox
			{
				id: 	alpha
				name: 	"alphaScale";
				label: 	qsTr("Cronbach's α");
			}

			CheckBox
			{
				id: 	lambda2
				name: 	"lambda2Scale";
				label: 	qsTr("Guttman's λ2");
			}

			CheckBox
			{
				id: 	lambda6
				name: 	"lambda6Scale";
				label: 	qsTr("Guttman's λ6");
			}

			CheckBox
			{
				id: 	glb
				name: 	"glbScale";
				label: 	qsTr("Greatest lower bound");
			}

			CheckBox { name: "averageInterItemCor";	label: qsTr("Average interitem correlation")}

			RowLayout {
				CheckBox { name: "meanScale";	label: qsTr("Mean");	id: mean}
				CheckBox { name: "sdScale";		label: qsTr("SD");		id: sd}

			}
			RadioButtonGroup
			{
				indent:		true
				enabled:	mean.checked || sd.checked
				title:		qsTr("")
				name:		"scoresMethod"

				RadioButton { value: "sumScores";	label: qsTr("of participants' sum scores"); checked: true}
				RadioButton { value: "meanScores";	label: qsTr("of participants' mean scores")}
			}

		}

		Group
		{
			title: qsTr("Individual Item Statistics")

			CheckBox
			{
				name: 		"omegaItem";
				label: 		qsTr("McDonald's ω  (if item dropped)");
				enabled: 	omega.checked
			}

			CheckBox
			{
				name: 		"alphaItem";
				label: 		qsTr("Cronbach's α (if item dropped)");
				enabled: 	alpha.checked
			}

			CheckBox
			{
				name: 		"lambda2Item";
				label: 		qsTr("Guttman's λ2 (if item dropped)");
				enabled: 	lambda2.checked
			}

			CheckBox
			{
				name: 		"lambda6Item";
				label: 		qsTr("Guttman's λ6 (if item dropped)");
				enabled: 	lambda6.checked
			}

			CheckBox
			{
				name: 		"glbItem";
				label: 		qsTr("Greatest lower bound (if item dropped)");
				enabled: 	glb.checked
			}

			CheckBox { name: "itemRestCor";	label: qsTr("Item-rest correlation")		}
			CheckBox { name: "meanItem";	label: qsTr("Mean")							}
			CheckBox { name: "sdItem";		label: qsTr("Standard deviation")			}
		}
	}

	Section
	{
		title: qsTr("Reverse-Scaled Items")

		VariablesForm
		{
			height: 150

			AvailableVariablesList 	{ name: "normalScaledItems"; 	title: qsTr("Normal-Scaled Items"); source: "variables" }
			AssignedVariablesList 	{ name: "reverseScaledItems"; 	title: qsTr("Reverse-Scaled Items") }
		}
	}

	Section
	{
		title: qsTr("Advanced Options")

		RadioButtonGroup
		{
				title: 	qsTr("Missing Values")
				name: 	"missingValues"

				RadioButton { value: "excludeCasesPairwise"; label: qsTr("Exclude cases pairwise"); checked: true}
				RadioButton { value: "excludeCasesListwise"; label: qsTr("Exclude cases listwise")}
		}

	Group
	{
			title: qsTr("Bootstrap")

			IntegerField
			{
				name: 			"noSamples"
				label: 			qsTr("No. of bootstrap samples")
				defaultValue: 	1000
				fieldWidth: 	50
				min: 			100
				max: 			1e7
			}

			RadioButtonGroup
			{
				title:		""
				name:		"bootType"
				enabled:	interval.checked

				RadioButton {value: "nonParametric"; label: qsTr("Non-parametric bootstrap"); checked: true}
				RadioButton {value: "parametric"; label: qsTr("Parametric bootstrap")}
			}

		}

		RadioButtonGroup
		{
			title: 		qsTr("McDonald's ω Estimation")
			name: 		"omegaMethod"
			enabled: 	omega.checked

			RadioButton
			{
				value: 		"cfa"
				label: 		qsTr("CFA")
				checked: 	true

					CheckBox
					{
						name: 		"fitMeasures"
						label: 		qsTr("Single factor model fit")
					}

					RadioButtonGroup
					{
						title:		qsTr("Interval")
						name:		"omegaInterval"
						enabled:	interval.checked

						RadioButton
						{

							value: 		"omegaAnalytic"
							label: 		qsTr("Analytic interval")
							checked: 	true
						}

						RadioButton
						{
							value: 	"omegaBoot"
							label: 	qsTr("Bootstrapped interval")
						}
					}
				}
				RadioButton { value: "pfa"; label: qsTr("PFA")}
		}

		Group
		{
			title: 		qsTr("Cronbach's α Estimation")
			enabled: 	alpha.checked

			RadioButtonGroup
			{
				title:	""
				name: 	"alphaMethod"

				RadioButton
				{
					value:		"alphaUnstand"
					label:		qsTr("Unstandardized")
					checked:	true
				}

				RadioButton
				{
					value:	"alphaStand"
					label:	qsTr("Standardized")
				}
			}

			RadioButtonGroup
			{
				title:		qsTr("Interval")
				name:		"alphaInterval"
				enabled:	interval.checked

				RadioButton
				{
					value:		"alphaAnalytic"
					label:		qsTr("Analytic interval")
					checked:	true
				}

				RadioButton
				{
					value: 	"alphaBoot"
					label: 	qsTr("Bootstrapped interval")
				}
			}
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
					helpPage:						"toolTip/sampleSavingFreq"
				}
			}


		}
	}
}
