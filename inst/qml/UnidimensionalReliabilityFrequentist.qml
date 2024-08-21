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
			allowedColumns: ["scale"]
		}
	}

	Section
	{
		title: 	qsTr("Analysis")

		Group
		{
			title: qsTr("Scale Statistics")


			CIField
			{
				name:			"ciLevel";
				label:			qsTr("Confidence Interval");
				defaultValue:	95;
			}
			

			CheckBox
			{
				id:			omega
				name:		"scaleOmega"
				label:		qsTr("Coefficient ω")
				checked:	true
			}

			CheckBox
			{
				id: 	alpha
				name: 	"scaleAlpha";
				label: 	qsTr("Coefficient α");
			}

			CheckBox
			{
				id: 	lambda2
				name: 	"scaleLambda2";
				label: 	qsTr("Guttman's λ2");
			}

			CheckBox
			{
				id: 	lambda6
				name: 	"scaleLambda6";
				label: 	qsTr("Guttman's λ6");
			}

			CheckBox
			{
				id: 	glb
				name: 	"scaleGreatestLowerBound";
				label: 	qsTr("Greatest lower bound");
			}

			CheckBox { name: "averageInterItemCorrelation";	label: qsTr("Average interitem correlation")}

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

			CheckBox
			{
				name: 		"itemDeletedOmega";
				label: 		qsTr("Coefficient ω  (if item dropped)");
				enabled: 	omega.checked
			}

			CheckBox
			{
				name: 		"itemDeletedAlpha";
				label: 		qsTr("Coefficient α (if item dropped)");
				enabled: 	alpha.checked
			}

			CheckBox
			{
				name: 		"itemDeletedLambda2";
				label: 		qsTr("Guttman's λ2 (if item dropped)");
				enabled: 	lambda2.checked
			}

			CheckBox
			{
				name: 		"itemDeletedLambda6";
				label: 		qsTr("Guttman's λ6 (if item dropped)");
				enabled: 	lambda6.checked
			}

			CheckBox
			{
				name: 		"itemDeletedGreatestLowerBound";
				label: 		qsTr("Greatest lower bound (if item dropped)");
				enabled: 	glb.checked
			}

			CheckBox { name: "itemRestCorrelation";	label: qsTr("Item-rest correlation")		}
			CheckBox { name: "itemMean";			label: qsTr("Mean")							}
			CheckBox { name: "itemSd";			label: qsTr("Standard deviation")			}
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
				name: 	"naAction"

				RadioButton { value: "pairwise"; label: qsTr("Exclude cases pairwise"); checked: true}
				RadioButton { value: "listwise"; label: qsTr("Exclude cases listwise")}
		}

	Group
	{
			title: qsTr("Bootstrap")

			IntegerField
			{
				name: 			"bootstrapSamples"
				label: 			qsTr("No. of bootstrap samples")
				defaultValue: 	1000
				fieldWidth: 	50
				min: 			100
				max: 			1e7
			}

			RadioButtonGroup
			{
				title:		""
				name:		"bootstrapType"

				RadioButton {value: "nonParametric"; label: qsTr("Non-parametric bootstrap"); checked: true}
				RadioButton {value: "parametric"; label: qsTr("Parametric bootstrap")}
			}

		}

		RadioButtonGroup
		{
			title: 		qsTr("Coefficient ω Estimation")
			name: 		"omegaEstimationMethod"
			enabled: 	omega.checked

			RadioButton
			{
				value: 		"cfa"
				label: 		qsTr("CFA")
				checked: 	true

				CheckBox
				{
					name: 		"omegaFitMeasures"
					label: 		qsTr("Single factor model fit")
				}

				RadioButtonGroup
				{
					title:		qsTr("Interval")
					name:		"omegaIntervalMethod"

					RadioButton
					{

						value: 		"analytic"
						label: 		qsTr("Analytic interval")
						checked: 	true
					}

					RadioButton
					{
						value: 	"bootstrapped"
						label: 	qsTr("Bootstrapped interval")
					}
				}
			}
			
			RadioButton
			{
				value: "pfa"
				label: qsTr("PFA")
			}

			CheckBox
			{
				name:		"standardizedLoadings"
				label:		qsTr("Standardized factor loadings");
			}
		}

		Group
		{
			title: 		qsTr("Coefficient α Estimation")
			enabled: 	alpha.checked

			RadioButtonGroup
			{
				title:	""
				name: 	"alphaType"

				RadioButton
				{
					value:		"unstandardized"
					label:		qsTr("Unstandardized")
					checked:	true
				}

				RadioButton
				{
					value:	"standardized"
					label:	qsTr("Standardized")
				}
			}

			RadioButtonGroup
			{
				title:		qsTr("Interval")
				name:		"alphaIntervalMethod"
				enabled:	interval.checked

				RadioButton
				{
					value:		"analytic"
					label:		qsTr("Analytic interval")
					checked:	true
				}

				RadioButton
				{
					value: 	"bootstrapped"
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
					name:				"samplesSavingDisabled"
					label:				qsTr("Do not save samples")
					checked:			false
				}
				HelpButton
				{
					toolTip: 			qsTr("Click to learn more about saving the samples.")
					helpPage:			"toolTip/sampleSavingFreq"
				}
			}


		}
	}
}
