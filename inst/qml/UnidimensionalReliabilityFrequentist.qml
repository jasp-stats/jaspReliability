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
import JASP.Controls

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
				label:			qsTr("Confidence interval");
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
				id: 	splithalf
				name: 	"scaleSplithalf";
				label: 	qsTr("Split-half coefficient");
			}

			CheckBox { name: "averageInterItemCorrelation";	label: qsTr("Average interitem correlation")}

			RowLayout {
				CheckBox { name: "scaleMean";		label: qsTr("Mean");			id: mean		}
				CheckBox { name: "scaleVar";		label: qsTr("Variance"); 	id: variance}
				CheckBox { name: "scaleSd";			label: qsTr("SD");				id: sd			}

			}
			RadioButtonGroup
			{
				indent:		true
				enabled:	mean.checked || sd.checked || variance.checked
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
				label: 			qsTr("Confidence interval");
				defaultValue: 	95
			}

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
				name: 		"itemDeletedSplithalf";
				label: 		qsTr("Split-half coefficient (if item dropped)");
				enabled: 	splithalf.checked
			}
			
			CheckBox { name: "itemRestCorrelation";	label: qsTr("Item-rest correlation")		}
			RowLayout {
				CheckBox { name: "itemMean";			label: qsTr("Mean")							}
				CheckBox { name: "itemVar";			label: qsTr("Variance")			}
				CheckBox { name: "itemSd";			label: qsTr("SD")			}
			}

		}
	}

	Section
	{
		title: qsTr("Reverse-Scaled Items")

		VariablesForm
		{
			height: 150

			AvailableVariablesList 	{ name: "normalScaledItems"; 	title: qsTr("Normal-Scaled Items"); source: "variables" }
			AssignedVariablesList 	{ 
				name: "reverseScaledItems"; 	
				title: qsTr("Reverse-Scaled Items")
				 allowedColumns: ["scale"]
			}
		}
	}

	Section
	{
		title: qsTr("Advanced Options")
		Group
		{
			Layout.rowSpan: 2
			title: qsTr("Confidence intervals")
			RadioButtonGroup
			{
				title: qsTr("Reliability coefficients")
				name:		"intervalMethod"
				id: intervalMethod
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
			}
			RadioButtonGroup
			{
				title : qsTr("Variance and SD")
				name:		"intervalMethodVar"
				RadioButton
				{
					value: 		"chisq"
					label: 		qsTr("Chisq-based")
					checked: 	true
				}
				RadioButton
				{
					value: 	"twostep"
					label: 	qsTr("Non-parametric")

				}
			}
		}

		RadioButtonGroup
		{
			title: 		qsTr("Coefficient ω Estimation")
			name: 		"omegaEstimationMethod"
			enabled: 	omega.checked
			id: omegaEst

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
			title: 		qsTr("Coefficients")

			RadioButtonGroup
 			{
 				title: ""
 				name: "coefficientType"

 				RadioButton{ value: "unstandardized"; label: qsTr("Unstandardized"); checked: true }
 				RadioButton{ value: "standardized";	label: qsTr("Standardized")}
			}
		}

		RadioButtonGroup
		{
				title: 	qsTr("Missing Values")
				name: 	"naAction"

				RadioButton { value: "pairwise"; label: qsTr("Pairwise"); checked: true}
				RadioButton { value: "listwise"; label: qsTr("Delete listwise")}
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

		IntegerField
		{
			name: 					"hiddenScaleThreshold"
			defaultValue: 	preferencesModel.thresholdScale
			visible: 				false
		}
	}
}
