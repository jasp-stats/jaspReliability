//
// Copyright (C) 2013-2025 University of Amsterdam
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

Form
{
	FactorsForm
	{
		id:						factors
		name:					"factors"
		initNumberFactors:		2
		allowedColumns:			["scale"]
		keepAvailableVariables:	true
		info:					qsTr("Assign the items to their factors. The analysis requires at least two factors with at least two items each. An item may be assigned to more than one factor (cross-loading); the bi-factor model does not support cross-loadings.")
	}

	DropDown
	{
		id:				modelType
		name:			"modelType"
		label:			qsTr("Model")
		values:
		[
			{ label: qsTr("Second-order"),			value: "secondOrder"	},
			{ label: qsTr("Bi-factor"),				value: "biFactor"		},
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
				name:			"ciLevel"
				label:			qsTr("Confidence interval")
				defaultValue:	95
				info:			qsTr("Coverage of the confidence interval for the reliability coefficients.")
			}

			RadioButtonGroup
			{
				name:	"meanSdScoresMethod"
				title:	qsTr("Mean and standard deviation of")
				info:	qsTr("Whether the mean and standard deviation in the scale table are based on participants' sum scores or mean scores across items.")

				RadioButton { value: "sumScores";	label: qsTr("participants' sum scores"); checked: true	}
				RadioButton { value: "meanScores";	label: qsTr("participants' mean scores")				}
			}
		}

		Group
		{
			title: qsTr("Item statistics")

			CheckBox
			{
				name:	"itemDeletedOmegaT"
				label:	qsTr("McDonald's ω_t (if item dropped)")
				info:	qsTr("McDonald's omega_t for the remaining items when this item is dropped. The model is refit once per item.")
			}

			CheckBox
			{
				name:		"itemDeletedOmegaH"
				label:		qsTr("McDonald's ω_h (if item dropped)")
				enabled:	modelType.currentValue !== "correlated"
				info:		qsTr("McDonald's omega_h for the remaining items when this item is dropped. Not available for the correlated-factors model.")
			}

			CheckBox
			{
				name:	"itemRestCorrelation"
				label:	qsTr("Item-rest correlation")
				info:	qsTr("Pearson correlation between each item and the sum of the remaining items.")
			}
		}

		Group
		{
			title: qsTr("Model fit")

			CheckBox
			{
				name:	"fitMeasures"
				label:	qsTr("Fit measures")
				info:	qsTr("Fit measures of the factor model: chi-square, CFI, TLI, RMSEA, SRMR, AIC and BIC.")
			}

			CheckBox
			{
				name:	"standardizedLoadings"
				label:	qsTr("Standardized factor loadings")
				info:	qsTr("Table of standardized loadings of the items on their group factors. For the bi-factor model the loadings on the general factor are added as a column; for the second-order model the loadings of the group factors on the general factor are shown in a separate table.")
			}
		}
	}

	Section
	{
		title: qsTr("Reverse-Scaled Items")

		VariablesForm
		{
			height: 150

			AvailableVariablesList	{ name: "normalScaledItems"; title: qsTr("Normal-Scaled Items"); source: factors.name }
			AssignedVariablesList
			{
				name:			"reverseScaledItems"
				title:			qsTr("Reverse-Scaled Items")
				info:			qsTr("Items assigned here are recoded (reverse-scored) before the analysis.")
				allowedColumns:	["scale"]
			}
		}
	}

	Section
	{
		title: qsTr("Advanced Options")

		RadioButtonGroup
		{
			name:	"intervalMethod"
			title:	qsTr("Confidence intervals")
			info:	qsTr("Analytic intervals are the Wald-type intervals of the factor model. Because omega is a ratio, these intervals are symmetric and can extend beyond zero or one; the bootstrapped interval resamples the data and refits the model, which is slower but does not make that assumption.")

			RadioButton
			{
				value:		"analytic"
				label:		qsTr("Analytic interval")
				checked:	true
			}

			RadioButton
			{
				value:	"bootstrapped"
				label:	qsTr("Bootstrapped interval")

				IntegerField
				{
					name:			"bootstrapSamples"
					label:			qsTr("No. of bootstrap samples")
					defaultValue:	1000
					fieldWidth:		50
					min:			100
					max:			1e5
					info:			qsTr("Number of bootstrap replications. The factor model is refit for every replication, so higher values take proportionally longer.")
				}
			}
		}

		RadioButtonGroup
		{
			title:	qsTr("Missing Values")
			name:	"naAction"
			info:	qsTr("Full information maximum likelihood uses all available observations; listwise deletion removes any row with a missing value.")

			RadioButton { value: "fiml";		label: qsTr("Full information maximum likelihood"); checked: true	}
			RadioButton { value: "listwise";	label: qsTr("Exclude cases listwise")								}
		}

		Group
		{
			title:	qsTr("Repeatability")
			info:	qsTr("Set a random seed to reproduce the same bootstrap results across runs.")

			CheckBox
			{
				name:				"setSeed"
				label:				qsTr("Set seed")
				childrenOnSameRow:	true

				IntegerField
				{
					name:			"seed"
					label:			""
					defaultValue:	1234
					fieldWidth:		100
					min:			1
					max:			1e9
				}
			}
		}

		Group
		{
			title: qsTr("Samples")

			CheckBox
			{
				name:		"samplesSavingDisabled"
				label:		qsTr("Disable saving samples")
				checked:	false
				info:		qsTr("When checked, bootstrap samples are not stored in the output file. Reduces file size but changing the confidence level requires resampling.")
			}
		}
	}
}
