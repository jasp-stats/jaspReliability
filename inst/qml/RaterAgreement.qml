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

Form
{
	VariablesForm
	{
		height: 300

		AvailableVariablesList { name: "allVariablesList"; info: qsTr("All variables available in the dataset.") }

		AssignedVariablesList
		{
			name: 			"variables"
			title: 			qsTr("Variables")
			allowedColumns: ["nominal", "ordinal", "scale"]
			info:			qsTr("Rating variables to include. Whether a variable represents a rater or a subject/item depends on the data structure setting.")
		}

		RadioButtonGroup
		{
			name: "dataStructure"
			title: qsTr("Data Structure")
			columns: 2
			info: qsTr("Specify whether raters are arranged in columns (default) or rows in the dataset.")
			RadioButton
			{
				name:    "ratersInColumns"
				label:   qsTr("Raters are in columns")
				checked: true
				info:    qsTr("Each column is one rater and each row is one subject or item being rated.")
			}

			RadioButton
			{
				name:  "ratersInRows"
				label: qsTr("Raters are in rows")
				info:  qsTr("Each row is one rater and each column is one subject or item being rated.")
			}
		}
	}




	CheckBox
	{
		name:   "cohensKappa"
		label:  qsTr("Cohen's kappa")
		info:   qsTr("Measures agreement between exactly two raters. When more than two raters are entered, all pairwise combinations are computed.")

		RadioButtonGroup
		{
			name: "cohensKappaType"
			info: qsTr("Unweighted kappa treats all disagreements equally. Weighted kappa accounts for the degree of disagreement and requires ordinal ratings.")
			RadioButton
			{
				value:   "unweighted"
				label:   qsTr("Unweighted")
				checked: true
				info:    qsTr("All disagreements are treated as equal regardless of their magnitude.")
			}

			RadioButton
			{
				value: "weighted"
				label: qsTr("Weighted")
				info:  qsTr("Disagreements are penalised according to their magnitude. Requires ordinal ratings.")

				RadioButtonGroup
				{
					name: "weightType"
					info: qsTr("Weighting scheme applied to disagreements between ordinal categories.")
					RadioButton { value: "quadratic"; label: qsTr("Quadratic weights"); checked: true; info: qsTr("Penalises larger disagreements quadratically; sensitive to large discrepancies.") }
					RadioButton { value: "linear";    label: qsTr("Linear weights");    info: qsTr("Penalises disagreements proportionally to their size.") }
				}
			}
		}
	}

Group
{
	CheckBox
	{
		id:    krippAlphaOpt
		name:  "krippendorffsAlpha"
		label: qsTr("Krippendorff's alpha")
		info:  qsTr("Measures agreement among two or more raters. Applicable to nominal, ordinal, interval, or ratio data.")

		DropDown
		{
			name:   "krippendorffsAlphaMethod"
			label:  qsTr("Method")
			info:   qsTr("Level of measurement determines how disagreements are quantified in the alpha calculation.")
			values:
			[
				{ label: qsTr("Nominal"),  value: "nominal"  },
				{ label: qsTr("Ordinal"),  value: "ordinal"  },
				{ label: qsTr("Interval"), value: "interval" },
				{ label: qsTr("Ratio"),    value: "ratio"    }
			]
		}
	}

	CheckBox
	{
		id:    kendallWOpt
		name:  "kendallW"
		label: qsTr("Kendall's W")
		info:  qsTr("Measures concordance among multiple raters on ordinal rankings. Ranges from 0 (no agreement) to 1 (perfect agreement).")

		CheckBox
		{
			name:    "correctForTies"
			label:   qsTr("Correct for ties")
			checked: true
			info:    qsTr("Apply a correction to Kendall's W when tied ranks are present in the data. Recommended whenever ties can occur.")
		}
	}
}

	CheckBox
	{
		name:  "fleissKappa"
		label: qsTr("Fleiss' kappa")
		info:  qsTr("Measures agreement among two or more raters on nominal categories. Generalises Cohen's kappa to multiple raters.")
	}






	Section
	{
		title: qsTr("Advanced Options")
		info:  qsTr("Options for confidence intervals and bootstrap settings.")

		CheckBox
		{
			name:              "ci"
			label:             qsTr("Confidence interval")
			checked:           true
			childrenOnSameRow: true
			id: ciOpt
			info:              qsTr("Report a confidence interval for each agreement coefficient. Cohen's kappa and Fleiss' kappa use asymptotic CIs; Krippendorff's alpha and Kendall's W use bootstrap CIs.")

			CIField
			{
				name:         "ciLevel"
				label:        ""
				defaultValue: 95
				enabled:      ciOpt.checked
				info:         qsTr("Width of the confidence interval.")
			}
		}

		IntegerField
		{
			name:         "bootstrapSamples"
			label:        qsTr("No. of bootstrap samples for CI")
			defaultValue: 1000
			fieldWidth:   50
			min:          100
			max:          10000000
			enabled:      ciOpt.checked && (krippAlphaOpt.checked || kendallWOpt.checked)
			info:         qsTr("Number of bootstrap replications used to compute confidence intervals. Higher values give more stable estimates. Only applies to Krippendorff's alpha and Kendall's W, whose CIs are bootstrap-based.")
		}

		SetSeed {
			enabled: ciOpt.checked && (krippAlphaOpt.checked || kendallWOpt.checked)
		}
	}
}
