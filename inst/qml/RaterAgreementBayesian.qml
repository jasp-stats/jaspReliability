//
// Copyright (C) 2013-2026 University of Amsterdam
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

		AvailableVariablesList { name: "allVariablesList"; info: qsTr("All variables available in the dataset.") }

		AssignedVariablesList
		{
			name: 			"variables"
			title: 			qsTr("Variables")
			allowedColumns: ["nominal", "ordinal"]
			info:			qsTr("Rating variables to include. Whether a variable represents a rater or a subject/item depends on the data structure setting.")
		}
	}

	RadioButtonGroup
	{
		Layout.columnSpan: 2
		name: "dataStructure"
		title: qsTr("Data Structure")
		columns: 2
		info: qsTr("Specify whether raters are arranged in columns (default) or rows in the dataset.")
		RadioButton
		{
			value:   "ratersInColumns"
			label:   qsTr("Raters are in columns")
			checked: true
			info:    qsTr("Each column is one rater and each row is one subject or item being rated.")
		}

		RadioButton
		{
			value: "ratersInRows"
			label: qsTr("Raters are in rows")
			info:  qsTr("Each row is one rater and each column is one subject or item being rated.")
		}
	}

	Group
	{
		title: qsTr("Coefficients")

		CheckBox
		{
			name:   "cohensKappa"
			label:  qsTr("Cohen's kappa")
			info:   qsTr("Chance-corrected agreement between two raters, with each rater's own category proportions. Use it when the same raters rated every subject. With more than two raters, every rater pair gets its own posterior.")

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

		CheckBox
		{
			name:  "fleissKappa"
			label: qsTr("Fleiss' kappa")
			info:  qsTr("Chance-corrected agreement with the category proportions pooled over both raters of a pair, so raters are treated as interchangeable. Use it when different raters rated different subjects. With two raters it equals Fleiss' kappa and Scott's pi; with more than two raters, every rater pair gets its own posterior.")
		}
	}

	Group
	{
		title: qsTr("Output")

		CheckBox
		{
			name:              "ci"
			label:             qsTr("Credible interval")
			checked:           true
			childrenOnSameRow: true
			info:              qsTr("Report the highest posterior density interval for each agreement coefficient.")

			CIField
			{
				name:         "ciLevel"
				label:        ""
				defaultValue: 95
				info:         qsTr("Width of the credible interval.")
			}
		}

		CheckBox
		{
			name:  "observedAndChanceAgreement"
			label: qsTr("Observed and chance agreement")
			info:  qsTr("Report the posterior means of the observed agreement and of the agreement expected by chance, from which kappa is computed.")
		}

		CheckBox
		{
			name:  "posteriorPlot"
			label: qsTr("Posterior plots")
			info:  qsTr("Plot the posterior density of each coefficient for every rater pair, with the credible interval shaded.")
		}
	}

	Section
	{
		title: qsTr("Prior")
		info:  qsTr("Prior on the cell probabilities of each rater pair's agreement table.")

		DoubleField
		{
			name:         "dirichletPriorConcentration"
			label:        qsTr("Dirichlet prior concentration")
			defaultValue: 0.001
			min:          0
			inclusive:    JASP.MaxOnly
			max:          100
			decimals:     3
			fieldWidth:   60
			info:         qsTr("Concentration of the Dirichlet prior on each cell of a rater pair's agreement table, that is, the number of prior pseudo-observations per cell. Small values let the data dominate. Larger values pull kappa towards 0, especially with many categories, because the pseudo-observations also fall into the disagreement cells.")
		}
	}

	Section
	{
		title: qsTr("Advanced Options")
		info:  qsTr("Options for the posterior sampling.")

		IntegerField
		{
			name:         "samples"
			label:        qsTr("No. of posterior samples")
			defaultValue: 5000
			fieldWidth:   60
			min:          100
			max:          10000000
			info:         qsTr("Number of draws from the posterior distribution of the cell probabilities. Higher values give more stable estimates.")
		}

		SetSeed {}
	}
}
