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
			allowedColumns: ["nominal", "ordinal"]
			allowTypeChange: true
			maxLevels: 50
			minNumericLevels: 2
			id: variables
			info: qsTr("Items to include in the SEM analysis. Must be ordinally or nominally scaled (dichotomous or polytomous).")
			onCountChanged:
			{
				var newValues = []
				for (var i = 2; i <= variables.count; i++)
				{
					if (variables.count % i == 0)
					{
						newValues.push(i)
					}
				}
				feldtNumberOfSplits.values = newValues
				mollenkopfFeldtNumberOfSplits.values = newValues

				var newValuesLord = []
				for (var ii = 2; ii < variables.count; ii++)
				{
					if (variables.count % ii == 0)
					{
						newValuesLord.push(ii)
					}
				}

				lord2NumberOfSplits.values = newValuesLord
			}

		}
	}

	Group
	{
		Layout.rowSpan: 2
		title: qsTr("Split-Test Methods")
		info:  qsTr("CTT-based methods that split the test into parts to estimate the conditional error variance per sum score.")

		CheckBox
		{
			name:   	"thorndike"
			label:  	qsTr("Thorndike")
			id: thorndike
			info:		qsTr("Splits items into odd- and even-numbered halves. Reorder variables in the list to change the split.")
		}
		CheckBox
		{
			name:   	"feldt"
			label:  	qsTr("Feldt")
			id: feldt
			info:		qsTr("Splits items into multiple equal parts to estimate error variance per score group.")
			DropDown
			{
				name: "feldtNumberOfSplits"
				label: qsTr("Number of splits")
				values: [2, 3, 4, 5, 6, 7, 8, 9, 10]
				id: feldtNumberOfSplits
				info: qsTr("How many equal parts to split the test into. Must be a divisor of the number of items.")
			}
		}
		CheckBox
		{
			name:   	"mollenkopfFeldt"
			label:  	qsTr("Mollenkopf-Feldt")
			id: mollenkopfFeldt
			info:		qsTr("Splits items into parts and models score differences with polynomial regression to estimate the conditional SEM.")
			DropDown
			{
				name: "mollenkopfFeldtNumberOfSplits"
				label: qsTr("Number of splits")
				values: [2, 3, 4, 5, 6, 7, 8, 9, 10]
				id: mollenkopfFeldtNumberOfSplits
				info: qsTr("How many equal parts to split the test into. Must be a divisor of the number of items.")
			}
			IntegerField
			{
				name: "mollenkopfFeldtPolyDegree"
				label: qsTr("Degree of polynomial")
				min: 2
				max: 8
				defaultValue: 2
				info: qsTr("Degree of the polynomial used to predict score differences, e.g., 3 fits Y = X + X² + X³.")
			}

		}
	}

	Group
	{
		CheckBox
		{
			name:   	"anova"
			label:  	qsTr("ANOVA")
			info:		qsTr("Estimates conditional SEM from a repeated-measures ANOVA using the ICC(3,k) approach (Emons, 2023).")
		}

		CheckBox
		{
			name:   	"irt"
			label:  	qsTr("IRT")
			info:		qsTr("IRT-based SEM using the 2PLM for dichotomous items or the graded response model for polytomous items. Assumes a single latent variable.")
		}
	}

	Group
	{
		title: qsTr("Binomial Methods")
		info:  qsTr("Methods based on the binomial error model. Only available for dichotomous (non-ordinal) items.")

		CheckBox
		{
			enabled: !variables.columnsTypes.includes("ordinal")
			name:   	"lord"
			label:  	qsTr("Lord")
			info:		qsTr("Lord's binomial method based on the number of correct and incorrect responses.")
		}
		CheckBox
		{
			enabled: !variables.columnsTypes.includes("ordinal")
			name:   	"keats"
			label:  	qsTr("Keats")
			info:		qsTr("Keats' correction of Lord's method; uses a reliability coefficient to reduce bias.")
		}
		CheckBox
		{
			enabled: !variables.columnsTypes.includes("ordinal")
			name:   	"lord2"
			label:  	qsTr("Lord generalized")
			id: lord2
			info:		qsTr("Lord's method generalised to multiple test parts.")
			DropDown
			{
				name: "lord2NumberOfSplits"
				label: qsTr("Number of splits")
				values: [2, 3, 4, 5, 6, 7, 8, 9, 10]
				id: lord2NumberOfSplits
				info: qsTr("How many parts to split the test into. Must be a divisor of the number of items (excluding the last value).")
			}
		}
	}

	Section
	{
		title: qsTr("Options")

		CheckBox
		{
			name: "sumScoreCiTable"
			label: qsTr("Sum score table")
			childrenOnSameRow: true
			info: qsTr("Table of sum scores with normal-theory confidence intervals computed from the estimated SEM.")
			CIField
			{
				name:			"ciLevelTable";
				label:			qsTr("CI:");
				defaultValue:	95;
			}
		}

		CheckBox
		{
			name: "userReliability"
			label: qsTr("User defined reliability")
			childrenOnSameRow: true
			info: qsTr("Override the reliability estimate used for the unconditional SEM and the Keats method.")

			DoubleField
			{
				name: "reliabilityValue"
				label: ""
				max: 1
				defaultValue: .5
			}
		}

		IntegerField
		{
			name: "minimumGroupSize"
			label: qsTr("Minimum number of observations per score group")
			min: 1
			defaultValue: 20
			info: qsTr("Score groups with fewer observations are merged with adjacent groups before estimating the SEM (default: 20). Applies to the Thorndike, Feldt, ANOVA, and Lord generalized methods.")
		}

		CheckBox
		{
			name: "hideTable"
			label: qsTr("Hide SEM table")
			info: qsTr("Hide the main SEM table showing conditional error estimates per sum score and method.")
		}

	}

	Section
	{
		title: qsTr("Plots")

		CheckBox
		{
			name: "histogramCounts"
			label: qsTr("Histogram of sum score counts")
			info: qsTr("Histogram of the number of respondents per sum score group.")
		}

		CheckBox
		{
			name: "pointPlots"
			label: qsTr("Plot per method")
			info: qsTr("Separate point plot per method showing the conditional SEM values across sum scores.")
		}

		CheckBox
		{
			name: "combinedPointPlot"
			label: qsTr("Plot all SEMs")
			info: qsTr("Combined plot displaying all methods' SEM values in a single panel for easy comparison.")
		}

		Group {
			CheckBox
			{
				name: "sumScoreCiPlots"
				label: qsTr("Sum score plots")
				id: sumScoreCiPlots
				childrenOnSameRow: true
				info: qsTr("Plot sum scores with confidence interval bands per method. Intervals are normal-theory CIs using the estimated SEM.")
				CIField
				{
					name:			"ciLevelPlots";
					label:			qsTr("CI");
					defaultValue:	95;
				}
			}
			RowLayout
			{
				Label
				{
						text: ""
						Layout.leftMargin: 2.5 // Adjust the value as needed for the desired indentation
				}
				CheckBox
				{
					name: "sumScoreCiPlotsCutoff"
					label: qsTr("Display cutoff score")
					enabled: sumScoreCiPlots.checked
					childrenOnSameRow: true
					info: qsTr("Add a horizontal reference line at the specified cut score.")
					DoubleField
					{
						name: "sumScoreCiPlotsCutoffValue"
					}
				}
			}
		}


	}


}
