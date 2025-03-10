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

		CheckBox
		{
			name:   	"thorndike"
			label:  	qsTr("Thorndike")
			id: thorndike
		}
		CheckBox
		{
			name:   	"feldt"
			label:  	qsTr("Feldt")
			id: feldt
			DropDown 
			{
				name: "feldtNumberOfSplits"
				label: qsTr("Number of splits")
				values: []
				id: feldtNumberOfSplits
			}
		}
		CheckBox
		{
			name:   	"mollenkopfFeldt"
			label:  	qsTr("Mollenkopf-Feldt")
			id: mollenkopfFeldt
			DropDown 
			{
				name: "mollenkopfFeldtNumberOfSplits"
				label: qsTr("Number of splits")
				values: []
				id: mollenkopfFeldtNumberOfSplits
			}
			IntegerField
			{
				name: "mollenkopfFeldtPolyDegree"
				label: qsTr("Degree of polynomial")
				min: 2
				max: 8
				defaultValue: 2
			}

		}
	}

	Group
	{
		CheckBox
		{
			name:   	"anova"
			label:  	qsTr("ANOVA")
		}
		
		CheckBox
		{
			name:   	"irt"
			label:  	qsTr("IRT")
		}
	}

	Group
	{
		title: qsTr("Binomial Methods")

		CheckBox
		{
			enabled: !variables.columnsTypes.includes("ordinal")
			name:   	"lord"
			label:  	qsTr("Lord")
		}
		CheckBox
		{
			enabled: !variables.columnsTypes.includes("ordinal")
			name:   	"keats"
			label:  	qsTr("Keats")
		}
		CheckBox
		{
			enabled: !variables.columnsTypes.includes("ordinal")
			name:   	"lord2"
			label:  	qsTr("Lord's compound")
			id: lord2
			DropDown 
			{
				name: "lord2NumberOfSplits"
				label: qsTr("Number of splits")
				values: []
				id: lord2NumberOfSplits
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
		}

		CheckBox 
		{
			name: "hideTable"
			label: qsTr("Hide SEM table")
		}
		
	}

	Section
	{
		title: qsTr("Plots")

		CheckBox
		{
			name: "histogramCounts"
			label: qsTr("Histogram of sum score counts")
		}	

		CheckBox
		{
			name: "pointPlots"
			label: qsTr("Plot per method")
		}

		CheckBox
		{
			name: "combinedPointPlot"
			label: qsTr("Plot all SEMs")
		}

		Group {
			CheckBox
			{
				name: "sumScoreCiPlots"
				label: qsTr("Sum score plots")
				id: sumScoreCiPlots
				childrenOnSameRow: true
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
					DoubleField
					{
						name: "sumScoreCiPlotsCutoffValue"
					}
				}
			}
		}


	}


}
