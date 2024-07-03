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
import QtQuick.Layouts 1.3
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
			allowedColumns: ["nominal", "ordinal"]
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
				lord2NumberOfSplits.values = newValues

			}
		}
	}

	Group
	{
		Layout.rowSpan: 2
		title: qsTr("Split-test methods")

		CheckBox
		{
			name:   	"thorndike"
			label:  	qsTr("Thorndike")
			
		}
		CheckBox
		{
			name:   	"feldt"
			label:  	qsTr("Feldt")
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
				max: 10
				defaultValue: 3
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
		title: qsTr("Binomial methods")

		CheckBox
		{
			name:   	"lord"
			label:  	qsTr("Lord")
		}
		CheckBox
		{
			name:   	"keats"
			label:  	qsTr("Keats")
		}
		CheckBox
		{
			name:   	"lord2"
			label:  	qsTr("Lord's compound")
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
			label: qsTr("Minimum observations per group") 
			min: 1
			defaultValue: 10
		}

		CheckBox 
		{
			name: "hideTable"
			label: qsTr("Hide table")
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
			label: qsTr("Combined plot")
		}	
	}


}
