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
			allowedColumns: ["nominal", "nominalText", "ordinal"]
			id: variables
		}
	}

	IntegerField
	{
		name: "responseCategories"
		label: qsTr("<b>Number of response categories</b>")
		defaultValue: 2
		id:responseCategories
	}

	Group{}
	Group
	{
		Layout.rowSpan: 2
		title: qsTr("Split-test methods")

		CheckBox
		{
			name:   	"thorndike"
			label:  	qsTr("Thorndike method")
			
		}
		CheckBox
		{
			name:   	"feldt"
			label:  	qsTr("Feldt method")
			IntegerField
			{
				name: "feldtNumberOfSplits"
				label: qsTr("Number of splits")
				min: 2
				max: variables.count > 1 ? variables.count : 2
				defaultValue: 2
			}
		}
		CheckBox
		{
			name:   	"mollenkopfFeldt"
			label:  	qsTr("Mollenkopf-Feldt method")
			IntegerField
			{
				name: "mollenkopfFeldtNumberOfSplits"
				label: qsTr("Number of splits")
				min: 2
				max: variables.count > 1 ? variables.count : 2
				defaultValue: 2
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
			label:  	qsTr("ANOVA method")
		}
		
		CheckBox
		{
			name:   	"irt"
			label:  	qsTr("IRT method")
		}
	}

	Group
	{
		title: qsTr("Binomial methods")

		CheckBox
		{
			name:   	"lord"
			label:  	qsTr("Lord method")
			enabled: responseCategories.value == 2
		}
		CheckBox
		{
			name:   	"keats"
			label:  	qsTr("Keats method")
			enabled: responseCategories.value == 2
		}
		CheckBox
		{
			name:   	"lord2"
			label:  	qsTr("Lord-2 method")
			enabled: responseCategories.value == 2
			IntegerField
			{
				name: "lord2NumberOfSplits"
				label: qsTr("Number of splits")
				min: 2
				max: variables.count > 1 ? (variables.count/2) : 2
				defaultValue: 2
			}
		}
	}

	Section
	{
		title: qsTr("Options")

		CheckBox 
		{
			name: "userReliability"
			label: qsTr("Use own reliability")
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
		
	}

	Section
	{
		title: qsTr("Plots")

		CheckBox
		{
			name: "histogramCounts"
			label: qsTr("Histogram of score counts")
		}	

		CheckBox
		{
			name: "pointPlots"
			label: qsTr("Points plot per method")
		}

		CheckBox
		{
			name: "combinedPointPlot"
			label: qsTr("Combined point plot")
		}	
	}


}
