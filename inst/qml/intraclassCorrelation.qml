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

	Group
	{

		RadioButtonGroup
		{
			title: qsTr("Each subject is rated by...")
			name: "iccType"

			RadioButton
			{
				value:   	"icc1"
				label:  	qsTr("a different rater (randomly selected)")
				checked: 	true
			}

			RadioButton
			{
				value:   	"icc2"
				label:  	qsTr("the same set of randomly selected raters/tests")
				checked: 	false
			}

			RadioButton
			{
				value:   	"icc3"
				label:  	qsTr("the same fixed set of raters/tests")
				checked: 	false
			}
		}

		CheckBox
		{
			name:   	"iccRatingAverage"
			label:  	qsTr("Ratings are averaged")
			checked: 	false
		}
	}

	CheckBox 
	{
		name: 				"intervalOn"
		label:				qsTr("Confidence Interval")
		checked: 			true
		childrenOnSameRow: 	true

		CIField 
		{      
			name: 		"confidenceIntervalValue";   
			label: 		"";
			defaultValue: 95;
		}
	}
	
	Section
	{
		title: qsTr("Plots")

		VariablesForm
		{
			preferredHeight: jaspTheme.smallDefaultVariablesFormHeight

			AvailableVariablesList
			{
				name:	"availablePlottingComponents"
				title:	qsTr("Plot components")
				source:	"variables"
			}

			AssignedPairsVariablesList
			{
				id:					pairs
				name:				"pairs"
				title:				qsTr("Measurement Pairs")
				suggestedColumns: 	["scale", "ordinal"]
			}
		}
		
		CheckBox
		{
			name:   	"descriptivesBlandAltman"
			label:  	qsTr("Bland-Altman plot")
			
			CheckBox 
			{
				name: 				"ciDisplay"
				label:				qsTr("Confidence Interval")
				childrenOnSameRow: 	true

				CIField 
				{      
					name: 		"ciValue";   
					label: 		"";
					defaultValue: 95;
				}
				
				CheckBox 
				{
					name: 		"ciShading"
					label:		qsTr("Shading")
				}
			}
			
			CheckBox 
			{
					name: 		"descriptivesBlandAltmanTable"
					label:		qsTr("Bland-Altman table")
			}
		}
	}
}
