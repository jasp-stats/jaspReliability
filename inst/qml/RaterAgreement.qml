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

		AvailableVariablesList { name: "allVariablesList" }

		AssignedVariablesList
		{
			name: 			"variables"
			title: 			qsTr("Variables")
			allowedColumns: ["nominal"]
		}
		
		RadioButtonGroup
		{
			name: "dataStructure"
			title: qsTr("Data Structure")
			columns: 2
			RadioButton
			{
			name:		"ratersInColumns"
			label:		qsTr("Raters are in columns")
			checked:	true
			}

			RadioButton
			{
			name:	"ratersInRows"
			label:	qsTr("Raters are in rows")
			}
		}
	}


	Group
	{
		CheckBox
		{
			name:   	"cohensKappa"
			label:  	qsTr("Cohen's kappa")
			checked: 	true
			
			RadioButtonGroup
			{
				name:		"cohensKappaType"
				RadioButton
				{
					value:		"unweighted"
					label:		qsTr("Unweighted")
					checked:	true
				}

				RadioButton
				{
					value:	"weighted"
					label:	qsTr("Weighted")

					RadioButtonGroup
					{
						name: "weightType"
						RadioButton { value: "quadratic"; label: qsTr("Quadratic weights"); checked:	true}
						RadioButton { value: "linear"; label: qsTr("Linear weights")}
					}
				}
			}
		}
	
		CheckBox
		{
			name:   	"fleissKappa"
			label:  	qsTr("Fleiss' kappa")
			checked: 	true
		}
	
		CheckBox
		{
			name:   	"krippendorffsAlpha"
			label:  	qsTr("Krippendorff's alpha")
			checked: 	true
			
			DropDown
			{
				name:	"krippendorffsAlphaMethod"
				label:	qsTr("Method")
				values:
				[
					{label: qsTr("Nominal"),		value: "nominal"},
					{label: qsTr("Ordinal"),		value: "ordinal"},
					{label: qsTr("Interval"),		value: "interval"},
					{label: qsTr("Ratio"),			value: "ratio"}
				]
			}
			IntegerField
			{
				name: 			"krippendorffsAlphaBootstrapSamplesForCI"
				label: 			qsTr("No. of bootstrap samples for CI")
				defaultValue: 	1000
				fieldWidth: 	50
				min: 			100
				max: 			1e7
			}
			SetSeed{}
		}
	}

	CheckBox 
	{
		name: 				"ci"
		label:				qsTr("Confidence Interval")
		checked: 			true
		childrenOnSameRow: 	true

		CIField 
		{      
			name:			"ciLevel"
			label:			""
			defaultValue:	95
		}
	}
}
