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
			allowedColumns: ["scale"]
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
			name:   	"averagedRating"
			label:  	qsTr("Ratings are averaged")
			checked: 	false
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
