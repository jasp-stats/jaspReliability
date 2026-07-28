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

		AssignedPairsVariablesList
		{
			id:					pairs
			name:				"pairs"
			title:				qsTr("Measurement Pairs")
			allowedColumns: 	["scale"]
			info:				qsTr("Pairs of scale variables to compare. Each pair produces a Bland-Altman plot of the mean versus the difference between measurements.")
		}
	}

	CheckBox
	{
		name: 				"ci"
		label:				qsTr("Confidence Interval")
		info:				qsTr("Display confidence intervals around the mean difference and its upper and lower limits of agreement in the plot.")

		CIField
		{
			name:			"ciLevel";
			label:			"";
			defaultValue:	95;
		}

		CheckBox
		{
			name: 		"ciShading"
			label:		qsTr("Shading")
			info:		qsTr("Shade the confidence interval regions around the mean difference and limits of agreement.")

			CheckBox
			{
				name: 	"ciShadingWithColour"
				label:	qsTr("Use colour")
				info:	qsTr("Use colour (rather than grey) to fill the shaded confidence regions.")
			}
		}


	}

	CheckBox
	{
			name: 		"blandAltmanTable"
			label:		qsTr("Bland-Altman table")
			info:		qsTr("Table showing the mean difference, the upper and lower limits of agreement, and their confidence intervals.")
	}
}
