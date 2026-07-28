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
			info:			qsTr("Rating variables to include. Each variable is one rater; each row is a subject being rated.")
		}
	}

	Group
	{

		RadioButtonGroup
		{
			title: qsTr("Each subject is rated by...")
			name: "iccType"
			info:  qsTr("Determines the ICC model. ICC(1): each subject rated by a different random rater. ICC(2): all subjects rated by the same random sample of raters. ICC(3): all subjects rated by the same fixed raters. See Shrout & Fleiss (1979).")

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
			info:		qsTr("When ratings are averaged across raters, the k-rater version of the selected ICC is computed (e.g. ICC(2,k) instead of ICC(2,1)) per Shrout & Fleiss (1979). Averaging generally raises the reliability estimate.")
		}
	}

	CheckBox
	{
		name: 				"ci"
		label:				qsTr("Confidence Interval")
		checked: 			true
		childrenOnSameRow: 	true
		info:				qsTr("Report a confidence interval for the ICC estimate based on the F-distribution.")
		CIField
		{
			name:			"ciLevel"
			label:			""
			defaultValue:	95
		}
	}
}
