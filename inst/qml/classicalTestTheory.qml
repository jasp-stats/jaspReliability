//
// Copyright (C) 2013-2018 University of Amsterdam
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
import JASP
import JASP.Controls

Form
{
	info: qsTr("Classical Test Theory (CTT) is a traditional framework for analyzing test scores. It focuses on observed scores as the primary measure of an individual's ability and assumes that these scores consist of a true score and measurement error. CTT provides simple reliability indices, such as Cronbach's alpha, to assess the consistency of test scores but doesn't account for item characteristics or the latent ability of individuals, limiting its ability to address more complex assessments.")

	VariablesForm
	{
		AvailableVariablesList
		{
			name: 				"variablesList"
			info:				qsTr("Select the variables to be analyzed.")
		}

		AssignedVariablesList
		{
			name: 				"items"
			title: 				qsTr("Item Scores")
			allowedColumns: 	["scale"]
			info:				qsTr("Select the variables representing item scores for the analysis.")
		}
	}

	Group
	{
		title:					qsTr("Data")

		CheckBox
		{
			name:				"customMaxScore"
			text:				qsTr("First row is maximum item score")
			info:				qsTr("Check this if the first row in your data represents maximum item scores.")
		}
	}

	Column
	{
		spacing:				20 * preferencesModel.uiScale

		Group
		{
			title:				qsTr("Display")

			CheckBox
			{
				name:			"explanatoryText"
				text:			qsTr("Explanatory text")
				info:			qsTr("Include explanatory text in the report.")
			}
		}

		Group
		{
			title:				qsTr("Tables")

			CheckBox
			{
				name:			"tableDescriptives"
				text:			qsTr("Descriptive statistics")
				info:			qsTr("Generate a table of descriptive statistics for the test (i.e., sum) scores.")
			}

			CheckBox
			{
				name:			"tableCronbachsAlpha"
				text:			qsTr("Test reliability")
				info:			qsTr("Calculate test reliability using Cronbach's Alpha.")

				CIField
				{
					name:		"tableCronbachsAlphaCI"
					text:		qsTr("Confidence interval")
					info:		qsTr("Specify the confidence interval for Cronbach's Alpha.")
				}
			}

			CheckBox
			{
				name:			"tableItemStatistics"
				text:			qsTr("Item information")
				info:			qsTr("Generate a table of item-level statistics.")
			}
		}

		Group
		{
			title:				qsTr("Plots")

			CheckBox
			{
				name:			"plotHistogram"
				text:			qsTr("Histogram of sum scores")
				info:			qsTr("Create a histogram of the test (i.e., sum) scores.")
			}

			CheckBox
			{
				name:			"plotItems"
				text:			qsTr("Item difficulty and discrimination")
				info:			qsTr("Generate a plot showing the item difficulty and discrimination parameters.")
			}

			CheckBox
			{
				name:			"plotCorrelationHeatmap"
				text:			qsTr("Item correlations")
				info:			qsTr("Generate a heatmap of correlations between item scores.")

				CheckBox
				{
					name:		"plotCorrelationHeatmapShowValues"
					text:		qsTr("Display values")
					info:		qsTr("Display correlation values in the heatmap.")
				}
			}
		}
	}
}
