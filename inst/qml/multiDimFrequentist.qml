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
import QtQuick.Layouts  1.3
import JASP.Controls 	1.0
import JASP.Theme		1.0
import JASP.Widgets 	1.0

Form
{
	FactorsForm
	{
		id: factors
		name: "factors"
		initNumberFactors: 2
		allowAll: true
	}

	Section
	{
		title: 	qsTr("Analysis")
		Group
		{
			title: qsTr("Scale Statistics")

			CheckBox
			{
				name: 				"intervalOn"
				label:				qsTr("Confidence Interval")
				checked: 			true
				childrenOnSameRow: 	true
				id:					interval

				CIField
				{
					name:			"confidenceIntervalValue";
					label:			"";
					defaultValue:	95;
				}
			}

			CheckBox
			{
				id:			omegah
				name:		"omegaHScale"
				label:		qsTr("McDonald's ω_h")
			}
			CheckBox
			{
				id:			omegat
				name:		"omegaTScale"
				label:		qsTr("McDonald's ω_t")
			}
			CheckBox { name: "averageInterItemCor";	label: qsTr("Average interitem correlation")}

			RowLayout {
				CheckBox { name: "meanScale";	label: qsTr("Mean");	id: mean}
				CheckBox { name: "sdScale";		label: qsTr("SD");		id: sd}

			}
			RadioButtonGroup
			{
				indent:		true
				enabled:	mean.checked || sd.checked
				title:		qsTr("")
				name:		"scoresMethod"

				RadioButton { value: "sumScores";	label: qsTr("of participants' sum scores"); checked: true}
				RadioButton { value: "meanScores";	label: qsTr("of participants' mean scores")}
			}
		}

		Group
		{
			title: qsTr("Individual Item Statistics")
			CheckBox { name: "itemRestCor";	label: qsTr("Item-rest correlation")		}
			CheckBox { name: "meanItem";	label: qsTr("Mean")							}
			CheckBox { name: "sdItem";		label: qsTr("Standard deviation")			}
		}
	}
	Section
	{
		title: qsTr("Reverse-Scaled Items")

		VariablesForm
		{
			height: 150

			AvailableVariablesList 	{ name: "normalScaledItems"; 	title: qsTr("Normal-Scaled Items"); source: factors.name }
			AssignedVariablesList 	{ name: "reverseScaledItems"; 	title: qsTr("Reverse-Scaled Items") }
		}
	}

	Section
	{
		title: qsTr("Advanced Options")

		RadioButtonGroup
		{
				title: 	qsTr("Missing Values")
				name: 	"missingValues"

				RadioButton { value: "excludeCasesPairwise"; label: qsTr("Exclude cases pairwise"); checked: true}
				RadioButton { value: "excludeCasesListwise"; label: qsTr("Exclude cases listwise")}
		}


		RadioButtonGroup
		{
			title:	qsTr("Factor Method")
			name:	"factorMethod"
			enabled: interval.checked
			RadioButton{
				value:		"cfa"
				label:		qsTr("CFA")
				checked:	true
				CheckBox
				{
					name: "fitmeasures"
					label: qsTr("Fit measures")
				}
			}
			RadioButton{ value: "efa"; label: qsTr("EFA"); checked: false}
		}

		RadioButtonGroup
		{
			title:	qsTr("Interval")
			name:	"intervalType"
			enabled: interval.checked
			RadioButton{
				value: "wald"
				label: qsTr("Wald-type")
				checked: true
			}
			RadioButton
			{
				value:	"boot"
				label:	qsTr("Bootstrap:")

				IntegerField
				{
					name: "bootNo"
					label: qsTr("Samples:")
					defaultValue: 1000
					fieldWidth: 50
					min: 100
					max: 1e6
				}
				CheckBox
				{
					name: 				"setSeed"
					label: 				qsTr("Set seed")
					childrenOnSameRow: 	true

					IntegerField
					{
						name: 			"seed"
						label: 			""
						defaultValue: 	1234
						fieldWidth: 	50
						min: 			1
						max: 			1e9
					}
				}
			}
		}
		Group
		{
			title: qsTr("Samples")

			RowLayout
			{
				CheckBox
				{
					name:		"disableSampleSave"
					label:		qsTr("Disable saving samples")
					checked:	false
				}
				HelpButton
				{
					toolTip: 	qsTr("Click to learn more about saving the samples.")
					helpPage:	"toolTip/sampleSavingFreq"
				}
			}
		}
	}
}


