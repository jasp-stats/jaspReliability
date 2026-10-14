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
			allowedColumns: ["scale"]
			info:			qsTr("Items/variables to include in the reliability analysis. Must be scale variables.")
		}
	}

	Section
	{
		title: 	qsTr("Analysis")

		Group
		{
			title: qsTr("Scale Statistics")
			info:  qsTr("Reliability coefficients and descriptive statistics computed for the full scale.")

			CIField
			{
				name:			"ciLevel";
				label:			qsTr("Confidence interval");
				defaultValue:	95;
				info:			qsTr("Coverage of the confidence intervals for scale reliability statistics.")
			}
			

			CheckBox
			{
				id:			omega
				name:		"scaleOmega"
				label:		qsTr("McDonald's ω")
				checked:	true
				info:		qsTr("McDonald's omega for unidimensional data based on the single-factor model. The denominator uses model-implied total variance.")
			}

			CheckBox
			{
				id: 	alpha
				name: 	"scaleAlpha";
				label: 	qsTr("Cronbach's α");
				info:	qsTr("Cronbach's alpha. For binary items this equals KR-20.")
			}

			CheckBox
			{
				id: 	lambda2
				name: 	"scaleLambda2";
				label: 	qsTr("Guttman's λ2");
				info:	qsTr("Guttman's lambda 2, a lower bound for reliability.")
			}

			CheckBox
			{
				id: 	splithalf
				name: 	"scaleSplithalf";
				label: 	qsTr("Split-half coefficient");
				info:	qsTr("Splits items into two halves (odd/even by default). Unstandardized: Flanagan-Rulon coefficient; Standardized: Spearman-Brown coefficient.")
			}

			CheckBox { name: "averageInterItemCorrelation";	label: qsTr("Average interitem correlation"); info: qsTr("Mean of all pairwise Pearson correlations between items.")}

			RowLayout {
				CheckBox { name: "scaleMean";		label: qsTr("Mean");			id: mean		}
				CheckBox { name: "scaleVar";		label: qsTr("Variance"); 	id: variance}
				CheckBox { name: "scaleSd";			label: qsTr("SD");				id: sd			}

			}
			RadioButtonGroup
			{
				indent:		true
				enabled:	mean.checked || sd.checked || variance.checked
				name:		"meanSdScoresMethod"
				info:		qsTr("Whether the mean, variance, and SD are based on sum scores or mean scores across items.")

				RadioButton { value: "sumScores";	label: qsTr("of participants' sum scores"); checked: true}
				RadioButton { value: "meanScores";	label: qsTr("of participants' mean scores")}
			}

		}

		Group
		{
			title: qsTr("Individual Item Statistics")
			info:  qsTr("Reliability coefficients and descriptive statistics per item, including what-if statistics when that item is removed.")

			CIField
			{
				name: 			"itemCiLevel";
				label: 			qsTr("Confidence interval");
				defaultValue: 	95
				info:			qsTr("Coverage of the confidence intervals for item-level statistics.")
			}

			CheckBox
			{
				name: 		"itemDeletedOmega";
				label: 		qsTr("McDonald's ω  (if item dropped)");
				enabled: 	omega.checked
				info:		qsTr("Omega of the remaining items when this item is removed from the scale.")
			}

			CheckBox
			{
				name: 		"itemDeletedAlpha";
				label: 		qsTr("Cronbach's α (if item dropped)");
				enabled: 	alpha.checked
				info:		qsTr("Alpha of the remaining items when this item is removed from the scale.")
			}

			CheckBox
			{
				name: 		"itemDeletedLambda2";
				label: 		qsTr("Guttman's λ2 (if item dropped)");
				enabled: 	lambda2.checked
				info:		qsTr("Lambda 2 of the remaining items when this item is removed from the scale.")
			}

			CheckBox
			{
				name: 		"itemDeletedSplithalf";
				label: 		qsTr("Split-half coefficient (if item dropped)");
				enabled: 	splithalf.checked
				info:		qsTr("Split-half coefficient of the remaining items when this item is removed from the scale.")
			}
			
			CheckBox { name: "itemRestCorrelation";	label: qsTr("Item-rest correlation"); info: qsTr("Correlation of each item with the sum of the remaining items.")}
			RowLayout {
				CheckBox { name: "itemMean";			label: qsTr("Mean")							}
				CheckBox { name: "itemVar";			label: qsTr("Variance")			}
				CheckBox { name: "itemSd";			label: qsTr("SD")			}
			}

		}
	}

	Section
	{
		title: qsTr("Reverse-Scaled Items")

		VariablesForm
		{
			height: 150

			AvailableVariablesList 	{ name: "normalScaledItems"; 	title: qsTr("Normal-Scaled Items"); source: "variables" }
			AssignedVariablesList 	{ 
				name: "reverseScaledItems"; 	
				title: qsTr("Reverse-Scaled Items")
				 allowedColumns: ["scale"]
			}
		}
	}

	Section
	{
		title: qsTr("Advanced Options")
		Group
		{
			Layout.rowSpan: 2
			title: qsTr("Confidence intervals")
			info:  qsTr("Control the method used to construct confidence intervals for reliability coefficients, variance, and SD.")
			RadioButtonGroup
			{
				title: qsTr("Reliability coefficients")
				name:		"intervalMethod"
				id: intervalMethod
				info:	qsTr("Analytic intervals use normal-theory standard errors (van der Ark, 2024). Bootstrapped intervals use percentile resampling.")
				RadioButton
				{
					value: 		"analytic"
					label: 		qsTr("Analytic interval")
					checked: 	true
				}
				RadioButton
				{
					value: 	"bootstrapped"
					label: 	qsTr("Bootstrapped interval")

					IntegerField
					{
						name: 			"bootstrapSamples"
						label: 			qsTr("No. of bootstrap samples")
						defaultValue: 	1000
						fieldWidth: 	50
						min: 			100
						max: 			1e7
						info:			qsTr("Number of bootstrap replications. Higher values yield more stable interval estimates.")
					}

					RadioButtonGroup
					{
						title:		""
						name:		"bootstrapType"
						info:		qsTr("Non-parametric bootstrap resamples the data; parametric bootstrap samples from a multivariate normal with the estimated parameters.")

						RadioButton {value: "nonParametric"; label: qsTr("Non-parametric bootstrap"); checked: true}
						RadioButton {value: "parametric"; label: qsTr("Parametric bootstrap")}
					}
				}
			}
			RadioButtonGroup
			{
				title : qsTr("Variance and SD")
				name:		"intervalMethodVar"
				info:		qsTr("Chi-square-based intervals assume normality; non-parametric intervals use a two-step bootstrap procedure (van der Ark, 2024).")
				RadioButton
				{
					value: 		"chisq"
					label: 		qsTr("Chisq-based")
					checked: 	true
				}
				RadioButton
				{
					value: 	"twostep"
					label: 	qsTr("Non-parametric")

				}
			}
		}

		RadioButtonGroup
		{
			title: 		qsTr("McDonald's ω Estimation")
			name: 		"omegaEstimationMethod"
			enabled: 	omega.checked
			id: omegaEst
			info:		qsTr("CFA fits the single-factor model via confirmatory factor analysis; PFA uses principal factor analysis.")

			RadioButton
			{
				value: 		"cfa"
				label: 		qsTr("CFA")
				checked: 	true

				CheckBox
				{
					name: 		"omegaFitMeasures"
					label: 		qsTr("Single factor model fit")
					info:		qsTr("Chi-square test, RMSEA, SRMR and other fit indices for the single-factor model underlying omega.")
				}
			}
			
			RadioButton
			{
				value: "pfa"
				label: qsTr("PFA")
			}

			CheckBox
			{
				name:		"standardizedLoadings"
				label:		qsTr("Standardized factor loadings");
				info:		qsTr("Table of standardized loadings from the single-factor model.")
			}
		}

		Group
		{
			title: 		qsTr("Coefficients")
			info:		qsTr("Unstandardized coefficients use the raw covariance matrix; standardized use the correlation matrix.")

			RadioButtonGroup
 			{
 				title: ""
 				name: "coefficientType"

 				RadioButton{ value: "unstandardized"; label: qsTr("Unstandardized"); checked: true }
 				RadioButton{ value: "standardized";	label: qsTr("Standardized")}
			}
		}

		RadioButtonGroup
		{
				title: 	qsTr("Missing Values")
				name: 	"naAction"
				info:	qsTr("Pairwise uses all available observations per covariance pair; listwise deletes any row with a missing value.")

				RadioButton { value: "pairwise"; label: qsTr("Pairwise"); checked: true}
				RadioButton { value: "listwise"; label: qsTr("Delete listwise")}
		}

		Group
		{
			title: qsTr("Repeatability")
			info:  qsTr("Set a random seed to reproduce the same bootstrap results across runs.")

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
					fieldWidth: 	100
					min: 			1
					max: 			1e9
				}
			}
		}

		Group
		{
			title: qsTr("Samples")

			CheckBox
			{
				name:				"samplesSavingDisabled"
				label:				qsTr("Do not save samples")
				checked:			false
				info:				qsTr("When checked, bootstrap samples are not stored in the output file. This reduces file size but may slow down re-running the analysis, because samples are precomputed and cached for speed.")
			}
		}

		IntegerField
		{
			name: 					"hiddenScaleThreshold"
			defaultValue: 	preferencesModel.thresholdScale
			visible: 				false
		}
	}
}
