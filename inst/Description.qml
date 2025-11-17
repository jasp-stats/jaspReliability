import QtQuick
import JASP.Module

Description
{
	title		: qsTr("Reliability")
	description	: qsTr("Quantify the reliability of test scores")
	icon		: "reliability_icon_classic.svg"
	hasWrappers:  false
	
	GroupTitle
	{
		title:	qsTr("Classical")
		icon:	"reliability_icon_classic.svg"
	}

	Analysis
	{
		title:	qsTr("Reliability")
		qml: 	"UnidimensionalReliabilityFrequentist.qml"
		func: 	"unidimensionalReliabilityFrequentist"
	}

	Analysis
	{
		title:	qsTr("Intraclass Correlation")
		qml: 	"IntraclassCorrelation.qml"
		func: 	"intraclassCorrelation"
	}

	Analysis
	{
		title:	qsTr("Rater Agreement")
		qml: 	"RaterAgreement.qml"
		func: 	"raterAgreement"
	}
	Analysis
	{
		title:	qsTr("Bland-Altman Plots")
		qml: 	"BlandAltman.qml"
		func: 	"blandAltman"
	}

	Analysis
	{
		title:	qsTr("Standard Error of Measurement")
		qml: 	"StandardErrorOfMeasurement.qml"
		func: 	"standardErrorOfMeasurement"
	}

	Separator {}

	GroupTitle
	{
		title:	qsTr("Bayesian")
		icon:	"reliability_icon_bayesian.svg"
	}
	Analysis
	{
		menu: 	qsTr("Reliability")
		title: 	qsTr("Bayesian Reliability")
		qml: 	"UnidimensionalReliabilityBayesian.qml"
		func: 	"unidimensionalReliabilityBayesian"
	}

}
