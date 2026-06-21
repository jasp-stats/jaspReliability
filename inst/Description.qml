import QtQuick
import JASP.Module

Description
{
	title		: qsTr("Reliability")
	description	: qsTr("Quantify the reliability of test scores")
	icon		: "reliability_icon_classic.svg"
	hasWrappers:  false
	preloadData: true
	
	GroupTitle
	{
		title:	qsTr("Classical")
		icon:	"reliability_icon_classic.svg"
	}

	Analysis
	{
		title:	qsTr("Reliability")
		qml: 	"ReliabilityUnidimensionalFrequentist.qml"
		func: 	"reliabilityUnidimensionalFrequentist"
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
		title: 	qsTr("Bayesian Reliability")
		menu: 	qsTr("Reliability")
		qml: 	"ReliabilityUnidimensionalBayesian.qml"
		func: 	"reliabilityUnidimensionalBayesian"
	}
	Analysis
	{
		title: 	qsTr("Bayesian Multidimensional Reliability")
		menu: 	qsTr("Multidimensional Reliability")
		qml: 	"ReliabilityMultidimensionalBayesian.qml"
		func: 	"reliabilityMultidimensionalBayesian"
	}

}
