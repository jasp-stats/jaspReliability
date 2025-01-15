import QtQuick
import JASP.Module

Description
{
	name		: "jaspReliability"
	title		: qsTr("Reliability")
	description	: qsTr("Quantify the reliability of test scores")
	version			: "0.19.2"
	author		: "JASP Team"
	maintainer	: "Julius M. Pfadt <julius.pfadt@gmail.com>"
	website		: "https://github.com/jasp-stats/jaspReliability"
	license		: "GPL (>= 2)"
	icon		: "reliability_icon_classic.svg"
	preloadData: true


	GroupTitle
	{
		title:	qsTr("Classical")
		icon:	"reliability_icon_classic.svg"
	}

	Analysis
	{
		title:	qsTr("Unidimensional Reliability")
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
		preloadData: false
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
		menu: 	qsTr("Unidimensional Reliability")
		title: 	qsTr("Bayesian Unidimensional Reliability")
		qml: 	"UnidimensionalReliabilityBayesian.qml"
		func: 	"unidimensionalReliabilityBayesian"
	}

}
