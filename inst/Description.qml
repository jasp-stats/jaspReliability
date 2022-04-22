import QtQuick 		2.12
import JASP.Module 	1.0

Description
{
	name		: "jaspReliability"
	title		: qsTr("Reliability")
	description	: qsTr("Quantify the reliability of test scores")
	version		: "0.15"
	author		: "Julius M. Pfadt, Don van den Bergh & Eric-Jan Wagenmakers"
	maintainer	: "Julius M. Pfadt <julius.pfadt@gmail.com>"
	website		: "https://github.com/jasp-stats/jaspReliability"
	license		: "GPL (>= 2)"
	icon		: "reliability_icon_classic.svg"


	GroupTitle
	{
		title:	qsTr("Classical")
		icon:	"reliability_icon_classic.svg"
	}

	Analysis
	{
		title:	qsTr("Unidimensional Reliability")
		qml: 	"uniDimFrequentist.qml"
		func: 	"reliabilityUniDimFrequentist"
	}

	Analysis
	{
		title:	qsTr("Intraclass Correlation")
		qml: 	"intraclassCorrelation.qml"
		func: 	"intraclassCorrelation"
	}
	
		Analysis
	{
		title:	qsTr("Rater Agreement")
		qml: 	"raterAgreement.qml"
		func: 	"raterAgreement"
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
		qml: 	"uniDimBayesian.qml"
		func: 	"reliabilityUniDimBayesian"
	}

}
