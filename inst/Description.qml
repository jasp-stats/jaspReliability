import QtQuick 		2.12
import JASP.Module 	1.0

Description
{
	name		: "Reliability"
	title		: qsTr("Reliability")
	description	: qsTr("This module offers the standard Frequentist and Bayesian Reliability analyses.")
	icon		: "reliability_icon_classic.svg"
	version		: "0.13.1"
	author		: "Julius M. Pfadt, Don van den Bergh & Eric-Jan Wagenmakers"
	maintainer	: "Julius M. Pfadt <julius.pfadt@gmail.com>"
	website		: "jasp-stats.org"
	license		: "GPL (>= 2)"
	
	Package	{ name: "Rdpack" 						}
	Package	{ name: "Bayesrel" 						}
	Package	{ name: "Rcsdp";	version: "0.1.57"; 	}

	GroupTitle
	{
		title:	qsTr("Classical")
		icon:	"reliability_icon_classic.svg"
	}

	Analysis
	{
		title	: qsTr("Single-Test Reliability Analysis")
		qml     : "ReliabilityFrequentist.qml"
		func	: "reliabilityFrequentist"
	}

	Separator {}

	GroupTitle
	{
		title:	qsTr("Bayesian")
		icon:	"reliability_icon_bayesian.svg"
	}

	Analysis
	{
		menu : qsTr("Single-Test Reliability Analysis")
		title: qsTr("Bayesian Single-Test Reliability Analysis")
		qml  : "ReliabilityBayesian.qml"
		func : "reliabilityBayesian"
	}
}
