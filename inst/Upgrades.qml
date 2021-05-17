import QtQuick 		2.12
import JASP.Module 	1.0

Upgrades
{
	Upgrade
	{
			functionName: 		"reliabilityFrequentist"
			newFunctionName:	"reliabilityUniDimFrequentist"
			fromVersion:		"0.14.3"
			toVersion:			"0.15"

			ChangeRename
			{
				from:	"mcDonaldScale"
				to:		"omegaScale"
			}
			ChangeRename
			{
				from:	"guttman2Scale"
				to:		"lambda2Scale"
			}
			ChangeRename
			{
				from:	"guttman6Scale"
				to:		"lambda6Scale"
			}
			ChangeRename
			{
				from:	"mcDonaldItem"
				to:		"omegaItem"
			}
			ChangeRename
			{
				from:	"guttman2Item"
				to:		"lambda2Item"
			}
			ChangeRename
			{
				from:	"guttman6Item"
				to:		"lambda6Item"
			}
			ChangeJS
			{
				name:		"bootType"
				jsFunction:	function(options) { return options["bootType"] === "bootNonpara" ? "nonParametric" : "parametric"; }
			}
			ChangeRename
			{
				from:	"omegaEst"
				to:		"omegaMethod"
			}
			ChangeRename
			{
				from:	"seedValue"
				to:		"seed"
			}
		}

	Upgrade
		{
			functionName: 		"reliabilityBayesian"
			newFunctionName:	"reliabilityUniDimBayesian"
			fromVersion:		"0.14.3"
			toVersion:			"0.15"

			ChangeRename
			{
				from:	"mcDonaldScale"
				to:		"omegaScale"
			}
			ChangeRename
			{
				from:	"guttman2Scale"
				to:		"lambda2Scale"
			}
			ChangeRename
			{
				from:	"guttman6Scale"
				to:		"lambda6Scale"
			}
			ChangeRename
			{
				from:	"mcDonaldItem"
				to:		"omegaItem"
			}
			ChangeRename
			{
				from:	"guttman2Item"
				to:		"lambda2Item"
			}
			ChangeRename
			{
				from:	"guttman6Item"
				to:		"lambda6Item"
			}
			ChangeRename
			{
				from:	"seedValue"
				to:		"seed"
			}


		}
}
