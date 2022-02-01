import QtQuick 		2.12
import JASP.Module 	1.0

Upgrades
{
	Upgrade
	{
		functionName: 		"IntraclassCorrelation"
		newFunctionName:	"intraclassCorrelation"
		fromVersion:		"0.16"
		toVersion:			"0.16.1"
	}

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
			condition:	function(options) { return options["guttman2Scale"] !== undefined }
			from:		"guttman2Scale"
			to:			"lambda2Scale"
		} // for 0.15 the name of the coefficient changed from guttman2 to lambda2
		ChangeSetValue
		{
			condition:	function(options) { return options["guttman2Scale"] === undefined }
			name:		"lambda2Scale"
			jsonValue:	false
		} // in jasp versions before 0.13 the coefficient name did not exist, hence this!

		ChangeRename
		{
			condition:	function(options) { return options["guttmanScale"] !== undefined }
			from:	"guttmanScale"
			to:		"lambda6Scale"
		} // in jasp versions before 0.13 the name of the coefficient was guttmanScale
		ChangeRename
		{
			condition:	function(options) { return options["guttman6Scale"] !== undefined }
			from:	"guttman6Scale"
			to:		"lambda6Scale"
		} // for 0.15 the coefficient name changed again.

		ChangeSetValue
		{
			condition:	function(options) { return options["scoresMethod"] === undefined }
			name:		"scoresMethod"
			jsonValue:	"meanScores"
		} // in older jasp versions there were no mean and sd options (rowSums and rowMeans),
		// just the colMeans. This way, the refreshed result is closer to the old versions

		ChangeRename
		{
			from:	"mcDonaldItem"
			to:		"omegaItem"
		}

		ChangeRename
		{
			condition:	function(options) { return options["guttman2Item"] !== undefined }
			from:		"guttman2Item"
			to:			"lambda2Item"
		} // for 0.15 the name of the coefficient changed from guttman2 to lambda2
		ChangeSetValue
		{
			condition:	function(options) { return options["guttman2Item"] === undefined }
			name:		"lambda2Item"
			jsonValue:	false
		} // in jasp versions before 0.13 the coefficient name did not exist, hence this!

		ChangeRename
		{
			condition:	function(options) { return options["guttmanItem"] !== undefined }
			from:	"guttmanItem"
			to:		"lambda6Item"
		} // in jasp versions before 0.13 the name of the coefficient was guttmanItem
		ChangeRename
		{
			condition:	function(options) { return options["guttman6Item"] !== undefined }
			from:	"guttman6Item"
			to:		"lambda6Item"
		} // for 0.15 the coefficient name changed again.

		ChangeJS
		{
			condition:	function(options) { return options["bootType"] !== undefined }
			name:		"bootType"
			jsFunction:	function(options) { return options["bootType"] === "bootNonpara" ? "nonParametric" : "parametric"; }
		}
		ChangeSetValue
		{
			condition:	function(options) { return options["bootType"] === undefined }
			name:		"bootType"
			jsonValue:	"nonParametric"
		}

		// before 0.13 omegaEst and seedValue were not defined
		ChangeRename
		{
			condition:	function(options) { return options["omegaEst"] !== undefined }
			from:	"omegaEst"
			to:		"omegaMethod"
		}
		ChangeRename
		{
			condition:	function(options) { return options["seedValue"] !== undefined }
			from:	"seedValue"
			to:		"seed"
		}

		// before 0.13 the radiobuttongroup was named alphaScaleStandardized
		ChangeRename
		{
			condition:	function(options) { return options["alphaScaleStandardized"] !== undefined }
			from:		"alphaScaleStandardized"
			to:			"alphaMethod"
		}
		ChangeJS // this also changes the values (names) of the radiobuttons
		{
			condition:	function(options) { return options["alphaScaleStandardized"] !== undefined }
			name:		"alphaMethod"
			jsFunction:	function(options)
			{
				switch(options["alphaMethod"])
				{
					case "_1unstandardized":	return "alphaUnstand";
					case "_2standardized":	return "alphaStand";
				}
			}
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

		ChangeSetValue
		{
			condition:	function(options) { return options["scoresMethod"] === undefined }
			name:		"scoresMethod"
			jsonValue:	"meanScores"
		} // in older jasp versions there were no mean and sd options (rowSums and rowMeans),
		// just the colMeans. This way, the refreshed result is closer to the old versions

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
