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

// ------------------ stuff for R syntax ----------------------
	Upgrade
	{
		functionName: 		"raterAgreement"
		fromVersion:		"0.16.4"
		toVersion:			"0.17.0"

		ChangeRename
		{
			from:	"cohensWeightedOrNot"
			to:		"cohensKappaType"
		}
		ChangeJS
		{
			name:		"cohensKappaType"
			jsFunction:	function(options)
			{
				switch(options["cohensKappaType"])
				{
					case "cohensUnweighted":	return "unweighted";
					case "cohensWeighted":		return "weighted";
				}
			}
		}
		ChangeRename
		{
			from:	"alphaMethod"
			to:		"krippendorffsAlphaMethod"
		}
		ChangeRename
		{
			from:	"kappaIntervalOn"
			to:		"ci"
		}
		ChangeRename
		{
			from:	"kappaConfidenceIntervalValue"
			to:		"ciLevel"
		}
	}

	Upgrade
	{
		functionName: 		"intraclassCorrelation"
		fromVersion:		"0.16.4"
		toVersion:			"0.17.0"

		ChangeRename
		{
			from:	"iccRatingAverage"
			to:		"averagedRating"
		}
		ChangeRename
		{
			from:	"intervalOn"
			to:		"ci"
		}
		ChangeRename
		{
			from:	"confidenceIntervalValue"
			to:		"ciLevel"
		}
	}

	Upgrade
	{
		functionName: 		"reliabilityUniDimFrequentist"
		newFunctionName:	"unidimensionalReliabilityFrequentist"
		fromVersion:		"0.16.4"
		toVersion:			"0.17.0"

		ChangeRename
		{
			from:	"intervalOn"
			to:		"ci"
		}
		ChangeRename
		{
			from:	"confidenceIntervalValue"
			to:		"ciLevel"
		}
		ChangeRename
		{
			from:	"omegaScale"
			to:		"scaleOmega"
		}
		ChangeRename
		{
			from:	"alphaScale"
			to:		"scaleAlpha"
		}
		ChangeRename
		{
			from:	"lambda2Scale"
			to:		"scaleLambda2"
		}
		ChangeRename
		{
			from:	"lambda6Scale"
			to:		"scaleLambda6"
		}
		ChangeRename
		{
			from:	"glbScale"
			to:		"scaleGreatestLowerBound"
		}
		ChangeRename
		{
			from:	"averageInterItemCor"
			to:		"averageInterItemCorrelation"
		}
		ChangeRename
		{
			from:	"meanScale"
			to:		"scaleMean"
		}
		ChangeRename
		{
			from:	"sdScale"
			to:		"scaleSd"
		}
		ChangeRename
		{
			from:	"scoresMethod"
			to:		"meanSdScoresMethod"
		}
		ChangeRename
		{
			from:	"omegaItem"
			to:		"itemDeletedOmega"
		}
		ChangeRename
		{
			from:	"alphaItem"
			to:		"itemDeletedAlpha"
		}
		ChangeRename
		{
			from:	"lambda2Item"
			to:		"itemDeletedLambda2"
		}
		ChangeRename
		{
			from:	"lambda6Item"
			to:		"itemDeletedLambda6"
		}
		ChangeRename
		{
			from:	"glbItem"
			to:		"itemDeletedGreatestLowerBound"
		}
		ChangeRename
		{
			from:	"itemRestCor"
			to:		"itemRestCorrelation"
		}
		ChangeRename
		{
			from:	"meanItem"
			to:		"itemMean"
		}
		ChangeRename
		{
			from:	"sdItem"
			to:		"itemSd"
		}
		ChangeRename
		{
			from:	"missingValues"
			to:		"naAction"
		}
		ChangeJS
		{
			name:		"naAction"
			jsFunction:	function(options)
			{
				switch(options["naAction"])
				{
					case "excludeCasesPairwise":	return "pairwise";
					case "excludeCasesListwise":	return "listwise";
				}
			}
		}
		ChangeRename
		{
			from:	"noSamples"
			to:		"bootstrapSamples"
		}
		ChangeRename
		{
			from:	"bootType"
			to:		"bootstrapType"
		}
		ChangeRename
		{
			from:	"omegaMethod"
			to:		"omegaEstimationMethod"
		}
		ChangeRename
		{
			from:	"fitMeasures"
			to:		"omegaFitMeasures"
		}
		ChangeRename
		{
			from:	"omegaInterval"
			to:		"omegaIntervalMethod"
		}
		ChangeJS
		{
			name:		"omegaIntervalMethod"
			jsFunction:	function(options)
			{
				switch(options["omegaIntervalMethod"])
				{
					case "omegaAnalytic":	return "analytic";
					case "omegaBoot":		return "bootstrapped";
				}
			}
		}
		ChangeRename
		{
			from:	"alphaMethod"
			to:		"alphaType"
		}
		ChangeJS
		{
			name:		"alphaType"
			jsFunction:	function(options)
			{
				switch(options["alphaType"])
				{
					case "alphaUnstand":	return "unstandardized";
					case "alphaStand":		return "standardized";
				}
			}
		}
		ChangeRename
		{
			from:	"alphaInterval"
			to:		"alphaIntervalMethod"
		}
		ChangeJS
		{
			name:		"alphaIntervalMethod"
			jsFunction:	function(options)
			{
				switch(options["alphaIntervalMethod"])
				{
					case "alphaAnalytic":	return "analytic";
					case "alphaBoot":		return "bootstrapped";
				}
			}
		}
		ChangeRename
		{
			from:	"disableSampleSave"
			to:		"samplesSavingDisabled"
		}

	}

	Upgrade
	{
		functionName: 		"reliabilityUniDimBayesian"
		newFunctionName:	"unidimensionalReliabilityBayesian"
		fromVersion:		"0.16.4"
		toVersion:			"0.17.0"

		ChangeRename
		{
			from:	"credibleIntervalValueScale"
			to:		"scaleCiLevel"
		}
		ChangeRename
		{
			from:	"omegaScale"
			to:		"scaleOmega"
		}
		ChangeRename
		{
			from:	"alphaScale"
			to:		"scaleAlpha"
		}
		ChangeRename
		{
			from:	"lambda2Scale"
			to:		"scaleLambda2"
		}
		ChangeRename
		{
			from:	"lambda6Scale"
			to:		"scaleLambda6"
		}
		ChangeRename
		{
			from:	"glbScale"
			to:		"scaleGreatestLowerBound"
		}
		ChangeRename
		{
			from:	"averageInterItemCor"
			to:		"averageInterItemCorrelation"
		}
		ChangeRename
		{
			from:	"meanScale"
			to:		"scaleMean"
		}
		ChangeRename
		{
			from:	"sdScale"
			to:		"scaleSd"
		}
		ChangeRename
		{
			from:	"scoresMethod"
			to:		"meanSdScoresMethod"
		}
		ChangeRename
		{
			from:	"credibleIntervalValueItem"
			to:		"itemCiLevel"
		}
		ChangeRename
		{
			from:	"omegaItem"
			to:		"itemDeletedOmega"
		}
		ChangeRename
		{
			from:	"alphaItem"
			to:		"itemDeletedAlpha"
		}
		ChangeRename
		{
			from:	"lambda2Item"
			to:		"itemDeletedLambda2"
		}
		ChangeRename
		{
			from:	"lambda6Item"
			to:		"itemDeletedLambda6"
		}
		ChangeRename
		{
			from:	"glbItem"
			to:		"itemDeletedGreatestLowerBound"
		}
		ChangeRename
		{
			from:	"plotItem"
			to:		"itemDeletedPlot"
		}
		ChangeRename
		{
			from:	"orderItem"
			to:		"itemDeletedPlotOrdered"
		}
		ChangeRename
		{
			from:	"orderType"
			to:		"itemDeletedPlotOrderedType"
		}
		ChangeJS
		{
			name:		"itemDeletedPlotOrderedType"
			jsFunction:	function(options)
			{
				switch(options["itemDeletedPlotOrderedType"])
				{
					case "orderItemMean":	return "mean"
					case "orderItemKL":		return "kullbackLeibler"
					case "orderItemKS":		return "kolmogorovSmirnov"
				}
			}
		}
		ChangeRename
		{
			from:	"itemRestCor"
			to:		"itemRestCorrelation"
		}
		ChangeRename
		{
			from:	"meanItem"
			to:		"itemMean"
		}
		ChangeRename
		{
			from:	"sdItem"
			to:		"itemSd"
		}
		ChangeRename
		{
			from:	"plotPosterior"
			to:		"posteriorPlot"
		}
		ChangeRename
		{
			from:	"fixXRange"
			to:		"posteriorPlotFixedRange"
		}
		ChangeRename
		{
			from:	"dispPrior"
			to:		"posteriorPlotPriorDisplayed"
		}
		ChangeRename
		{
			from:	"probTable"
			to:		"probabilityTable"
		}
		ChangeRename
		{
			from:	"probTableValueLow"
			to:		"probabilityTableLowerBound"
		}
		ChangeRename
		{
			from:	"probTableValueHigh"
			to:		"probabilityTableUpperBound"
		}
		ChangeRename
		{
			from:	"shadePlots"
			to:		"posteriorPlotShaded"
		}
		ChangeRename
		{
			from:	"noSamples"
			to:		"samples"
		}
		ChangeRename
		{
			from:	"noBurnin"
			to:		"burnin"
		}
		ChangeRename
		{
			from:	"noThin"
			to:		"thinning"
		}
		ChangeRename
		{
			from:	"noChains"
			to:		"chains"
		}
		ChangeRename
		{
			from:	"disableSampleSave"
			to:		"samplesSavingDisabled"
		}
		ChangeRename
		{
			from:	"iwScale"
			to:		"inverseWishartPriorScale"
		}
		ChangeRename
		{
			from:	"iwDf"
			to:		"inverseWishartPriorDf"
		}
		ChangeRename
		{
			from:	"igShape"
			to:		"inverseGammaPriorShape"
		}
		ChangeRename
		{
			from:	"igScale"
			to:		"inverseGammaPriorScale"
		}
		ChangeRename
		{
			from:	"loadMean"
			to:		"normalPriorMean"
		}
		ChangeRename
		{
			from:	"missingValues"
			to:		"naAction"
		}
		ChangeJS
		{
			name:		"naAction"
			jsFunction:	function(options)
			{
				switch(options["naAction"])
				{
					case "excludeCasesPairwise":	return "imputation";
					case "excludeCasesListwise":	return "listwise";
				}
			}
		}
		ChangeRename
		{
			from:	"dispPPC"
			to:		"omegaPosteriorPredictiveCheck"
		}
		ChangeRename
		{
			from:	"fitMeasures"
			to:		"omegaFitMeasures"
		}
		ChangeRename
		{
			from:	"credibleIntervalValueFitMeasures"
			to:		"omegaFitMeasuresCiLevel"
		}
		ChangeRename
		{
			from:	"fitCutoffSat"
			to:		"omegaFitMeasuresCutoffRmsea"
		}
		ChangeRename
		{
			from:	"fitCutoffNull"
			to:		"omegaFitMeasuresCutoffCfiTli"
		}
		ChangeRename
		{
			from:	"dispLoadings"
			to:		"standardizedLoadings"
		}
		ChangeRename
		{
			from:	"stdCoeffs"
			to:		"coefficientType"
		}
		ChangeJS
		{
			name:		"coefficientType"
			jsFunction:	function(options)
			{
				switch(options["coefficientType"])
				{
					case "unstand":	return "undstandardized";
					case "stand":	return "standardized";
				}
			}
		}
		ChangeRename
		{
			from:	"pointEst"
			to:		"pointEstimate"
		}

	}
}
