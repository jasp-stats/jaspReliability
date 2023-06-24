#
# Copyright (C) 2021 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' @export
raterAgreement <- function(jaspResults, dataset, options) {

  ready <- length(options[["variables"]]) > 1

  dataset <- .readDataCohensFleissKappa(dataset, options)

  if (options[["cohensKappa"]])
    jaspResults[["cohensKappa"]] <- .computeCohensKappaTable(dataset, options, ready)
  if (options[["fleissKappa"]])
    jaspResults[["fleissKappa"]] <- .computeFleissKappaTable(dataset, options, ready)
  if (options[["krippendorffsAlpha"]]) {
    if (options[["ci"]])
      .kripAlphaBoot(jaspResults, dataset, options, ready)
    jaspResults[["krippendorffsAlpha"]] <- .computeKrippendorffsAlphaTable(jaspResults, dataset, options, ready)
  }

  return()
}

# Read in the dataset (copied from .reliabilityReadData)
.readDataCohensFleissKappa <- function(dataset, options) {
  variables <- unlist(options[["variables"]])
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.factor = variables)
  }
  return(dataset)
}

.computeCohensKappaTable <- function(dataset, options, ready) {

  weighted <- options[["cohensKappaType"]] == "weighted"

  weightedString <- ifelse(weighted, "Weighted", "Unweighted")

  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Cohen's %s kappa", weightedString))
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "cKappa", title = gettextf("%s kappa", weightedString), type = "number")
  jaspTable$position <- 1

  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "cohensKappa",
      "cohensKappaType",
      "ci",
      "ciLevel"
    )
  )


  formattedCIPercent <- format(
    100 * options[["ciLevel"]],
    digits = 3,
    drop0trailing = TRUE
  )


  if (ready) {
    #calculate Cohen's Kappas
    nPairs <- ncol(dataset) * (ncol(dataset) - 1) / 2
    out_kappa <- psych::cohen.kappa(dataset, alpha = 1 - options[["ciLevel"]])
    if (nPairs == 1) {
      allKappaData <- list(out_kappa)
      allPairStrings <- paste(options[["variables"]], collapse = " - ")
    } else {
      allKappaData <- out_kappa[2:(nPairs + 1)]
      allPairStrings <- sub(" ", " - ", names(out_kappa[2:(nPairs + 1)]))
    }

    #Extract Kappas and CIs
    allKappas <- c()
    allSE <- c()
    allLowerBounds <- c()
    allUpperBounds <- c()

    for (i in allKappaData) {
      kappaData <- i$confid
      k <- ifelse(weighted, 2, 1)
      allKappas <- c(allKappas, kappaData[k, 2])
      allSE <- c(allSE, sqrt(ifelse(weighted, i$var.weighted, i$var.kappa)))
      allLowerBounds <- c(allLowerBounds, kappaData[k, 1])
      allUpperBounds <- c(allUpperBounds, kappaData[k, 3])
    }

    averageKappa <- mean(allKappas)

    tableData <- list("ratings" = c("Average kappa", allPairStrings),
                      "cKappa" = c(averageKappa, allKappas))
    footnote <- gettextf('%i subjects/items and %i raters/measurements.', nrow(dataset), ncol(dataset))

    if (options[["ci"]]) {
      jaspTable$addColumnInfo(name = "SE", title = gettext("SE"), type = "number")
      jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      tableData[["SE"]] <- c(NA, allSE)
      tableData[["CIL"]] <- c(NA, allLowerBounds)
      tableData[["CIU"]] <- c(NA, allUpperBounds)
      footnote <- paste(footnote, gettext('Confidence intervals are asymptotic.'))
    }


    #if weighted kappa option is on but data only has 2 levels
    if (weighted && length(levels(unlist(dataset))) < 3)
      footnote <- paste(footnote, gettext('If there are only 2 levels, weighted kappa is equal to unweighted kappa.'))

    jaspTable$setData(tableData)
    jaspTable$addFootnote(footnote)
  }
  return(jaspTable)
}

.computeFleissKappaTable <- function(dataset, options, ready) {

  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Fleiss' kappa"))
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "fKappa", title = gettext("Fleiss' kappa"), type = "number")
  jaspTable$position <- 2

  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "fleissKappa",
      "ci",
      "ciLevel"
    )
  )

  formattedCIPercent <- format(
    100 * options[["ciLevel"]],
    digits = 3,
    drop0trailing = TRUE
  )

  if (ready) {
    #calculate Fleiss' Kappa
    allKappaData <- .fleissKappaMod(dataset, detail = TRUE)
    overallKappa <- allKappaData$value
    categoryKappas <- allKappaData$detail[, 1]
    overallSE <- allKappaData$se
    categorySE <- allKappaData$se_cat
    alpha <- 1 - options[["ciLevel"]]

    # for nominal text data we want the rating text to be displayed:
    # if the transformation to numeric goes wrong:
    if (anyNA(as.numeric(as.character(unique(unlist(dataset)))))) {
      categories <- sort(unique(unlist(dataset)))
    } else { # if data is ordinal (can be transformed to numeric)
      categories <- sort(as.numeric(as.character(unique(unlist(dataset)))))
    }
    ratings <- c("Overall", as.character(categories))
    # when the data was read with columns as factors the kappa function
    # would not order the categories numerically, this is a fix
    categoryKappas <- categoryKappas[as.character(categories)]

    tableData <- list("ratings" = ratings,
                      "fKappa"  = c(overallKappa, categoryKappas))
    footnote <- gettextf('%i subjects/items and %i raters/measurements.', nrow(dataset), ncol(dataset))

    if (options[["ci"]]) {
      nCategories <- length(categories)
      SE <- c(overallSE, rep(categorySE, nCategories))
      overallCI <- overallKappa + c(-1, 1) * qnorm(1 - alpha / 2) * overallSE
      categoryCIL <- categoryKappas - qnorm(1 - alpha / 2) * categorySE
      categoryCIU <- categoryKappas + qnorm(1 - alpha / 2) * categorySE
      jaspTable$addColumnInfo(name = "SE", title = gettext("SE"), type = "number")
      jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      tableData[["SE"]] <- SE
      tableData[["CIL"]] <- c(overallCI[1], categoryCIL)
      tableData[["CIU"]] <- c(overallCI[2], categoryCIU)
      footnote <- paste(footnote, gettext('Confidence intervals are asymptotic.'))
    }

    jaspTable$setData(tableData)
    jaspTable$addFootnote(footnote)
  }
  return(jaspTable)
}

.computeKrippendorffsAlphaTable <- function(jaspResults, dataset, options, ready) {
  # Create the JASP Table
  jaspTable <- createJaspTable(title = "Krippendorff's alpha")
  jaspTable$addColumnInfo(name = "method", title = gettext("Method"), type = "string")
  jaspTable$addColumnInfo(name = "kAlpha", title = "Krippendorff's alpha", type = "number")
  jaspTable$position <- 2

  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "krippendorffsAlpha",
      "ci",
      "ciLevel",
      "krippenDorffsAlphaDataStructure"
    )
  )

  formattedCIPercent <- format(
    100 * options[["ciLevel"]],
    digits = 3,
    drop0trailing = TRUE
  )

  if (ready) {
    #calculate Krippendorff's alpha
    if (options[["krippendorffsAlphaDataStructure"]] == "ratersInColumns") {
      kAlphaData <- t(as.matrix(dataset))
    } else { # raters in rows
      kAlphaData <- as.matrix(dataset)
    }
    method <- options[["krippendorffsAlphaMethod"]]
    kAlpha <- irr::kripp.alpha(kAlphaData, method)

    tableData <- list("method" = paste0(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method))),
                      "kAlpha" = kAlpha$value)
    footnote <- gettextf('%i subjects/items and %i raters/measurements.', kAlpha$subjects, kAlpha$raters)

    if (options[["ci"]]) {
      alphas <- jaspResults[["bootstrapSamples"]]$object
      conf <- options[["ciLevel"]]
      confs <- (1 + c(-conf, conf)) / 2
      CIs <- quantile(alphas, probs = confs, na.rm = TRUE)

      jaspTable$addColumnInfo(name = "SE", title = gettext("SE"), type = "number")
      jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      tableData[["SE"]] <- sd(alphas)
      tableData[["CIL"]] <- CIs[1]
      tableData[["CIU"]] <- CIs[2]
    }
    jaspTable$setData(tableData)
    jaspTable$addFootnote(footnote)
  }

  return(jaspTable)
}

.kripAlphaBoot <- function(jaspResults, dataset, options, ready) {
   if (!ready || !is.null(jaspResults[["bootstrapSamples"]]$object))
    return()

  bootstrapSamples <- createJaspState()
  method <- options[["krippendorffsAlphaMethod"]]
  n <- nrow(dataset)
  alphas <- numeric(options[["krippendorffsAlphaBootstrapSamplesForCI"]])

  jaspBase::.setSeedJASP(options)

  for (i in seq_len(options[["krippendorffsAlphaBootstrapSamplesForCI"]])) {
    bootData <- as.matrix(dataset[sample.int(n, size = n, replace = TRUE), ])
    alphas[i] <- irr::kripp.alpha(t(bootData), method = method)$value
  }
  bootstrapSamples$object <- alphas
  jaspResults[["bootstrapSamples"]] <- bootstrapSamples
  jaspResults[["bootstrapSamples"]]$dependOn(options = c(
    "variables",
    "krippendorffsAlpha",
    "ci",
    "krippendorffsAlphaBootstrapSamplesForCI",
    "setSeed", "seed"))
  return()
}



#####################################################
#copied from IRR package, modified to output Std. Err.
######################################################

.fleissKappaMod <- function(ratings, exact = FALSE, detail = FALSE)
{
  ratings <- as.matrix(na.omit(ratings))
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  lev <- levels(as.factor(ratings))
  for (i in 1:ns) {
    frow <- factor(ratings[i, ], levels = lev)
    if (i == 1)
      ttab <- as.numeric(table(frow))
    else ttab <- rbind(ttab, as.numeric(table(frow)))
  }
  ttab <- matrix(ttab, nrow = ns)
  agreeP <- sum((apply(ttab^2, 1, sum) - nr)/(nr * (nr - 1))/ns)
  if (!exact) {
    method <- "Fleiss' Kappa for m Raters"
    chanceP <- sum(apply(ttab, 2, sum)^2)/(ns * nr)^2
  }
  else {
    method <- "Fleiss' Kappa for m Raters (exact value)"
    for (i in 1:nr) {
      rcol <- factor(ratings[, i], levels = lev)
      if (i == 1)
        rtab <- as.numeric(table(rcol))
      else rtab <- rbind(rtab, as.numeric(table(rcol)))
    }
    rtab <- rtab/ns
    chanceP <- sum(apply(ttab, 2, sum)^2)/(ns * nr)^2 -
      sum(apply(rtab, 2, var) * (nr - 1)/nr)/(nr - 1)
  }
  value <- (agreeP - chanceP)/(1 - chanceP)
  if (!exact) {
    pj <- apply(ttab, 2, sum)/(ns * nr)
    qj <- 1 - pj
    varkappa <- (2/(sum(pj * qj)^2 * (ns * nr * (nr - 1)))) *
      (sum(pj * qj)^2 - sum(pj * qj * (qj - pj)))
    SEkappa <- sqrt(varkappa)
    u <- value/SEkappa
    p.value <- 2 * (1 - pnorm(abs(u)))
    if (detail) {
      pj <- apply(ttab, 2, sum)/(ns * nr)
      pjk <- (apply(ttab^2, 2, sum) - ns * nr * pj)/(ns *
                                                       nr * (nr - 1) * pj)
      kappaK <- (pjk - pj)/(1 - pj)
      varkappaK <- 2/(ns * nr * (nr - 1))
      SEkappaK <- sqrt(varkappaK)
      uK <- kappaK/SEkappaK
      p.valueK <- 2 * (1 - pnorm(abs(uK)))
      tableK <- as.table(round(cbind(kappaK, uK, p.valueK),
                               digits = 3))
      rownames(tableK) <- lev
      colnames(tableK) <- c("Kappa", "z", "p.value")
    }
  }
  if (!exact) {
    if (!detail) {
      rval <- list(method = method, subjects = ns, raters = nr,
                   irr.name = "Kappa", value = value)
    }
    else {
      ###############
      #BEGIN CHANGES
      ###############
      rval <- list(method = method, subjects = ns, raters = nr,
                   irr.name = "Kappa", value = value, detail = tableK, se = SEkappa, se_cat = SEkappaK)
      ############
      #END CHANGES
      #############
    }
    rval <- c(rval, stat.name = "z", statistic = u, p.value = p.value)
  }
  else {
    rval <- list(method = method, subjects = ns, raters = nr,
                 irr.name = "Kappa", value = value)
  }
  class(rval) <- "irrlist"
  return(rval)
}

