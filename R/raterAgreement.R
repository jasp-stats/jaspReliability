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

  dataset <- .raterAgreementHandleData(dataset, options)

  anyCoefficient <- options[["cohensKappa"]] || options[["fleissKappa"]] ||
                    options[["krippendorffsAlpha"]] || options[["kendallW"]]

  if (!anyCoefficient)
    .raterAgreementPlaceholderTable(jaspResults, options, ready)

  if (options[["cohensKappa"]])
    jaspResults[["cohensKappa"]] <- .computeCohensKappaTable(dataset, options, ready)
  if (options[["fleissKappa"]])
    jaspResults[["fleissKappa"]] <- .computeFleissKappaTable(dataset, options, ready)
  if (options[["krippendorffsAlpha"]]) {
    if (options[["ci"]])
      .kripAlphaBoot(jaspResults, dataset, options, ready)
    jaspResults[["krippendorffsAlpha"]] <- .computeKrippendorffsAlphaTable(jaspResults, dataset, options, ready)
  }
  if (options[["kendallW"]]) {
    if (options[["ci"]])
      .kendallWBootRA(jaspResults, dataset, options, ready)
    jaspResults[["kendallW"]] <- .computeKendallWTable(jaspResults, dataset, options, ready)
  }

  return()
}

.raterAgreementPlaceholderTable <- function(jaspResults, options, ready) {
  if (!is.null(jaspResults[["placeholder"]]))
    return()

  jaspTable <- createJaspTable(title = gettext("Agreement Coefficient"))
  jaspTable$info <- gettext("Overview of all selected agreement coefficients with standard errors and confidence intervals.")
  jaspTable$addColumnInfo(name = "coefficient", title = gettext("Coefficient"), type = "string")
  jaspTable$addColumnInfo(name = "estimate",    title = gettext("Estimate"),    type = "number")
  jaspTable$addColumnInfo(name = "SE",          title = gettext("SE"),          type = "number")
  jaspTable$addColumnInfo(name = "CIL",         title = gettext("Lower"),       type = "number")
  jaspTable$addColumnInfo(name = "CIU",         title = gettext("Upper"),       type = "number")
  if (ready)
    jaspTable$addFootnote(gettext("Check one of the coefficients to start the analysis."))
  jaspTable$dependOn(options = c("cohensKappa", "fleissKappa", "krippendorffsAlpha", "kendallW", "variables"))
  jaspResults[["placeholder"]] <- jaspTable
}

.raterAgreementHandleData <- function(dataset, options) {

  if (options[["dataStructure"]] == "ratersInColumns") {
    return(dataset)
  }

  # raters in rows: transpose so that rows are subjects and columns are raters
  dataset <- as.data.frame(t(dataset))

  return(dataset)
}

# complete-case numeric ratings for Kendall's W. Kendall only uses ranks WITHIN each
# rater, so (ordered) factors can safely become their level codes -- as.matrix() would
# turn them into character labels that rank alphabetically instead of by declared order.
# Character data (e.g. after the raters-in-rows transpose) is parsed as numbers,
# unparseable entries become NA.
.kendallWRatings <- function(dataset) {
  cols <- lapply(dataset, function(x) if (is.factor(x)) as.numeric(x) else x)
  mat  <- do.call(cbind, cols)
  colnames(mat) <- colnames(dataset)
  if (is.character(mat))
    mat <- matrix(suppressWarnings(as.numeric(mat)), nrow = nrow(mat), dimnames = dimnames(mat))
  return(mat[stats::complete.cases(mat), , drop = FALSE])
}

# Krippendorff's alpha compares ratings ACROSS raters, so values must stay aligned
# between columns. as.matrix() would convert (ordered) factors to their labels, which
# breaks the ordinal/interval/ratio metrics (labels order alphabetically or coerce to
# NA). Instead map all columns onto the union of their levels, preserving each column's
# declared level order, so identical labels get identical codes across raters.
.krippAlphaRatings <- function(dataset) {
  isDiscrete <- vapply(dataset, function(x) is.factor(x) || is.character(x), logical(1L))
  if (!any(isDiscrete))
    return(as.matrix(dataset))
  allLevels <- unique(unlist(lapply(dataset, function(x) {
    if (is.factor(x)) levels(x) else unique(as.character(x[!is.na(x)]))
  })))
  cols <- lapply(dataset, function(x) as.numeric(factor(as.character(x), levels = allLevels)))
  mat  <- do.call(cbind, cols)
  colnames(mat) <- colnames(dataset)
  return(mat)
}

.computeCohensKappaTable <- function(dataset, options, ready) {

  weighted <- options[["cohensKappaType"]] == "weighted"

  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettext("Cohen's kappa"))
  jaspTable$info <- gettext("Cohen's kappa: chance-corrected agreement between exactly two raters. Ranges from -1 (worse than chance) to 1 (perfect agreement).")
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "cKappa", title = gettextf("kappa"), type = "number")
  jaspTable$position <- 1

  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "cohensKappa",
      "cohensKappaType",
      "ci",
      "ciLevel",
      "weightType",
      "dataStructure"
    )
  )


  formattedCIPercent <- format(
    100 * options[["ciLevel"]],
    digits = 3,
    drop0trailing = TRUE
  )


  if (ready) {

    if (any(options[["variables.types"]] == "scale")) {
      jaspTable$setError(gettext("Cohen's kappa requires nominal or ordinal variables. Remove scale variables or change their type."))
      return(jaspTable)
    }

    if (weighted && any(options[["variables.types"]] == "nominal")) {
      jaspTable$setError(gettext("Weighted Cohen's kappa requires ordinal variables. Remove nominal variables or change their type."))
      return(jaspTable)
    }

    if (nrow(dataset) > 2) { # psych gives an error when there are not at least 3 subjects rated
      #calculate Cohen's Kappas
      possiblePairs <- combn(ncol(dataset), 2)
      nPairs <- ncol(possiblePairs)

      out_kappa <- psych::cohen.kappa(dataset, alpha = 1 - options[["ciLevel"]],
                                      w.exp = ifelse(options[["weightType"]] == "quadratic", 2, 1))
      # if weightType = linear, the exponent should be 1

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

      # because cohens kappa uses pairwise rater agreements the number of subjects overall is not listwise deleted
      if (options[["dataStructure"]] == "ratersInColumns") {
        dataCount <- dataset[rowSums(is.na(dataset)) < ncol(dataset), ]
      } else {
        dataCount <- dataset[, colSums(is.na(dataset)) < nrow(dataset)]
      }

      footnote <- gettextf('%1$i subjects/items and %2$i raters/measurements.', nrow(dataCount), ncol(dataCount))
      if (anyNA(dataset))
        footnote <- gettextf('%1$s Based on pairwise complete cases.', footnote)

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
    } else {
      jaspTable$setError(gettext("Cohen's kappa requires at least 3 subjects/items."))
    }


  }

  return(jaspTable)
}

.computeFleissKappaTable <- function(dataset, options, ready) {

  # Create the JASP Table
  jaspTable <- createJaspTable(title = gettextf("Fleiss' kappa"))
  jaspTable$info <- gettext("Fleiss' kappa: generalization of Cohen's kappa for three or more raters assigning subjects to nominal categories.")
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "fKappa", title = gettext("Fleiss' kappa"), type = "number")
  jaspTable$position <- 2

  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "fleissKappa",
      "ci",
      "ciLevel",
      "dataStructure"
    )
  )

  formattedCIPercent <- format(
    100 * options[["ciLevel"]],
    digits = 3,
    drop0trailing = TRUE
  )

  if (ready) {

    if (any(options[["variables.types"]] == "scale")) {
      jaspTable$setError(gettext("Fleiss' kappa requires nominal or ordinal variables. Remove scale variables or change their type."))
      return(jaspTable)
    }

    if (nrow(dataset) < 3) {
      jaspTable$setError(gettext(
        "Fleiss' kappa requires at least 3 subjects/items. Check whether raters are in columns or rows."
      ))
      return(jaspTable)
    }

    if (ncol(dataset) < 2) {
      jaspTable$setError(gettext(
        "Fleiss' kappa requires at least 2 raters/measurements."
      ))
      return(jaspTable)
    }

    #calculate Fleiss' Kappa
    allKappaData <- irr::kappam.fleiss(dataset)
    overallKappa <- allKappaData$value
    alpha        <- 1 - options[["ciLevel"]]

    ns <- allKappaData$subjects
    nr <- allKappaData$raters

    # categories must come from the analyzed (listwise-complete) data: a category that only
    # occurs in an incomplete row is not part of the results
    completeRatings <- as.matrix(stats::na.omit(dataset)) # same complete cases as irr uses
    present         <- unique(as.character(as.vector(completeRatings)))

    # display order: numeric when all categories are numbers, otherwise factor-level/alphabetical
    rawCategories <- unique(unlist(dataset))
    rawCategories <- rawCategories[!is.na(rawCategories)]
    if (anyNA(suppressWarnings(as.numeric(as.character(rawCategories))))) {
      categories <- as.character(sort(rawCategories))
    } else {
      categories <- as.character(rawCategories)[order(as.numeric(as.character(rawCategories)))]
    }
    categories <- categories[categories %in% present]
    ratings    <- c("Overall", categories)

    # compute the per-category kappas and all SEs directly (Fleiss, Nee & Landis, 1979) --
    # irr's detail table rounds kappa and z to 3 decimals, so reconstructing the SE as
    # kappa/z is inaccurate and breaks down entirely (0/0) for zero kappa
    counts <- t(apply(completeRatings, 1, function(row) table(factor(as.character(row), levels = categories))))
    pj     <- colSums(counts) / (ns * nr)
    qj     <- 1 - pj
    pjk    <- (colSums(counts^2) - ns * nr * pj) / (ns * nr * (nr - 1) * pj)

    categoryKappas <- (pjk - pj) / (1 - pj)
    categorySE     <- rep(sqrt(2 / (ns * nr * (nr - 1))), length(categories))
    overallSE      <- sqrt((2 / (sum(pj * qj)^2 * (ns * nr * (nr - 1)))) *
                             (sum(pj * qj)^2 - sum(pj * qj * (qj - pj))))

    tableData <- list("ratings" = ratings,
                      "fKappa"  = c(overallKappa, categoryKappas))

    footnote <- gettextf('%1$i subjects/items and %2$i raters/measurements.', allKappaData$subjects, allKappaData$raters)
    if (anyNA(dataset))
      footnote <- gettextf('%1$s Based on listwise complete cases.', footnote)


    if (options[["ci"]]) {
      nCategories <- length(categories)
      SE <- c(overallSE, categorySE)
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
  jaspTable <- createJaspTable(title = gettext("Krippendorff's alpha"))
  jaspTable$info <- gettext("Krippendorff's alpha: reliability coefficient applicable to any number of raters, any scale level (nominal/ordinal/interval/ratio), and incomplete data.")
  jaspTable$addColumnInfo(name = "method", title = gettext("Method"), type = "string")
  jaspTable$addColumnInfo(name = "kAlpha", title = gettext("Krippendorff's alpha"), type = "number")
  jaspTable$position <- 2

  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "krippendorffsAlpha",
      "krippendorffsAlphaMethod",
      "ci",
      "ciLevel",
      "dataStructure",
      "bootstrapSamples"
    )
  )

  formattedCIPercent <- format(
    100 * options[["ciLevel"]],
    digits = 3,
    drop0trailing = TRUE
  )

  if (ready) {
    #calculate Krippendorff's alpha
    method  <- options[["krippendorffsAlphaMethod"]]
    ratings <- .krippAlphaRatings(dataset)
    kAlpha  <- irr::kripp.alpha(t(ratings), method) # the irr-package expects raters to be in rows.

    tableData <- list("method" = paste0(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method))),
                      "kAlpha" = kAlpha$value)

    footnote <- gettextf('%1$i subjects/items and %2$i raters/measurements.', kAlpha$subjects, kAlpha$raters)
    if (anyNA(dataset))
      footnote <- gettextf('%1$s Based on pairwise complete cases.', footnote)

    if (options[["ci"]] && !is.null(jaspResults[["bootstrapSamples"]])) {
      alphas <- jaspResults[["bootstrapSamples"]]$object
      conf <- options[["ciLevel"]]
      confs <- (1 + c(-conf, conf)) / 2
      CIs <- quantile(alphas, probs = confs, na.rm = TRUE)

      jaspTable$addColumnInfo(name = "SE", title = gettext("SE"), type = "number")
      jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI", formattedCIPercent))
      tableData[["SE"]] <- sd(alphas, na.rm = TRUE)
      tableData[["CIL"]] <- CIs[1]
      tableData[["CIU"]] <- CIs[2]

      nFailed <- sum(is.na(alphas))
      if (nFailed > 0)
        footnote <- paste(footnote, gettextf("%1$i of %2$i bootstrap samples could not be computed and were excluded from the CI.", nFailed, length(alphas)))
    }
    jaspTable$setData(tableData)
    jaspTable$addFootnote(footnote)
  }

  return(jaspTable)
}

.kripAlphaBoot <- function(jaspResults, dataset, options, ready) {
   if (!ready || !is.null(jaspResults[["bootstrapSamples"]]))
    return()

  bootstrapSamples <- createJaspState()
  bootstrapSamples$dependOn(options = c(
    "variables",
    "krippendorffsAlpha",
    "krippendorffsAlphaMethod",
    "ci",
    "bootstrapSamples",
    "dataStructure",
    "setSeed", "seed"))
  jaspResults[["bootstrapSamples"]] <- bootstrapSamples

  method  <- options[["krippendorffsAlphaMethod"]]
  ratings <- .krippAlphaRatings(dataset)
  alphas  <- rep(NA_real_, options[["bootstrapSamples"]])
  n       <- nrow(ratings)

  jaspBase::.setSeedJASP(options)

  for (i in seq_len(options[["bootstrapSamples"]])) {
    bootData <- ratings[sample.int(n, size = n, replace = TRUE), , drop = FALSE]
    alpha    <- try(irr::kripp.alpha(t(bootData), method = method)$value, silent = TRUE)
    if (!jaspBase::isTryError(alpha))
      alphas[i] <- alpha
  }
  bootstrapSamples$object <- alphas
  return()
}

.kendallWBootRA <- function(jaspResults, dataset, options, ready) {
  if (!ready || !is.null(jaspResults[["kendallWBootstrapSamples"]]))
    return()

  if (any(options[["variables.types"]] == "nominal"))
    return()

  # validation and listwise deletion must precede the bootstrap: resampling raw rows can
  # produce replicates with too few complete cases, erroring deep inside irr::kendall()
  ratings <- .kendallWRatings(dataset)
  if (nrow(ratings) < 2 || any(is.infinite(ratings)))
    return() # the table shows the validation error

  bootstrapSamples <- createJaspState()
  bootstrapSamples$dependOn(options = c(
    "variables", "kendallW", "ci", "bootstrapSamples",
    "correctForTies", "dataStructure", "setSeed", "seed"
  ))
  jaspResults[["kendallWBootstrapSamples"]] <- bootstrapSamples

  n       <- nrow(ratings)
  correct <- options[["correctForTies"]]
  ws      <- rep(NA_real_, options[["bootstrapSamples"]])

  jaspBase::.setSeedJASP(options)

  for (i in seq_len(options[["bootstrapSamples"]])) {
    bootData <- ratings[sample.int(n, size = n, replace = TRUE), , drop = FALSE]
    w        <- try(irr::kendall(bootData, correct = correct)$value, silent = TRUE)
    if (!jaspBase::isTryError(w))
      ws[i] <- w
  }

  bootstrapSamples$object <- ws
  return()
}

.computeKendallWTable <- function(jaspResults, dataset, options, ready) {
  formattedCIPercent <- format(
    100 * options[["ciLevel"]],
    digits        = 3,
    drop0trailing = TRUE
  )

  jaspTable <- createJaspTable(title = gettext("Kendall's W"))
  jaspTable$info <- gettext("Kendall's coefficient of concordance W: measures agreement of rankings across multiple raters. Ranges from 0 (no agreement) to 1 (perfect concordance).")
  jaspTable$addColumnInfo(name = "W",     title = "W",                   type = "number")
  jaspTable$addColumnInfo(name = "chisq", title = gettext("Chi-square"), type = "number")
  jaspTable$addColumnInfo(name = "df",    title = "df",                  type = "integer")
  jaspTable$addColumnInfo(name = "p",     title = "p",                   type = "pvalue")
  jaspTable$position <- 3
  jaspTable$dependOn(options = c(
    "variables", "kendallW", "correctForTies", "ci", "ciLevel",
    "bootstrapSamples", "dataStructure", "setSeed", "seed"
  ))

  if (!ready)
    return(jaspTable)

  if (any(options[["variables.types"]] == "nominal")) {
    jaspTable$setError(gettext(
      "Kendall's W requires ordinal or scale variables. Remove nominal variables."
    ))
    return(jaspTable)
  }

  ratings <- .kendallWRatings(dataset)

  if (nrow(ratings) < 2) {
    jaspTable$setError(gettext(
      "Kendall's W requires at least 2 complete subjects/items (rows without missing ratings)."
    ))
    return(jaspTable)
  }

  if (any(is.infinite(ratings))) {
    jaspTable$setError(gettext("Kendall's W cannot be computed: the data contain infinite values."))
    return(jaspTable)
  }

  result <- irr::kendall(ratings, correct = options[["correctForTies"]])

  tableData <- list(
    W     = result$value,
    chisq = result$statistic,
    df    = result$subjects - 1L,
    p     = result$p.value
  )

  footnote <- gettextf("%1$i subjects/items and %2$i raters/measurements.", result$subjects, result$raters)
  if (anyNA(dataset))
    footnote <- gettextf("%1$s Based on listwise complete cases.", footnote)

  if (options[["ci"]] && !is.null(jaspResults[["kendallWBootstrapSamples"]])) {
    ws    <- jaspResults[["kendallWBootstrapSamples"]]$object
    conf  <- options[["ciLevel"]]
    probs <- (1 + c(-conf, conf)) / 2
    CIs   <- quantile(ws, probs = probs, na.rm = TRUE)

    jaspTable$addColumnInfo(name = "SE",  title = gettext("SE"),    type = "number")
    jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number",
                            overtitle = gettextf("%s%% CI", formattedCIPercent))
    jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number",
                            overtitle = gettextf("%s%% CI", formattedCIPercent))
    tableData[["SE"]]  <- stats::sd(ws, na.rm = TRUE)
    tableData[["CIL"]] <- CIs[[1L]]
    tableData[["CIU"]] <- CIs[[2L]]
    footnote <- paste(footnote, gettext("Confidence intervals are based on bootstrap."))

    nFailed <- sum(is.na(ws))
    if (nFailed > 0)
      footnote <- paste(footnote, gettextf("%1$i of %2$i bootstrap samples could not be computed and were excluded from the CI.", nFailed, length(ws)))
  }

  jaspTable$addFootnote(footnote)
  jaspTable$addFootnote(gettext("Chi-square test is valid for large samples only."))
  if (!options[["correctForTies"]] && !is.null(result$error))
    jaspTable$addFootnote(gettext("Ties are present in the ratings, so the uncorrected coefficient may be inaccurate. Consider enabling the tie correction."), symbol = gettext("Warning:"))
  jaspTable$setData(tableData)
  return(jaspTable)
}
