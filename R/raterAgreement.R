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
  jaspTable$info <- gettext("Empty placeholder table shown while no agreement coefficient is selected; check a coefficient to obtain results.")
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

# union of the levels of all columns. Order matters: the codes derived from it act as
# ordinal positions. Numeric-looking labels are ordered by value; otherwise the sets are
# merged preserving each column's declared level order (raters may carry subsets of the
# levels, e.g. a rater who never used a middle category).
.raterAgreementUnionLevels <- function(dataset) {
  levelSets <- lapply(dataset, function(x) {
    if (is.factor(x)) levels(x) else unique(as.character(x[!is.na(x)]))
  })
  union <- unique(unlist(levelSets))

  if (!anyNA(suppressWarnings(as.numeric(union))))
    return(union[order(as.numeric(union))])

  merged <- levelSets[[which.max(lengths(levelSets))]]
  for (levelSet in levelSets) {
    pos <- 0L
    for (level in levelSet) {
      i <- match(level, merged)
      if (is.na(i)) {
        merged <- append(merged, level, after = pos)
        pos    <- pos + 1L
      } else {
        pos <- i
      }
    }
  }
  return(merged)
}

.raterAgreementIsDiscrete <- function(dataset) {
  vapply(dataset, function(x) is.factor(x) || is.character(x), logical(1L))
}

.raterAgreementAddCiColumns <- function(jaspTable, options) {
  ciPercent <- format(100 * options[["ciLevel"]], digits = 3, drop0trailing = TRUE)
  jaspTable$addColumnInfo(name = "SE",  title = gettext("SE"),    type = "number")
  jaspTable$addColumnInfo(name = "CIL", title = gettext("Lower"), type = "number", overtitle = gettextf("%s%% CI", ciPercent))
  jaspTable$addColumnInfo(name = "CIU", title = gettext("Upper"), type = "number", overtitle = gettextf("%s%% CI", ciPercent))
}

# every column mapped onto integer codes over the shared union of levels, so identical
# labels get identical codes across raters and codes respect the declared level order
.raterAgreementUnionCodes <- function(dataset, allLevels = .raterAgreementUnionLevels(dataset)) {
  mat <- vapply(dataset, function(x) as.numeric(factor(as.character(x), levels = allLevels)),
                numeric(nrow(dataset)))
  return(matrix(mat, nrow = nrow(dataset), dimnames = list(NULL, colnames(dataset))))
}

.raterAgreementHandleData <- function(dataset, options) {

  if (options[["dataStructure"]] == "ratersInColumns") {
    return(dataset)
  }

  # raters in rows: transpose so that rows are subjects and columns are raters.
  # t() would strip factor levels, so encode discrete columns onto the union of their
  # levels FIRST and rebuild factors (union levels, declared order) after transposing.
  isDiscrete <- .raterAgreementIsDiscrete(dataset)
  if (any(isDiscrete) && !all(isDiscrete))
    .quitAnalysis(gettext("With raters in rows, all subject/item columns must have the same measurement type. The data mix categorical and continuous columns; change the column types so they match."))

  if (any(isDiscrete)) {
    allLevels <- .raterAgreementUnionLevels(dataset)
    tCodes    <- t(.raterAgreementUnionCodes(dataset, allLevels))
    dataset   <- as.data.frame(lapply(seq_len(ncol(tCodes)),
                                      function(j) factor(allLevels[tCodes[, j]], levels = allLevels)))
  } else {
    dataset <- as.data.frame(t(as.matrix(dataset)))
  }
  colnames(dataset) <- paste0("Rater", seq_len(ncol(dataset)))

  return(dataset)
}

# complete-case numeric ratings for Kendall's W. Kendall uses ranks within each rater,
# so discrete data become union-level codes -- as.matrix() would turn factors into
# character labels that rank alphabetically instead of by declared order.
.kendallWRatings <- function(dataset) {
  if (any(.raterAgreementIsDiscrete(dataset))) {
    mat <- .raterAgreementUnionCodes(dataset)
  } else {
    mat <- as.matrix(dataset)
  }
  return(mat[stats::complete.cases(mat), , drop = FALSE])
}

# Krippendorff's alpha compares ratings ACROSS raters, so values must stay aligned
# between columns. as.matrix() would convert (ordered) factors to their labels, which
# breaks the ordinal/interval/ratio metrics (labels order alphabetically or coerce to
# NA). Nominal/ordinal: map all columns onto the union of their levels (declared order),
# so identical labels get identical codes. Interval/ratio: the metric needs the actual
# numeric distances, so labels are parsed as numbers (unparseable entries become NA and
# are caught by .krippAlphaValidate).
.krippAlphaRatings <- function(dataset, method) {
  if (!any(.raterAgreementIsDiscrete(dataset)))
    return(as.matrix(dataset))
  if (method %in% c("interval", "ratio")) {
    cols <- lapply(dataset, function(x) suppressWarnings(as.numeric(as.character(x))))
    mat  <- do.call(cbind, cols)
    colnames(mat) <- colnames(dataset)
    return(mat)
  }
  return(.raterAgreementUnionCodes(dataset))
}

# returns NULL if alpha can be computed, otherwise a user-facing error message
.krippAlphaValidate <- function(ratings, dataset, method) {
  if (ncol(ratings) < 2)
    return(gettext("Krippendorff's alpha requires at least 2 raters/measurements."))
  if (method %in% c("interval", "ratio") && any(is.na(ratings) & !is.na(as.matrix(dataset))))
    return(gettext("Krippendorff's alpha with an interval or ratio method requires numeric ratings. Use the nominal or ordinal method, or recode the data."))
  if (method == "ratio" && any(ratings < 0, na.rm = TRUE))
    return(gettext("Krippendorff's alpha with the ratio method requires non-negative ratings."))
  vals <- ratings[!is.na(ratings)]
  if (any(is.infinite(vals)))
    return(gettext("Krippendorff's alpha cannot be computed: the data contain infinite values."))
  if (length(vals) == 0L || !any(rowSums(!is.na(ratings)) >= 2))
    return(gettext("Krippendorff's alpha cannot be computed: no subject/item was rated by at least 2 raters (no pairable observations)."))
  if (length(unique(vals)) < 2)
    return(gettext("Krippendorff's alpha is not estimable: the ratings do not vary."))
  return(NULL)
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
      "variables.types",
      "cohensKappa",
      "cohensKappaType",
      "ci",
      "ciLevel",
      "weightType",
      "dataStructure"
    )
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

    if (ncol(dataset) < 2) {
      jaspTable$setError(gettext("Cohen's kappa requires at least 2 raters/measurements."))
      return(jaspTable)
    }

    if (nrow(dataset) > 2) { # psych gives an error when there are not at least 3 subjects rated
      #calculate Cohen's Kappas
      possiblePairs <- combn(ncol(dataset), 2)
      nPairs <- ncol(possiblePairs)

      # psych::cohen.kappa derives category order alphabetically from labels, which breaks
      # the (weighted) distance between declared ordinal levels -- feed it the level codes
      # from the union of all columns' levels instead
      cohenData <- dataset
      if (any(.raterAgreementIsDiscrete(cohenData)))
        cohenData <- as.data.frame(.raterAgreementUnionCodes(cohenData))

      out_kappa <- try(psych::cohen.kappa(cohenData, alpha = 1 - options[["ciLevel"]],
                                          w.exp = ifelse(options[["weightType"]] == "quadratic", 2, 1)),
                       silent = TRUE)
      # if weightType = linear, the exponent should be 1

      if (jaspBase::isTryError(out_kappa)) {
        jaspTable$setError(gettext("Cohen's kappa could not be computed. Ensure that each pair of raters has jointly rated subjects/items."))
        return(jaspTable)
      }

      if (nPairs == 1) {
        allKappaData <- list(out_kappa)
        allPairStrings <- paste(colnames(dataset), collapse = " - ")
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

      # subjects contribute pairwise: count rows with at least 2 ratings
      nSubjects <- sum(rowSums(!is.na(dataset)) >= 2)

      footnote <- gettextf('%1$i subjects/items and %2$i raters/measurements.', nSubjects, ncol(dataset))
      if (anyNA(dataset)) {
        footnote <- gettextf('%1$s Based on pairwise complete cases.', footnote)
        # pairwise sample sizes can differ, so report n per pair
        pairN <- apply(possiblePairs, 2, function(pair) sum(stats::complete.cases(dataset[, pair])))
        jaspTable$addColumnInfo(name = "n", title = gettext("n"), type = "integer")
        tableData[["n"]] <- c(NA, pairN)
      }

      if (options[["ci"]]) {
        .raterAgreementAddCiColumns(jaspTable, options)
        tableData[["SE"]] <- c(NA, allSE)
        tableData[["CIL"]] <- c(NA, allLowerBounds)
        tableData[["CIU"]] <- c(NA, allUpperBounds)
        footnote <- paste(footnote, gettext('Confidence intervals are asymptotic.'))
      }


      #if weighted kappa option is on but data only has 2 levels
      if (weighted && length(unique(stats::na.omit(unlist(lapply(dataset, as.character))))) < 3)
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
  jaspTable$info <- gettext("Fleiss' kappa: generalization of Cohen's kappa for two or more raters assigning subjects to nominal categories.")
  jaspTable$addColumnInfo(name = "ratings", title = gettext("Ratings"), type = "string")
  jaspTable$addColumnInfo(name = "fKappa", title = gettext("Fleiss' kappa"), type = "number")
  jaspTable$position <- 2

  #dependencies
  jaspTable$dependOn(
    options = c(
      "variables",
      "variables.types",
      "fleissKappa",
      "ci",
      "ciLevel",
      "dataStructure"
    )
  )

  if (ready) {

    if (any(options[["variables.types"]] == "scale")) {
      jaspTable$setError(gettext("Fleiss' kappa requires nominal or ordinal variables. Remove scale variables or change their type."))
      return(jaspTable)
    }

    if (ncol(dataset) < 2) {
      jaspTable$setError(gettext(
        "Fleiss' kappa requires at least 2 raters/measurements."
      ))
      return(jaspTable)
    }

    # everything label-based from here on: as.matrix() on mixed factor/numeric columns
    # would apply format() and pad numbers (" 7" vs "7"), breaking label matching, so
    # convert columns to character explicitly (irr gets the same clean labels below)
    fleissData      <- as.data.frame(lapply(dataset, as.character))

    # validate the ANALYZED (listwise-complete) data, not the raw rows: irr deletes
    # incomplete rows before computing
    completeRatings <- as.matrix(stats::na.omit(fleissData)) # same complete cases as irr uses
    if (nrow(completeRatings) < 3) {
      jaspTable$setError(gettext(
        "Fleiss' kappa requires at least 3 complete subjects/items (rows without missing ratings). Check whether raters are in columns or rows."
      ))
      return(jaspTable)
    }

    present <- unique(as.character(as.vector(completeRatings)))

    # categories are the labels present in the analyzed data, in numeric order when all
    # labels are numbers and in declared level order otherwise (unlist(dataset) must NOT
    # be used here: it returns factor CODES when factor and numeric columns are mixed)
    if (!anyNA(suppressWarnings(as.numeric(present)))) {
      categories <- present[order(as.numeric(present))]
    } else {
      allLevels  <- .raterAgreementUnionLevels(dataset)
      categories <- allLevels[allLevels %in% present]
    }

    if (length(categories) < 2) {
      jaspTable$setError(gettext(
        "Fleiss' kappa is not estimable: all ratings fall into a single category."
      ))
      return(jaspTable)
    }

    #calculate Fleiss' Kappa
    allKappaData <- irr::kappam.fleiss(fleissData)
    overallKappa <- allKappaData$value
    alpha        <- 1 - options[["ciLevel"]]

    ns <- allKappaData$subjects
    nr <- allKappaData$raters

    ratings <- c("Overall", categories)

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
      SE <- c(overallSE, categorySE)
      overallCI <- overallKappa + c(-1, 1) * qnorm(1 - alpha / 2) * overallSE
      categoryCIL <- categoryKappas - qnorm(1 - alpha / 2) * categorySE
      categoryCIU <- categoryKappas + qnorm(1 - alpha / 2) * categorySE
      .raterAgreementAddCiColumns(jaspTable, options)
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
      "variables.types",
      "krippendorffsAlpha",
      "krippendorffsAlphaMethod",
      "ci",
      "ciLevel",
      "dataStructure",
      "bootstrapSamples",
      "setSeed",
      "seed"
    )
  )

  if (ready) {
    #calculate Krippendorff's alpha
    method  <- options[["krippendorffsAlphaMethod"]]
    ratings <- .krippAlphaRatings(dataset, method)

    validationError <- .krippAlphaValidate(ratings, dataset, method)
    if (!is.null(validationError)) {
      jaspTable$setError(validationError)
      return(jaspTable)
    }

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

      .raterAgreementAddCiColumns(jaspTable, options)
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

  method  <- options[["krippendorffsAlphaMethod"]]
  ratings <- .krippAlphaRatings(dataset, method)

  if (!is.null(.krippAlphaValidate(ratings, dataset, method)))
    return() # the table shows the validation error

  bootstrapSamples <- createJaspState()
  bootstrapSamples$dependOn(options = c(
    "variables",
    "variables.types",
    "krippendorffsAlpha",
    "krippendorffsAlphaMethod",
    "ci",
    "bootstrapSamples",
    "dataStructure",
    "setSeed", "seed"))
  jaspResults[["bootstrapSamples"]] <- bootstrapSamples

  alphas   <- rep(NA_real_, options[["bootstrapSamples"]])
  n        <- nrow(ratings)
  tRatings <- t(ratings) # the irr-package expects raters to be in rows; transpose once, resample columns

  jaspBase::.setSeedJASP(options)

  for (i in seq_len(options[["bootstrapSamples"]])) {
    bootData <- tRatings[, sample.int(n, size = n, replace = TRUE), drop = FALSE]
    alpha    <- try(irr::kripp.alpha(bootData, method = method)$value, silent = TRUE)
    if (!jaspBase::isTryError(alpha))
      alphas[i] <- alpha
  }
  bootstrapSamples$object <- alphas
  return()
}

.kendallWBootRA <- function(jaspResults, dataset, options, ready) {
  if (!ready || !is.null(jaspResults[["kendallWBootstrapSamples"]]))
    return()

  # bootstrap CIs are only offered for the tie-corrected coefficient: resampling with
  # replacement introduces ties, which distort the uncorrected W and would put the
  # uncorrected point estimate on a different scale than its own CI
  if (!options[["correctForTies"]])
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
    "variables", "variables.types", "kendallW", "ci", "bootstrapSamples",
    "dataStructure", "setSeed", "seed"
  ))
  jaspResults[["kendallWBootstrapSamples"]] <- bootstrapSamples

  n  <- nrow(ratings)
  ws <- rep(NA_real_, options[["bootstrapSamples"]])

  jaspBase::.setSeedJASP(options)

  for (i in seq_len(options[["bootstrapSamples"]])) {
    bootData <- ratings[sample.int(n, size = n, replace = TRUE), , drop = FALSE]
    w        <- try(irr::kendall(bootData, correct = TRUE)$value, silent = TRUE)
    if (!jaspBase::isTryError(w))
      ws[i] <- w
  }

  bootstrapSamples$object <- ws
  return()
}

.computeKendallWTable <- function(jaspResults, dataset, options, ready) {
  jaspTable <- createJaspTable(title = gettext("Kendall's W"))
  jaspTable$info <- gettext("Kendall's coefficient of concordance W: measures agreement of rankings across multiple raters. Ranges from 0 (no agreement) to 1 (perfect concordance). Significance is assessed with the large-sample chi-square test and the F test, which performs better for small samples.")
  jaspTable$addColumnInfo(name = "W",     title = "W",        type = "number")
  jaspTable$addColumnInfo(name = "chisq", title = "\u03C7\u00B2", type = "number", overtitle = gettext("Chi-square test"))
  jaspTable$addColumnInfo(name = "df",    title = "df",       type = "integer", overtitle = gettext("Chi-square test"))
  jaspTable$addColumnInfo(name = "p",     title = "p",        type = "pvalue",  overtitle = gettext("Chi-square test"))
  jaspTable$addColumnInfo(name = "F",     title = "F",        type = "number",  overtitle = gettext("F test"))
  jaspTable$addColumnInfo(name = "df1",   title = "df1",      type = "number",  overtitle = gettext("F test"))
  jaspTable$addColumnInfo(name = "df2",   title = "df2",      type = "number",  overtitle = gettext("F test"))
  jaspTable$addColumnInfo(name = "pF",    title = "p",        type = "pvalue",  overtitle = gettext("F test"))
  jaspTable$position <- 3
  jaspTable$dependOn(options = c(
    "variables", "variables.types", "kendallW", "correctForTies", "ci", "ciLevel",
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

  if (ncol(dataset) < 2) {
    jaspTable$setError(gettext("Kendall's W requires at least 2 raters/measurements."))
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

  # constant ratings: the uncorrected W is a meaningless 0 and the corrected W is NaN
  if (all(apply(ratings, 2, stats::var) == 0)) {
    jaspTable$setError(gettext("Kendall's W is not estimable: the rankings do not vary."))
    return(jaspTable)
  }

  result <- irr::kendall(ratings, correct = options[["correctForTies"]])

  if (!is.finite(result$value)) {
    jaspTable$setError(gettext("Kendall's W is not estimable: the rankings do not vary."))
    return(jaspTable)
  }

  # F approximation (Kendall & Babington Smith, 1939; Legendre, 2005); the chi-square
  # test is unreliable especially for small numbers of raters (Marozzi, 2014 -- note that
  # Legendre's simulations did not examine this F test)
  W     <- result$value
  m     <- result$raters
  n     <- result$subjects
  Fstat <- (m - 1) * W / (1 - W)
  df1   <- n - 1 - 2 / m
  df2   <- df1 * (m - 1)

  fTestValid <- df1 > 0 && W < 1 # F is undefined at perfect concordance

  tableData <- list(
    W     = W,
    chisq = result$statistic,
    df    = result$subjects - 1L,
    p     = result$p.value,
    F     = if (fTestValid) Fstat else NA,
    df1   = if (fTestValid) df1   else NA,
    df2   = if (fTestValid) df2   else NA,
    pF    = if (fTestValid) stats::pf(Fstat, df1, df2, lower.tail = FALSE) else NA
  )

  footnote <- gettextf("%1$i subjects/items and %2$i raters/measurements.", result$subjects, result$raters)
  if (anyNA(dataset))
    footnote <- gettextf("%1$s Based on listwise complete cases.", footnote)

  if (options[["ci"]] && !options[["correctForTies"]]) {
    footnote <- paste(footnote, gettext("Bootstrap CIs are only available with the tie correction enabled: resampling with replacement introduces ties, which distort the uncorrected coefficient."))
  } else if (options[["ci"]] && !is.null(jaspResults[["kendallWBootstrapSamples"]])) {
    ws    <- jaspResults[["kendallWBootstrapSamples"]]$object
    conf  <- options[["ciLevel"]]
    probs <- (1 + c(-conf, conf)) / 2
    CIs   <- quantile(ws, probs = probs, na.rm = TRUE)

    .raterAgreementAddCiColumns(jaspTable, options)
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
  if (!fTestValid && W == 1)
    jaspTable$addFootnote(gettext("The F test is undefined for perfect concordance (W = 1)."))
  else if (!fTestValid)
    jaspTable$addFootnote(gettext("The F test is not available for this design: it requires n - 1 - 2/m > 0 (n = subjects/items, m = raters)."))
  if (!options[["correctForTies"]] && !is.null(result$error))
    jaspTable$addFootnote(gettext("Ties are present in the ratings, so the uncorrected coefficient may be inaccurate. Consider enabling the tie correction."), symbol = gettext("Warning:"))
  jaspTable$setData(tableData)
  return(jaspTable)
}
