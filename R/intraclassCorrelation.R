IntraclassCorrelation <- function(jaspResults, dataset, options) {

  dataset <- .reliabilityReadData(dataset, options)

  if (ncol(dataset) > 2) {
    jaspResults[["table"]] <- handleIntraclassCorrelation(dataset, options)
  } else {
    # TODO: message to select data?
  }

  return()
}

handleIntraclassCorrelation <- function (dataset, options) {

  # Get the ICC type e.g. "ICC1k"
  type <- toupper(options["iccType"])
  if (options[["iccRatingAverage"]]) {
    type <- paste0(type, "k")
  }

  # Compute the ICC
  full_results <- psych::ICC(
    dataset,
    alpha = 1 - options[["confidenceIntervalValue"]]
  )
  icc_results <- full_results$results
  
  # Select correct ICC
  icc <- icc_results[icc_results$type == type, ]
  rownames(icc) <- NULL
  
  # Round all numeric columns
  numeric_columns <- unlist(lapply(icc, is.numeric))
  icc[ , numeric_columns] <- round(icc[ , numeric_columns], digits = 3)

  # Only report relevant columns
  cols <- c("type", "ICC")
  # Report CI
  if (options[["intervalOn"]]) {
    cols <- c(cols, "lower bound", "upper bound")
  }
  icc <- icc[, cols]

  jaspTable <- createJaspTable(title = gettext("Intraclass Correlation"))
  jaspTable$dependOn(options = c("variables", "intervalOn", "confidenceIntervalValue", "iccType", "iccRatingAverage"))
  # jaspTable$addColumnInfo(name = "estimate", title = gettext("Estimate"), type = "string")
  jaspTable$setData(icc)
  jaspTable$addFootnote(gettextf("%s subjects and %s judges. ICC type as referenced by Shrout & Fleiss (1979).", full_results$n.obs, full_results$n.judge))

  return(jaspTable)
}
