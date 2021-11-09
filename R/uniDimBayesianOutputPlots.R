

.BayesianPosteriorPlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainer"]]$object))
    return()

  plotContainer <- createJaspContainer(gettext("Posterior Plots"))
  plotContainer$dependOn(options = c("plotPosterior", "shadePlots", "probTable", "probTableValueLow",
                                     "probTableValueHigh", "fixXRange", "dispPrior", "credibleIntervalValueScale",
                                     "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale"))

  derivedOptions <- model[["derivedOptions"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  idxSelected   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

  if (options[["shadePlots"]] && options[["probTable"]]) {
    if (options[["probTableValueLow"]] > options[["probTableValueHigh"]]) {
      low <- options[["probTableValueHigh"]]
      high <- options[["probTableValueLow"]]
    } else {
      low <- options[["probTableValueLow"]]
      high <- options[["probTableValueHigh"]]
    }
    shadePlots <- c(low, high)
  } else {
    shadePlots <- NULL
  }

  if (options[["plotPosterior"]] && is.null(model[["empty"]])) {
    n.item <- model[["k"]]

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (is.null(plotContainer[[nmsObjsNoGreek[i]]])) {
        if (options[["dispPrior"]]) {
          if (nm == "omegaScale") {
            startProgressbar(2e3)
          } else {
            startProgressbar(4e3)
          }
          prior <- .samplePrior(n.item, nm, progressbarTick, options[["iwScale"]], options[["iwDf"]],
                                options[["igShape"]], options[["igScale"]])
        } else {
          prior <- NULL
        }

        p <- .makeSinglePosteriorPlot(model[[nm]], model[["scaleResults"]][["cred"]][[nm]], nmsLabs[[i]],
                                      options[["fixXRange"]], shadePlots,
                                      options[["dispPrior"]], prior)
        plotObj <- createJaspPlot(plot = p, title = nmsObjs[i])
        plotObj$dependOn(options = names(idxSelected[i]))
        plotObj$position <- i
        plotContainer[[nmsObjsNoGreek[i]]] <- plotObj

      }
    }
    plotContainer$position <- 4
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["plotContainer"]] <- plotContainer
  }

  return()
}

.makeSinglePosteriorPlot <- function(coefList, coefResults, nms, fixXRange, shade = NULL, priorTrue, priorSample) {

  # TODO: consider precomputing all densities (maybe with kernsmooth?) and reducing memory that way

  samp_tmp <- as.vector(coefList[["samp"]])
  if (fixXRange) {
    d <- stats::density(samp_tmp, from = 0, to = 1, n = 2^10)
  } else {
    d <- stats::density(samp_tmp, n = 2^10)
  }
  datDens <- data.frame(x = d$x, y = d$y)
  datPrior <- data.frame(x = priorSample$x, y = priorSample$y)


  xBreaks <- jaspGraphs::getPrettyAxisBreaks(datDens$x)

  # max height posterior is at 90% of plot area; remainder is for credible interval
  ymax <- max(d$y) / .9
  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, ymax))
  ymax <- max(yBreaks)
  scaleCriRound <- round(coefResults, 3)
  datCri <- data.frame(xmin = scaleCriRound[1L], xmax = scaleCriRound[2L], y = .925 * ymax)
  height <- (ymax - .925 * ymax) / 2
  if (fixXRange) {
    if (datCri$xmin == datCri$xmax) { # if two zeros, the interval is merged together
      datTxt <- data.frame(x = c(datCri$xmin, datCri$xmax),
                           y = 0.985 * ymax,
                           label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
                           stringsAsFactors = FALSE)
    } else {
      datTxt <- data.frame(x = c(datCri$xmin - .08, datCri$xmax + .08),
                           y = 0.985 * ymax,
                           label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
                           stringsAsFactors = FALSE)
    }

  } else {
    datTxt <- data.frame(x = c(datCri$xmin, datCri$xmax),
                         y = 0.985 * ymax,
                         label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
                         stringsAsFactors = FALSE)
  }


  if (datCri$xmin[1L] < 0) {
    datCri$xmin[1L] <- 0
    datTxt$x[1L] <- 0
    datTxt$label[1L] <- "< 0"
  } else if (datCri$xmin[1L] == 0) {
    datCri$xmin[1L] <- 0
    datTxt$x[1L] <- 0
    datTxt$label[1L] <- "0"
  }

  # if the bounds are less than 0.05 away from 0 or 1, expand the axis by 0.1 so the credible interval text does not
  # get chopped off.
  xExpand <- .1 * ((c(0, 1) - datTxt$x) <= 0.05)
  if (fixXRange && max(datDens$y, na.rm = T) >= 1000) xExpand <- c(xExpand[1], .05)
  # with large numebrs on the y-axis, the x-axis labels to the right get cut off sometimes, when the range is fixed

  g <- ggplot2::ggplot(data = datDens, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(size = .85) +
    ggplot2::geom_errorbarh(data = datCri, mapping = ggplot2::aes(xmin = xmin, xmax = xmax, y = y),
                            height = height, inherit.aes = FALSE) +
    ggplot2::geom_text(data = datTxt, mapping = ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
                       size = 5) +
    ggplot2::scale_y_continuous(name = gettext("Density"), breaks = yBreaks, limits = range(yBreaks)) +
    ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, expand = xExpand, limits = range(xBreaks))

  if (!is.null(shade)) {
    datFilter <- datDens[datDens[["x"]] >= shade[1] & datDens[["x"]] <= shade[2], ]
    if (length(datFilter$x) == 0)
      datFilter <- data.frame(x = 0, y = 0)
    g <- g + ggplot2::geom_ribbon(data = datFilter, mapping = ggplot2::aes(ymin = 0, ymax = y),
                                  fill = "grey", alpha = 0.95) +
      ggplot2::geom_line(size = .85)
  }

  if (priorTrue) {
    g <- g + ggplot2::geom_line(data = datPrior, mapping = ggplot2::aes(x = x, y = y),
                                linetype = "dashed", size = .85) +
      ggplot2::scale_x_continuous(name = nms, breaks = xBreaks, limits = range(xBreaks),
                                  expand = xExpand)

  }

  return(jaspGraphs::themeJasp(g))

}



.BayesianIfItemPlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainerItem"]]$object))
    return()

  plotContainerItem <- createJaspContainer(gettext("If Item Dropped Posterior Plots"))
  plotContainerItem$dependOn(options = c("variables", "plotItem",
                                         "credibleIntervalValueItem", "orderType", "orderItem",
                                         "omegaItem", "alphaItem", "lambda2Item", "lambda6Item", "glbItem"))

  derivedOptions <- model[["derivedOptions"]]

  selected <- derivedOptions[["itemDroppedSelectedItem"]]
  idxSelected   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables_item"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

  if (options[["orderItem"]]) {
    ordering <- options[["orderType"]]
  } else {
    ordering <- NULL
  }

  if (is.null(model[["empty"]]) && options[["plotItem"]]) {
    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (is.null(plotContainerItem[[nmsObjsNoGreek[i]]])) {

        if (length(options[["variables"]]) < 3) {
          plotObjItem <- createJaspPlot(plot = NULL, title = nmsObjs[i], width = 400)
          plotObjItem$setError(gettext("Please enter at least 3 variables for the if item dropped plot"))
        } else {
          # use the coefficient item object in the model list, and the object directly before it,
          # which is always the corresponding scale object
          prevNumber <- which(names(model) == nm) - 1
          name <- unlist(strsplit(nm, "Item"))
          coefPos <- grep(name, names(model[["scaleResults"]][["est"]]))
          p <- .makeIfItemPlot(model[[nm]], model[[prevNumber]],
                               model[["itemResults"]][["est"]][[nm]],
                               model[["scaleResults"]][["est"]][[coefPos]],
                               nmsLabs[[i]],
                               options[["credibleIntervalValueItem"]],
                               ordering = ordering, model[["itemsDropped"]])
          plotObjItem <- createJaspPlot(plot = p, title = nmsObjs[i], width = 400)
          plotObjItem$dependOn(options = names(idxSelected[i]))
          plotObjItem$position <- i
          if (is.null(p)) {
            plotObjItem$setError(gettext("KLD ordering failed because two variables have perfect correlation"))
          }
        }
        plotContainerItem[[nmsObjsNoGreek[i]]] <- plotObjItem

      }
    }
    plotContainerItem$position <- 5
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["plotContainerItem"]] <- plotContainerItem
  }

  return()
}

.makeIfItemPlot <- function(coefItem, coefScale, coefItemEst, coefScaleEst, nms, int, ordering, variables) {
  n_row <- length(variables)
  lower <- (1 - int) / 2
  upper <- int + (1 - int) / 2

  samp_tmp <- as.vector(coefScale[["samp"]])
  dat <- data.frame(as.matrix(samp_tmp), row.names =  NULL)
  names(dat) <- "value"
  dat$colos <- "1"
  dat$var <- "original"

  dat_del <- t(as.matrix(as.data.frame(coefItem[["itemSamp"]])))
  names <- decodeColNames(variables)

  for (i in n_row:1) {
    tmp <- as.data.frame(dat_del[i, ])
    colnames(tmp) <- "value"
    tmp$var <- names[i]
    tmp$colos <- "2"
    dat <- rbind(dat, tmp)
  }
  dat$var <- factor(dat$var, levels = unique(dat$var))

  if (!is.null(ordering)) {
    est <- as.data.frame(coefItemEst)
    est[n_row + 1, ] <- 1
    colnames(est) <- "value"
    est$name <- c(names, "original")

    if (ordering == "orderItemMean") {
      dists <- abs(coefScaleEst - coefItemEst)
      dists[length(dists) + 1] <- 0
      est <- est[order(dists, decreasing = FALSE), ]
      dat$var <- factor(dat$var, levels = c(est$name))

    } else if (ordering == "orderItemKL") {
      samps <- coefItem[["itemSamp"]]
      og_samp <- samp_tmp

      dists <- try(apply(samps, 2, .KLD.statistic, y = og_samp)) # kl divergence
      ### when there are only three variables and two of them have almost perfect correlation, KLD fails
      if (any(round(cor(samps)[lower.tri(cor(samps))], 3) == 1) && inherits(dists, "try-error")) {
        return(NULL)
      } else {
        dists[length(dists) + 1] <- 0
        est <- est[order(dists), ]
        dat$var <- factor(dat$var, levels = c(est$name))
      }

    } else if (ordering == "orderItemKS") {
      samps <- coefItem[["itemSamp"]]
      og_samp <- samp_tmp
      dists <- apply(samps, 2, .ks.test.statistic, y = og_samp) # ks distance
      dists[length(dists) + 1] <- 0
      est <- est[order(dists), ]
      dat$var <- factor(dat$var, levels = c(est$name))
    }
  }

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = value, y = var, fill = colos)) +
    ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = c(lower, 0.5, upper),
                                  alpha = .85, show.legend = FALSE, scale = 1) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(colour = "black")) +
    ggplot2::ylab(gettext("Item Dropped")) +
    ggplot2::xlab(nms) +
    ggplot2::scale_fill_grey() +
    ggplot2::scale_y_discrete(expand = ggplot2::expand_scale(mult = c(0.1, 0.25)))

  return(jaspGraphs::themeJasp(g))

}



.omegaPosteriorPredictive <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["omegaPPC"]]$object))
    return()

  if (options[["dispPPC"]] && options[["omegaScale"]] && is.null(model[["empty"]])) {

    ll <- model[["omegaScale"]][["loadings"]]
    rr <- model[["omegaScale"]][["residuals"]]
    cobs <- model[["data_cov"]]

    k <- ncol(cobs)
    nsamp <- nrow(ll)
    ee_impl <- matrix(0, nsamp, k)
    for (i in seq_len(nsamp)) {
      ctmp <- ll[i, ] %*% t(ll[i, ]) + diag(rr[i, ])
      dtmp <- MASS::mvrnorm(model[["n"]], rep(0, k), ctmp)
      ee_impl[i, ] <- eigen(cov(dtmp), only.values = TRUE)$values
    }

    eframe <- data.frame(number = seq(1, k), eigen_value = eigen(cobs)$values)
    eframe$eigen_sim_low <- apply(ee_impl, 2, quantile, prob = .025)
    eframe$eigen_sim_up <- apply(ee_impl, 2, quantile, prob = .975)

    yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(eframe$eigen_sim_up)))

    g <- ggplot2::ggplot(eframe, mapping = ggplot2::aes(x = number, y = eigen_value)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = eigen_sim_low, ymax = eigen_sim_up), color = "grey55", width = 0.2,
                             size = 1) +
      ggplot2::geom_point(size = 2.25) +
      ggplot2::xlim(c(1, k)) +
      ggplot2::scale_y_continuous(name = gettext("Eigenvalue"), breaks = yBreaks, limits = range(yBreaks)) +
      ggplot2::scale_x_continuous(name = gettext("Eigenvalue No."),
                                  breaks = seq(1, k),
                                  expand = ggplot2::expand_scale(mult = c(.1, .1)))

    g <- jaspGraphs::themeJasp(g)

    plot <- createJaspPlot(plot = g, title = "Posterior Predictive Check Omega", width = 350)
    plot$dependOn(options = c("dispPPC", "omegaScale"))

    plot$position <- 6
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["omegaPPC"]] <- plot
  }

  return()
}



.BayesianTracePlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainerTP"]]$object))
    return()

  if (is.null(model[["empty"]]) && options[["tracePlot"]]) {

    plotContainerTP <- createJaspContainer(gettext("Convergence Traceplot"))
    plotContainerTP$dependOn(options = c("tracePlot", "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale",
                                         "glbScale"))

    derivedOptions <- model[["derivedOptions"]]
    selected <- derivedOptions[["selectedEstimatorsPlots"]]
    idxSelected   <- which(selected)
    nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
    nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]]
    nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]

    for (j in seq_along(idxSelected)) {
      i <- idxSelected[j]
      nm <- names(idxSelected[j])

      if (is.null(plotContainerTP[[nmsObjsNoGreek[i]]])) {

        p <- .makeTracePlot(model[[nm]], nmsLabs[[i]])
        plotObjTP <- createJaspPlot(plot = p, title = nmsObjs[i], width = 400)
        plotObjTP$dependOn(options = names(idxSelected[i]))
        plotObjTP$position <- i
        plotContainerTP[[nmsObjsNoGreek[i]]] <- plotObjTP

      }
    }

    plotContainerTP$position <- 7
    stateContainer <- .getStateContainerB(jaspResults)
    stateContainer[["plotContainerTP"]] <- plotContainerTP

  }

  return()
}


.makeTracePlot <- function(coefList, nms) {

  dd <- coefList[["samp"]]
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, length(dd[1, ])))

  dv <- cbind(dd[1, ], 1, seq(1, ncol(dd)))
  for (j in 2:nrow(dd)) {
    dv <- rbind(dv, cbind(dd[j, ], j, seq(1, ncol(dd))))
  }
  dat <- data.frame(dv)
  colnames(dat) <- c("Value", "chain", "Iterations")
  dat$chain <- as.factor(dat$chain)

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = Iterations, y = Value, colour = chain)) +
    ggplot2::geom_line(size = .3) +
    ggplot2::ylab(nms) +
    ggplot2::scale_x_continuous(name = gettext("Iterations"), breaks = xBreaks,
                                limits = range(xBreaks),
                                expand = ggplot2::expand_scale(mult = c(0.05, 0.1)))

  return(jaspGraphs::themeJasp(g))

}
