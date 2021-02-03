

.BayesianPosteriorPlot <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["plotContainer"]]$object))
    return()

  plotContainer <- createJaspContainer(gettext("Posterior Plots"))
  plotContainer$dependOn(options = c("plotPosterior", "shadePlots", "probTable", "probTableValueLow",
                                     "probTableValueHigh", "fixXRange", "dispPrior", "credibleIntervalValueScale",
                                     "alphaScale", "omegaScale", "lambda2Scale", "lambda6Scale", "glbScale"))

  derivedOptions <- model[["derivedOptions"]]
  selected <- derivedOptions[["selectedEstimatorsPlots"]]
  indices   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]
  idMatchedNames <- which(!is.na(charmatch(names(model), names(selected[indices]))))

  if (options[["shadePlots"]] && options[["probTable"]])
    shadePlots <- c(options[["probTableValueLow"]], options[["probTableValueHigh"]])
  else
    shadePlots <- NULL

  if (options[["plotPosterior"]] && is.null(model[["empty"]])) {
    n.item <- model[["k"]]
    priors <- Bayesrel:::priors[[as.character(n.item)]]
    # the prior names dont match the model names, thus rename the priors in their original order
    names(priors) <- c("alphaScale", "lambda2Scale", "lambda6Scale", "glbScale", "omegaScale")

    z <- 1
    for (i in indices) {
      if (is.null(plotContainer[[nmsObjsNoGreek[i]]])) {
        prior <- priors[[grep(names(model)[idMatchedNames[z]], names(priors))]]
        p <- .makeSinglePosteriorPlot(model[[idMatchedNames[z]]], nmsLabs[[i]],
                                                         options[["fixXRange"]], shadePlots,
                                                         options[["dispPrior"]], prior)
        plotObj <- createJaspPlot(plot = p, title = nmsObjs[i])
        plotObj$dependOn(options = names(indices[i]))
        plotObj$position <- i
        plotContainer[[nmsObjsNoGreek[i]]] <- plotObj

      }
      z <-  z+1
    }
    plotContainer$position <- 4
    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["plotContainer"]] <- plotContainer
  }

  return()
}

.makeSinglePosteriorPlot <- function(coefList, nms, fixXRange, shade = NULL, priorTrue, priorSample) {

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
  scaleCriRound <- round(coefList[["cred"]], 3)
  datCri <- data.frame(xmin = scaleCriRound[1L], xmax = scaleCriRound[2L], y = .925 * ymax)
  height <- (ymax - .925 * ymax) / 2
  if (fixXRange) {
    if (datCri$xmin == datCri$xmax) { # if two zeros, the interval is merged together
      datTxt <- data.frame(x = c(datCri$xmin, datCri$xmax),
                           y = 0.985 * ymax,
                           label = sapply(c(datCri$xmin, datCri$xmax), format, digits = 3, scientific = -1),
                           stringsAsFactors = FALSE)
    } else {
      datTxt <- data.frame(x = c(datCri$xmin -.08, datCri$xmax + .08),
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
  if (fixXRange && max(datDens$y, na.rm=T) >= 1000) xExpand <- c(xExpand[1], .05)
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
      datFilter <- data.frame(x=0, y=0)
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
  # fixes issue that unchecking the scale coefficient box, does not uncheck the item-dropped coefficient box:
  for (i in 1:5) {
    if (!derivedOptions[["selectedEstimators"]][i]) {
      derivedOptions[["itemDroppedSelectedItem"]][i] <- derivedOptions[["selectedEstimators"]][i]
    }
  }

  selected <- derivedOptions[["itemDroppedSelectedItem"]]
  indices   <- which(selected)
  nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
  nmsObjs   <- derivedOptions[["namesEstimators"]][["tables_item"]]
  nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]
  idMatchedNames <- which(!is.na(charmatch(names(model), names(selected[indices]))))


  if (options[["orderItem"]]) {
    ordering <- options[["orderType"]]
  } else {
    ordering <- NULL
  }


  if (is.null(model[["empty"]]) && options[["plotItem"]]) {
    z <- 1
    for (i in indices) {
      if (is.null(plotContainerItem[[nmsObjsNoGreek[i]]])) {

        # use the coefficient item object in the model list, and the object directly before it,
        # which is always the corresponding scale object
        p <- .makeIfItemPlot(model[[idMatchedNames[z]]], model[[idMatchedNames[z]-1]], nmsLabs[[i]],
                                                options[["credibleIntervalValueItem"]],
                                                ordering = ordering, model[["itemsDropped"]])
        plotObjItem <- createJaspPlot(plot = p, title = nmsObjs[i], width = 400)
        plotObjItem$dependOn(options = names(indices[i]))
        plotObjItem$position <- i
        plotContainerItem[[nmsObjsNoGreek[i]]] <- plotObjItem

      }
      z <- z+1
    }
    plotContainerItem$position <- 5
    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["plotContainerItem"]] <- plotContainerItem
  }

  return()
}

.makeIfItemPlot <- function(coefItem, coefScale, nms, int, ordering, variables) {
  n_row <- length(variables)
  lower <- (1-int)/2
  upper <- int + (1-int)/2

  samp_tmp <- as.vector(coefScale[["samp"]])
  dat <- data.frame(as.matrix(samp_tmp), row.names =  NULL)
  names(dat) <- "value"
  dat$colos <- "1"
  dat$var <- "original"

  item_tmp <- apply(coefItem[["itemSamp"]], 3, as.vector)
  dat_del <- t(as.matrix(as.data.frame(item_tmp)))
  names <- decodeColNames(variables)

  for (i in n_row:1){
    tmp <- as.data.frame(dat_del[i, ])
    colnames(tmp) <- "value"
    tmp$var <- names[i]
    tmp$colos <- "2"
    dat <- rbind(dat, tmp)
  }
  dat$var <- factor(dat$var, levels = unique(dat$var))

  if (!is.null(ordering)) {
    est <- as.data.frame(coefItem[["itemEst"]])
    est[n_row + 1, ] <- 1
    colnames(est) <- "value"
    est$name <- c(names, "original")

    if (ordering == "orderItemMean") {
      dists <- abs(coefScale[["est"]] - coefItem[["itemEst"]])
      dists[length(dists)+1] <- 0
      est <- est[order(dists, decreasing = F), ]
      dat$var <- factor(dat$var, levels = c(est$name))

    } else if (ordering == "orderItemKL") {
      samps <- item_tmp
      og_samp <- samp_tmp
      dists <- apply(samps, 2, .KLD.statistic, y = og_samp) # kl divergence
      dists[length(dists)+1] <- 0
      est <- est[order(dists), ]
      dat$var <- factor(dat$var, levels = c(est$name))

    } else if (ordering == "orderItemKS") {
      samps <- item_tmp
      og_samp <- samp_tmp
      dists <- apply(samps, 2, .ks.test.statistic, y = og_samp) # ks distance
      dists[length(dists)+1] <- 0
      est <- est[order(dists), ]
      dat$var <- factor(dat$var, levels = c(est$name))
    }
  }

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = value, y = var, fill = colos)) +
    ggridges::stat_density_ridges(quantile_lines = T, quantiles = c(lower, 0.5, upper),
                                  alpha = .85, show.legend = F, scale = 1) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(colour = "black")) +
    ggplot2::ylab(gettext("Item Dropped")) +
    ggplot2::xlab(nms) +
    ggplot2::scale_fill_grey() +
    ggplot2::scale_y_discrete(expand = ggplot2::expand_scale(mult = c(0.1, 0.25)))

  return(jaspGraphs::themeJasp(g))

}



.BayesianPosteriorPredictive <- function(jaspResults, model, options) {

  if (!is.null(.getStateContainerB(jaspResults)[["omegaPPC"]]$object))
    return()

  if (!is.null(model[["omegaScale"]]) && options[["dispPPC"]] && options[["omegaScale"]]) {

    ll <- model[["omegaScale"]][["loadings"]]
    rr <- model[["omegaScale"]][["residuals"]]
    cimpl <- ll %*% t(ll) + diag(rr)
    cobs <- model[["data_cov"]]
    k <- ncol(cobs)
    eframe <- data.frame(number = seq(1, k), eigen_value = eigen(cobs)$values)
    ee_impl <- matrix(0, 1e3, k)
    for (i in 1:1e3) {
      dtmp <- MASS::mvrnorm(model[["n"]], rep(0, k), cimpl)
      ee_impl[i, ] <- eigen(cov(dtmp))$values
    }
    eframe$eigen_sim_low <- apply(ee_impl, 2, quantile, prob = .025)
    eframe$eigen_sim_up<- apply(ee_impl, 2, quantile, prob = .975)
    leg_pos <- (max(eframe$eigen_value) + min(eframe$eigen_value)) * .75
    yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, max(eframe$eigen_sim_up)))

    g <- ggplot2::ggplot(eframe, mapping = ggplot2::aes(x = number, y = eigen_value)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = eigen_sim_low, ymax = eigen_sim_up), color = "grey60", width = 0.2) +
      ggplot2::geom_point() +
      ggplot2::xlim(c(1, k)) +
      ggplot2::scale_y_continuous(name = gettext("Eigenvalue"), breaks = yBreaks, limits = range(yBreaks)) +
      ggplot2::xlab(gettext("Eigenvalue No.")) +
      ggplot2::scale_x_continuous(expand = ggplot2::expand_scale(mult = c(.1, .1)))

    g <- jaspGraphs::themeJasp(g)

    plot <- createJaspPlot(plot = g, title = "Posterior Predictive Check Omega", width = 350)
    plot$dependOn(options = c("dispPPC", "omegaScale"))

    plot$position <- 6
    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["omegaPPC"]] <- plot
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
    indices   <- which(selected)
    nmsLabs   <- derivedOptions[["namesEstimators"]][["plots"]]
    nmsObjs   <- derivedOptions[["namesEstimators"]][["tables"]]
    nmsObjsNoGreek   <- derivedOptions[["namesEstimators"]][["plotsNoGreek"]]
    idMatchedNames <- which(!is.na(charmatch(names(model), names(selected[indices]))))

    xlim <- (options[["noSamples"]] - options[["noBurnin"]]) / options[["noThin"]]
    z <- 1
    for (i in indices) {
      if (is.null(plotContainerTP[[nmsObjsNoGreek[i]]])) {

        p <- .makeTracePlot(model[[idMatchedNames[z]]], nmsLabs[[i]])
        plotObjTP <- createJaspPlot(plot = p, title = nmsObjs[i], width = 400)
        plotObjTP$dependOn(options = names(indices[i]))
        plotObjTP$position <- i
        plotContainerTP[[nmsObjsNoGreek[i]]] <- plotObjTP

      }
      z <- z+1
    }

    plotContainerTP$position <- 7
    stateContainerB <- .getStateContainerB(jaspResults)
    stateContainerB[["plotContainerTP"]] <- plotContainerTP

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

