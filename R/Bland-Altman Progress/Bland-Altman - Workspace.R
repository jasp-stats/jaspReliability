setwd('/Users/lucat/OneDrive/Dokumente/Uni Amsterdam/Internship/Example Data/')
write.csv(df,"BlandAtlman.csv", row.names = FALSE)

df <- read.csv("LikertData_7Lvl.csv", header = T)
df <- df[,-c(1,2)]
for (i in 1:ncol(df)) {
  df[, i] <- factor(df[, i], levels = 1:nLevels)
}


library(blandr)
?rnorm
set.seed(2)
measurement1 <- rnorm(100)
measurement2 <- rnorm(100)

set.seed(2)
statistics.results <- blandr.statistics(rnorm(100) , rnorm(100))
blandr.plot.ggplot(statistics.results)

View(blandr.statistics)
View(blandr.plot.ggplot)



library(BlandAltmanLeh)

bland.altman.plot(measurement1, measurement2, xlab="mean measurement",
                  ylab="differences", conf.int= .95, graph.sys = "ggplot2")
View(bland.altman.plot)
debug(bland.altman.plot)


View(bland.altman.ggplot2)

xx <- bland.altman.stats(measurement1, measurement2)
values <- data.frame(m = xx$means, d = xx$diffs)


# Function Content:

# Errors: (1) different length of variables & (2) obs below 2 not working

place <- function(measure1, measure2, ci = 0.95, ciDisplay = TRUE, ciShading = TRUE){

  df <- data.frame(group1 = measure1, group2 = measure2)
  df <- na.omit(df)

  # Main components
  diffs <- df[[1]] - df[[2]]
  means <- (df[[1]] + df[[2]])/2
  values <- data.frame(m = means, d = diffs)
  obs <- length(df[[1]])
  meanDiffs <- mean(diffs)
  criticalDiff <- 1.96 * sd(diffs) #1.96 used instead of suggested 2 by Bland & Altman
  lowerLimit <- meanDiffs - criticalDiff
  upperLimit <- meanDiffs + criticalDiff
  lines <- c(lowerLimit, meanDiffs, upperLimit)
  t1 <- qt((1 - ci)/2, df = obs - 1)
  t2 <- qt((ci + 1)/2, df = obs - 1)

  # CI calculations based on Bland & Altman, 1986
  lowerLimitCiLower <- lowerLimit + t1 * sqrt(sd(diffs)^2 * 3/obs)
  lowerLimitCiUpper <- lowerLimit + t2 * sqrt(sd(diffs)^2 * 3/obs)
  meanDiffCiLower <- meanDiffs + t1 * sd(diffs)/sqrt(obs)
  meanDiffCiUpper <- meanDiffs + t2 * sd(diffs)/sqrt(obs)
  upperLimitCiLower <- upperLimit + t1 * sqrt(sd(diffs)^2 * 3/obs)
  upperLimitCiUpper <- upperLimit + t2 * sqrt(sd(diffs)^2 * 3/obs)
  CiLines <- c(lowerLimitCiLower, lowerLimitCiUpper, meanDiffCiLower, meanDiffCiUpper, upperLimitCiLower, upperLimitCiUpper)

  # Bland-Altman Plot
  if(max(CiLines) > max(diffs) && min(CiLines) < min(diffs)){
    yBreaks <- jaspGraphs::getPrettyAxisBreaks(CiLines, n = 10)
  } else if(max(CiLines) < max(diffs) && min(CiLines) > min(diffs)){
    yBreaks <- jaspGraphs::getPrettyAxisBreaks(diffs, n = 10)
  } else if(max(CiLines) < max(diffs) && min(CiLines) < min(diffs)){
    yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(min(CiLines),diffs), n = 10)
  } else if(max(CiLines) > max(diffs) && min(CiLines) > min(diffs)){
    yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(min(diffs),CiLines), n = 10)
  }
  yBreaks <- c(floor(min(yBreaks)), yBreaks[!yBreaks%%1], ceiling(max(yBreaks)))
  xBreaks <- jaspGraphs::getPrettyAxisBreaks(means, n = 10)
  xBreaks <- c(floor(min(xBreaks)), xBreaks[!xBreaks%%1], ceiling(max(xBreaks)))

  p <- ggplot2::ggplot(values, ggplot2::aes(x = m, y = d)) +
    ggplot2::geom_point(size = 3) + ggplot2::geom_hline(yintercept = lines, linetype = 2, size = 1) +
    ggplot2::xlab("Mean of Measurements") +
    ggplot2::ylab("Difference of Measurements")

  if(ciDisplay){
    p <- p + ggplot2::geom_hline(yintercept = CiLines, linetype = 2, size = 0.5)

    if (ciShading){
      p <- p + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                                 ymin = CiLines[3],
                                 ymax = CiLines[4],
                                 fill = "blue", alpha = 0.3) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = CiLines[5],
                          ymax = CiLines[6],
                          fill = "green", alpha = 0.3) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = CiLines[1],
                          ymax = CiLines[2],
                          fill = "red", alpha = 0.3)
    }
  }

  p <- p + ggplot2::scale_y_continuous(breaks = yBreaks, limits = range(yBreaks),
                                       sec.axis = ggplot2::dup_axis(~., breaks = yBreaks, labels = NULL, name = "")) +
    ggplot2::theme(axis.ticks.y.right = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks = xBreaks, limits = range(xBreaks))+
    jaspGraphs::geom_rangeframe(sides = "rbl") +
    jaspGraphs::themeJaspRaw() +
    ggplot2::theme(plot.margin = ggplot2::margin(5))

  return(p)
}

place(measurement1, measurement2, ciDisplay = F)

set.seed(12)
place(rnorm(100), rnorm(100))

undebug(place)



yBreaks <- jaspGraphs::getPrettyAxisBreaks(diffs, n = 10)
yBreaks[!yBreaks%%1]
c(floor(min(yBreaks)), yBreaks[!yBreaks%%1], ceiling(max(yBreaks)))
jaspGraphs::getPrettyAxisBreaks(diffs, n = 10)
max(CiLines) > max(diffs)
min(CiLines) > min(diffs)
