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
