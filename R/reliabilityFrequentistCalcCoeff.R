

# .reliabilityFrequentistAlpha <- function(jaspResults, dataset, options, model) {
#   out <- list()
#   # is alpha even checked?
#   if (options[["alphaScale"]]) {
#
#     if (!is.null(.getStateContainer(jaspResults)[["alphaObj"]]$object))
#       return(.getStateContainer(jaspResults)[["alphaObj"]]$object)
#
#     # alpha unstandardized
#     if (options[["alphaMethod"]] == "alphaUnstand") {
#
#       out[["est"]] <- Bayesrel:::applyalpha(model[["data_cov"]])
#
#       # do we need an interval estimate?
#       if (options[["intervalOn"]]) {
#         ciValue <- options[["confidenceIntervalValue"]]
#
#         # should the interval be analytic
#         if (options[["alphaInterval"]] == "alphaAnalytic") {
#           out[["conf"]] <- Bayesrel:::ciAlpha(1 - ciValue, nrow(dataset), model[["data_cov"]])
#
#         } else {
#           alpha_samp <- model[["sampCov"]]
#           if (is.null(alpha_samp)) {
#             alpha_samp <- apply(model[["bootsamp"]], 1, Bayesrel:::applyalpha)
#             out[["sampCov"]] <- alpha_samp
#           }
#
#           out[["conf"]] <- quantile(alpha_samp, probs = c((1-ciValue)/2, 1-(1-ciValue)/2))
#         }
#       }
#
#       # do we have to compute item dropped values
#       if (options[["alphaItem"]]) {
#         out[["itemDropped"]] <- apply(model[["itemDroppedCovs"]], 1, Bayesrel:::applyalpha)
#       }
#
#
#     } else { # alpha standardized
#       ccor <- cov2cor(model[["data_cov"]])
#       out[["est"]] <- Bayesrel:::applyalpha(ccor)
#
#       # do we need an interval estimate?
#       if (options[["intervalOn"]]) {
#         ciValue <- options[["confidenceIntervalValue"]]
#
#         # should the interval be analytic
#         if (options[["alphaInterval"]] == "alphaAnalytic") {
#           out[["conf"]] <- Bayesrel:::ciAlpha(1 - ciValue, nrow(dataset), ccor)
#
#         } else {
#           alpha_samp <- model[["sampCor"]]
#
#           if (is.null(alpha_samp)) {
#             alpha_samp <- numeric(options[["noSamples"]])
#
#             for (i in 1:options[["noSamples"]]) {
#               alpha_samp[i] <- Bayesrel:::applyalpha(cov2cor(model[["bootsamp"]][i, ,]))
#             }
#             out[["sampCor"]] <- alpha_samp
#
#           }
#           out[["conf"]] <- quantile(alpha_samp, probs = c((1-ciValue)/2, 1-(1-ciValue)/2), na.rm=T)
#
#         }
#       }
#
#       # do we have to compute item dropped values
#       if (options[["alphaItem"]]) {
#         out[["itemDropped"]] <- numeric(ncol(dataset))
#         for (i in 1:ncol(dataset)){
#           out[["itemDropped"]][i] <- Bayesrel:::applyalpha(ccor[-i, -i])
#         }
#       }
#     }
#
#     if (!is.null(.getStateContainer(jaspResults))) {
#       alphaObject <- createJaspState(out)
#       # alphaObject$dependOn(...)
#       .getStateContainer(jaspResults)[["alphaObj"]] <- alphaObject
#     }
#   }
#
#   return(out)
# }


.reliabilityFrequentistGuttman2 <- function(jaspResults, dataset, options, model) {
  
  if (!is.null(.getStateContainer(jaspResults)[["guttman2Obj"]]$object))
      return(.getStateContainer(jaspResults)[["guttman2Obj"]]$object)
  
  out <- model[["guttman2"]]
  if (is.null(out))
    out <- list()
  # is coefficient even checked?
  if (options[["guttman2Scale"]]) {

    out[["est"]] <- Bayesrel:::applyalpha(model[["data_cov"]])

    # do we need an interval estimate?
    if (options[["intervalOn"]]) {

      ciValue <- options[["confidenceIntervalValue"]]
      guttman2_samp <- apply(model[["bootsamp"]], 1, Bayesrel:::applylambda2)
      out[["sampCov"]] <- guttman2_samp
      out[["conf"]] <- quantile(guttman2_samp, probs = c((1-ciValue)/2, 1-(1-ciValue)/2))

    }

    # do we have to compute item dropped values
    if (options[["guttman2Item"]]) {
      out[["itemDropped"]] <- apply(model[["itemDroppedCovs"]], 1, Bayesrel:::applylambda2)
    }

    stateContainer <- .getStateContainer(jaspResults)
    stateContainer[["guttman2Obj"]] <- createJaspState(out)
  }
  return(out)
}
