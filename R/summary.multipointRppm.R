#' generic summary of multipointRppm objects
#'
#' @param object `ppm`; a model fit with `spatstat.model::ppm`
#' or the wrapper `multipointR::fitModel`
#' @param quick `logical` whether or not to calculate standard errors
#' @param fast `logical`; whether or not to use the default `summary.ppm`/ 
#' `summary.kppm` functions or make use of the sparse matrix calculation of 
#' `vcov.kppm`. Has no effect for non `kppm` objects and defaults to the 
#' normal `summary`
#' @param eps.rmax `numerical` from `vcov.kppm`, check directly
#' @param ... other parameters passed on to `summary`
#'
#' @returns A list with summary information, see 
#' \code{\link[spatstat.model]{summary.kppm}()} or 
#' \code{\link[spatstat.model]{summary.ppm}()}
#'
#' @export
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#'
#' mdl <- fitModel(
#'     spe = speSub,
#'     marks = "cellType",
#'     formula = as.formula(
#'         "Keratin_Tumour ~ spatstat.geom::distfun(CD8_T_cell)"
#'     ),
#'     method.kppm = "mincon" #to reduce compute time
#' )
#' 
#' summary(mdl)
#' @author Martin Emons, adapted from
#' \code{\link[spatstat.model]{summary.kppm}()} written by Adrian Baddeley, 
#' Rolf Turner and Ege Rubak
summary.multipointRppm <- function(
  object,
  quick=FALSE,
  fast=TRUE,
  eps.rmax=0.01,
  ...
){
  class(object) <- setdiff(class(object), "multipointRppm")

  # only use the fast path for kppm objects
  if (!(inherits(object, "kppm"))) {
      return(summary(object, quick = quick, ...))
  }
  
  if(fast && isFALSE(quick)){
    #written by claude Opus 5.5 and adapted code from spatstat.model::summary.kppm
    #licensed under GPL-3 to Adrian Baddeley, Rolf Turner and Ege Rubak
    result <- summary(object, quick=TRUE, ...)
    theta <- stats::coef(object)
    if (length(theta) > 0) {
      vc <- stats::vcov(object, fast=fast, eps.rmax = eps.rmax)
    }
    if (!is.null(vc)) {
                se <- if (is.matrix(vc)) 
                  sqrt(diag(vc))
                else if (length(vc) == 1) 
                  sqrt(vc)
                else NULL
            }
    if (!is.null(se)) {
                two <- stats::qnorm(0.975)
                lo <- theta - two * se
                hi <- theta + two * se
                zval <- theta/se
                pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)
                psig <- cut(pval, c(0, 0.001, 0.01, 0.05, 1), 
                  labels = c("***", "**", "*", "  "), include.lowest = TRUE)
                result$coefs.SE.CI <- data.frame(Estimate = theta, 
                  S.E. = se, CI95.lo = lo, CI95.hi = hi, Ztest = psig, 
                  Zval = zval)
                }
  }else{
    result <- summary(object, quick=quick, ...)
  }
  return(result)
}
