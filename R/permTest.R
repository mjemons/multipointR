#' Permutation test for a single point pattern LRT between the model and a null
#'
#' @param mdl `ppm` the original model
#' @param mdl0 `ppm` the null model
#' @param nsim `integer` the number of simulations
#' @param studentised `logical` whether or not to perform a studentisation; a 
#' division of the test statistic by an estimate of the  standard errors. This
#' makes the permutation distribution more stable to non-equal variances.
#'
#' @returns named `list` with the permutation distribution and the p-value
#' @export
#' 
#' @examples
#' 
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#'
#' mdl <- fitModel(
#'    spe = speSub,
#'    marks = "cellType",
#'    formula = as.formula(
#'        "Keratin_Tumour ~ spatstat.geom::distfun(CD8_T_cell)"
#'    )
#' )
#' mdl0 <- fitModel(
#'    spe = speSub,
#'    marks = "cellType",
#'    formula = as.formula(
#'        "Keratin_Tumour ~ 1"
#'    )
#' )
#'
#' permTest(mdl, mdl0, nsim = 10)

permTest <- function(mdl, 
  mdl0, 
  nsim = 99,
  studentised = TRUE){
  #obtain the original data
  #pp <- spatstat.model::data.ppm(mdl)
  pp <- speToPPP(spe = mdl$spe, marks = mdl$marks)
  #recover the response 
  call <- stats::getCall(mdl)$Q
  response <- as.character(formula.tools::lhs(call))
  #perform the permutations
  sims <- spatstat.random::rlabel(pp, nsim = nsim)

  #written with the help of Opus 5
  Dobs <- 2 * (as.numeric(stats::logLik(mdl)) - as.numeric(stats::logLik(mdl0)))
  Dsim <- sapply(sims, function(Y) {
    Xp <- spatstat.geom::unmark(Y[Y$marks %in% response, drop = TRUE])
    f1 <- stats::update(mdl, Xp)
    f0 <- stats::update(mdl0, Xp)
    Dexp <- 2 * (as.numeric(stats::logLik(f1)) - as.numeric(stats::logLik(f0)))
    return(Dexp)
  })

  pVal <- (1 + sum(Dsim >= Dobs, na.rm = TRUE)) / (1 + nsim)
  return(list(Dobs = Dobs, Dsim = Dsim, pVal = pVal))
}