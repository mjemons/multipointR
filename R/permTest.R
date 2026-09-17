#' Permutation test for a single point pattern LRT between the model and a null
#'
#' @param mdl `ppm` the original model
#' @param coefficient `character` the coefficient(s) to test
#' @param nsim `integer` the number of simulations
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
#'
#' permTest(mdl, coefficient = "spatstat.geom::distfun(CD8_T_cell)", nsim = 10)

permTest <- function(mdl, 
  coefficient,
  nsim = 99){
  #obtain the original data
  #pp <- spatstat.model::data.ppm(mdl)
  pp <- speToPPP(spe = mdl$spe, marks = mdl$marks)
  #recover the response 
  call <- stats::getCall(mdl)$Q
  response <- as.character(formula.tools::lhs(call))
  #perform the permutations
  sims <- spatstat.random::rshift(pp, nsim = nsim)

  #need to rewrite the coefficient to correspond to `multipointR` internals
  coefficient <- gsub("[()]|::", ".", coefficient)
  fm0 <- stats::as.formula(paste("~ . -", coefficient))
  mdlP <-  stats::update(mdl, interaction = NULL)
  mdl0 <- stats::update(mdl, fm0, interaction = NULL)

  #written with the help of Opus 5
  Dobs <- 2 * (as.numeric(stats::logLik(mdlP)) - as.numeric(stats::logLik(mdl0)))
  Dsim <- sapply(sims, function(Y) {
    Xp <- spatstat.geom::unmark(Y[Y$marks %in% response, drop = TRUE])
    f1 <- stats::update(mdlP, Xp)
    f0 <- stats::update(mdl0, Xp)
    Dexp <- 2 * (as.numeric(stats::logLik(f1)) - as.numeric(stats::logLik(f0)))
    return(Dexp)
  })

  pVal <- (1 + sum(Dsim >= Dobs, na.rm = TRUE)) / (1 + nsim)
  return(list(Dobs = Dobs, Dsim = Dsim, pVal = pVal))
}