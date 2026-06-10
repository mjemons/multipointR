#' define the Interaction Model 
#'
#' @param interaction `character`; Formula specifying whether to fit a `Hardcore`,
#' `Strauss`, `Fiksel` or `StraussHard` process to the data
#' @param cellspacing `numeric` how much spacing should be accounted for in the 
#' interaction process due to the cell body. If this is not provided, the cell spacing
#' parameter is estimated from the data as the minimum nearest neighbour distance
#' divided by $n(n+1)$ as done in `spatstat.model::Hardcore`. In the case of `StraussHard`,
#' only the Strauss interaction radius can be user-provided, the Hardcore interaction is 
#' always estimated from the data.
#' @param response `character`; the response for the process, the lhs of the formula object
#' @param data `list` the interaction model for the Gibbs process
#'
#' @returns object of class `interact`
#'
#' @export
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#' formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)")
#' #define the response
#' response <- as.character(formula.tools::lhs(formula))
#' 
#' #deparse the Formula and extract the data
#' out <- deparseFormula(spe = speSub,
#'                       response = response,
#'                       formula = formula,
#'                       marks = "cellType", 
#'                       lambda = NULL,
#'                       threshold = 10)
#' 
#' interactionModel <- defineInteractionModel(interaction = "Hardcore",
#'                                            cellspacing = NA,
#'                                            response = response,
#'                                            data = out[["data"]])
#' 
defineInteractionModel <- function(interaction,
                                    cellspacing,
                                    response,
                                    data){
  #calculate the minimum nearest neighbour distance if `is.null(cellspacing)`
  ### adapted from spatstat.model::Hardcore GPL-2 licensed
  if(length(cellspacing)>0 || is.na(cellspacing)){
    minNnDist <- minnndist(data[[response]])
    nX <- npoints(data[[response]])
    cellspacing <- minNnDist * nX/(nX+1)
  }
  ### end of directly adapted code ### 
  if(is.null(interaction)){
    interactionModel <- interaction
    message(paste0("You have specified a Poisson process - make sure that the assumption 
    of physical overlap is justified in your data"))
  }
  else if(interaction == "Strauss"){
    interactionModel <- spatstat.model::Strauss(r = cellspacing)
  }else if(interaction == "Hardcore"){
    interactionModel <- spatstat.model::Hardcore(hc = cellspacing)
  }else if(interaction == "StraussHard"){
    #optimise the model parameters
    rs <- expand.grid(r=seq(cellspacing+0.1, cellspacing + 5, by=0.5),
                    hc=cellspacing)
    pg <- spatstat.model::profilepl(rs, spatstat.model::StraussHard, data[[response]], verbose = FALSE, fast = TRUE)
    interactionModel <- pg$fit$interaction
  }else if(interaction == "Fiksel"){
    #optimise the model parameters
    rs <- expand.grid(r=seq(cellspacing+0.1, cellspacing + 5, by=0.5),
                      hc=cellspacing,
                      kappa=seq(0.5,2, by=0.5))
    pg <- spatstat.model::profilepl(rs, spatstat.model::Fiksel, data[[response]], verbose = FALSE, fast = TRUE)
    interactionModel <- pg$fit$interaction
  }else{
    warning(paste0("Interaction model ", interaction, " is not implemented"))
  }
  return(interactionModel)
}