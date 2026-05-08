#' Fit point process models to one FOV of a `SpatialExperiment`
#' object
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#' @param model `character`; the `spatstat.model` function to use for computation.
#' Typically one of `ppm`, `kppm` or `dppm`.
#' @param marks `character`; the column with the labels e.g. cell types
#' @param formula `formula`; the formula to pass to the `ppm` function
#' @param family `detpointprocfamily`; Family to use in the point process model.
#' One of `dppGauss`, `dppMatern`, `dppCauchy`, `dppBessel` or `dppPowerExp`
#' @param threshold `numeric`; a threshold to apply on the minimum number of
#' points a point pattern needs to have to fit a `ppm` model to it.
#' @param interaction `character`; Formula specifying whether to fit a `Hardcore`,
#' `Strauss`, `Fiksel` or `StraussHard` process to the data
#' @param lambda `im` or `NULL`; offset of the intensity to include in the model to account for
#' the underlying inhomogeneity. If this is not user provided, will be estimated
#'  as inhomogeneous intensity via diggle correction
#' @param cellspacing `numeric` how much spacing should be accounted for in the 
#' interaction process due to the cell body. If this is not provided, the cell spacing
#' parameter is estimated from the data as the minimum nearest neighbour distance
#' divided by $n(n+1)$ as done in `spatstat.model::Hardcore`. In the case of `StraussHard`,
#' only the Strauss interaction radius can be user-provided, 
#' the Hardcore interaction is always estimated from the data.
#' @param ... other parameters passed on to `ppm` model from `spatstat.model`
#'
#' @returns `list`; result from a `dppm` model in `spatstat.model`
#' @export
#'
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#'
#' mdl <- fitModel(spe = speSub,
#'                 marks = "cellType",
#'                 formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)")
#' )
#' 
#' @importFrom mgcv s
#' @importFrom spatstat.model dppm
#' @importFrom spatstat.model kppm
#' @importFrom spatstat.model ppm
#' @import spatstat.geom
fitModel <- function(spe,
                     model = "ppm",
                     marks,
                     formula,
                     family = spatstat.model::dppGauss(),
                     threshold = NULL,
                     interaction = "StraussHard",
                     lambda = NULL, 
                     cellspacing = NA,
                     ...){
  #some type assertions
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(is(marks, "character"))

  #define the response
  response <- as.character(formula.tools::lhs(formula))
 
  #deparse the Formula and extract the data
  out <- deparseFormula(spe = spe,
                        response = response,
                        formula = formula,
                        marks = marks, 
                        lambda = lambda,
                        threshold = threshold)
  
  #if the output of the deparsing is NULL, return a NULL model
  if(is.null(out)){
    return(NULL)
  }
  
  formula <- out[["formula"]]
  data <- out[["data"]]
  
  #parametrise the interaction model
  interactionModel <- defineInteractionModel(interaction = interaction,
                                            cellspacing = cellspacing,
                                            response = response,
                                            data = data)

  mdl <- do.call(model, 
    args = list(Q = formula,
                data = data,
                interaction = interactionModel,
                ...)
  )
  #add covariates to the mdl list
  mdl$colData <- colData(spe)
  return(mdl)
}
