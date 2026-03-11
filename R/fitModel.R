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
#' @param cellspacing `numeric` how much spacing should be accounted for in the 
#' Hardcore process due to the cell body
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
                     cellspacing = NULL,
                     ...){
  #some type assertions
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(is(marks, "character"))
 

  #convert spe to ppp object
  pp <- speToPPP(spe, marks = marks)
  response <- as.character(formula.tools::lhs(formula))

  #subset the ppp object to the response mark
  ppResponse <- spatstat.geom::unmark(pp[pp$marks %in% response, drop = TRUE])
  if(!is.null(threshold) && spatstat.geom::npoints(ppResponse) <= threshold){
    message(paste0("There were less than ",  threshold, " points to compute an intensity on"))
        return(NULL)
  }
  #rhs variables of the formula
  rhs <- formula.tools::rhs(formula)
  #separate variables from terms (variable plus functions)
  rhsVars <- all.vars(rhs)

  #extract and store response
  data <- list(x = SpatialExperiment::spatialCoords(spe)[,SpatialExperiment::spatialCoordsNames(spe)[1]],
    y = SpatialExperiment::spatialCoords(spe)[,SpatialExperiment::spatialCoordsNames(spe)[2]])
  
  data[[response]] = ppResponse

  for(var in levels(pp$marks)){
    if(var %in% rhsVars){
      #convert to pp object
      pp = spatstat.geom::unmark(pp[pp$marks %in% var, drop = TRUE])
      #get position in the rhs list
      position <- which(all.names(rhs) == var)
      #get the function that is being applied to the object
      fun <- all.names(rhs)[position-1]
      #apply the function and store in data
      data[[paste0(fun,".",var,".")]] = do.call(fun, 
        args = list(X=pp,
                    x = pp))
      #rename the function accordingly
      deparsed <- formula.tools::rhs.vars(formula)
      #get the position in the function
      formula_position <- which(grepl(var, deparsed))
      #replace function name
      newVar <- gsub(paste("[()]"), ".", deparsed[formula_position])
      #get the rhs formula vars
      formulaVars <- formula.tools::rhs.vars(formula)
      #overwrite the variable at the changed position
      formulaVars[formula_position] <- newVar
      #new formula
      formula <- stats::as.formula(paste(formula.tools::lhs.vars(formula)," ~ ", paste(formulaVars, collapse = "+")))
    }
  }
  #fit the model. If there is a cellspacing value indicated, this will be a
  #Gibbs point process with a Hardcore spacing between points.
  if(!is.null(cellspacing)){
    mdl <- do.call(model, 
      args = list(Q = formula,
                  family = family,
                  data = data,
                  interaction = spatstat.model::AreaInter(cellspacing),
                  use.gam = TRUE,
                  ...)
    )
  }else{
    #else Poisson process without spacing
    mdl <- do.call(model, 
      args = list(Q = formula,
                  family = family,
                  data = data,
                  use.gam = TRUE,
                  ...)
    )
  }
  #add covariates to the mdl list
  mdl$colData <- colData(spe)
  return(mdl)
}
