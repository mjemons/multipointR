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
#' `Strauss`, or `StraussHard` process to the data
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
 

  #convert spe to ppp object
  pp <- speToPPP(spe, marks = marks)
  response <- as.character(formula.tools::lhs(formula))

  #subset the ppp object to the response mark
  ppResponse <- spatstat.geom::unmark(pp[pp$marks %in% response, drop = TRUE])
  #calculate the minimum nearest neighbour distance if `is.null(cellspacing)`
  ### adapted from spatstat.model::Hardcore GPL-2 licensed
  if(length(cellspacing)>0 || is.na(cellspacing)){
    minNnDist <- minnndist(ppResponse)
    nX <- npoints(ppResponse)
    cellspacing <- minNnDist * nX/(nX+1)
  }
  ### end of directly adapted code ### 
  
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

  #calculate inhomogeneous intensity offset if not provided
  #if this is not provided, calculate it, else take the user
  #provided offset
  if(is.null(lambda)){
    lambda  <- stats::density(ppResponse, positive = TRUE)
    data[["lambda"]] <- lambda
  }else{
    data[["lambda"]] <- lambda
  }

  for(var in levels(pp$marks)){
    if(var %in% rhsVars){
      #convert to pp object
      ppVar <- spatstat.geom::unmark(pp[pp$marks %in% var, drop = TRUE])
      #get position in the rhs list
      position <- which(all.names(rhs) == var)
      #get the function that is being applied to the object
      fun <- all.names(rhs)[position-1]
      #apply the function and store in data
      data[[paste0(fun,".",var,".")]] = do.call(fun, 
        args = list(X=ppVar,
                    x = ppVar))
      #rename the function accordingly
      deparsed <- formula.tools::rhs.vars(formula)
      #get the position in the function
      formula_position <- which(grepl(var, deparsed))
      #string split the variable 
      splitVar <- strsplit(deparsed[formula_position], "[()]")[[1]]
      #remove empty strings
      splitVar <- splitVar[nzchar(splitVar)]
      #remove the already transformed part
      if(length(splitVar)>2){
        transform <- splitVar[-c(length(splitVar)-1, length(splitVar))]
        stopifnot("Only one composite function can be passed in this function" = length(transform) == 1)
        #replace function name
        newVar <- paste0(transform,"(",fun,".",var,".",")")
      }else{
        newVar <- paste0(fun,".",var,".")
      }
      #get the rhs formula vars
      formulaVars <- formula.tools::rhs.vars(formula)
      #overwrite the variable at the changed position
      formulaVars[formula_position] <- newVar
      #new formula
      formula <- stats::as.formula(paste(formula.tools::lhs.vars(formula)," ~ ", attr(stats::terms(formula), "variables")[attr(stats::terms(formula), "offset")+1] ,"+", paste(formulaVars, collapse = "+")))
    }
  }
  #extract model variables that are in the point pattern marks but have not
  #been added as a covariate to the data
  missingVars <- rhsVars[rhsVars %in% levels(pp$marks) == FALSE & 
                          rhsVars %in% names(data) == FALSE]
  if(length(missingVars)>=1){
    message(paste0("The covariate(s) ", missingVars, " is missing"))
        return(NULL)
  }
  #fit the model. If there is a cellspacing value indicated, this will be the
  #Gibbs point process with a Hardcore spacing between points. Else, this
  #value is estimated from the data
  #TODO: Make a try catch in case there is a fit failure to return NULL instead
  #of breaking the entire process.
  if(interaction == "Strauss"){
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
  mdl <- do.call(model, 
    args = list(Q = formula,
                family = family,
                data = data,
                interaction = interactionModel,
                ...)
  )
  #add covariates to the mdl list
  mdl$colData <- colData(spe)
  return(mdl)
}
