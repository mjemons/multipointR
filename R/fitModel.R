#' Fit a determinantal point process model to one FOV of a `SpatialExperiment`
#' object
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#' @param fun `character`; the `spatstat.model` function to use for computation.
#' Typically one of `ppm`, `kppm` or `dppm`.
#' @param marks `character`; the column with the labels e.g. cell types
#' @param response `character`; The mark whos intensity is modelled as response
#' @param distanceTo `character` | `owin`; optional, an character specifyng
#' the mark of the ppp to which the distance `distfun` from `spatstat.geom`
#' shall be computed. Alternatively, this can be a segmented object passed as
#' owin
#' @param polygon `owin`; optional, a polygon to add as a covariate to the data
#' @param inhomogeneous `logical`; Whether a correction for inhomogeneous
#' distribution of cell types with a B-spline should be perfomred.
#' @param family `detpointprocfamily`; Family to use in the point process model.
#' One of `dppGauss`, `dppMatern`, `dppCauchy`, `dppBessel` or `dppPowerExp`
#' @param threshold `numeric`; a threshold to apply on the minimum number of
#' points a point pattern needs to have to fit a `dppm` model to it.
#' @param ... other parameters passed on to `dppm` model from `spatstat.model`
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
#'                 response = "Keratin_Tumour",
#'                 distanceTo = "CD8_T_cell")
#' @importFrom splines bs
#' @importFrom spatstat.model dppm
#' @importFrom spatstat.model kppm
#' @importFrom spatstat.model ppm
fitModel <- function(spe,
                     fun = "dppm",
                     marks,
                     response,
                     distanceTo = NULL,
                     polygon = NULL,
                     inhomogeneous = FALSE,
                     family = spatstat.model::dppGauss(),
                     threshold = NULL,
                     ...){
  #some type assertions
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(is(marks, "character"))
  stopifnot(is(response, "character"))
  #stopifnot(is(distanceTo, "character" || "owin"))
  stopifnot(is(inhomogeneous, "logical"))

  #convert spe to ppp object
  pp <- speToPPP(spe, marks = marks)
  #subset the ppp object to the response mark
  ppResponse <- spatstat.geom::unmark(pp[pp$marks %in% response, drop = TRUE])
  if(!is.null(threshold) && spatstat.geom::npoints(ppResponse) <= threshold){
    return(paste0("There were less than ",  threshold, "points to fit"))
  }
  data <- list(ppResponse = ppResponse)

  #in case that distanceTo is not null compute it
  if(!is.null(distanceTo)){
    #if it is a character, it is a mark in the point pattern
    if(is(distanceTo, "character")){
      object <- pp[pp$marks %in% distanceTo, drop = TRUE]
      if(!is.null(threshold) && spatstat.geom::npoints(object) <= threshold){
        return(paste0("There were less than ",  threshold, "points to compute
                      distances too"))
      }
    }
    #else it is a shape which
    else if(is(distanceTo, "owin")){
      object <- distanceTo
    }
    #if it is neither, quite with an error
    else{
      print("Error in computing distfun to an object that is neither a mark in
            the ppp an owin object")
      return(NULL)
    }
    distanceFun <- spatstat.geom::distfun(object)
    data <- c(data, distanceFun = distanceFun)
  }
  if(length(data)>1){
    #create formula object
    formula <- stats::as.formula(paste("ppResponse ~ 1 + ", paste(names(data)[c(-1)],
                                                                  collapse="+")))
  }
  else{
    #create formula object
    formula <- stats::as.formula(paste("ppResponse ~ 1" ))
  }
  #correct for spatial inhomogeneity with a spline basis for x and y
  if(inhomogeneous){
    formula <- stats::as.formula(paste("ppResponse ~ 1 + ", paste(names(data)[c(-1)],
                                                             collapse="+"),
                                       "+bs(x,5)+bs(y,5)"))
  }
  #fit the model
  mdl <- do.call(fun,
                 args = list(formula = formula,
                             family = family,
                             data = data,
                             ...)
  )
  return(mdl)
}
