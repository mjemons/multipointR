#' Fit determinantal point process model across all images in the
#' `SpatialExperiment` object
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#' @param imageId `character`; column in the colData of the `SpatialExperiment`
#' object specifying the image ID
#' @param imageLs `list`; a list specifying the subset of all images. If `NULL`
#' then all images are considered
#' @param model `character`; the `spatstat.model` function to use for computation.
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
#' @param ncores `numeric`; the number of cores to used for parallel processing
#' @param ... other parameters passed on to `dppm` model from `spatstat.model`
#'
#' @returns `list`; result from a `dppm` model in `spatstat.model`
#' @export
#'
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#'
#' mdlLs <- fitModelAcrossImages(spe = spe,
#'                 imageId = "imageID",
#'                 imageLs = list("1", "2"),
#'                 marks = "cellType",
#'                 response = "Keratin_Tumour",
#'                 distanceTo = "CD8_T_cell",
#'                 threshold = 10)
#' @importFrom splines bs
#' @importFrom spatstat.model dppm
#' @importFrom spatstat.model kppm
fitModelAcrossImages <- function(spe,
                                 model = "ppm",
                                 imageId,
                                 imageLs = NULL,
                                 marks,
                                 response,
                                 distanceTo = NULL,
                                 polygon = NULL,
                                 inhomogeneous = FALSE,
                                 family = spatstat.model::dppGauss(),
                                 threshold = NULL,
                                 ncores = 1,
                                 ...){
  if(is.null(imageLs)){
    imageLs <- spe[[imageId]] |> unique() |> as.factor()
  }

  mdlLs <- parallel::mclapply(imageLs, function(image){
    speSub <- spe[, colData(spe)[[imageId]] == image]
    mdl <- fitModel(spe = speSub,
                    model = model,
                    marks = marks,
                    response = response,
                    distanceTo = distanceTo,
                    polygon = polygon,
                    inhomogeneous = inhomogeneous,
                    family = family,
                    threshold = threshold,
                    ...)
    return(mdl)
  }, mc.cores = ncores)
  return(mdlLs)
}
