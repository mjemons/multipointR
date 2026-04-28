#' Fit point process models across all images in the
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
#' @param formula `formula`; the formula to pass to the `ppm` function
#' @param family `detpointprocfamily`; Family to use in the point process model.
#' One of `dppGauss`, `dppMatern`, `dppCauchy`, `dppBessel` or `dppPowerExp`
#' @param threshold `numeric`; a threshold to apply on the minimum number of
#' points a point pattern needs to have to fit a `ppm` model to it.
#' @param cellspacing `numeric` how much spacing should be accounted for in the 
#' Hardcore process due to the cell body. If this is not provided, the cell spacing
#' parameter is estimated from the data
#' @param ncores `numeric`; the number of cores to used for parallel processing
#' @param verbose `logical`; whether to print informations on the fitting
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
#'                 formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)"),
#'                 threshold = 10)
#' 
#' @importFrom mgcv s
#' @importFrom spatstat.model dppm
#' @importFrom spatstat.model kppm
fitModelAcrossImages <- function(spe,
                                model = "ppm",
                                imageId,
                                imageLs = NULL,
                                marks,
                                formula,
                                family = spatstat.model::dppGauss(),
                                threshold = NULL,
                                cellspacing = NA,
                                ncores = 1,
                                verbose = TRUE,
                                ...){
  if(is.null(imageLs)){
    imageLs <- spe[[imageId]] |> unique() |> as.factor()
  }

  mdlLs <- parallel::mclapply(imageLs, function(image){
    speSub <- spe[, colData(spe)[[imageId]] == image]
    if(verbose){
      message(paste0("Fitting ", model, " to image ", image))
    }
    mdl <- fitModel(spe = speSub,
                    model = model,
                    marks = marks,
                    formula = formula,
                    family = family,
                    threshold = threshold,
                    cellspacing = cellspacing,
                    ...)
    return(mdl)
  }, mc.cores = ncores)
  return(mdlLs)
}
