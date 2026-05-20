#' fit $n$ univariate Models for all $n$ images
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
#' @param lambda `im` or `NULL`; offset of the intensity to include in the model to account for
#' the underlying inhomogeneity. If this is not user provided, will be estimated
#'  as inhomogeneous intensity via diggle correction
#' @param interaction `character`; Formula specifying whether to fit a `Hardcore`,
#' `Strauss`, `StraussHard` or `Fiksel` process to the data
#' @param cellspacing `numeric` how much spacing should be accounted for in the 
#' Hardcore process due to the cell body. If this is not provided, the cell spacing
#' parameter is estimated from the data
#' @param ncores `numeric`; the number of cores to used for parallel processing
#' @param verbose `logical`; whether to print informations on the fitting
#' @param ... other parameters passed on to `dppm` model from `spatstat.model`
#'
#' @returns `mppm` object of the shared fit across all images
#'
#' @export
#' @examples 
#' spe <- SpatialDatasets::spe_Keren_2018()
#' 
#' out <- .fitSingelModelsPerImage(spe = spe,
#'                 imageId = "imageID",
#'                 imageLs = list("1", "2"),
#'                 marks = "cellType",
#'                 formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)"),
#'                 threshold = 10)
#' 
.fitSingelModelsPerImage <- function(spe,
                                model = "ppm",
                                imageId,
                                imageLs = NULL,
                                marks,
                                formula,
                                family = spatstat.model::dppGauss(),
                                threshold = NULL,
                                interaction = "StraussHard",
                                cellspacing = NA,
                                lambda = NULL,
                                ncores = 1,
                                verbose = TRUE,
                                ...){
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
                    interaction = interaction,
                    lambda = lambda,
                    cellspacing = cellspacing,
                    ncores = ncores,
                    verbose = verbose,
                    ...)
    return(mdl)
  }, mc.cores = ncores)
  return(mdlLs)
}

#' fit a shared univariate Models for all $n$ images simultaneously
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
#' @param interaction `character`; Formula specifying whether to fit a `Hardcore`,
#' `Strauss`, `StraussHard` or `Fiksel` process to the data
#' @param lambda `im` or `NULL`; offset of the intensity to include in the model to account for
#' the underlying inhomogeneity. If this is not user provided, will be estimated
#'  as inhomogeneous intensity via diggle correction 
#' @param cellspacing `numeric` how much spacing should be accounted for in the 
#' Hardcore process due to the cell body. If this is not provided, the cell spacing
#' parameter is estimated from the data
#' @param ncores `numeric`; the number of cores to used for parallel processing
#' @param verbose `logical`; whether to print informations on the fitting
#' @param ... other parameters passed on to `mppm` model from `spatstat.model`
#'
#' @returns list`; result from a `ppm` model in `spatstat.model`
#'
#' @export
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' 
#' out <- .fitSharedModelAcrossImages(spe = spe,
#'                 imageId = "imageID",
#'                 imageLs = list("1", "2"),
#'                 marks = "cellType",
#'                 formula = as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)"),
#'                 threshold = 10)
#' 
.fitSharedModelAcrossImages <- function(spe,
                                model = "ppm",
                                imageId,
                                imageLs = NULL,
                                marks,
                                formula,
                                family = spatstat.model::dppGauss(),
                                threshold = NULL,
                                interaction = "Fiksel",
                                cellspacing = NA,
                                lambda = NULL,
                                ncores = 1,
                                verbose = TRUE,
                                ...){
  #we do not need the assays anymore, therefore we set them to NULL
  SummarizedExperiment::assays(spe) <- list()
  #for computational reasons, remove the rowData as we don't need them
  SummarizedExperiment::rowData(spe) <- S4Vectors::DataFrame(row.names = rownames(spe))
  #create a hyperframe object
  list_of_lists <- lapply(imageLs, function(image){
    speSub <- spe[, colData(spe)[[imageId]] == image]
    response <- as.character(formula.tools::lhs(formula))

    #deparse the Formula and extract the data
    out <- deparseFormula(spe = speSub,
                          response = response,
                          formula = formula,
                          marks = marks,
                          lambda = lambda, 
                          threshold = threshold)
    #if the output of the deparsing is NULL, return a NULL model
    if(is.null(out)){
      return(NULL)
    }
    #parametrise the interaction model
    interactionModel <- defineInteractionModel(interaction = interaction,
                                              cellspacing = cellspacing,
                                              response = response,
                                              data = out[["data"]])
    
    return(c(out[["data"]], list(formula = out[["formula"]], interact = interactionModel)))
  })
  ### code from Claude.ai
  # Transpose: list of rows -> list of columns
  cols <- do.call(Map, c(list(list), list_of_lists))

  #the columns for anything else than `ppp` or `interaction` or `im` objects have to be flat
  #therefore, unlist them 
  cols <- lapply(cols, function(col) {
    if (all(sapply(col, function(x) length(x) == 1 && (is.factor(x) || is.character(x) || is.numeric(x))))) {
      return(unlist(col))
    } else {
      return(col)
    }
  })

  # Build hyperframe
  hf <- do.call(spatstat.geom::hyperframe, cols)
  # remove the coordinates as these are in the `ppp` object already
  hf[["x"]] <- NULL
  hf[["y"]] <- NULL
  ### end of code from Claude.ai
  formula <- hf[["formula"]] |> unique()
  #the formula is nested, take it apart
  formula <- formula[[1]]
  #extract the lhs of the formula as response
  response <- as.character(formula.tools::lhs(formula))
  #for mppm we need to separate fixed from random effects
  fixedEffects <- reformulas::nobars(formula)
  randomEffects <- reformulas::findbars(formula)
  if(is.null(randomEffects)){
    mdl <- spatstat.model::mppm(formula=formula, 
      data=hf, 
      interaction = as.hyperframe(Interaction = hf[["interact"]]),
      ...)
  }else{
    #build the two formulae for fixed and random effects
    feFormula <- stats::as.formula(paste(deparse(fixedEffects)), env = baseenv())
    reFormula <- stats::as.formula(
                paste("~", paste0(sapply(randomEffects, deparse), collapse = " + ")), env = baseenv()
    )
    mdl <- spatstat.model::mppm(formula=feFormula, 
      random = reFormula, data=hf,
      interaction = as.hyperframe(hf[["interact"]]),
      ...)
  }
  return(mdl)
}

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
#' @param lambda `im` or `NULL`; offset of the intensity to include in the model to account for
#' the underlying inhomogeneity. If this is not user provided, will be estimated
#'  as inhomogeneous intensity via diggle correction
#' @param interaction `character`; Formula specifying whether to fit a `Hardcore`,
#' `Strauss`, `StraussHard` or `Fiksel` process to the data
#' @param cellspacing `numeric` how much spacing should be accounted for in the 
#' Hardcore process due to the cell body. If this is not provided, the cell spacing
#' parameter is estimated from the data
#' @param ncores `numeric`; the number of cores to used for parallel processing
#' @param lambda `im` or `NULL`; offset of the intensity to include in the model to account for
#' the underlying inhomogeneity. If this is not user provided, will be estimated
#'  as inhomogeneous intensity via diggle correction
#' @param sharedModel `logical`; whether or not to fit one univariate model
#' per image ($n$ models in total) or estimate one shared univariate model across all images.
#' @param verbose `logical`; whether to print informations on the fitting
#' @param ... other parameters passed on to `dppm` model from `spatstat.model`
#'
#' @returns `list`; result from a `ppm` model in `spatstat.model` or a `mppm` model
#' @export
#'
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#'
#' mdl <- fitModelAcrossImages(spe = spe,
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
                                interaction = "Fiksel",
                                cellspacing = NA,
                                lambda = NULL,
                                sharedModel = TRUE,
                                ncores = 1,
                                verbose = TRUE,
                                ...){
  if(is.null(imageLs)){
    imageLs <- spe[[imageId]] |> unique() |> as.factor()
  }

  if(sharedModel){
    out <- .fitSharedModelAcrossImages(spe,
                                      model = "mppm",
                                      imageId = imageId,
                                      imageLs = imageLs,
                                      marks = marks,
                                      formula = formula,
                                      family = family,
                                      threshold = threshold,
                                      interaction = interaction,
                                      lambda = lambda,
                                      cellspacing = cellspacing,
                                      ncores = ncores,
                                      verbose = verbose,
                                      ...

    )
  }else{
    out <- .fitSingelModelsPerImage(spe,
                                    model = model,
                                    imageId = imageId,
                                    imageLs = imageLs,
                                    marks = marks,
                                    formula = formula,
                                    family = family,
                                    threshold = threshold,
                                    interaction = interaction,
                                    lambda = lambda,
                                    cellspacing = cellspacing,
                                    ncores = ncores,
                                    verbose = verbose,
                                    ...
    )
  }
  return(out)
}
