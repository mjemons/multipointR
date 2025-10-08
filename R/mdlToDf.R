
#' Convert the list of spatial models to a dataframe
#'
#' @param mdlLs `list`; A list of all the models fit with e.g. `ppm`
#' @param imageCovariates `list`; A list of all the covariates on the
#' image level to add to the model dataframe
#'
#' @returns DataFrame of the model coefficients
#'
#' @export
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
#' mdlDf <- mdlToDf(mdlLs = mdlLs,
#'                  imageCovariates = c("imageID",
#'                                      "tumour_type"))
#' 
mdlToDf <- function(mdlLs, 
                    imageCovariates = c("imageID")){
  dfTotal <- data.frame()
  dfTotal <- lapply(mdlLs, function(mdl){
    if(is.null(mdl)){
      return(NULL)
    }else{
      dfCoef <- stats::coef(summary(mdl))
      dfCoef$covariate <- rownames(dfCoef)
      imageCovariateDf <- mdl$colData |> subset(,colnames(mdl$colData) %in% imageCovariates) |> unique()
      dfCoef <- cbind(dfCoef, imageCovariateDf)
      return(dfCoef)
    }
  }) |> dplyr::bind_rows()
  return(dfTotal)
}