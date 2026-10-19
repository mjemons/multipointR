#' Convert the list of spatial models to a dataframe
#'
#' @param mdlLs `list`; A list of all the models fit with e.g. `ppm`
#' @param imageCovariates `list`; A list of all the covariates on the
#' image level to add to the model dataframe
#' @param fine `logical` Logical flag passed to `summary.ppm` controlling
#' the coarse or finer estimate of the covariance matrix
#' @param fast `logical` Logical flag passed to `summary.multipointRppm` 
#' controlling whether tapering (using sparse matrices from Matrix) 
#' should be used to speed up calculations.
#' Warning: This is expected to underestimate the true asymptotic 
#' variances/covariances.
#'
#' @returns DataFrame of the model coefficients
#'
#' @export
#' @examples
#' \dontrun{
#' spe <- SpatialDatasets::spe_Keren_2018()
#'
#' mdlLs <- fitModelAcrossImages(
#'     spe = spe,
#'     imageId = "imageID",
#'     imageLs = list("1"),
#'     marks = "cellType",
#'     sharedModel = FALSE,
#'     formula = as.formula(
#'         "Keratin_Tumour ~ spatstat.geom::distfun(CD8_T_cell)"
#'     ),
#'     threshold = 10,
#'     method.kppm = "mincon" #to reduce compute time
#' )
#'
#' mdlDf <- mdlToDf(
#'     mdlLs = mdlLs,
#'     imageCovariates = c(
#'         "imageID",
#'         "tumour_type"
#'     )
#' )
#'}
#' @importFrom dplyr bind_rows
mdlToDf <- function(mdlLs,
    imageCovariates = c("imageID"),
    fine=FALSE,
    fast=TRUE) {
    dfTotal <- data.frame()
    dfTotal <- lapply(mdlLs, function(mdl) {
        if (is.null(mdl)) {
            return(NULL)
        } else {
            dfCoef <- stats::coef(summary(mdl, fine=fine,fast=fast))
            dfCoef$covariate <- rownames(dfCoef)
            imageCovariateDf <- SummarizedExperiment::colData(mdl$spe) |>
                subset(, colnames(SummarizedExperiment::colData(mdl$spe)) 
                    %in% imageCovariates) |>
                unique()
            dfCoef <- cbind(dfCoef, imageCovariateDf)
            return(dfCoef)
        }
    }) |> dplyr::bind_rows()
    return(dfTotal)
}
