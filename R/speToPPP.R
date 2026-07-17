### Code from spatialFDA written by Martin Emons, Samuel Gunz, Fabian Scheipl
### and Mark Robinson licensed under GPLv3 ###

#' Convert SpatialExperiment object to ppp object
#'
#' @param df `DataFrame`; x and y coordinates from the corresponding
#' SpatialExperiment and the ColData
#' @param marks `character`; the column with the labels e.g. cell types
#' @param continuous `logical`; indicating whether the marks are continuous
#' defaults to FALSE
#' @param window `owin`; An observation window of the point pattern of class.
#' @return `ppp`; object for use with `spatstat` functions
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
.dfToppp <- function(df, marks = NULL, continuous = FALSE, window = NULL) {
    # type checking
    stopifnot(is(df, "data.frame"))
    # this definition of the window is quite conservative
    # - can be set explicitly
    pp <- spatstat.geom::as.ppp(data.frame(x = df$x, y = df$y),
        W = spatstat.geom::owin(
            c(
                base::min(df$x) - 1.49e-8,
                base::max(df$x) + 1.49e-8
            ),
            c(
                base::min(df$y) - 1.49e-8,
                base::max(df$y) + 1.49e-8
            )
        )
    )
    # set the marks
    if (!continuous) {
        spatstat.geom::marks(pp) <- factor(df[[marks]])
    } else {
        spatstat.geom::marks(pp) <- base::subset(df,
            select =
                names(df) %in% marks
        )
    }
    # if window exist, set is as new window and potentially exclude some points
    if (!is.null(window)) {
        pp <- spatstat.geom::as.ppp(spatstat.geom::superimpose(pp, W = window))
    }

    return(pp)
}

#' Transform a SpatialExperiment into a dataframe
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#'
#' @return `DataFrame`; x and y coordinates from the corresponding
#' SpatialExperiment and the colData
#' @importFrom methods is
.speToDf <- function(spe) {
    stopifnot(is(spe, "SpatialExperiment"))
    df <- data.frame(
        x = SpatialExperiment::spatialCoords(spe)[, 1],
        y = SpatialExperiment::spatialCoords(spe)[, 2]
    )
    df <- cbind(df, colData(spe))
}

### new code ###

#' convert a SpatialExperiment object into a point pattern
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#' @param marks `character`; the column with the labels e.g. cell types
#' @param continuous `logical`; indicating whether the marks are continuous
#' defaults to FALSE
#' @param window `owin`; An observation window of the point pattern of class.
#'
#' @return A ppp object for use with any `spatstat` package
#' @export
#'
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "6")
#' pp <- speToPPP(speSub, mark = "cellType")
#'
speToPPP <- function(spe, marks, continuous = FALSE, window = NULL) {
    df <- .speToDf(spe)
    pp <- .dfToppp(df, continuous = continuous, window = window, marks = marks)
}
