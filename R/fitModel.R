#' Fit point process models to one FOV of a `SpatialExperiment`
#' object
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#' @param model `character`; the `spatstat.model` function to
#' use for computation. Typically one of `ppm`, `kppm` or `dppm`.
#' @param marks `character`; the column with the labels e.g. cell types
#' @param formula `formula`; the formula to pass to the `ppm` function
#' @param family `detpointprocfamily`; Family to use in the point process model.
#' One of `dppGauss`, `dppMatern`, `dppCauchy`, `dppBessel` or `dppPowerExp`
#' @param threshold `numeric`; a threshold to apply on the minimum number of
#' points a point pattern needs to have to fit a `ppm` model to it.
#' @param interaction `character`; Formula specifying whether to fit
#'  a `Hardcore`, `Strauss`, `Fiksel` or `StraussHard` process to the data
#' @param lambda `im` or `NULL`; offset of the intensity to include in the
#' model to account for the underlying inhomogeneity.
#' If this is not user provided, will be estimated as inhomogeneous
#' intensity via diggle correction
#' @param cellspacing `numeric` how much spacing should be accounted for in the
#' interaction process due to the cell body. If this is not provided,
#' the cell spacing parameter is estimated from the data as the minimum
#' nearest neighbour distance divided by $n(n+1)$ as done in
#' `spatstat.model::Hardcore`. In the case of `StraussHard`,
#' only the Strauss interaction radius can be user-provided, the Hardcore
#' interaction is always estimated from the data.
#' @param improve.type `character` the improve.type argument 
#' from spatstat.model::ppm.ppp, passed on directly
#' @param relaxed `logical` whether or not to perform a relaxed fit with only
#' the non zero coefficients from `glmnet` improvement
#' @param selectionExclude `character`; Formula term specified to be removed
#' from the selection of the Lasso. This is necessary with 
#' `improve.type = "enet"`, because else the $p$-values will be not valid due
#' to obvious double dipping. The inference then becomes an LRT with and without
#' the variable indicated in `selectionExclude`.
#' @param ... other parameters passed on to `ppm` model from `spatstat.model`
#'
#' @returns `list`; result from a `dppm` model in `spatstat.model`
#' @export
#'
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#'
#' mdl <- fitModel(
#'     spe = speSub,
#'     marks = "cellType",
#'     formula = as.formula(
#'         "Keratin_Tumour ~ spatstat.geom::distfun(CD8_T_cell)"
#'     )
#' )
#'
#' @importFrom mgcv s
#' @importFrom spatstat.model dppm
#' @importFrom spatstat.model kppm
#' @importFrom spatstat.model ppm
fitModel <- function(
    spe,
    model = "ppm",
    marks,
    formula,
    family = spatstat.model::dppGauss(),
    threshold = NULL,
    interaction = "Fiksel",
    lambda = NULL,
    cellspacing = NA,
    improve.type = NULL,
    relaxed = FALSE,
    selectionExclude = NULL,
    ...
) {
    # some type assertions
    stopifnot(
        is(spe, "SpatialExperiment"),
        is.character(marks),
        is(formula, "formula"),
        is.null(threshold) || is.numeric(threshold),
        is.character(interaction) || is.null(interaction),
        is.null(lambda) || is(lambda, "im"),
        is.na(cellspacing) || is.numeric(cellspacing),
        is.null(selectionExclude) || is.character(selectionExclude)
    )

    # we do not need the assays anymore, therefore we set them to NULL
    SummarizedExperiment::assays(spe) <- list()
    # for computational reasons, remove the rowData as we don't need them
    SummarizedExperiment::rowData(spe) <- S4Vectors::DataFrame(
        row.names = rownames(spe)
    )

    if (identical(improve.type, "enet")) {
        requireNamespace("glmnet", quietly = TRUE)
    }
    # define the response
    response <- as.character(formula.tools::lhs(formula))

    # if we do selection of variables with `improve.type = "enet"` we need
    # to exclude the inferential variable because else the p-values will be
    # not valid
    if(improve.type == "enet" && !is.null(selectionExclude)){
        termLabels <- base::attr(stats::terms(formula), "term.labels")
        
        missingTerms <- setdiff(selectionExclude, termLabels)

        if (length(missingTerms) > 0) {
            stop(
                "The following `selectionExclude` terms are not present ",
                "in the model formula: ",
                paste(missingTerms, collapse = ", "),
                ". Supply multiple terms as a character vector, e.g. ",
                '`c("log(x)", "density(Endothelial)")`.'
            )
        }
        selectionTerms <- setdiff(termLabels, selectionExclude)
        #keep the full formula for later
        fullFormula <- formula
        formula <- stats::reformulate(selectionTerms, 
            response = response
        )

        # deparse the full Formula and extract the data - we need this to deparse
        # also the excluded term
        outFull <- deparseFormula(
            spe = spe,
            response = response,
            formula = fullFormula,
            marks = marks,
            lambda = lambda,
            threshold = threshold
        )
        fullFormula <- outFull[["formula"]]
        fullData <- outFull[["data"]]
    }else{
        fullFormula <- NULL
        fullData <- NULL
    }
    # deparse the Formula and extract the data
    out <- deparseFormula(
        spe = spe,
        response = response,
        formula = formula,
        marks = marks,
        lambda = lambda,
        threshold = threshold
    )

    # if the output of the deparsing is NULL, return a NULL model
    if (is.null(out)) {
        return(NULL)
    }

    formula <- out[["formula"]]
    data <- out[["data"]]

    # small hack to overwrite the reduced dataframe with the full dataframe
    # to be able to refit later
    if(!is.null(fullData)){
        data <- fullData
    }

    # parametrise the interaction model
    interactionModel <- defineInteractionModel(
        interaction = interaction,
        cellspacing = cellspacing,
        response = response,
        data = data
    )

    # fix for `sf` object evaluation as suggested by Adrian Baddeley
    transformSf <- function(z) {
        lapply(z, function(x) {
            if (inherits(x, "sf")) spatstat.geom::as.owin(x) else x
        })
    }

    mdl <- do.call(model,
        args = list(
            Q = formula,
            data = transformSf(data),
            interaction = interactionModel,
            improve.type = improve.type,
            ...
        )
    )
    #if we fit an elastic net, some coefficients can be zeroed out
    #in that case it can be advantageous to refit the model with only
    #the non-zero coefficients
    #coded with claude.ai
    if(improve.type == "enet" && relaxed == TRUE){
        allCoefs <- stats::coef(mdl)
        selectionFormula <- mdl$trend 
        # extract the model matrix and the terms from the formula
        mm <- stats::model.matrix(mdl)
        termLabels <- base::attr(stats::terms(selectionFormula), "term.labels")
        assignVec  <- base::attr(mm, "assign")  

        # keep any enet fit coefficient which is greater zero -> if one 
        # spline bases is >0 then the entire spline basis will be kept
        keepTerms <- c()
        for (i in base::seq_along(termLabels)) {
            cols <- base::which(assignVec == i)
            blockCoefs <- allCoefs[base::names(allCoefs) 
                %in% base::colnames(mm)[cols]]
            if (length(blockCoefs) > 0 && any(blockCoefs != 0)) {
            keepTerms <- c(keepTerms, termLabels[i])
            }
        }
        
        # problem with pre-evaluated spatstat function handling
        # improved by GPT 5.6
        if(!is.null(fullFormula)){
            fullTerms <- base::attr(stats::terms(fullFormula), "term.labels")
            excludedTerms <- setdiff(fullTerms, termLabels)
        }else{
            excludedTerms <- character()
        }
        # add the selectionExclude argument back
        refitTerms <- unique(c(
            keepTerms,
            excludedTerms
        ))

        # build the formula from the intact terms not from the model matrix
        # due to the spline bases
        refitFormula <- stats::reformulate(refitTerms, 
            response = response
        )

        # refit unpenalised with only the non-zero coefficients
        mdl <- stats::update(mdl, Q = refitFormula, improve.type = "none")
        message("Model was fit with relaxed enet. The number of coefficients can
        therefore be different than an unregularised fit. The post-selection
        p-values are only approximate")
    }
    # add covariates to the mdl list
    mdl$colData <- colData(spe)
    class(mdl) <- c("multipointRppm", class(mdl))
    return(mdl)
}


#' Plot mulitpointR `ppm` objects
#'
#' @param x `ppm`; a model fit with `spatstat.model::ppm`
#' or the wrapper `multipointR::fitModel`
#' @param type `character`; the type to plot, one of "intensity" or "trend"
#' @param ... further arguments passed to geom_raster
#'
#' @returns a `ggplot2` object of the model trend/intensity surfance
#'
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#'
#' mdl <- fitModel(
#'     spe = speSub,
#'     marks = "cellType",
#'     formula = as.formula(
#'         "Keratin_Tumour ~ spatstat.geom::distfun(CD8_T_cell)"
#'     )
#' )
#' plot(mdl)
#' @export
#' @method plot multipointRppm
#' @importFrom rlang .data
plot.multipointRppm <- function(x, type = "trend", ...) {
    ### coded with the help of claude.ai ###
    stopifnot(spatstat.geom::verifyclass(x, "ppm"))
    # extract the response `ppp` object
    pp_df <- as.data.frame(x$Q$data)
    # extract the trend image
    mdl_img <- stats::predict(x, type = type)
    # convert the image to a dataframe
    mdl_df <- as.data.frame((mdl_img))

    p <- ggplot2::ggplot(
        mdl_df,
        ggplot2::aes(x = .data[["x"]], y = .data[["y"]])
    ) +
        ggplot2::geom_raster(ggplot2::aes(fill = .data[["value"]])) +
        ggplot2::scale_fill_viridis_c(option = "magma", name = type) +
        ggplot2::geom_point(
            data = pp_df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
            shape = 1,
            size = 1.5,
            color = "white",
            stroke = 0.15
        ) +
        ggplot2::coord_equal() +
        ggplot2::theme_light() +
        ggplot2::labs(
            title =
                paste0("Fitted ", type, " surface"), x = "x", y = "y"
        )

    return(p)
}
