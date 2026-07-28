#' deparse the Formula for spatstat and extract the relevant data
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#' @param formula `formula`; the formula to pass to the `ppm` function
#' @param response `character`; the response for the process,
#' the lhs of the formula object
#' @param marks `character`; the column with the labels e.g. cell types
#' @param lambda `im` or `NULL`; offset of the intensity to include in the
#' model to account for the underlying inhomogeneity.
#' If this is not user provided, will be estimated as inhomogeneous
#' intensity via diggle correction
#' @param threshold `numeric`; a threshold to apply on the minimum number of
#' points a point pattern needs to have to fit a `ppm` model to it.
#'
#' @returns named `list` with both the updated formula and the data for fitting
#'
#' @export
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#' formula <- stats::as.formula(
#'     "Keratin_Tumour ~ spatstat.geom::distfun(CD8_T_cell)"
#' )
#' # define the response
#' response <- as.character(formula.tools::lhs(formula))
#'
#' # deparse the Formula and extract the data
#' out <- deparseFormula(
#'     spe = speSub,
#'     response = response,
#'     formula = formula,
#'     marks = "cellType",
#'     lambda = NULL,
#'     threshold = 10
#' )
#'
deparseFormula <- function(
    spe,
    formula,
    response,
    marks,
    lambda,
    threshold
) {
    # convert spe to ppp object
    pp <- speToPPP(spe, marks = marks)
    # subset the ppp object to the response mark
    ppResponse <- spatstat.geom::unmark(pp[pp$marks %in% response, drop = TRUE])

    if (!is.null(threshold) &&
        spatstat.geom::npoints(ppResponse) <= threshold) {
        message(
            "There were less than ", threshold,
            " points to compute an intensity on"
        )
        return(NULL)
    }
    # separate random from fixed effects for formula deparsing
    randomEffects <- reformulas::findbars(formula)
    formula <- reformulas::nobars(formula)
    # rhs variables of the formula
    rhs <- formula.tools::rhs(formula)
    # separate variables from terms (variable plus functions)
    rhsVars <- all.vars(rhs)
    # extract and store response
    data <- list(
        x = SpatialExperiment::spatialCoords(spe)[
            ,
            SpatialExperiment::spatialCoordsNames(spe)[1]
        ],
        y = SpatialExperiment::spatialCoords(spe)[
            ,
            SpatialExperiment::spatialCoordsNames(spe)[2]
        ]
    )

    data[[response]] <- ppResponse

    # calculate inhomogeneous intensity offset if not provided
    # if this is not provided, calculate it, else take the user
    # provided offset
    if (is.null(lambda)) {
        lambda <- spatstat.explore::density.ppp(ppResponse,
            sigma = spatstat.explore::bw.ppl(ppResponse),
            positive = TRUE,
            diggle = TRUE,
            edge = TRUE,
            kernel = "gaussian"
        )
        data[["lambda"]] <- lambda
    } else {
        data[["lambda"]] <- lambda
    }
    Vars <- c()
    for (var in levels(pp$marks)) {
        if (var %in% rhsVars) {
            ppVar <- spatstat.geom::unmark(pp[pp$marks %in% var, drop = TRUE])

            funExpr <- get_fun_call(rhs, var)
            funObj <- eval(funExpr, envir = parent.frame())
            funLabel <- gsub("::", ".", deparse1(funExpr))

            data[[paste0(funLabel, ".", var, ".")]] <- do.call(
                funObj,
                args = list(X = ppVar, x = ppVar)
            )

            Vars <- c(Vars, var)
        }
    }
    # deparse the formula
    deparsed <- formula |>
        formula.tools::rhs() |>
        deparse1()
    # split the terms of the formula
    splitTerms <- strsplit(deparsed, "\\+")[[1]] |> trimws()
    for (var in Vars) {
        # get all the instances where the evaluated spatstat function was called
        mask <- grepl(var, splitTerms)
        # split function from variable
        term <- splitTerms[mask]
        # split interactions
        # line optimised by claude.ai
        terms <- strsplit(term, "\\+|\\*|/|\\\\|(?<!:):(?!:)",
            perl = TRUE
        )[[1]] |> trimws()
        # check again for var
        maskTerms <- grepl(var, terms)
        # extract the innermost function call
        # -> this will be the spatstat function on the
        # ppp object
        inner <- get_innermost_call(terms[maskTerms])
        # exchange () and :: with "." to signal that this was already evaluated
        inner_modified <- gsub("[()]|::", ".", inner)
        # put this back in the original function
        splitTerms[mask] <- sub(inner, inner_modified, term, fixed = TRUE)
    }
    formula <- stats::as.formula(
        paste(
            formula.tools::lhs(formula),
            "~",
            paste(splitTerms, collapse = " + ")
        ),
        env = parent.frame()
    )
    # extract model variables that are not in the point pattern
    # marks and have not been added as a covariate to the data
    missingVars <- rhsVars[rhsVars %in% levels(pp$marks) == FALSE &
        rhsVars %in% names(data) == FALSE]

    # check wether the missingVars are in the colData of the `spe`
    for (missingVar in missingVars) {
        # check wether the missingVar is in the colData of the `spe`
        if (missingVar %in% colnames(colData(spe))) {
            data[[missingVar]] <- unique(colData(spe)[[missingVar]])
            # if the provided object is a `sfe` check whether
            # the missing variable is in the annotation geometries
        } else if (is(spe, "SpatialFeatureExperiment") &&
            missingVar %in%
                names(SpatialFeatureExperiment::annotGeometries(spe))) {
            data[[missingVar]] <- SpatialFeatureExperiment::annotGeometry(
                spe, missingVar
            )
        } else {
            message("The covariate(s) ", missingVar, " is missing")
            return(NULL)
        }
    }

    ### code optimised with claude.ai
    # Handle the random effect grouping variable separately
    if (!is.null(randomEffects)) {
        groupVar <- trimws(gsub(".*\\|", "", deparse1(randomEffects[[1]])))
        if (groupVar %in% colnames(colData(spe))) {
            data[[groupVar]] <- unique(colData(spe)[[groupVar]])
        } else {
            message(
                "Random effect grouping variable '", groupVar,
                "' not found in colData"
            )
            return(NULL)
        }
        reFormula <- stats::as.formula(
            paste(". ~ . +", paste0("(", vapply(
                randomEffects,
                deparse, character(1)
            ), ")", collapse = " + ")),
            env = parent.frame()
        )

        formula <- stats::update(formula, reFormula)
    }
    return(list(formula = formula, data = data))
}
