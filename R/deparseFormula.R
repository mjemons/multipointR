#' deparse the Formula for spatstat and extract the relevant data
#'
#' @param spe `SpatialExperiment`; object subset to a single image
#' @param formula `formula`; the formula to pass to the `ppm` function
#' @param response `character`; the response for the process, the lhs of the formula object
#' @param marks `character`; the column with the labels e.g. cell types
#' @param lambda `im` or `NULL`; offset of the intensity to include in the model to account for
#' the underlying inhomogeneity. If this is not user provided, will be estimated
#'  as inhomogeneous intensity via diggle correction
#' @param threshold `numeric`; a threshold to apply on the minimum number of
#' points a point pattern needs to have to fit a `ppm` model to it.
#'
#' @returns named `list` with both the updated formula and the data for fitting
#'
#' @export
#' @examples
#' spe <- SpatialDatasets::spe_Keren_2018()
#' speSub <- subset(spe, , imageID == "5")
#' formula = stats::as.formula("Keratin_Tumour ~ distfun(CD8_T_cell)")
#' #define the response
#' response <- as.character(formula.tools::lhs(formula))
#' 
#' #deparse the Formula and extract the data
#' out <- deparseFormula(spe = speSub,
#'                       response = response,
#'                       formula = formula,
#'                       marks = "cellType", 
#'                       lambda = NULL,
#'                       threshold = 10)
#' 
deparseFormula <- function(spe,
                          formula,
                          response,
                          marks,
                          lambda,
                          threshold){
  #convert spe to ppp object
  pp <- speToPPP(spe, marks = marks)
  #subset the ppp object to the response mark
  ppResponse <- spatstat.geom::unmark(pp[pp$marks %in% response, drop = TRUE])
  
  if(!is.null(threshold) && spatstat.geom::npoints(ppResponse) <= threshold){
    message(paste0("There were less than ",  threshold, " points to compute an intensity on"))
        return(NULL)
  }
  #separate random from fixed effects for formula deparsing
  randomEffects <- reformulas::findbars(formula)
  formula     <- reformulas::nobars(formula)
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
    lambda  <- stats::density(ppResponse, sigma = spatstat.explore::bw.ppl(ppResponse), positive = TRUE)
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
      # formula <- stats::as.formula(paste(formula.tools::lhs.vars(formula)," ~ ", 
      # attr(stats::terms(formula), "variables")[attr(stats::terms(formula), "offset")+1] ,
      # "+", paste(formulaVars, collapse = "+")), env = baseenv())
      ### code optimised by claude.ai
      offset_idx <- attr(stats::terms(formula), "offset")
      offsetTerm <- if (!is.null(offset_idx)) {
        deparse(attr(stats::terms(formula), "variables")[[offset_idx + 1]])
        } else {
        NULL
        }

      rhs <- paste(c(offsetTerm, formulaVars), collapse = " + ")
        
      formula <- stats::as.formula(
        paste(formula.tools::lhs.vars(formula), "~", rhs),
        env = parent.frame()
      )
    }
  }
  #extract model variables that are in the point pattern marks but have not
  #been added as a covariate to the data
  missingVars <- rhsVars[rhsVars %in% levels(pp$marks) == FALSE &
                          rhsVars %in% names(data) == FALSE]

  if (length(missingVars) >= 1) {
    message(paste0("The covariate(s) ", paste(missingVars, collapse=", "), " is missing"))
    return(NULL)
  }
  ### code optimised with claude.ai
  # Handle the random effect grouping variable separately
  if (!is.null(randomEffects)) {
    groupVar <- trimws(gsub(".*\\|", "", deparse(randomEffects[[1]])))
    if (groupVar %in% colnames(colData(spe))) {
      data[[groupVar]] <- unique(colData(spe)[[groupVar]])
    } else {
      message(paste0("Random effect grouping variable '", groupVar, "' not found in colData"))
      return(NULL)
    }
    reFormula <- stats::as.formula(
                paste(". ~ . +", paste0("(", sapply(randomEffects, deparse), ")", collapse = " + ")), env = parent.frame()
    )
    
    formula <- stats::update(formula, reFormula)
  }
  return(list(formula = formula, data = data))
}