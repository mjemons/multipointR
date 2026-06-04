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
    lambda  <- spatstat.explore::density.ppp(ppResponse, 
      sigma = spatstat.explore::bw.ppl(ppResponse), 
      positive = TRUE, 
      diggle = TRUE,
      edge = TRUE,
      kernel = "gaussian")
    data[["lambda"]] <- lambda
  }else{
    data[["lambda"]] <- lambda
  }
  Vars <- c()
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
      ### code by claude.ai
      collect_terms <- function(expr) {
        if (!is.call(expr)) {
          # Base case: bare variable
          return(deparse1(expr))
        }
        
        op <- deparse1(expr[[1]])
        
        if (op %in% c("+", "-", "*", "/", ":")) {
          # Binary operator — recurse into both sides
          unlist(lapply(as.list(expr[-1]), collect_terms))
        } else {
          # Function call like distfun(...) — return as-is
          deparse1(expr)
        }
      }
      ### end code by claude.ai
      deparsed <- unique(collect_terms(rhs))
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
      Vars <- c(Vars, var)
    }
  }
  ### written by claude.ai
  get_innermost_call <- function(term) {
    expr <- parse(text = trimws(term))[[1]]
    
    # recurse until no more calls
    while (is.call(expr)) {
      # find the first argument that is itself a call
      inner <- Filter(is.call, as.list(expr[-1]))
      if (length(inner) == 0) break
      expr <- inner[[1]]
    }
    return(deparse1(expr))
  }
  ### end code by claude.ai
  #deparse the formula
  deparsed <- formula |> formula.tools::rhs() |> deparse()
  #split the terms of the formula
  splitTerms <- strsplit(deparsed, "\\+")[[1]] |> trimws()
  for(var in Vars){
    #get all the instances where the evaluated spatstat function was called
    mask <- grepl(var,splitTerms)
    #split function from variable
    term <- splitTerms[mask]
    #split interactions
    terms <- strsplit(term, "[+:*/\\\\]")[[1]] |> trimws()
    #check again for var
    maskTerms <- grepl(var,terms)
    #extract the innermost function call -> this will be the spatstat function on the
    #ppp object
    inner <- get_innermost_call(terms[maskTerms])
    #exchange () with "." to signal that this was already evaluated
    inner_modified <- gsub("[()]", ".", inner)
    #put this back in the original function
    splitTerms[mask] <- sub(inner, inner_modified, term, fixed = TRUE)
  }
  formula <- stats::as.formula(paste(formula.tools::lhs(formula), 
                        "~", 
                        paste(splitTerms, collapse = " + ")),
                        env = parent.frame())
  #extract model variables that are not in the point pattern marks and have not
  #been added as a covariate to the data
  missingVars <- rhsVars[rhsVars %in% levels(pp$marks) == FALSE &
                          rhsVars %in% names(data) == FALSE]

  #check wether the missingVars are in the colData of the `spe` 
  #TODO: add an option to check here for `annotGeometries` in 
  #an `sfe` object -> could be a way to store segmented regions
  for (missingVar in missingVars) {
    #check wether the missingVar is in the colData of the `spe` 
    if (missingVar %in% colnames(colData(spe))) {
      data[[missingVar]] <- unique(colData(spe)[[missingVar]])
    } else {
      message(paste0("The covariate(s) ", missingVar, " is missing"))
      return(NULL)
    }
  }

  ### code optimised with claude.ai
  # Handle the random effect grouping variable separately
  if (!is.null(randomEffects)){
    groupVar <- trimws(gsub(".*\\|", "", deparse(randomEffects[[1]])))
    if (groupVar %in% colnames(colData(spe))) {
      data[[groupVar]] <- unique(colData(spe)[[groupVar]])
    }else{
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