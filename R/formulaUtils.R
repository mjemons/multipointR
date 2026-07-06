### written by claude.ai
#' Function to obtain the innermost call of a 
#' formula expression
#'
#' @param term the term of the function to de-nest
#'
#' @returns the innermost call of a nested funciton
#'
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

### written by claude.ai

#' Title
#'
#' @param expr part of a formula object
#' @param var variable to extract
#'
#' @returns the extracted formula call of the bespoke variable
#'
get_fun_call <- function(expr, var) {
  if (!is.call(expr)) return(NULL)
  args <- as.list(expr)[-1]
  # does this call directly contain var as one of its arguments?
  for (a in args) {
    if (identical(a, as.symbol(var))) {
      return(expr[[1]])   # the function part, e.g. quote(spatstat.geom::distfun)
    }
  }
  # otherwise recurse into the arguments
  for (a in args) {
    res <- get_fun_call(a, var)
    if (!is.null(res)) return(res)
  }
  NULL
}
### end code by claude.ai