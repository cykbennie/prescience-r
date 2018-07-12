#' @title PRESCIENCE estimates for coefficients
#'
#' @description
#' \code{coef.select} produces the estimated coefficients of the PRESCIENCE
#' produced by the \code{\link{select}} function.
#'
#' @param object an object of class \code{select}.
#' @param ... additional parameters.
#' @method coef select
#' @export
#' @return the input object is returned silently.
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @examples
#' results <- select(auto ~ dcost + cars + dovtt + divtt,
#' data = transportation, nfoc = 1, q = 1, bound = 10)
#'
#' coef(results)


coef.select <- function(object, ...){
	if (!inherits(object, "select"))
		stop("Object must be of class 'select'")

    bhat <- object$bhat
    return(bhat)
}
