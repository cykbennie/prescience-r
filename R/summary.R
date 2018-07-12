#' @title Summarize the PRESCIENCE results
#'
#' @description
#' \code{summary.select} summarizes the results of PRESCIENCE
#' produced by the \code{\link{select}} function.
#'
#' @param object an object of class \code{select}.
#' @param ... additional parameters.
#' @method summary select
#' @export
#' @return the input object is returned silently.
#' @author Yankang (Bennie) Chen <yankang.chen@@columbia.edu>
#' @examples
#' results <- select(auto ~ dcost + cars + dovtt + divtt,
#' data = transportation, nfoc = 1, q = 1, bound = 10)
#'
#' summary(results)


summary.select <- function(object, ...){
	if (!inherits(object, "select"))
		stop("Object must be of class 'select'")

	status <- object$status
	tol <- object$tolerance
    bhat <- object$bhat
    score <- object$score
    gap <- object$gap
    rtime <- object$rtime
    ncount <- object$ncount

	cat("Optimization returned status:", status, "\n")
	cat("Tolerance level:", tol, "\n\n")
	cat("Parameter estimates:\n")
	print(bhat)
	cat("\nGurobi score:", score, "\n")
	cat("MIO gap:", gap, "\n")
	cat("CPU time (in seconds):", rtime, "\n")
	cat("Branch-and-bound nodes:", ncount, "\n")
}
