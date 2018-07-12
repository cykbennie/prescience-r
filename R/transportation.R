#' @title Transportation Mode Choice by Horowitz (1993)
#' 
#' @description A dataset of the transportation mode choice containing
#' 842 observations sampled randomly from Washington, D.C., area
#' transportation study.
#' 
#' @docType data
#' @keywords datasets
#' @name transportation
#' @usage transportation
#' @format A data frame with 842 rows and 5 variables. The variables
#' are as follows:
#' \describe{
#'   \item{auto}{Transportation mode choice with 1 for automobile and 0 otherwise}
#'   \item{dcost}{Transit fare - automobile travel cost (in dollars)} 
#'   \item{cars}{Number of cars owned by the traveler's household}
#'   \item{dovtt}{Transit out-of-vehicle travel time - automobile out-of-vehicle travel time (in minutes)}
#'   \item{divtt}{Transit in-vehicle travel time - automobile in-vehicle travel time (in minutes)}
#' }
#' @source The \code{transportation} data were obtained from
#' \emph{Semiparametric estimation of a work-trip mode choice model}
#' by Joel L.Horowitz (1993).
#' \url{https://www.sciencedirect.com/science/article/pii/030440769390113J}.
NULL