% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.select}
\alias{summary.select}
\title{Summarize the PRESCIENCE results}
\usage{
\method{summary}{select}(object, ...)
}
\arguments{
\item{object}{an object of class \code{select}.}

\item{...}{additional parameters.}
}
\value{
the input object is returned silently.
}
\description{
\code{summary.select} summarizes the results of PRESCIENCE
produced by the \code{\link{select}} function.
}
\examples{
results <- select(auto ~ dcost + cars + dovtt + divtt,
data = transportation, nfoc = 1, q = 1, bound = 10)

summary(results)
}
\author{
Yankang (Bennie) Chen <yankang.chen@columbia.edu>
}
