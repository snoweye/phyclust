\name{.boundary.method}
\alias{.boundary.method}
\docType{data}
\title{Boundary Methods for Population Proportions}
\description{
  Methods used in EM Algorithms to deal with boundary problems of
  population proportions, \eqn{\eta_k}{eta_k}.
  The first element is the default value.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.boundary.method}
\format{
  A character vector contains implemented boundary methods in C.
}
\details{
  The boundary value 0 of the population proportions makes the log
  likelihood as -Inf. Since degeneracy of subpopulations can affect
  the maximizing processes in EM steps. This problem is usually caused
  by bad initializations, and may suggest that number of cluster
  \eqn{K}{K} may be too large.

  Two methods have been implemented when any \eqn{\eta_k}{eta_k} less than
  the lower bound (\eqn{1/N}{1/N} or \eqn{1e-16}{1e-16}).
  The \code{ADJUST} (default) will adjust the \eqn{\eta_k}{eta_k} to the
  lower bound, and the \code{IGNORE} will stop the iterations and return
  errors.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.init.procedure}},
  \code{\link{.init.method}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.boundary.method
}
}
\keyword{data}
