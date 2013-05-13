\name{.identifier}
\alias{.identifier}
\docType{data}
\title{Identifiers for Evolution Models}
\description{
  Identifiers for evolution models identify the \eqn{Q_k}{Q_k} matrix and
  evolution time \eqn{t_k}{t_k} for subpopulations.
  The first element is the default value.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.identifier}
\format{
  A character vector contains implemented identifiers in C.
}
\details{
  Four major identifiers are implemented in C, \code{EE}, \code{EV},
  \code{VE}, and \code{VV}. The first letter indicates the structure
  for \eqn{Q_k}{Q_k} matrix, and the second letter indicates the
  evolution time \eqn{t_k}{t_k} for subpopulations. \code{E} and
  \code{V} indicate equal and varied for all subpopulations.
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
  \code{\link{.EMC}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.identifier
}
}
\keyword{data}
