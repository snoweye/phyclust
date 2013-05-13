\name{.label.method}
\alias{.label.method}
\docType{data}
\title{Label Method}
\description{
  An object stores label method for un-, semi-, and general semi-
  supervised clustering..
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.label.method}
\format{
  A character vector contains implemented evolution distances in C.
}
\details{
  This vector stores possible label methods implemented in \pkg{phyclust}.
  The default value is \code{NONE} for unsupervised clustering.
  \code{SEMI} is for semi-supervised clustering, and
  \code{GENERAL} is for general semi-supervised clustering.
  Only un- and semi-supervised clustering are implemented.
}

\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.label.method
}
}
\keyword{data}
