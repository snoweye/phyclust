\name{.EMC}
\alias{.EMC}
\docType{data}
\title{EM Control}
\description{
  An default template object stores controlling options for \code{phyclust}
  to perform EM algorithms. This object combines all other read-only
  objects and more required options for EM algorithms. This is essential
  for \code{phyclust} and other related functions.
}
\usage{.EMC}
\format{
  A list contains all controlling options
}
\details{
  A list created by \code{.EMControl} contains all controlling options for
  EM algorithms. This list will be directly passed to C codes and control
  the every things of EM algorithms.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.EMC
}
}
\keyword{datasets}
