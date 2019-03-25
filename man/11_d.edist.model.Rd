\name{.edist.model}
\alias{.edist.model}
\docType{data}
\title{Evolution Distance Model}
\description{
  Evolution distances based on certain evolution models as in \pkg{ape}. 
  The implemented model is used in \code{phyclust.edist} and initializations
  of EM algorithms. The first element is the default value.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.edist.model}
\format{
  A character vector contains implemented evolution distances in C.
}
\details{
  This vector stores possible evolution distances implemented in \pkg{phyclust}.
  The default value is \code{D_JC69} computed form the probability of
  \emph{JC69} model.
}
\references{
  Phylogenetic Clustering Website:
  \url{https://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{phyclust.edist}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.edist.model
}
}
\keyword{data}
