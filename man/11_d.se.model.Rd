\name{.se.model}
\alias{.se.model}
\docType{data}
\title{Sequencing Error Model}
\description{
  An object stores sequencing error models.
}
\usage{.se.model}
\format{
  A character vector contains all possible sequencing models. 
}
\details{
  Currently, only a \code{CONVOLUTION} model is implemented.
}
\references{
  Phylogenetic Clustering Website:
  \url{https://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.show.option}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.se.model
}
}
\keyword{datasets}
