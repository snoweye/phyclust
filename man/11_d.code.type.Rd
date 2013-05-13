\name{.code.type}
\alias{.code.type}
\docType{data}
\title{Code Types of Dataset and Substitution Models}
\description{
  Indicate the types of codes for datasets and functions.
  The first element is the default value.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.code.type}
\format{
  A character vector contains implemented code types in C.
}
\details{
  Two possible types are implemented,
  "NUCLEOTIDE" (default) and "SNP", used in
  data transfers and indicating substitution models.
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
  \code{\link{.substitution.model}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.code.type
}
}
\keyword{data}
