\name{.show.option}
\alias{.show.option}
\title{Show Available Options}
\description{
  This function show available options for functions in \pkg{phyclust}.
}
\usage{
  .show.option()
}
\details{
  This function show some available options for functions in \pkg{phyclust}.
  They are used in \code{.EMControl}, \code{phyclust}, ... etc, and
  options are stored in several objects separately. They will be passed 
  into C, so the elemental order are important. Basically, they are all
  read-only objects.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.boundary.method}},
  \code{\link{.code.type}},
  \code{\link{.edist.model}},
  \code{\link{.em.method}},
  \code{\link{.EMC}},
  \code{\link{.EMControl}},
  \code{\link{.identifier}},
  \code{\link{.init.method}},
  \code{\link{.init.procedure}},
  \code{\link{.nucleotide}},
  \code{\link{.snp}},
  \code{\link{.substitution.model}},
}
\examples{
\dontrun{
.show.option()
}
}
\keyword{programming}
