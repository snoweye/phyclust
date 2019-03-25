\name{gen.star.tree}
\alias{gen.star.tree}
\title{Generate a Rooted Tree with a Star Shape}
\description{
  Generate a rooted tree with a star shape based on a sequence calls
  of several functions.
}
\usage{
gen.star.tree(N, total.height = 1)
}
\arguments{
  \item{N}{number of leaves.}
  \item{total.height}{total tree height.}
}
\details{
  A tree with a star shape means that all internal branches are 0 and
  all leaf branches are equal.

  This function combining with \code{gen.seq.HKY} or \code{gen.seq.SNP}
  is used in simulation studies and bootstrap tree samples 
}
\value{
  Return a rooted tree in \code{Class phylo} with a star shape.
}
\references{
  Phylogenetic Clustering Website:
  \url{https://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{ms}}.
  \code{\link{as.star.tree}},
  \code{\link{get.rooted.tree.height}}.
  \code{\link{rescale.rooted.tree}}.
  \code{\link[ape]{as.phylo}},
  \code{\link[ape]{plot.phylo}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

ret.star <- gen.star.tree(5)
plot(ret.star, type = "u")
}
}
\keyword{simulation}
