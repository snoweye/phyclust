\name{rescale.rooted.tree}
\alias{rescale.rooted.tree}
\title{Rescale a Rooted Tree's Height}
\description{
  This function rescaled the input rooted tree height by a scale.height.
}
\usage{
rescale.rooted.tree(rooted.tree, scale.height = 1)
}
\arguments{
  \item{rooted.tree}{a rooted tree in \code{Class phylo}.}
  \item{scale.height}{a scale to all branch lengths.}
}
\details{
  The \code{rooted.tree} should be in a \code{phylo} class of \pkg{ape},
  and may be created by \code{ms}.

  \code{scale.height} is a positive number multiplying on the lengths of
  all branches of the rooted tree.
}
\value{
  Return a rooted tree in \code{Class phylo} with scaled branches.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{ms}},
  \code{\link[ape]{read.tree}},
  \code{\link[ape]{as.phylo}},
  \code{\link{get.rooted.tree.height}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
ret.ms <- ms(5, 1, opts = paste("-T", sep = " "))
tree.ms <- read.tree(text = ret.ms[3])
get.rooted.tree.height(tree.ms)

tree.scaled <- rescale.rooted.tree(tree.ms, scale.height = 2)
get.rooted.tree.height(tree.scaled)
}
}
\keyword{simulation}
