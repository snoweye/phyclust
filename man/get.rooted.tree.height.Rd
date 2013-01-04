\name{get.rooted.tree.height}
\alias{get.rooted.tree.height}
\title{Get a Rooted Tree Height}
\description{
  This function gets a rooted tree height, and only meaningful for
  a \strong{ultrametric} tree which has the equal height from the root to
  all leaves.
}
\usage{
get.rooted.tree.height(rooted.tree,
    tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{rooted.tree}{a rooted tree in \code{Class phylo}.}
  \item{tol}{for \code{is.ultrametric} of \pkg{ape}.}
}
\details{
  The \code{rooted.tree} should be in a \code{phylo} class of \pkg{ape},
  and should be ultrametric that may be created by \code{ms}.
}
\value{
  Return the rooted tree height.
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
  \code{\link[ape]{is.ultrametric}},
  \code{\link{rescale.rooted.tree}}.
}
\examples{
\dontrun{
set.seed(1234)
ret.ms <- ms(5, 1, opts = paste("-T", sep = " "))
tree.ms <- read.tree(text = ret.ms[3])
is.ultrametric(tree.ms)
get.rooted.tree.height(tree.ms)
}
}
\keyword{programming}
