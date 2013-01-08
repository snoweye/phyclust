\name{as.star.tree}
\alias{as.star.tree}
\alias{Class phylo}
\title{Coerce a Rooted Tree to a Star Tree in Class phylo}
\description{
  Coerce a rooted tree generating by \code{ms} to a star tree and
  maintain a bifurcation structure.
}
\usage{
as.star.tree(rooted.tree, keep.bifurcation = TRUE)
}
\arguments{
  \item{rooted.tree}{a rooted tree in \code{Class phylo}.}
  \item{keep.bifurcation}{keep a bifurcation structure.}
}
\details{
  A tree with a star shape means that all internal branches are 0 and
  all leaf branches are equal.

  The \code{rooted.tree} should be in a \code{phylo} class of \pkg{ape},
  and may be created by \code{ms}.

  Basically, it is a list with an attribute that the class is phylo, and
  the other elements are:
  \itemize{
    \item{'edge'}{edge ids.}
    \item{'Nnode'}{number of internal nodes.}
    \item{'tip.lab'}{number of tips (leaves).}
    \item{'edge.length'}{length of edges.}
  }

  If \code{keep.bifurcation} is TRUE, then internal branches are set to be 0
  and leaves branches are set to the original tree height. Otherwise,
  the internal branches will be dropped from rooted.tree.
}
\value{
  Return a rooted tree in \code{Class phylo} with a star shape.
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
  \code{\link[ape]{plot.phylo}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
ret.ms <- ms(5, 1, opts = paste("-T", sep = " "))
tree.ms <- read.tree(text = ret.ms[3])
str(tree.ms)
(tree.star <- as.star.tree(tree.ms))

# Plot results
par(mfrow = c(1, 2))
plot(tree.ms, type = "u", main = "original tree")
plot(tree.star, type = "u", main = "as star tree")
}
}
\keyword{programming}
