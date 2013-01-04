\name{plotnj}
\alias{plotnj}
\title{Plot an Unrooted Trees.}
\description{
  This is a enhanced version of \code{plot.phylo} in \pkg{ape} which
  can plot trees in \code{Class phylo} including neighbor-joining trees,
  unrooted trees, trees with star shapes, ... etc.
}
\usage{
plotnj(unrooted.tree, X.class = NULL, type = "u", main = NULL,
   show.tip.label = FALSE, show.node.label = FALSE,
   edge.width = 1, edge.width.class = edge.width, ...)
}
\arguments{
  \item{unrooted.tree}{an unrooted tree in \code{Class phylo}.}
  \item{X.class}{class ids indicated for all tips.}
  \item{type}{plot types, see \code{plot.phylo} in \pkg{ape} for details.}
  \item{main}{main label.}
  \item{show.tip.label}{show tip label if available.}
  \item{show.node.label}{show node label if available.}
  \item{edge.width}{edge width for all internal branches if \code{X.class} is set.}
  \item{edge.width.class}{edge width for tip branches if \code{X.class} is set.}
  \item{\dots}{other options passed to \code{plot.phylo}.}
}
\details{
  This function is built to plot unrooted trees, but it may also apply
  for other trees in \code{Class phylo}.

  \code{type} can be "u", "p", "c", "f", "r" as in \code{plot.phylo}.

  If \code{X.class} is set, then the tip branches will be drawn with
  colors by class ids, and the colors are controlled by \code{.color}.
  The width of branches is controlled by \code{edge.width} for all internal
  branches and by \code{edge.width.class} for tip branches.
}
\value{
  Return a tree plot.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link[ape]{plot.phylo}},
  \code{\link{.Color}}.
}
\examples{
\dontrun{
set.seed(1234)
ret.ms <- ms(nsam = 24, opts = "-T -G 0.5")
tree.anc <- read.tree(text = ret.ms[3])

is.rooted(tree.anc)
tree.new <- as.star.tree(tree.anc)
X.class <- rep(1:6, each = 4)

par(mfrow = c(2, 2))
plotnj(tree.anc, X.class, type = "u", edge.width.class = 2,
       main = "unrooted tree")
plotnj(tree.new, X.class, type = "u", edge.width.class = 2,
       main = "star tree")
plotnj(tree.anc, X.class, type = "c", edge.width.class = 2,
       main = "unrooted tree in cladogram")
plotnj(tree.new, X.class, type = "r", edge.width.class = 2,
       main = "star tree in radial")
}
}
\keyword{plot}
