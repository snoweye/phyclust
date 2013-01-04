\name{gen.equal.star.anc.dec}
\alias{gen.equal.star.anc.dec}
\title{Generate Comprehensive Trees.}
\description{
  Generate comprehensive trees for simulation studies.
}
\usage{
gen.equal.star.anc.dec(K, N.K, rate.f = 0.5)
}
\arguments{
  \item{K}{number of clusters, \eqn{K}{K}.}
  \item{N.K}{number of sequences for each cluster, a vector with length \code{K}.}
  \item{rate.f}{\eqn{r_f}{r_f}, growth rate ratio of ancestral and descendent trees.}
}
\details{
  These functions generates an ancestral tree in K tips and
  generates descendent trees according to \code{N.K} tips.
  All trees, ancestral and descendent, are coerced to star shapes
  and scaled their heights to fit the ratio \code{rate.f}, and
  the final tree has total height 1.
  The returns are stored in a list, and the final tree is stored
  with a name \code{equal.star}.
}
\value{
  A list contains all information of generation and results including:
  \item{'K'}{number of clusters.}
  \item{'N.K}{number of sequences for each cluster.}
  \item{'rate.f'}{\eqn{r_f}{r_f}, growth rate ratio of ancestral and descendent trees.}
  \item{'anc'}{an ancestral tree.}
  \item{'dec'}{all descendent trees.}
  \item{'equalstar'}{a tree that descendants are equal star trees.}
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{gen.unit.K}}.
}
\examples{
\dontrun{
set.seed(1234)
tree.K <- gen.equal.star.anc.dec(6, rep(3:5, 2),
                                 rate.f = 0.7)
X.class <- as.numeric(gsub("d(.).(.)", "\\\\1",
                           tree.K$equal.star$tip.label))

# Plot results
plotnj(tree.K$equal.star, X.class, type = "p",
       edge.width.class = 2, main = "equal.star")
axis(1)
}
}
\keyword{simulation}
