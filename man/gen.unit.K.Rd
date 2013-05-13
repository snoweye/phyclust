\name{gen.unit.K}
\alias{gen.unit.K}
\title{Generate Comprehensive Trees.}
\description{
  Generate comprehensive trees for simulation studies.
}
\usage{
gen.unit.K(K, N.K, rate.anc = 10, rate.dec = 10)
}
\arguments{
  \item{K}{number of clusters, \eqn{K}{K}.}
  \item{N.K}{number of sequences for each cluster, a vector with length \code{K}.}
  \item{rate.anc}{\eqn{r_a}{r_a}, growth rate of ancestral tree.}
  \item{rate.dec}{\eqn{r_d}{r_d}, growth rate of descendent tree.}
}
\details{
  These functions generates an ancestral tree in K tips and
  generates descendent trees according to \code{N.K} tips,
  and returns several types of trees, \code{org}, \code{equal},
  \code{max}, and \code{star}, as the following:
  \itemize{
    \item 'org': original tree, adjacent the descendent trees to the ancestral
                 tree.
    \item 'equal': descendent trees are scaled to the average height and
                   attached to the ancestral tree, then scale the total
                   height to be 1.
    \item 'max': descendent trees are attached to the ancestral tree,
                 then scale the maximum height to be 1.
    \item 'star': descendent trees are applied \code{as.star.tree} and
                  attached to the ancestral tree,
                  then scale the maximum height to be 1.
  }
}
\value{
  A list contains all information of generation and results including:
  \item{'K'}{number of clusters.}
  \item{'N.K}{number of sequences for each cluster.}
  \item{'rate.anc'}{\eqn{r_a}{r_a}, growth rate of ancestral tree.}
  \item{'rate.dec'}{\eqn{r_d}{r_d}, growth rate of descendent tree.}
  \item{'height.anc'}{height of ancestral tree.}
  \item{'height.dec'}{height of all descendent trees.}
  \item{'anc'}{an ancestral tree.}
  \item{'dec'}{all descendent trees.}
  \item{'org'}{an original tree.}
  \item{'equal'}{a three that descendants are all equal height.}
  \item{'max'}{a tree that descendants are scaled by the maximum height.}
  \item{'star'}{a tree that descendants are star trees.}
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{gen.equal.star.anc.dec}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

# For gen.unit.K()
set.seed(1234)
tree.K <- gen.unit.K(6, rep(3:5, 2),
                     rate.anc = 0.7, rate.dec = 1.1)
X.class <- as.numeric(gsub("d(.)(.*)", "\\\\1",
                           tree.K$org$tip.label))

# Plot results
par(mfrow = c(2, 2))
plotnj(tree.K$org, X.class, type = "p",
       edge.width.class = 2, main = "org")
axis(1)
plotnj(tree.K$equal, X.class, type = "p",
       edge.width.class = 2, main = "equal")
axis(1)
plotnj(tree.K$max, X.class, type = "p",
       edge.width.class = 2, main = "max")
axis(1)
plotnj(tree.K$star, X.class, type = "p",
       edge.width.class = 2, main = "star")
axis(1)
}
}
\keyword{simulation}
