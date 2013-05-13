\name{bootstrap.star.trees}
\alias{bootstrap.star.trees}
\title{Bootstrap a Star Tree from a Fitted Model.}
\description{
  This function bootstraps a star tree from a model fitted by \code{phyclust}.
  Each cluster corresponds to a star tree and a center sequence where
  sequences will evolve from. This function is called by
  \code{bootstrap.star.trees.seq} to generate sequences.
  The fitted model can be varied in \code{.identifier}.
}
\usage{
bootstrap.star.trees(ret.phyclust, min.n.class = 1)
}
\arguments{
  \item{ret.phyclust}{a phyclust object in \code{Class phyclust}.}
  \item{min.n.class}{minimum number of sequences for a cluster.}
}
\details{
  \code{ret.phyclust} is a phyclust object in \code{Class phyclust} which is usually
  fitted by \code{phyclust}, or returned by \code{phyclust.m.step}.

  \code{min.n.class} is a boundary condition to avoid degenerate clusters
  when some population proportions, \eqn{\eta_k}{eta_k}, are small in the
  fitted model.
}
\value{
  Return a list containing \eqn{K}{K} star trees according to \code{ret.phyclust}.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{phyclust}},
  \code{\link{bootstrap.seq}},
  \code{\link{bootstrap.star.trees.seq}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)

ret.1 <- phyclust(seq.data.toy$org, 2, EMC = EMC.1)
ret.trees <- bootstrap.star.trees(ret.1)
}
}
\keyword{bootstrap}
