\name{bootstrap.star.trees.seq}
\alias{bootstrap.star.trees.seq}
\title{Bootstrap Sequences from a Fitted Model.}
\description{
  This function bootstraps sequences from a model fitted by \code{phyclust}
  by combining two functions \code{bootstrap.star.trees} and
  \code{bootstrap.seq}.
  The fitted model can be varied in \code{.identifier}.
}
\usage{
bootstrap.star.trees.seq(ret.phyclust, min.n.class = 1)
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
  Return a list containing two elements, and both are corresponding to the
  model of \code{ret.phyclust}, including:

  \item{trees}{a list, \eqn{K}{K} star trees according to \code{ret.phyclust}}
  \item{seq}{a list, sequences in \eqn{K}{K} clusters}
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
  \code{\link{bootstrap.star.trees}},
  \code{\link{bootstrap.seq}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)

ret.1 <- phyclust(seq.data.toy$org, 2, EMC = EMC.1)
ret.all <- bootstrap.star.trees.seq(ret.1)
}
}
\keyword{bootstrap}
