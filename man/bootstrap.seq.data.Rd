\name{bootstrap.seq.data}
\alias{bootstrap.seq.data}
\title{Bootstrap a seq.data from a Fitted Model.}
\description{
  This function simplifies the bootstrap function
  \code{bootstrap.star.trees.seq()}, and only
  return a list object with class \code{seq.data}.
}
\usage{
bootstrap.seq.data(ret.phyclust, min.n.class = 1)
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
  Return an object in \code{Class seq.data} as
  the result from \code{read.*()}.
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
  \code{\link{Class seq.data}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)

ret.1 <- phyclust(seq.data.toy$org, 2, EMC = EMC.1)
(ret.all <- bootstrap.seq.data(ret.1))
}
}
\keyword{bootstrap}
