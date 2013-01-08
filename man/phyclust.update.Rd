\name{phyclust.update}
\alias{phyclust.update}
\title{Update phyclust Results}
\description{
  This function will run the EM algorithm on initial parameters
  specified by users or from other initial procedures. All parameters
  (Eta, Mu, Q, ...) in this function will be updated.
}
\usage{
phyclust.update(X, EMC = .EMC, ret.phyclust = NULL, K = NULL,
    Eta = NULL, Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    label = NULL, byrow = TRUE)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{EMC}{EM control.}
  \item{ret.phyclust}{an object with the class \code{phyclust}.}
  \item{K}{number of clusters.}
  \item{Eta}{proportion of subpopulations, \eqn{\eta_k}{eta_k}, length = \code{K}, sum to 1.}
  \item{Mu}{centers of subpopulations, dim = \eqn{K\times L}{K*L}, each row is a center.}
  \item{pi}{equilibrium probabilities, each row sums to 1.}
  \item{kappa}{transition and transversion bias.}
  \item{Tt}{total evolution time, \eqn{t}{t}.}
  \item{label}{label of sequences for semi-supervised clustering.}
  \item{byrow}{advanced option for \code{X}, default = TRUE.}
}
\details{
  This function is equivalent to run \code{exhaustEM} on one specified
  initial parameters, and no initial procedure is involved. While this
  function is a little bit different to run \code{phyclust} with
  \code{manual.id} where \code{Mu} will be reestimated as the new initials.
  Simply speaking, this function only runs the EM algorithm given the
  initial parameters.

  All the input arguments are the same as the inputs of the functions
  \code{phyclust} and \code{phyclust.em.step}.
}
\value{
  This function returns an object with class \code{phyclust}.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{phyclust}},
  \code{\link{find.best}},
  \code{\link{phyclust.se}},
  \code{\link{phyclust.se.update}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)
X <- seq.data.toy$org

(ret.1 <- phyclust(X, 2, EMC = EMC.1))
(ret.2 <- phyclust.update(X, ret.phyclust = ret.1))
}
}
\keyword{phyclust}

