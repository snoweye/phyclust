\name{phyclust.logL}
\alias{phyclust.logL}
\title{Log-Likelihood of phyclust}
\description{
  This computes a log-likelihood value of \code{phyclust}.
}
\usage{
phyclust.logL(X, ret.phyclust = NULL, K = NULL, Eta = NULL,
    Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    label = NULL)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{ret.phyclust}{an object with the class \code{phyclust}.}
  \item{K}{number of clusters.}
  \item{Eta}{proportion of subpopulations, \eqn{\eta_k}{eta_k}, length = \code{K}, sum to 1.}
  \item{Mu}{centers of subpopulations, dim = \eqn{K\times L}{K*L}, each row is a center.}
  \item{pi}{equilibrium probabilities, each row sums to 1.}
  \item{kappa}{transition and transversion bias.}
  \item{Tt}{total evolution time, \eqn{t}{t}.}
  \item{substitution.model}{substitution model.}
  \item{identifier}{identifier.}
  \item{code.type}{code type.}
  \item{label}{label of sequences for semi-supervised clustering.}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.

  Either input \code{ret.phyclust} or all other arguments for this function.
  \code{ret.phyclust} can be obtain either from an EM iteration of
  \code{phyclust} or from a M step of \code{phyclust.m.step}.

  If \code{label} is inputted, the label information will be used to
  calculate log likelihood (complete-data), even the \code{ret.phyclust}
  is the result of unsupervised clustering.
}
\value{
  This function returns a log-likelihood value of \code{phyclust}.
}
\references{
  Phylogenetic Clustering Website:
  \url{https://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{phyclust}},
  \code{\link{phyclust.em.step}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)
X <- seq.data.toy$org

ret.1 <- phyclust(X, 2, EMC = EMC.1)
phyclust.logL(X, ret.phyclust = ret.1)

# For semi-supervised clustering.
semi.label <- rep(0, nrow(X))
semi.label[1:3] <- 1
phyclust.logL(X, ret.phyclust = ret.1, label = semi.label)
}
}
\keyword{phyclust stepwised}
