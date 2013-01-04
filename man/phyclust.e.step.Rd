\name{phyclust.e.step}
\alias{phyclust.e.step}
\title{One E-Step of phyclust}
\description{
  This is a single E-step of \code{phyclust}, usually following or followed
  by the other M-step.
}
\usage{
phyclust.e.step(X, ret.phyclust = NULL, K = NULL, Eta = NULL,
    Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    Z.state = TRUE, label = NULL)
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
  \item{Z.state}{see \sQuote{Details}.}
  \item{label}{label of sequences for semi-supervised clustering.}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.

  Either input \code{ret.phyclust} or all other arguments for this function
  except \code{Z.state}. \code{ret.phyclust} can be obtain either from an
  EM iteration of \code{phyclust} or from a M step of \code{phyclust.m.step}.

  \code{Z.state} indicates the return values of \eqn{Z_{nk}}{Z_nk}. If
  TRUE, the \code{Z.normalized} returned by this function will be
  posterior probabilities. Otherwise, it will be \code{logPt}, log of
  transition probabilities, \eqn{\log(\phi(\cdots))}{log(phi(...))}.

  If \code{label} is inputted, the label information will be used
  the E-step, even the \code{ret.phyclust} is the result of unsupervised
  clustering.
}
\value{
  This function returns a \eqn{Z_{nk}}{Z_nk} matrix with dimension =
  \eqn{N\times K}{N*K}. The values is dependent on \code{Z.state}, and
  they are either posterior probabilities if TRUE or transition
  probabilities otherwise.
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
  \code{\link{phyclust.em.step}},
  \code{\link{phyclust.m.step}}.
}
\examples{
\dontrun{
set.seed(1234)
EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)
X <- seq.data.toy$org

ret.1 <- phyclust(X, 2, EMC = EMC.1)
ret.2 <- phyclust.e.step(X, ret.phyclust = ret.1)
str(ret.2)

# For semi-supervised clustering.
semi.label <- rep(0, nrow(X))
semi.label[1:3] <- 1
ret.3 <- phyclust.e.step(X, ret.phyclust = ret.1, label = semi.label)
}
}
\keyword{phyclust stepwised}
