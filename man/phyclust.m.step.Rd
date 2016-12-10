\name{phyclust.m.step}
\alias{phyclust.m.step}
\title{One M-Step of phyclust}
\description{
  This is a single M-step of \code{phyclust}, usually following or followed
  by the other E-step.
}
\usage{
phyclust.m.step(X, ret.phyclust = NULL, K = NULL,
    pi = NULL, kappa = NULL, Tt = NULL, Z.normalized = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    label = NULL)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{ret.phyclust}{an object with the class \code{phyclust}.}
  \item{K}{number of clusters.}
  \item{pi}{equilibrium probabilities, each row sums to 1.}
  \item{kappa}{transition and transversion bias.}
  \item{Tt}{total evolution time, \eqn{t}{t}.}
  \item{Z.normalized}{posterior probabilities obtained from an E-step.}
  \item{substitution.model}{substitution model.}
  \item{identifier}{identifier.}
  \item{code.type}{code type.}
  \item{label}{label of sequences for semi-supervised clustering.}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.

  Either input \code{ret.phyclust} or all other arguments for this function.
  \code{ret.phyclust} can be obtained either from an EM iteration of
  \code{phyclust} or from a E step of \code{phyclust.e.step}.

  If \code{label} is inputted, the label information will be used
  the M-step and \code{Z.normalized} will be replaced, even the
  \code{ret.phyclust} is the result of unsupervised clustering.
}
\value{
  This function returns an object with class \code{phyclust}.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{phyclust}},
  \code{\link{phyclust.em.step}},
  \code{\link{phyclust.e.step}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
EMC.1 <- .EMC
EMC.1$short.iter <- 1
EMC.1$EM.iter <- 1

# Test with phyclust.
X <- seq.data.toy$org
ret.1 <- phyclust(X, 2, EMC = EMC.1)

# Test with an em step.
ret.em <- phyclust.em.step(X, ret.1)

# Test with an E- and M-step.
ret.1$Z.normalized <- phyclust.e.step(X, ret.phyclust = ret.1)
ret.m <- phyclust.m.step(X, ret.phyclust = ret.1)

# Test with 2 em steps.
set.seed(1234)
EMC.2 <- EMC.1
EMC.2$EM.iter <- 2
ret.2 <- phyclust(X, 2, EMC = EMC.2)

# Check logL.
phyclust.logL(X, ret.1)
phyclust.logL(X, ret.em)
phyclust.logL(X, ret.m)
phyclust.logL(X, ret.2)

# For semi-supervised.
semi.label <- rep(0, nrow(X))
semi.label[1:3] <- 1
ret.m.1 <- phyclust.m.step(X, ret.phyclust = ret.1, label = semi.label)
}
}
\keyword{phyclust stepwised}
