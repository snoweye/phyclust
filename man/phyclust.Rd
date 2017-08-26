\name{phyclust}
\alias{phyclust}
\alias{Class phyclust}
\title{The Main Function of phyclust}
\description{
  The main function of \pkg{phyclust} implements finite mixture
  models for sequence data that the mutation processes are modeled by
  evolution processes based on Continuous Time Markov Chain theory.
}
\usage{
phyclust(X, K, EMC = .EMC, manual.id = NULL, label = NULL, byrow = TRUE)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{K}{number of clusters.}
  \item{EMC}{EM control.}
  \item{manual.id}{manually input class ids.}
  \item{label}{label of sequences for semi-supervised clustering}
  \item{byrow}{advanced option for \code{X}, default = TRUE.}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.

  \code{EMC} contains all options used for EM algorithms.

  \code{manual.id} manually input class ids as an initialization only for
  the initialization method, 'manualMu'.

  \code{label} indicates the known clusters for labeled sequences which is a
  vector with length \code{N} and has values from \code{0} to \code{K}.
  \code{0} indicates clusters are unknown. \code{label = NULL} is for
  unsupervised clustering. Only un- and semi-supervised clustering are
  implemented.

  \code{byrow} used in bootstraps to avoid transposing matrix 'X'. If
  FALSE, then the 'X' should be have the dimension \eqn{L\times K}{L*K}.
}
\value{
  A list with class \code{phyclust} will be returned containing
  several elements as the following:

  \item{'N.X.org'}{number of sequences in the \code{X} matrix.}
  \item{'N.X.unique'}{number of unique sequences in the \code{X} matrix.}
  \item{'L'}{number of sites, length of sequences, number of column of the \code{X} matrix.}
  \item{'K'}{number of clusters.}
  \item{'Eta'}{proportion of subpopulations, \eqn{\eta_k}{eta_k}, length = \code{K}, sum to 1.}
  \item{'Z.normalized'}{posterior probabilities, \eqn{Z_{nk}}{Z_nk}, each row sums to 1.}
  \item{'Mu'}{centers of subpopulations, dim = \eqn{K\times L}{K*L}, each row is a center.}
  \item{'QA'}{Q matrix array, information for the evolution model, a list contains:
              \tabular{ll}{
              'pi'         \tab equilibrium probabilities, each row sums to 1. \cr
              'kappa'      \tab transition and transversion bias. \cr
              'Tt'         \tab total evolution time, \eqn{t}{t}. \cr
              'identifier' \tab identifier for QA. \cr
              }
  }
  \item{'logL'}{log likelihood values.}
  \item{'p'}{number of parameters.}
  \item{'bic'}{BIC, \eqn{-2\log L + p \log N}{-2logL + plogN}.}
  \item{'aic'}{AIC, \eqn{-2\log L + 2p}{-2logL + 2p}.}
  \item{'N.seq.site'}{number of segregating sites.}
  \item{'class.id'}{class id for each sequences based on the maximum posterior.}
  \item{'n.class'}{number of sequences in each cluster.}
  \item{'conv'}{convergence information, a list contains:
                \tabular{ll}{
                  'eps'         \tab relative error. \cr
                  'error'       \tab error if the likelihood decreased. \cr
                  'flag'        \tab convergence state. \cr
                  'iter'        \tab convergence iterations. \cr
                  'inner.iter'  \tab convergence of inner iterations other than EM. \cr
                  'cm.iter'     \tab convergence of CM iterations. \cr
                  'check.param' \tab parameter states. \cr
                }
  }
  \item{'init.procedure'}{initialization procedure.}
  \item{'init.method'}{initialization method.}
  \item{'substitution.model'}{substitution model.}
  \item{'edist.model'}{evolution distance model.}
  \item{'code.type'}{code type.}
  \item{'em.method'}{EM algorithm.}
  \item{'boundary.method'}{boundary method.}
  \item{'label.method'}{label method.}
}
\section{ToDo(s)}{
  \itemize{
    \item make a general class for \code{Q} and \code{QA}.
  }
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.EMC}},
  \code{\link{.EMControl}},
  \code{\link{find.best}},
  \code{\link{phyclust.se}}.
  \code{\link{phyclust.se.update}}.
}
\examples{
%\dontrun{
library(phyclust, quiet = TRUE)

X <- seq.data.toy$org

set.seed(1234)
(ret.1 <- phyclust(X, 3))

EMC.2 <- .EMC
EMC.2$substitution.model <- "HKY85"
# the same as EMC.2 <- .EMControl(substitution.model = "HKY85")

(ret.2 <- phyclust(X, 3, EMC = EMC.2))

# for semi-supervised clustering
semi.label <- rep(0, nrow(X))
semi.label[1:3] <- 1
(ret.3 <- phyclust(X, 3, EMC = EMC.2, label = semi.label))
%}
}
\keyword{phyclust}
