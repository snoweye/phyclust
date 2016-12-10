\name{find.best}
\alias{find.best}
\title{Find the Best Solution of phyclust}
\description{
  Based on input initialization procedures and methods,
  this function tries to find the best solution in terms of the highest
  log-likelihood value.
}
\usage{
find.best(X, K, EMC = .EMC, manual.id = NULL, byrow = TRUE,
    init.procedure = .init.procedure, init.method = .init.method,
    file.tmp = NULL, visible = FALSE, save.all = FALSE)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{K}{number of clusters.}
  \item{EMC}{EM control.}
  \item{manual.id}{manually input class ids.}
  \item{byrow}{advanced option for \code{X}, default = TRUE.}
  \item{init.procedure}{customized initialization procedures.}
  \item{init.method}{customized initialization methods.}
  \item{file.tmp}{a file for saving temporary results.}
  \item{visible}{TRUE for reporting iterations.}
  \item{save.all}{TRUE for saving all results.}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.

  \strong{Note}: gaps \code{-} are not supported yet, drop them from data.

  \code{EMC} contains all options used for EM algorithms.

  \code{manual.id} manually input class ids as an initialization only for
  the initialization method, 'manualMu'.

  \code{byrow} used in bootstraps to avoid transposing matrix 'X'. If
  FALSE, then the 'X' should be have the dimension \eqn{L\times K}{L*K}.

  \code{init.procedure} and \code{init.method} are methods for searching
  the best result. This function will try all combinations of these two
  options.

  \code{file.tmp} is used to save temporary results due to long computing.
  If \code{NULL}, there will no saving in each combinations.
}
\value{
  An list with class \code{phyclust} will be returned containing
  several elements, see \code{phyclust} for detail. 
}
\section{ToDo(s)}{
  \itemize{
    \item implement codes for gaps \code{-}.
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
  \code{\link{.EMControl}}.
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
EMC.1 <- .EMControl(exhaust.iter = 1, short.iter = 5, EM.iter = 5)
(ret.1 <- find.best(seq.data.toy$org, 2, EMC = EMC.1))
}
}
\keyword{programming}
