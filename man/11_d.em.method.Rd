\name{.em.method}
\alias{.em.method}
\docType{data}
\title{EM Methods and Algorithms}
\description{
  The varied EM algorithms are implemented in C.
  The first element is the default value.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.em.method}
\format{
  A character vector contains implemented EM algorithms in C.
}
\details{
  \code{EM} (default) stands for the standard EM algorithm, \code{ECM}
  stands for Expectation/Conditional Maximization algorithm, and \code{AECM}
  stands for Alternating ECM algorithm.
  The performance is roughly about \code{AECM} > \code{EM} ~ \code{ECM} which
  are dependent on the separations of data set.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Dempster, A. and Laird, N. and Rubin, D. (1977)
  \dQuote{Maximum Likelihood Estimation from Incomplete Data
          via the EM Algorithm},
  \emph{Journal of the Royal Statistical Society Series B}, \bold{39:3}, 1-38.

  Meng, X.-L. and Rubin, D. (1993)
  \dQuote{Maximum likelihood estimation via the ECM algorithm: A general
          framework},
  \emph{Biometrika}, \bold{80:2}, 511-567.

  Meng, X.-L. and van Dyk, D. (1997)
  \dQuote{The EM Algorithm --- an Old Folk-song Sung to a Fast New Tune
          (with discussion)},
  \emph{Journal of the Royal Statistical Society Series B}, \bold{59}, 511-567.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.EMC}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.em.method
}
}
\keyword{data}
