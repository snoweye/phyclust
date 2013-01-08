\name{phyclust.se}
\alias{phyclust.se}
\title{The Main Function of phyclust for Sequencing Error Models}
\description{
  The \code{phyclust.se} is an advanced function of \code{phyclust}.
  The \code{phyclust.se} implements finite mixture models for sequence data
  with sequencing errors. The same as \code{phyclust},
  the mutation processes are modeled by evolution processes based on
  Continuous Time Markov Chain theory.
}
\usage{
phyclust.se(X, K, EMC = .EMC, manual.id = NULL, byrow = TRUE)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{K}{number of clusters.}
  \item{EMC}{EM control.}
  \item{manual.id}{manually input class ids.}
  \item{byrow}{advanced option for \code{X}, default = TRUE.}
}
\details{
  \code{phyclust.se} considers mutations with sequencing error, but only
  for NUCLEOTIDE data and only the EM algorithm is
  implemented. Currently, \code{phyclust.se} implements the following steps:
  \tabular{ll}{
  1 \tab assume non-sequencing errors as the \code{phyclust} does. \cr
  2 \tab use the initialization as the \code{phyclust} does. \cr
  3 \tab run the \code{phyclust} to find the non-sequencing error MLE's. \cr
  4 \tab initial by the results of \code{phyclust}. \cr
  5 \tab update all parameters including probabilities of sequencing errors. \cr
  }

  See the help page of \code{\link{phyclust}} for the explanations of
  input arguments since both functions share the same arguments. Note that
  the underling model assumptions are very different of both functions.
}
\value{
  A list with class \code{phyclust} will be returned containing
  several elements as described in \code{\link{phyclust}}. But, the
  \code{phyclust.se} returns extra parameters for sequencing errors, and
  they are shown in the following:
  \item{'SE'}{a list returning parameters of sequencing error models, including:
              \tabular{ll}{
              'type'     \tab 'TRUE' for modeling sequencing errors.          \cr
              'model'    \tab 'CONVOLUTION', the only model implemented.      \cr
              'constant' \tab the constrained constant for sequencing errors. \cr
              'f.err'    \tab probability matrix, each row sums to 1.         \cr
              }
  }
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.EMC}},
  \code{\link{.EMControl}},
  \code{\link{phyclust.se}},
  \code{\link{phyclust.se.update}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

X <- seq.data.toy$org

set.seed(1234)
(ret.1 <- phyclust.se(X, 3))
}
}
\keyword{phyclust}
