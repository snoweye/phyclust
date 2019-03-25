\name{.substitution.model}
\alias{.substitution.model}
\docType{data}
\title{Substitution Models for Mutation Processes}
\description{
  An object stores substitution models for mutation processes
  for Continuous Time Markov Chain theory.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.substitution.model}
\format{
  A data frame contains two character vectors, \code{mode} and \code{code.type}.
}
\details{
  This data frame indicates substitution models implemented in C.
  \itemize{
    \item 'model': names of substitution models for mutations.
    \item 'code.type': code types of substitution models, either NUCLEOTIDE or SNP.
  }

  The major models are:
  \tabular{lll}{
    Model \tab Author and Publication            \tab Parameter \cr
    JC69  \tab Jukes and Cantor 1969.            \tab \eqn{t}{t} \cr
    K80   \tab Kimura 1980.                      \tab \eqn{\kappa, t}{kappa, t} \cr
    F81   \tab Felsenstein 1981.                 \tab \eqn{\pi, t}{pi, t} \cr 
    HKY85 \tab Hasegawa, Kishino, and Yano 1985. \tab \eqn{\pi, \kappa, t}{kappa, pi, t} \cr
  }
  Other models starting with \code{E_} use empirical frequencies for
  equilibrium probabilities.
}

\references{
  Phylogenetic Clustering Website:
  \url{https://snoweye.github.io/phyclust/}

  Jukes, T. H. and Cantor, C. R. (1969)
  \dQuote{Evolution of Protein Molecules},
  \emph{Mammalian Protein Metabolism}, \bold{3}, 21-132

  Kimura, M. (1980)
  \dQuote{A Simple Method for Estimating Evolutionary Rates of Base
          Substitutions through Comparative Studies of Nucleotide Sequences},
  \emph{Journal of Molecular Evolution}, \bold{16}, 111-120

  Felsenstein, J. (1981)
  \dQuote{Evolutionary Trees from DNA Sequences: A Maximum Likelihood
          Approach},
  \emph{Journal of Molecular Evolution}, \bold{17}, 368-376

  Hasegawa, M. and Kishino, H. and Yano, T. (1985)
  \dQuote{Dating of the Human-Ape Splitting by a Molecular Clock of
          Mitochondrial DNA},
  \emph{Journal of Molecular Evolution}, \bold{22:2}, 160-174
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.code.type}},
  \code{\link{.identifier}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

.substitution.model
}
}
\keyword{data}
