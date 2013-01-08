\name{phyclust.edist}
\alias{phyclust.edist}
\title{Evolution Distance of Sequences}
\description{
  This computes pair wised evolution distance of sequences.
}
\usage{
phyclust.edist(X, edist.model = .edist.model[1])
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{edist.model}{evolution distance model.}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.
}
\value{
  This function returns a object with class \code{dist}.
}
\section{ToDo(s)}{
  \itemize{
    \item incorporate \code{dist.dna} of \pkg{ape}.
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
  \code{\link{.edist.model}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

X <- rbind(c(0, 2, 1, 3, 0, 2, 2, 0, 3, 2, 2),
           c(0, 0, 1, 3, 2, 2, 1, 0, 3, 1, 2),
           c(0, 2, 1, 1, 0, 2, 1, 3, 0, 0, 1),
           c(2, 2, 1, 1, 0, 0, 2, 3, 0, 2, 1),
           c(2, 2, 1, 1, 0, 0, 2, 3, 1, 2, 0))
(ret <- phyclust.edist(X, edist.model = "D_HAMMING"))
str(ret)
as.matrix(ret)
plot(nj(ret), type = "u", no.margin = TRUE)
}
}
\keyword{programming}
