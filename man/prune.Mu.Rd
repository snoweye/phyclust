\name{prune.Mu}
\alias{prune.Mu}
\title{Prune the Center Sequences Mu}
\description{
  This function prune the center sequences Mu where the sites will be reset
  as GAPs if all members within the same cluster are all GAPs.
}
\usage{
prune.Mu(X, X.class, Mu, code.type = .code.type[1])
}
\arguments{
  \item{X}{numerical data matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{X.class}{class ids indicated for all sequences.}
  \item{Mu}{a center sequence with length \eqn{L}{L}.}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
}
\details{
  For each cluster indicated by \code{X.class}, this function will prune
  \code{Mu} and reset the sites as GAPs if all members within cluster
  are all GAPs. \code{Mu} are usually the returning values of
  \code{phyclust()}.
}
\value{
  A pruned \code{Mu} will be returned.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

X <- seq.data.toy$org
X[, 5] <- .nucleotide$nid[.nucleotide$code == "-"]
ret <- phyclust(X, 2)
Mu.GAPs <- prune.Mu(X, ret$class.id, ret$Mu)

ret$Mu[, 5]
Mu.GAPs[, 5]  # Replace by GAPs.
}
}
\keyword{programming}
