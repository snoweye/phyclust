\name{plotgaps}
\alias{plotgaps}
\title{Gaps Plots of Sequences for Visual Comparisons}
\description{
  This function provides gaps plots of data set to identify regions where
  gaps enriched. The plot show the proportions of context by sites and the
  diverse may be caused by mutations, sequencing errors, or alignment
  errors.
}
\usage{
plotgaps(X, code.type = .code.type[1], main = "Gaps Plot",
    xlab = "Sites", ylab = "Proportion", ...)
}
\arguments{
  \item{X}{numerical data matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
  \item{main}{main label, default = "Gaps Plot".}
  \item{xlab}{x-axis label, default = "Sites". }
  \item{ylab}{y-axis label, default = "Proportion". }
  \item{\dots}{other options passed to \code{plot}.}
}
\details{
  Proportions of gaps will be drawn.
}
\value{
  A gaps plot will be drawn.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{plotdots}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

# For nucleotide
set.seed(1234)
X <- seq.data.toy$org
X[sample(c(T, F), length(X), replace = TRUE, prob = c(0.05, 0.95))] <-
  .nucleotide$nid[.nucleotide$code == "-"]
plotgaps(X)
}
}
\keyword{plot}
