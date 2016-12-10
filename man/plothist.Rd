\name{plothist}
\alias{plothist}
\title{Plot Histogram to Compare Number of Mutations.}
\description{
  Plot histogram to compare number of mutations.
}
\usage{
plothist(X, X.class = NULL, Mu = NULL, fill.color = .Color,
   draw.all = TRUE, main = "Mutation counts",
   xlab = "Difference", ylab = "Counts", append = FALSE)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{X.class}{class ids indicated for all sequences.}
  \item{Mu}{a central sequence with length \eqn{L}{L}.}
  \item{fill.color}{color to fill the histogram.}
  \item{draw.all}{draw a histogram use all sequences.}
  \item{main}{main label, default = "Mutation counts".}
  \item{xlab}{x-axis label, default = "Difference". }
  \item{ylab}{y-axis label, default = "Counts". }
  \item{append}{overwrite histograms.}
}
\details{
  If \code{X.class} is set, the histograms will be drawn by classes and
  all sequences will be compared to the central sequence \code{Mu}.
  Otherwise, all sequences will be used to count mutations.
  \code{draw.all} is not effect if \code{X.class} is not set.

  If \code{Mu} is set, it will be used to compare to all other sequences
  to count mutations. Otherwise, the first sequence of \code{X} will be
  used, and the first sequence in the first class will be used if
  \code{X.class} is set. If \code{Mu} is a matrix, the first row will be
  used as the central sequence.
}
\value{
  Histograms will be drawn to show the number of mutations away from
  the central sequence.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{seqgen}},
  \code{\link{plotdots}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

X <- seq.data.toy$org
plothist(X)

# With class ids
X.class <- as.numeric(gsub(".*-(.*)", "\\\\1", seq.data.toy$seqname))
plothist(X, X.class)
}
}
\keyword{plot}
