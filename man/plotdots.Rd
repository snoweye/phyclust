\name{plotdots}
\alias{plotdots}
\title{Dots Plots of Sequences for Visual Comparisons}
\description{
  This function provides dots plots of data set given an idea how diverse
  the sequences are by drawing dots with different colors for all mutations.
}
\usage{
plotdots(X, X.class = NULL, Mu = NULL, code.type = .code.type[1],
    diff.only = TRUE, fill = FALSE, label = TRUE, with.gap = FALSE,
    xlim = NULL, ylim = NULL, main = "Dots Plot", xlab = "Sites",
    ylab = "Sequences", missing.col = "gray95", ...)
}
\arguments{
  \item{X}{numerical data matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{X.class}{class ids indicated for all sequences.}
  \item{Mu}{a center sequence with length \eqn{L}{L}.}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
  \item{diff.only}{draw the segregating sites only, default = TRUE.}
  \item{fill}{fill in all dots, default = FALSE.}
  \item{label}{indicate segregating sites, default = TRUE.}
  \item{with.gap}{pass to \code{find.consensus} if \code{Mu} is \code{NULL}, default = FALSE}
  \item{xlim}{limit of x-axis.}
  \item{ylim}{limit of y-axis.}
  \item{main}{main label, default = "Dots Plot".}
  \item{xlab}{x-axis label, default = "Sites". }
  \item{ylab}{y-axis label, default = "Sequences". }
  \item{missing.col}{color for the missing allele, default = NA. }
  \item{\dots}{other options passed to \code{plot}.}
}
\details{
  The first rows in \code{Mu} will be drawn entirely on dots plots
  in colors which are "green3", "blue2", "\#CC00CC", "red2", "gray",
  and "white", according the ids + 1. If \code{fill} is FALSE, other sequences
  will be drawn by the mutation sites comparing to the first sequences.
  Otherwise, they be drawn entirely.

  If \code{X.class} is set, the sequences will be drawn in cluster order.

  If \code{Mu} is \code{NULL}, the consensus sequence of \code{X} will be drawn.

  If \code{label} is TRUE, the bottom row will be drawn in color "orange"
  to indicate segregating sites.

  \code{with.gap} is only used when \code{Mu} is \code{NULL}.
}
\value{
  A dots plot will be drawn.
}
\references{
  Phylogenetic Clustering Website:
  \url{https://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{seqgen}},
  \code{\link{plothist}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

# For nucleotide
X <- seq.data.toy$org
par(mfrow = c(2, 2))
plotdots(X)
plotdots(X, diff.only = FALSE)
plotdots(X, diff.only = FALSE, label = FALSE)
plotdots(X, fill = TRUE, diff.only = FALSE, label = FALSE)

# With class ids
X.class <- as.numeric(gsub(".*-(.*)", "\\\\1", seq.data.toy$seqname))
plotdots(X, X.class)

# For SNP
X.SNP <- nid2sid(X)
plotdots(X.SNP, X.class)
}
}
\keyword{plot}
