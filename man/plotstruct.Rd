\name{plotstruct}
\alias{plotstruct}
\title{Struct Plots of Observations Based on Posterior Probabilities}
\description{
  This function provides structure plots of data set given based on posterior
  probabilities, the \code{Z.normalized} matrix.
}
\usage{
plotstruct(Z, X.class = NULL, sort.inside.class = TRUE,
    direction = "h", main = "Structure Plot", xlab = "Observations",
    ylab = "Posterior Probabilities", ...)
}
\arguments{
  \item{Z}{a Z matrix as Z.normalized in \code{Class phyclust}.}
  \item{X.class}{class ids indicated for all observations}
  \item{sort.inside.class}{sort observations inside class by max posteriors.}
  \item{direction}{either "h" for horizontal or "v" for vertical.}
  \item{main}{main label, default = "Structure Plot".}
  \item{xlab}{x-axis label, default = "Observations". }
  \item{ylab}{y-axis label, default = "Posterior Probabilities". }
  \item{\dots}{other options passed to \code{plot}.}
}
\details{
  The posterior probabilities in \code{ret.phyclust$Z.normalized} will
  be drawn in colors.

  If \code{X.class} is submitted, the plot will draw in the order of
  class ids and the \code{sort.inside.class} will be skipped.
}
\value{
  A structure plot will be drawn.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{phyclust}},
  \code{\link{find.best}},
  \code{\link{plotdots}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

set.seed(1234)
ret.1 <- phyclust(seq.data.toy$org, 3)
plotstruct(ret.1$Z.normalized)
windows()
plotstruct(ret.1$Z.normalized, sort.inside.class = FALSE)

# With class ids
X.class <- as.numeric(gsub(".*-(.*)", "\\\\1", seq.data.toy$seqname))
windows()
plotstruct(ret.1$Z.normalized, X.class = X.class)
}
}
\keyword{plot}
