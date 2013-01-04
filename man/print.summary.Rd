\name{print.object}
\alias{print.baseml}
\alias{print.ms}
\alias{print.phyclust}
\alias{print.Pt}
\alias{print.RRand}
\alias{print.seq.data}
\alias{print.seqgen}
\alias{summary.phyclust}
\title{Functions for Printing or Summarizing Objects According to Classes}
\description{
  Several classes are declared in \pkg{phyclust}, and these are functions
  to print and summary objects.
}
\usage{
\method{print}{baseml}(x, ...)
\method{print}{ms}(x, ...)
\method{print}{phyclust}(x, digits = max(4, getOption("digits") - 3), ...)
\method{print}{Pt}(x, ...)
\method{print}{RRand}(x, digits = max(4, getOption("digits") - 3), ...)
\method{print}{seq.data}(x, ...)
\method{print}{seqgen}(x, ...)
\method{summary}{phyclust}(object, ...)
}
\arguments{
  \item{x}{an object with the class attributes.}
  \item{digits}{for printing out numbers.}
  \item{object}{an object with the class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  These are useful functions for summarizing and debugging.

  For \code{ms}, \code{seqgen}, and \code{paml.baseml}, it will show the
  result as standalone versions on the STDOUT out with line by line.

  For other functions, they only show summaries of objects. Use
  \code{names} or \code{str} to explore the details.
}
\value{
  The results will cat or print on the STDOUT by default.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{ms}},
  \code{\link{paml.baseml}},
  \code{\link{phyclust}},
  \code{\link{phyclust.Pt}},
  \code{\link{RRand}},
  \code{\link{seqgen}}.
}
\examples{
\dontrun{
# Functions applied by directly type the names of objects.
}
}
\keyword{programming}
