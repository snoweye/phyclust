\name{file.read}
\alias{read.fasta}
\alias{read.fasta.format}
\alias{read.phylip}
\alias{read.phylip.format}
\title{Read Data from Files by Formats and Return a seq.data Object}
\description{
  Two major file formats are supported in \pkg{phyclust},
  \code{Format phylip} and \code{Format fasta}.
  \strong{These functions only read files in basic syntax}, and
  return an object in \code{Class seq.data}.
}
\usage{
read.fasta(filename, byrow = TRUE, code.type = .code.type[1], aligned = TRUE,
           sep = "")
read.fasta.format(filename, byrow = TRUE, aligned = TRUE, sep = "")

read.phylip(filename, byrow = TRUE, code.type = .code.type[1], sep = "")
read.phylip.format(filename, byrow = TRUE, sep = "")
}
\arguments{
  \item{filename}{a file name where data is read from.}
  \item{byrow}{advanced option, default = TRUE.}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
  \item{aligned}{indicate aligned data.}
  \item{sep}{use to split sites, "" (default) and "," for "CODON".}
}
\details{
  For unaligned sequences, \code{read.fasta} returns a list storing data.
  \code{read.phylip} is only for aligned data and returns a matrix.

  \code{read.fasta.format} and \code{read.phylip.format} will read in
  original coding without any transformation as \code{code.type = NULL}
  in \code{write.fasta} and \code{write.phylip}. Suppose these functions
  return an object \code{ret}, one can write other functions \code{ret2aa()}
  to post transform the coding and replace \code{ret$org} by the results of
  \code{ret2aa(ret$org.code)}.

  \code{byrow} indicates the data will be store by row or not. Usually,
  the default is TRUE. The FALSE is only for advance users
  with careful manipulations and for speeding up the bootstraps.

  \code{sep} can specify a character which is used to split sites in file.
  By default, "" denote no character between sites. Only "CODON" id requires
  a character to avoid ambiguity
}
\value{
  Return an object in \code{Class seq.data}.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{write.fasta}},
  \code{\link{write.phylip}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

# PHYLIP
data.path <- paste(.libPaths()[1], "/phyclust/data/crohn.phy", sep = "")
(my.snp <- read.phylip(data.path, code.type = "SNP"))

# FASTA
data.path <- paste(.libPaths()[1], "/phyclust/data/pony625.fas", sep = "")
(my.pony <- read.fasta(data.path))
}
}
\keyword{file I/O}
