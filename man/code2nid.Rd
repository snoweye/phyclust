\name{code2nid}
\alias{code2nid}
\alias{nid2code}
\title{Transfer Codes (A, G, C, T, -) and nids (0, 1, 2, 3, 4)}
\description{
  Transfer nucleotide codes (A, G, C, T, -) and nucleotide ids (0, 1, 2, 3, 4).
}
\usage{
### S3 methods for a list, vector or matrix (default).
code2nid(codeseq)
nid2code(nidseq, lower.case = TRUE)
}
\arguments{
  \item{codeseq}{a character vector contains nucleotide codes, A, G, C, T, or -.}
  \item{nidseq}{a numerical vector contains nucleotide ids, 0, 1, 2, 3, or 4.}
  \item{lower.case}{transfer in lower cases.}
}
\details{
  These functions are based on the internal object \code{.nucleotide} to
  transfer codes and nids.
}
\value{
  \code{code2nid} returns a numerical vector containing nucleotide ids, and
  \code{nid2code} returns a character vector containing nucleotide codes. 
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.nucleotide}},
  \code{\link{snp2sid}},
  \code{\link{sid2snp}},
  \code{\link{code2snp}},
  \code{\link{snp2code}}.
}
\examples{
\dontrun{
a <- c("A", "C", "G", "-", "T")
code2nid(a)
nid2code(code2nid(a)) 
}
}
\keyword{programming}
