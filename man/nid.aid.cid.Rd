\name{nid.aid.cid}
\alias{nid2aid}
\alias{nid2cid}
\alias{cid2aid}
\alias{aid2acode}
\alias{acode2aid}
\title{Transfer nids (0, 1, ..., 4) , aids (0, 1, ..., 21) and cids (0, 1, ..., 64)}
\description{
  Transfer nids (0, 1, ..., 4) , aids (0, 1, ..., 21) and cids (0, 1, ..., 64).
}
\usage{
### S3 methods for a list, vector or matrix (default).
nid2aid(nidseq, start = 1, end = NULL, drop.gap = FALSE, byrow = TRUE)
nid2cid(nidseq, start = 1, end = NULL, drop.gap = FALSE, byrow = TRUE)
cid2aid(cidseq)
aid2acode(aidseq, lower.case = FALSE)
acode2aid(acodeseq)
}
\arguments{
  \item{nidseq}{a numerical vector contains nucleotide ids, 0, 1, 2, 3, or 4.}
  \item{cidseq}{a numerical vector contains codon ids, 0, 1, ..., or 64.}
  \item{aidseq}{a numerical vector contains amino acid ids, 0, 1, ..., or 21.}
  \item{acodeseq}{a character vector contains amino acid codes.}
  \item{start}{the start site to translate.}
  \item{end}{the end site to translate.}
  \item{drop.gap}{ignore gaps if TRUE.}
  \item{byrow}{advanced option, default = TRUE.}
  \item{lower.case}{transfer in lower cases.}
}
\details{
  These functions are based on the internal object \code{.nucleotide},
  \code{.codon}, \code{.amino.acid}, and \code{.genetic.code} to
  transfer sequences.
}
\value{
  \code{nid2aid} and \code{cid2aid} returns a numerical vector containing
  amino acid ids, and \code{nid2cid} returns a numerical vector containing
  codon ids, \code{aid2acode} returns a character vector containing
  amino acid codes, and \code{acode2aid} returns a numerical vector containing
  amino acid ids.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.nucleotide}},
  \code{\link{.amino.acid}},
  \code{\link{.codon}},
  \code{\link{.genetic.code}},
  \code{\link{code2nid}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

### Test S3 methods by a vector
a.vector <- c("A", "C", "G", "-", "T", "A")
code2nid(a.vector)
nid2cid(code2nid(a.vector))
cid2aid(nid2cid(code2nid(a.vector)))
nid2aid(code2nid(a.vector))
aid2acode(nid2aid(code2nid(a.vector)))
acode2aid(aid2acode(nid2aid(code2nid(a.vector))))

### Test S3 methods by a matrix
a.matrix <- rbind(a.vector, a.vector, a.vector)
code2nid(a.matrix)
nid2cid(code2nid(a.matrix))
cid2aid(nid2cid(code2nid(a.matrix)))
nid2aid(code2nid(a.matrix))
aid2acode(nid2aid(code2nid(a.matrix)))
acode2aid(aid2acode(nid2aid(code2nid(a.matrix))))

### Test S3 methods by a list
a.list <- list(a.vector, a.vector)
code2nid(a.list)
nid2cid(code2nid(a.list))
cid2aid(nid2cid(code2nid(a.list)))
nid2aid(code2nid(a.list))
aid2acode(nid2aid(code2nid(a.list)))
acode2aid(aid2acode(nid2aid(code2nid(a.list))))
}
}
\keyword{programming}
