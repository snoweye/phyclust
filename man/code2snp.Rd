\name{code2snp}
\alias{code2snp}
\alias{snp2code}
\alias{nid2sid}
\alias{sid2nid}
\title{Transfer Nucleotide Codes / nids and SNPs / sids}
\description{
  Transfer nucleotide codes (A, G, C, T, -) and SNPs (1, 2, -).
  Transfer nucleotide ids (0, 1, 2, 3, 4) and SNP ids (0, 1, 2).
}
\usage{
### S3 methods for a list, vector or matrix (default).
code2snp(codeseq)
snp2code(snpseq, half = TRUE)

nid2sid(nidseq)
sid2nid(sidseq, half = TRUE)
}
\arguments{
  \item{codeseq}{a character vector contains nucleotide codes, A, G, C, T, or -.}
  \item{snpseq}{a character vector contains SNPs, 1, 2, or -.}
  \item{half}{nucleotide codes will be half assigned, see the \sQuote{Details} for more information.}

  \item{nidseq}{a numerical vector contains nucleotide ids, 0, 1, 2, 3, or 4.}
  \item{sidseq}{a numerical vector contains SNP ids, 0, 1, or 2.}
}
\details{
  These functions are based on the internal object \code{.nucleotide} and 
  \code{.snp} to transfer nucleotide codes and SNPs.
  For \code{code2snp}, A, G are transfered to 1, and C, T are transfered to 2.
  For \code{snp2code}, 1 is transfered half to A and G, and 2 is transfered
  half to C and T if \code{half = TRUE}. Otherwise, 1 is all transfered to A,
  and 2 is all transfered to C.
}
\value{
  \code{code2nid} returns a character vector containing nucleotide ids, and
  \code{nid2code} returns a character vector containing nucleotide codes. 

  \code{nid2sid} returns a numerical vector containing SNP ids, and
  \code{sid2nid} returns a numerical vector containing nucleotide ids.
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
  \code{\link{.snp}},
  \code{\link{code2nid}},
  \code{\link{nid2code}},
  \code{\link{snp2sid}},
  \code{\link{sid2snp}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

# For codes
a.vector <- c("A", "C", "G", "-", "T")
code2snp(a.vector)
snp2code(code2snp(a.vector))
snp2code(code2snp(a.vector), half = FALSE)

# For ids
a.sid.vector <- c(0, 2, 1, 4, 3)
nid2sid(a.sid.vector)
sid2nid(nid2sid(a.sid.vector))
sid2nid(nid2sid(a.sid.vector), half = FALSE)

# Test list
a.list <- list(a, a)
code2snp(a.list)
snp2code(code2snp(a.list))
snp2code(code2snp(a.list), half = FALSE)

# Test matrix
a.matrix <- rbind(a, a)
code2snp(a.matrix)
snp2code(code2snp(a.matrix))
snp2code(code2snp(a.matrix), half = FALSE)
}
}
\keyword{programming}
