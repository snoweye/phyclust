\name{snp2sid}
\alias{snp2sid}
\alias{sid2snp}
\title{Transfer SNP codes (1, 2, -) and sids (0, 1, 2)}
\description{
  Transfer SNP codes (1, 2, -) and SNP ids (0, 1, 2).
}
\usage{
snp2sid(snpseq)
sid2snp(sidseq)
}
\arguments{
  \item{snpseq}{a character vector contains SNP codes, 1, 2, or -.}
  \item{sidseq}{a numerical vector contains SNP ids, 0, 1, or 2.}
}
\details{
  This function is based on the internal object \code{.snp} to
  transfer SNP codes and SNP ids.
}
\value{
  \code{snp2sid} returns a numerical vector containing SNP ids, and
  \code{sid2snp} returns a character vector containing SNP codes. 
}
\references{
  Phylogenetic Clustering Website:
  \url{https://snoweye.github.io/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.snp}},
  \code{\link{code2nid}},
  \code{\link{nid2code}},
  \code{\link{code2snp}},
  \code{\link{snp2code}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

a <- c("1", "2", "1", "-", "2")
snp2sid(a)
sid2snp(snp2sid(a)) 
}
}
\keyword{programming}
