\name{standard.code}
\alias{.nucleotide}
\alias{.snp}
\alias{.codon}
\alias{.amino.acid}
\alias{.genetic.code}
\alias{.missing.code}
\docType{data}
\title{Standard Codes and ids for Nucleotides, SNPs, Codon, Amino Acid and Genetic Code}
\description{
  Standard codes and ids for nucleotides, SNPs, codon, amino acid and genetic
  code.
  All objects are used to transfer data.
  \strong{These are read-only objects and the elemental order is followed in C.}
}
\usage{
.nucleotide
.snp
.codon
.amino.acid
.genetic.code
}
\format{
  All objects are data frames containing ids and codes.
}
\details{
  \strong{Note}: All ids are coding started from \strong{0}.

  Nucleotides, A, G, C, T, and - have codes 0, 1, 2, 3, and 4 where
  "-" is for gaps.
  SNPs, 1, 2, and - have code codes 0, 1, and 2.

  These are objects are in data frames unlike other internal objects due
  to heavily used in processing data. The original data should be transfered
  to numerical codes in order to be passed to C codes. In C codes, we use
  integers, 0, 1, 2, ..., for coding nucleotides or SNPs and so on.

  Now, models and methods are implemented only for \code{.nucleotide} and
  \code{.snp}. Other objects are leaved for further extension.

  Data frames use factor formats as the default, and \code{as.character} is
  the way to transfer to the characters.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{code2nid}},
  \code{\link{nid2code}},
  \code{\link{snp2sid}},
  \code{\link{sid2snp}}.
}
\examples{
\dontrun{
.nucleotide
.snp
.codon
.amino.acid
.genetic.code
.missing.code
}
}
\keyword{data}
