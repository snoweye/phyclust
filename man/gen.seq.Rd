\name{gen.seq}
\alias{gen.seq.HKY}
\alias{gen.seq.SNP}
\title{Generate Sequences Given a Rooted Tree.}
\description{
  These functions call \code{seqgen} to generate sequences by
  evolutions models based on a rooted tree. \code{gen.seq.HKY} is to
  generate nucleotide sequences, and \code{gen.seq.SNP} is to generate
  SNP sequences.
}
\usage{
gen.seq.HKY(rooted.tree, pi, kappa, L, rate.scale = 1,
            anc.seq = NULL)
gen.seq.SNP(rooted.tree, pi, L, rate.scale = 1, anc.seq = NULL)
}
\arguments{
  \item{rooted.tree}{a rooted tree in \code{Class phylo}.}
  \item{pi}{equilibrium probabilities, sums to 1.}
  \item{kappa}{transition and transversion bias.}
  \item{L}{number of sites.}
  \item{rate.scale}{a scale to all branch lengths.}
  \item{anc.seq}{an ancestral sequence either in nids or sids, length = \eqn{L}{L}.}
}
\details{
  The \code{rooted.tree} should be in a \code{phylo} class of \pkg{ape},
  and may be created by \code{ms}.

  The \code{pi} has length 4 for nucleotide sequences, and 2 for SNP sequences.

  The \code{rate.scale} is equivalent to rescale the height of
  \code{rooted.tree}.
}
\value{
  Return an object in \code{Class seqgen}.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{gen.star.tree}},
  \code{\link{seqgen}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

# Generate a tree
set.seed(1234)
ret.ms <- ms(nsam = 5, nreps = 1, opts = "-T")
tree.ms <- read.tree(text = ret.ms[3])

# Generate nucleotide sequences
anc.HKY <- rep(0:3, 3)
pi.HKY <- c(0.2, 0.2, 0.3, 0.3)
kappa <- 1.1
L <- length(anc.HKY)
set.seed(1234)
paste(nid2code(anc.HKY, lower.case = FALSE), collapse = "")
(HKY.1 <- gen.seq.HKY(tree.ms, pi.HKY, kappa, L, anc.seq = anc.HKY))

# evolve 5 times longer
(HKY.2 <- gen.seq.HKY(tree.ms, pi.HKY, kappa, L,
                      rate.scale = 5, anc.seq = anc.HKY))

# Generate SNP sequences
anc.SNP <- rep(0:1, 6)
pi.SNP <- c(0.4, 0.6)
L <- length(anc.SNP)
set.seed(1234)
paste(sid2snp(anc.SNP), collapse = "")
(SNP.1 <- gen.seq.SNP(tree.ms, pi.SNP, L, anc.seq = anc.SNP))

# evolve 5 times longer
(SNP.2 <- gen.seq.SNP(tree.ms, pi.SNP, L, rate.scale = 5,
                      anc.seq = anc.SNP))
}
}
\keyword{simulation}
