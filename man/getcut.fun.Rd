\name{getcut.fun}
\alias{getcut.fun}
\title{Tzeng's Method: Finding the Best Number of Clusters}
\description{
  For SNP sequences only, Tzeng's method (2005) uses an evolution approach
  to group haplotypes based on a deterministic transformation of haplotype
  frequency. This function find the best number of clusters based on
  Shannon information content.
}
\usage{
getcut.fun(pp.org, nn, plot = 0)
}
\arguments{
  \item{pp.org}{frequency of haplotypes, sorted in decreasing order.}
  \item{nn}{number of haplotypes.}
  \item{plot}{illustrated in a plot.}
}
\details{
  \code{pp.org} is summarized from \code{X} in \code{haplo.post.prob},
  \code{nn} is equal to the number of rows of \code{X}.

  This function is called by \code{haplo.post.prob} to determine
  the best guess of number of clusters.
  See Tzeng (2005) and Shannon (1948) for details.
}
\value{
  Return the best guess of number of clusters.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Tzeng, J.Y. (2005)
  \dQuote{Evolutionary-Based Grouping of Haplotypes in Association Analysis},
  \emph{Genetics Epidemiology}, \bold{28}, 220-231.
  \url{http://www4.stat.ncsu.edu/~jytzeng/software.php}

  Shannon, C.E. (1948)
  \dQuote{A mathematical theory of communication},
  \emph{Bell System Tech J}, \bold{27}, 379-423, 623-656.
}
\author{
  Jung-Ying Tzeng.

  Maintain: Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{haplo.post.prob}}.
}
\examples{
\dontrun{
data.path <- paste(.libPaths()[1], "/phyclust/data/crohn.phy", sep = "")
my.snp <- read.phylip(data.path, code.type = "SNP")
ret <- haplo.post.prob(my.snp$org, ploidy = 1)
getcut.fun(sort(ret$haplo$hap.prob, decreasing = TRUE),
           nn = my.snp$nseq, plot = 1)
}
}
\keyword{Haploclust}
