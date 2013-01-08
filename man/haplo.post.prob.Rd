\name{haplo.post.prob}
\alias{haplo.post.prob}
\title{Tzeng's Method: Haplotype Grouping for SNP Sequences}
\description{
  For SNP sequences only, Tzeng's method (2005) uses an evolution approach
  to group haplotypes based on a deterministic transformation of haplotype
  frequency. This is a modified version of original function,
  \code{haplo.score.RD.unphased.fun}.
}
\usage{
  haplo.post.prob(X, ploidy = 2, skip.haplo = 1e-07, K = NULL)
}
\arguments{
  \item{X}{sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{ploidy}{ploidy, no effect for phase known, keep consistence only.}
  \item{skip.haplo}{lower bound of haplotypes frequencies.}
  \item{K}{number of clusters.}
}
\details{
  \code{X} should be a phase known haplotype data. For phase unknown and
  Tzeng's method (2006) are not tested yet.

  If \code{K} is NULL, the result of \code{getcut.fun} will be used.
}
\value{
  See the original paper and source codes' documents for details.
  The function returns a list contains:
  \item{'haplo'}{summarized data set in a list contains:
    \tabular{ll}{
      'haplotype' \tab unique haplotypes, dim = \eqn{N_{unique} \times L}{N.X.unique * L}. \cr
      'hap.prob'  \tab frequency of haplotypes. \cr
      'post'      \tab posterior probabilities of phase unknown haplotypes. \cr
      'hap1code'  \tab unique ids of 'haplotype'. \cr
      'hap2code'  \tab unique ids of 'haplotype', no effect if ploidy = 2. \cr
      'indx.subj' \tab id of subjects. \cr
    }
  }
  \item{'FD.id'}{unique ids of 'haplotype' for full dimension analysis.}
  \item{'RD.id'}{unique ids of 'haplotype' for reduced dimension analysis.}
  \item{'FD.post'}{posterior probabilities for full dimension analysis.}
  \item{'RD.post'}{posterior probabilities for reduced dimension analysis.}
  \item{'g.truncate'}{number of clusters}
}
\section{ToDo(s)}{
  \itemize{
    \item test codes for phased unknown cases.
  }
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Tzeng, J.Y. (2005)
  \dQuote{Evolutionary-Based Grouping of Haplotypes in Association Analysis},
  \emph{Genetics Epidemiology}, \bold{28}, 220-231.
  \url{http://www4.stat.ncsu.edu/~jytzeng/software.php}
}
\author{
  Jung-Ying Tzeng.

  Maintain: Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{getcut.fun}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

data.path <- paste(.libPaths()[1], "/phyclust/data/crohn.phy", sep = "")
my.snp <- read.phylip(data.path, code.type = "SNP")
ret <- haplo.post.prob(my.snp$org, ploidy = 1)
str(ret)
}
}
\keyword{Haploclust}
