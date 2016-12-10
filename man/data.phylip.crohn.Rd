\name{data.phylip.crohn}
\alias{Data Crohn's Disease}
\docType{data}
\title{Crohn's Disease SNP Dataset in the phylip Format}
\description{
  Crohn's disease dataset is published by Hugot, et al. (2001).
}
\format{
  A text file in phylip format is stored in the data subdirectory.
}
\details{
  Crohn's disease dataset is used to perform haplotype grouping used
  in Tzeng's paper (2005).

  Totally, 1102 haplotypes/SNP sequences and 8 sites.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}

  Hugot, J.P., et al. (2001)
  \dQuote{Association of NOD2 Leucine-Rich Repeat Variants with
          Susceptibility to Crohn's Disease},
  \emph{Nature}, \bold{411}, 599-603.

  Tzeng, J.Y. (2005)
  \dQuote{Evolutionary-Based Grouping of Haplotypes in Association Analysis},
  \emph{Genetics Epidemiology}, \bold{28}, 220-231.
  \url{http://www4.stat.ncsu.edu/~jytzeng/software.php}
}
\author{
  Hugot, J.P., et al. (2001).
}
\seealso{
  \code{\link{read.phylip}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

data.path <- paste(.libPaths()[1], "/phyclust/data/crohn.phy", sep = "")
# edit(file = data.path)
my.snp <- read.phylip(data.path, code.type = "SNP")
str(my.snp)
}
}
\keyword{datasets}
