\name{data.phylip.pony}
\alias{Data Pony 524}
\alias{Format phylip}
\docType{data}
\title{Great Pony 524 EIAV rev Dataset in the phylip Format}
\description{
  Great pony 524 EIAV dataset is published by Baccam, P., et al. (2003),
  and they are also available on NCBI database.
  There is a follow-up study, \code{Data Pony 625}.
}
\format{
  A text file in phylip format is stored in the data subdirectory.
}
\details{
  EIAV rev dataset contains 146 nucleotide sequences and 405 sites.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Belshan, M., et al. (2001)
  \dQuote{Genetic and Biological Variation in Equine Infectious Anemia
          Virus Rev Correlates with Variable Stages of Clinical Disease
          in an Experimentally Infected Pony},
  \emph{Virology}, \bold{279}, 185-200.
}
\author{
  Belshan, M., et al. (2001).
}
\seealso{
  \code{\link{read.fasta}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

data.path <- paste(.libPaths()[1], "/phyclust/data/pony524.phy", sep = "")
# edit(file = data.path)
my.pony.524 <- read.phylip(data.path)
str(my.pony.524)
}
}
\keyword{datasets}
