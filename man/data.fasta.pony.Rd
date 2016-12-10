\name{data.fasta.pony}
\alias{Data Pony 625}
\alias{Format fasta}
\docType{data}
\title{Great Pony 625 EIAV rev Dataset in the Fasta Format}
\description{
  Great pony 625 EIAV dataset is published by Baccam, P., et al. (2003),
  and they are also available on NCBI database.
  This is a follow-up study of \code{Data Pony 618}.
}
\format{
  A text file in fasta format is stored in the data subdirectory.
}
\details{
  EIAV rev dataset contains 62 nucleotide sequences and 406 sites.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}

  Baccam, P., et al. (2003)
  \dQuote{Subpopulations of Equine Infectious Anemia Virus Rev Coexist
          In Vivo and Differ in Phenotype},
  \emph{Journal of Virology}, \bold{77}, 12122-12131.
}
\author{
  Baccam, P., et al. (2003).
}
\seealso{
  \code{\link{read.phylip}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

data.path <- paste(.libPaths()[1], "/phyclust/data/pony625.fas", sep = "")
# edit(file = data.path)
my.pony.625 <- read.fasta(data.path)
str(my.pony.625)
}
}
\keyword{datasets}
