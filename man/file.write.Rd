\name{file.write}
\alias{write.fasta}
\alias{write.fasta.format}
\alias{write.phylip}
\alias{write.phylip.format}
\alias{write.paml}
\alias{write.paml.format}
\title{Write Data to Files by Formats}
\description{
  Two major file formats are supported in \pkg{phyclust},
  \code{Format phylip} and \code{Format fasta}.
  \strong{These functions only write files in basic syntax.}
}
\usage{
write.fasta(seqdata, filename, classid = NULL, seqname = NULL,
    width.line = 60, lower.case = FALSE, code.type = .code.type[1],
    sep = "")
write.fasta.format(seqdata, filename, classid = NULL, seqname = NULL,
    width.line = 60, sep = "")

write.phylip(seqdata, filename, classid = NULL, seqname = NULL,
    width.seqname = 10, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1], sep = "")
%    code.type = .code.type[1], sep = "", add.space = 10)
write.phylip.format(seqdata, filename, classid = NULL, seqname = NULL,
    width.seqname = 10, width.line = 60, sep = "")
%    width.seqname = 10, width.line = 60, sep = "", add.space = 10)

write.paml(seqdata, filename, classid = NULL, seqname = NULL,
    width.seqname = 10, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1], sep = "")
write.paml.format(seqdata, filename, classid = NULL, seqname = NULL,
    width.seqname = 10, width.line = 60, sep = "")
}
\arguments{
  \item{seqdata}{a matrix contains sequence ids as \code{X} in \code{phyclust}.}
  \item{filename}{a file name where data is written to.}
  \item{classid}{class id of sequences.}
  \item{seqname}{sequence names.}
  \item{width.seqname}{number of characters of sequence names to be stored.}
  \item{width.line}{width of lines for breaking lines.}
  \item{lower.case}{use lower case of letters to write}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
  \item{sep}{a character to split sites, "" (default) and "," for "CODON".}
%  \item{add.space}{add a space for every few sites.}
}
\details{
  \code{write.fasta}, \code{write.phylip}, and \code{write.paml} are general
  functions call \code{write.fasta.format}, \code{write.phylip.format} and
  \code{write.paml.format}.

  \code{write.fasta.format}, \code{write.phylip.format}, and
  \code{wirte.paml.format} will not do any transformation for input sequences,
  but directly write them into the file as \code{code.type = NULL} in
  \code{write.fasta}, \code{write.phylip} and \code{write.paml}.

  Note that PAML uses one of PHYLIP format to deal with sequence files, so
  \code{write.paml.format} is to write files in a different format of
  \code{write.phylip.format}. The main purpose of \code{write.paml} and
  \code{write.paml.format} is to generate files for \code{pamle.baseml}.

  \code{sep} can specify a character which is used to split sites in file.
  By default, "" denote no character between sites. Only "CODON" id requires
  a character to avoid ambiguity.
}
\value{
  Save a text file.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{read.fasta}},
  \code{\link{read.phylip}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

# PHYLIP
data.path <- paste(.libPaths()[1], "/phyclust/data/crohn.phy", sep = "")
my.snp <- read.phylip(data.path, code.type = "SNP")
write.phylip(my.snp$org, "new.crohn.phy", code.type = "SNP")

# FASTA
data.path <- paste(.libPaths()[1], "/phyclust/data/pony625.fas", sep = "")
(my.pony <- read.fasta(data.path))
write.fasta(my.pony$org, "new.pony.fas")

# PAML
write.paml(my.pony$org, "new.pony.pam")

# Amino acid in PHYLIP
aa.aid <- nid2aid(my.pony$org)
aa.acode <- aid2acode(aa.aid)
write.phylip(aa.aid, "new.pony.aa.phy", code.type = "AMINO_ACID")
write.phylip.format(aa.aid, "new.pony.aa.aid.phy", sep = ",")
write.phylip.format(aa.acode, "new.pony.aa.acode.phy")

# Amino acid in FASTA
write.fasta(aa.aid, "new.pony.aa.phy", code.type = "AMINO_ACID")
write.fasta.format(aa.aid, "new.pony.aa.aid.phy", sep = ",")
write.fasta.format(aa.acode, "new.pony.aa.acode.fas")

# Amino acid in PAML
write.paml(aa.aid, "new.pony.aa.pam", code.type = "AMINO_ACID")
write.paml.format(aa.aid, "new.pony.aa.aid.pam", sep = ",")
write.paml.format(aa.acode, "new.pony.aa.acode.pam")
}
}
\keyword{file I/O}
