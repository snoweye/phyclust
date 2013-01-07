\name{seqgen}
\alias{seqgen}
\alias{Class seqgen}
\title{Seq-Gen}
\description{
  This function modifies the original standalone code of \code{seq-gen}
  developed by Rambaut, A. and Grassly, N.C. (1997).
}
\usage{
seqgen(opts = NULL, rooted.tree = NULL, newick.tree = NULL, input = NULL,
       temp.file = NULL)
}
\arguments{
  \item{opts}{options as the standalone version.}
  \item{rooted.tree}{a rooted tree which sequences are generated according to.}
  \item{newick.tree}{a NEWICK tree which sequences are generated according to.}
  \item{input}{optional inputs of seq-gen, e.g. ancestral sequences.}
  \item{temp.file}{temporary file for seqgen output.}
}
\details{
  This function directly reuses the C code of \code{seq-gen} by arguments
  as input from the STDIN. The options \code{opts} is followed from the
  original \code{seq-gen} except an input tree.

  Input either a \code{rooted.tree} or a \code{newick.tree}, and
  \code{rooted.tree} should have a \code{Class phylo}.

  For examples, options commonly used in \code{phyclust} are:
  \itemize{
    \item "-m": set an evolution model, e.g. "-mHKY".
    \item "-t": set transition/transversion ratio, e.g. "-t0.7".
    \item "-f": equilibrium probabilities of A, C, G, and T, e.g. "-f0.1,0.2,0.3,0.4".
    \item "-l": length of sequences, e.g. "-l10".
    \item "-s": scale rate for the total height of input tree, "-s0.2".
    \item "-k": index of ancestral sequence in input file, see \code{gen.seq.HKY}.
  }
  These will return sequences in \code{Format phylip} which can be
  read by \code{read.seqgen} and transfered into an object with
  \code{Class seq.data}.

  \strong{Note}:
  \itemize{
    \item \code{input} and \code{rooted.tree}/\code{newick.tree} can not be
          submitted at the same time.
    \item \code{seq-gen} use the order A, C, G, T.
    \item \code{-t} is ts/tv ratio which is not equal to \eqn{\kappa}{kappa}.
    \item See more examples in \code{gen.seq.HKY} and \code{gen.seq.SNP}.
  }

  \code{temp.file} allows users to specify seqgen output file themselves, but
  this file will not be deleted nor converted into R after the call to
  \code{seqgen}. Users should take care the readings. By default, \code{seqgen}
  uses a system temp file to store the output which is converted into R
  after the call and is deleted after converting.
}
\value{
  This function returns a vector, and each element stores one line of STDOUT
  of \code{seq-gen} separated by newline. The vector stores in a class
  \code{seqgen}. The details of output format can found on the website
  \url{http://tree.bio.ed.ac.uk/software/seqgen/} and its manual.
}
\section{Warning(s)}{
  Carefully read the \code{seq-gen}'s original document before using the
  \code{seqgen} function.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Rambaut, A. and Grassly, N.C. (1997)
  \dQuote{Seq-Gen: An Application for the Monte Carlo Simulation of DNA
          Sequence Evolution along Phylogenetic Trees},
  \emph{Computer Applications In The Biosciences}, \bold{13:3}, 235-238.
  \url{http://tree.bio.ed.ac.uk/software/seqgen/}
}
\author{
  Rambaut, A. and Grassly, N.C. (1997).

  Maintain: Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{print.seqgen}},
  \code{\link[ape]{read.tree}},
  \code{\link{ms}},
  \code{\link{gen.seq.HKY}},
  \code{\link{gen.seq.SNP}}.
}
\examples{
\dontrun{
set.seed(123)
ret.ms <- ms(nsam = 5, nreps = 1, opts = "-T")
seqgen(opts = "-mHKY -l40 -s0.2", newick.tree = ret.ms[3])

}
}
\keyword{Seq-Gen}
