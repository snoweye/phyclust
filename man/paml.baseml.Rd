\name{paml.baseml}
\alias{paml.baseml}
\alias{paml.baseml.control}
\alias{paml.baseml.show.default}
\alias{Class baseml}
\title{Phylogenetic Analysis by Maximum Likelihood for Nucleotide Sequences}
\description{
  This function modifies the original standalone code of \code{baseml}
  in PAML developed by Yang (1997) for phylogenetic analysis by maximum
  likelihood. This function provides a way to generate an ancestral tree
  for given central sequences clustered by \code{phyclust}.
}
\usage{
paml.baseml(X, seqname = NULL, opts = NULL, newick.trees = NULL)
paml.baseml.control(...)
paml.baseml.show.default()
}
\arguments{
  \item{X}{nid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{seqname}{sequence names.}
  \item{opts}{options as the standalone version, provided by \code{paml.baseml.contol}.}
  \item{newick.trees}{a vector/list contains NEWICK trees for \code{runmode = 2}.}
  \item{...}{for other possible opts and values. See PAML manual for details.}
  \item{show}{show opts and values.}
}
\details{
  The function \code{paml.baseml} directly reuses the C code of \code{baseml}
  of PAML, and the function \code{paml.baseml.control} is to generate controls
  for \code{paml.baseml} as the file \code{baseml.ctl} of PAML.

  The \code{seqname} should be consistent with \code{X}, and the leaf nodes
  of \code{newick.trees}.

  The options \code{opts} is followed from the original \code{baseml.ctl}
  except \code{seqfile}, \code{treefile} and \code{outputfile} will be omitted.

  \code{paml.baseml.control} generates default \code{opts}, and
  \code{paml.baseml.show.default} displays annotations for the default
  \code{opts}.
}
\value{
  This function returns a list, and each element stores one line of outputs
  of \code{baseml} separated by newline. The list stores in a class
  \code{baseml}. All the output of \code{baseml} of PAML will be saved in
  several files, and these will be read in by \code{scan}. Further post
  processing can be done by parsing the returning vector. The details of
  output format can found on the website
  \url{http://abacus.gene.ucl.ac.uk/software/paml.html} and its manual.

  Note that some functionalities of \code{baseml} of PAML are changed in
  \code{paml.baseml} due to the complexity of input and output. The changes
  include such as disable the option \code{G} and rename the file \code{2base.t}
  to \code{pairbase.t}.

  Typically, the list contains the original output of \code{baseml} including
  \code{pairbase.t}, \code{mlb}, \code{rst}, \code{rst1}, and \code{rub} if they
  are not empty. The best tree (unrooted by default) will be stored in
  \code{best.tree} parsed from \code{mlb} based on the highest log likelihood.
  All output to \code{STDOUT} are stored in \code{stdout}. No \code{STDIN} are
  allowed.

  Note that the print function for the class \code{baseml} will only show
  the \code{best.tree} only. Use \code{str} or \code{names} to see the whole
  returns of the list.
}
\section{Warning(s)}{
  Carefully read the PAML's original document before using the
  \code{paml.baseml} function, and \code{paml.baseml} may not be ported
  well from \code{baseml} of PAML. Please double check with the standalone
  version.

  \code{baseml} may not be a well designed program, and may run slowly.
  If it were stuck, temporary files would all store at a directory obtained
  by \code{tempfile("/paml.baseml.")}.

  \code{baseml} has its own options and settings which may be different
  than \pkg{phyclust} and \pkg{ape}. For example, the following is from
  the PAML's document, ``In PAML, a rooted tree has a bifurcation at the root,
  while an unrooted tree has a trifurcation or multifurcation at the root.''
  i.e. \code{paml.baseml} uses a rooted result for an unrooted tree, as well
  as for a rooted tree.

  \code{baseml} also needs a sequence file which is dumped from R (duplicated
  from memory) for \code{paml.baseml}, and this file can be very big if
  sequences are too long or number of sequences is too large. Also,
  \code{paml.baseml} may take long time to search the best tree if data are
  large or initial trees are not provided.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://snoweye.github.io/phyclust/}

  Yang, Z. (1997)
  \dQuote{PAML: a program package for phylogenetic analysis by maximum
          likelihood},
  \emph{Computer Applications in BioSciences}, \bold{13}, 555-556.

  Yang, Z. (2007)
  \dQuote{PAML 4: a program package for phylogenetic analysis by maximum
          likelihood},
  \emph{Molecular Biology and Evolution}, \bold{24}, 1586-1591.
  \url{http://abacus.gene.ucl.ac.uk/software/paml.html}
}
\author{
  Yang, Z. (1997) and Yang, Z. (2007)

  Maintain: Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{print.baseml}},
  \code{\link{write.paml}}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

paml.baseml.show.default()

### Generate data.
set.seed(123)
ret.ms <- ms(nsam = 5, nreps = 1, opts = "-T")
ret.seqgen <- seqgen(opts = "-mHKY -l40 -s0.2", newick.tree = ret.ms[3])
(ret.nucleotide <- read.seqgen(ret.seqgen))
X <- ret.nucleotide$org
seqname <- ret.nucleotide$seqname

### Run baseml.
opts <- paml.baseml.control(model = 4, clock = 1)
(ret.baseml <- paml.baseml(X, seqname = seqname, opts = opts))
(ret.baseml.init <- paml.baseml(X, seqname = seqname, opts = opts,
   newick.trees = ret.ms[3]))
ret.ms[3]

### Unrooted tree.
opts <- paml.baseml.control(model = 4)
(ret.baseml.unrooted <- paml.baseml(X, seqname = seqname, opts = opts))

### More information.
opts <- paml.baseml.control(noisy = 3, verbose = 1, model = 4, clock = 1)
ret.more <- paml.baseml(X, seqname = seqname, opts = opts)
# ret.more$stdout

### Plot trees
par(mfrow = c(2, 2))
plot(read.tree(text = ret.ms[3]), main = "true")
plot(read.tree(text = ret.baseml$best.tree), main = "baseml")
plot(read.tree(text = ret.baseml.init$best.tree), main = "baseml with initial")
plot(unroot(read.tree(text = ret.baseml.unrooted$best.tree)),
     main = "baseml unrooted")
}
}
\keyword{PAML}

