\name{ms}
\alias{ms}
\alias{Class ms}
\title{Generating Samples under a Wright-Fisher Neutral Model of Genetic Variation}
\description{
  This function modifies the original standalone code of \code{ms()} developed by
  Hudson (2002) for generating samples/coalescent trees under a Wright-Fisher
  neutral model.
}
\usage{
ms(nsam = NULL, nreps = 1, opts = NULL, temp.file = NULL,
   tbs.matrix = NULL)
}
\arguments{
  \item{nsam}{number of samples/coalescent trees, usually greater than 2.}
  \item{nreps}{number of replications.}
  \item{opts}{options as the standalone version.}
  \item{temp.file}{temporary file for ms output.}
  \item{tbs.matrix}{a matrix for 'tbs' options given in opts.}
}
\details{
  This function directly reuses the C code of \code{ms} by arguments
  as input from the \code{opts}. The options \code{opts} is followed from the
  original \code{ms} except \code{nsam} and \code{nreps}.
  Note that stdin, stdout, and pipe are all disable from \code{opts}.

  For examples, options commonly used in \pkg{phyclust} are:
  \itemize{
    \item "-T": generate trees in a neutral model.
    \item "-G": generate trees with a population growth rate, e.g. "-G 0.5".
  }
  These will return trees in a NEWICK format which can be read by the
  \code{read.tree()} of \pkg{ape} and passed to \code{seqgen()} to generate
  sequences. 

  \code{temp.file} allows users to specify ms output file themselves, but
  this file will not be deleted nor converted into R after the call to
  \code{ms()}. Users should take care the readings. By default, \code{ms()}
  uses a system temp file to store the output which is converted into R
  after the call and is deleted after converting.

  \code{tbs.matrix} is a matrix to specify the values of \code{tbs} given
  in \code{opts}. See \code{demo('simu_ms_tbs')} for an example how to
  use this option.
}
\value{
  This function returns a vector, and each element stores one line of STDOUT
  of \code{ms()} separated by newline. The vector stores in a class \code{ms}.
  The details of output format can found on the website
  \url{http://home.uchicago.edu/~rhudson1/source.html} and its manual.
}
\section{Warning(s)}{
  Carefully read the \code{ms}'s original document before using the \code{ms()}
  function.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Hudson, R.R. (2002)
  \dQuote{Generating Samples under a Wright-Fisher Neutral Model of Genetic
          Variation},
  \emph{Bioinformatics}, \bold{18}, 337-338.
  \url{http://home.uchicago.edu/~rhudson1/source.html}
}
\author{
  Hudson, R.R. (2002).

  Maintain: Wei-Chen Chen \email{wccsnow@gmail.com}
}
\seealso{
  \code{\link{print.ms}()},
  \code{\link[ape]{read.tree}()},
  \code{\link[ape]{bind.tree}()},
  \code{\link{seqgen}()}.
}
\examples{
\dontrun{
library(phyclust, quiet = TRUE)

ms()

# an ancestral tree
set.seed(1234)
(ret.ms <- ms(nsam = 3, opts = "-T -G 0.1"))
(tree.anc <- read.tree(text = ret.ms[3]))
tree.anc$tip.label <- paste("a", 1:K, sep = "")

# adjacent descendant trees to the ancestral tree
K <- 3
N <- 12
N.k <- c(3, 4, 5)
ms.dec <- NULL         # a list to store trees of ms
tree.dec <- NULL       # a list to store the trees in phylo class
tree.joint <- tree.anc
for(k in 1:K){
  ms.dec[[k]] <- ms(N.k[k], opts = "-T -G 1.0")
  tree.dec[[k]] <- read.tree(text = ms.dec[[k]][3])
  tree.dec[[k]]$tip.label <- paste("d", k, ".", 1:N.k[k], sep = "")
  tree.joint <- bind.tree(tree.joint, tree.dec[[k]],
                          where = which(tree.joint$tip.label ==
                                        paste("a", k, sep = "")))
}
str(tree.joint)

# plot trees
par(mfrow = c(2, 3))
plot(tree.anc, main = paste("anc (", K, ")", sep = ""))
axis(1)
for(k in 1:K){
  plot(tree.dec[[k]], main = paste("dec", k, " (", N.k[k], ")", sep = ""))
  axis(1)
}
plot(tree.joint, main = paste("joint (", N, ")", sep = ""))
axis(1)

# use tbs option (an example from msdoc.pdf by Hudson, R.R.)
tbs.matrix <- matrix(c(3.0, 3.5, 5.0, 8.5), nrow = 2)
ret <- ms(nsam = 5, nreps = 2, opts = "-t tbs -r tbs 1000",
          tbs.matrix = tbs.matrix)
print(ret)
}
}
\keyword{ms}
