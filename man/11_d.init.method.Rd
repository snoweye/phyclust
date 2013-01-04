\name{.init.method}
\alias{.init.method}
\docType{data}
\title{Initialization Methods for EM Algorithms}
\description{
  The varied initialization methods are implemented in C.
  The first element is the default value.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.init.method}
\format{
  A character vector contains implemented initialization methods in C.
}
\details{
  \code{randomMu}, \code{NJ}, \code{randomNJ}, \code{PAM},
  \code{K-Medoids} and \code{manualMu} are implemented where
  the codes for the NJ are modified from \pkg{ape}, and the codes
  for the PAM method are modified from \pkg{cluster}.
  These methods are only provide initializations for EM algorithms.

  \itemize{
  \item{'randomMu'}{randomly picks centers and assigns all sequences near by
                    the center according an evolution distance.}
  \item{'NJ'}{bases on a neighbor-joining tree and partitions the tree by
              the long branches into subtrees to form clusters.}
  \item{'randomNJ'}{randomly partition a neighbor-joining tree into subtrees
                    to form clusters.}
  \item{'PAM'}{uses the partition around medoids algorithm
               to locate the centers of dataset.}
  \item{'K-Medoids'}{performs K-Means types algorithms to randomly and roughly
                     locate centers and form clusters.}
  \item{'manualMu'}{requires a vector containing class ids for all sequences.}
  }
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Saitou, N. and Nei, M. (1987),
  \dQuote{The Neighbor-Joining Method: A New Method for Reconstructing
          Phylogenetic Trees},
  \emph{Molecular Biology and Evolution}, \bold{4:4}, 406-425.

  Kaufman, L. and Rousseeuw, P.J. (1990),
  \emph{Finding Groups in Data: An Introduction to Cluster Analysis},
  Wiley.

  Theodoridis, S. and Koutroumbas, K. (2006),
  \emph{Pattern Recognition}, 3rd ed., pp. 635.
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.init.procedure}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
.init.method
}
}
\keyword{data}
