\name{.init.procedure}
\alias{.init.procedure}
\docType{data}
\title{Initialization Procedures for EM Algorithms}
\description{
  The varied initialization procedures are implemented in C.
  The first element is the default value.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{.init.procedure}
\format{
  A character vector contains implemented initialization procedures in C.
}
\details{
  \code{exhaustEM}, \code{emEM}, \code{RndEM}, and \code{RndpEM} are
  implemented. Based on initialization states given by a initialization
  method, see \code{.init.method} for more information. These procedures
  will search a better starting states for final EM steps.

  \itemize{
  \item{'exhaustEM'}{runs each initialization with EM steps
                     until convergent, and pick the best one of the
                     convergence as the return result.}
  \item{'emEM'}{uses few short EM steps to improve initialization, then
                pick the best of initialization state for long EM steps,
                and returns the final result.}
  \item{'RandEM'}{bases on initialization methods to generate initialization
                  states, the number is equal to short EM steps,
                  then pick the best of initialization state for long EM steps,
                  and returns the final result.}
  \item{'RandEM'}{bases on initialization methods to generate initialization
                  states and run a fixed number of EM steps,
                  until total steps exhaust short EM steps,
                  then pick the best of initialization state for long EM steps,
                  and returns the final result.}
  }
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}

  Biernacki, C. and Celeux, G. and Govaert, G. (2003)
  \dQuote{Choosing starting values for the EM algorithm for getting the
          highest likelihood in multivariate Gaussian mixture models},
  \emph{Computational Statistics and Data Analysis}, \bold{41:3}, 561-575.

  Maitra, R. (2009)
  \dQuote{Initializing partition-optimization algorithms},
  \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  \bold{6:1}, 114-157.
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.init.method}},
  \code{\link{.EMControl}},
  \code{\link{phyclust}}.
}
\examples{
\dontrun{
.init.procedure
}
}
\keyword{data}
