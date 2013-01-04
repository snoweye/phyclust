\name{.EMControl}
\alias{.EMControl}
\title{EM Control Generator}
\description{
  Generate an EM control (\code{.EMC}) controlling the options,
  methods, conditions and models of EM algorithms.
  As \code{.EMC}, this function generate a default template.
  One can either modify \code{.EMC} or employ this function to
  control EM algorithms.
}
\usage{
.EMControl(exhaust.iter = 1, fixed.iter = 5,
    short.iter = 100, EM.iter = 1000,
    short.eps = 1e-2, EM.eps = 1e-6,
    cm.reltol = 1e-8, cm.maxit = 5000,
    nm.abstol.Mu.given.QA = 1e-8, nm.reltol.Mu.given.QA = 1e-8,
    nm.maxit.Mu.given.QA = 500,
    nm.abstol.QA.given.Mu = 1e-8, nm.reltol.QA.given.Mu = 1e-8,
    nm.maxit.QA.given.Mu = 5000,
    est.non.seg.site = FALSE, max.init.iter = 50,
    init.procedure = .init.procedure[1],
    init.method = .init.method[1],
    substitution.model = .substitution.model$model[1],
    edist.model = .edist.model[1], identifier = .identifier[1],
    code.type = .code.type[1], em.method = .em.method[1],
    boundary.method = .boundary.method[1], min.n.class = 1,
    se.type = FALSE, se.model = .se.model[1], se.constant = 1e-2)
}
\arguments{
  \item{exhaust.iter}{number of iterations for "exhaustEM", default = 1.}
  \item{fixed.iter}{number of iterations for "RndpEM", default = 5.}
  \item{short.iter}{number of short-EM steps, default = 100.}
  \item{EM.iter}{number of long-EM steps, default = 1000.}
  \item{short.eps}{tolerance of short-EM steps, default = 1e-2.}
  \item{EM.eps}{tolerance of long-EM steps, default = 1e-6.}

  \item{cm.reltol}{relative tolerance for a CM step, default = 1e-8}
  \item{cm.maxit}{maximum number iteration for a CM step, default = 5000.}

  \item{nm.abstol.Mu.given.QA}{see \sQuote{Details}, default = 1e-8}
  \item{nm.reltol.Mu.given.QA}{see \sQuote{Details}, default = 1e-8}
  \item{nm.maxit.Mu.given.QA}{see \sQuote{Details}, default = 500.}
  \item{nm.abstol.QA.given.Mu}{see \sQuote{Details}, default = 1e-8}
  \item{nm.reltol.QA.given.Mu}{see \sQuote{Details}, default = 1e-8}
  \item{nm.maxit.QA.given.Mu}{ see \sQuote{Details}, default = 5000.}
  \item{est.non.seg.site}{estimate non-segregation sites, default = FALSE.}
  \item{max.init.iter}{maximum number of initialization iteration, default = 50.}

  \item{init.procedure}{initialization procedure, default = "exhaustEM".}
  \item{init.method}{initialization method, default = "randomMu".}
  \item{substitution.model}{substitution model, default = "JC69".}
  \item{edist.model}{evolution distance, default = "D\_J69".}
  \item{identifier}{identifier, default = "EE".}
  \item{code.type}{code type, default = "NUCLEOTIDE".}
  \item{em.method}{EM method, default = "EM".}
  \item{boundary.method}{boundary method, default = \code{ADJUST}.}
  \item{min.n.class}{minimum number of sequences in a cluster, default = 1.}

  \item{se.type}{sequencing error type, default = FALSE.}
  \item{se.model}{sequencing error model, default = "CONVOLUTION".}
  \item{se.constant}{constrained constant, default = 1e-2.}
}
\details{
  \code{exhaust.iter}, \code{fixed.iter}, \code{short.iter}, and
  \code{short.eps} are used to control the iterations of initialization
  procedures and methods.

  \code{EM.iter} and \code{EM.eps} are used to control the EM iterations.

  \code{cm.reltol} and \code{cm.maxit} are used to control the ECM
  iterations.

  Arguments starting with \code{nm.} are options for the Nelder-Mead
  method as in \code{optim}. The C codes of Nelder-Mead are modified
  from the R math library and the options are all followed.
  \code{abstol} and \code{reltol} are for absolute and relative tolerances. 
  \code{Mu.given.QA} is for maximizing the profile function of
  \eqn{\mu_k}{mu_k} given \eqn{Q_k}{Q_k}, and
  \code{QA.given.Mu} is for maximizing the profile function of
  \eqn{Q_k}{Q_k} given \eqn{\mu_k}{mu_k}.

  \code{est.non.seg.site} indicates whether to estimate the states of
  center sequences. If FALSE, the states will be fixed as the
  non segregating sites. Usually, there is no need to estimate.

  \code{max.init.iter} is for certain initialization methods, e.g.
  \code{randomNJ} and \code{K-Medoids} need few tries to obtain an
  appropriate initial state.

  \code{init.procedure} and \code{init.method} are for initializations.

  \code{min.n.class} is the minimum number of sequences in a cluster to
  avoid bad initialization state and degenerated clusters.

  \code{se.type}, \code{se.model}, and \code{se.constant} which are used only
  for sequencing error models and only for nucleotide data without labels.
}
\value{
  This function returns a list as \code{.EMC}. 

  The sequencing error controls are stored in
  \code{se.type}, \code{se.model}, and \code{se.constant}, for
  sequencing error type, model, and constrained constant of errors,
  respectively.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{.EMC}},
  \code{\link{.boundary.method}},
  \code{\link{.code.type}},
  \code{\link{.edist.model}},
  \code{\link{.em.method}},
  \code{\link{.identifier}},
  \code{\link{.init.method}},
  \code{\link{.init.procedure}},
  \code{\link{.substitution.model}},
  \code{\link[stats]{optim}},
  \code{\link{phyclust}},
  \code{\link{phyclust.se}}.
}
\examples{
\dontrun{
# The same as .EMC
.EMControl()

# Except code.type, all others are the same as .EMC
.EMControl(code.type = "SNP")
.EMControl(code.type = .code.type[2])
}
}
\keyword{programming}
