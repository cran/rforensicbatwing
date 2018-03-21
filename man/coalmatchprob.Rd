\name{Inference from STR data using coalescent theory}
\alias{batwing}
\alias{coalmatchprob}
\alias{print.batwing}
\alias{plot.forensicbatwing}
\title{Inference from STR data using coalescent theory (e.g. calculating forensic trace-suspect match probability)}
\description{
Estimates the trace-suspect match probability for a lineage haplotype of STR markers using coalescent theory.
}
\usage{
coalmatchprob(database, haplotype, 
  reps = 10, burnin = 0, treebetN = 10, Nbetsamp = 10, 
  muprior = "constant(0.003)", Nprior = "lognormal(9, 1)", alphaprior = NULL, 
  progress = TRUE, trace = FALSE)
batwing(database, 
  reps = 10, burnin = 0, treebetN = 10, Nbetsamp = 10, 
  muprior = "constant(0.003)", Nprior = "lognormal(9, 1)", alphaprior = NULL, 
  progress = TRUE, trace = FALSE) 
\method{print}{batwing}(x, ...)
\method{plot}{forensicbatwing}(x, ...)
}
\arguments{
  \item{database}{Reference STR database.}
  \item{haplotype}{Haplotype of the suspect.}
  \item{reps}{Number of output lines.}
  \item{burnin}{Number of reps to take before starting recording data.}
  \item{treebetN}{The number of times that changes to the genealogical tree are attempted before any changes to the hyperparameters are attempted. Thus BATWING outputs are separated by treebetN * Nbetsamp attempted tree updates.}
  \item{Nbetsamp}{The number of times that changes to hyperparameters are attempted between outputs.}
  \item{muprior}{Either a single prior distribution for the mutation rate or a vector of prior distributions (one for each locus).
    If only one prior is supplied, the same mutation rate is used for all loci. If one prior per locus is supplied, each locus has its own chain of mutation rates.}
  \item{Nprior}{Prior distribution of the effective population size.}
  \item{alphaprior}{If NULL, there is no growth (constant population size). If a prior distribution is specified, this gives exponential growth at rate alpha at all times.}
  \item{progress}{Whether to print progress or not.}
  \item{trace}{Whether to print extra trace information or not.}
  \item{x}{A \code{batwing} or \code{forensicbatwing} object.}
  \item{...}{Not used}
}
\details{
Note that the \code{batwing} function runs a standard coalescent inference as described in I.J. Wilson (1999, 2003).

Note that, in contrast to the original BATWING program, migration is not supported. Neither is BATWING's \code{sizemodel=2} (constant-sized population up to a time from where there is exponential growth).

Valid prior distributions:
\itemize{
\item{\code{uniform(v1, v2)} }{uniform on the interval \code{(v1, v2)}.}
\item{\code{constant(v1)} }{constant value \code{v1}.}
\item{\code{normal(v1, v2)} }{Normal distribution with mean = \code{v1} and sd = \code{v2}.}
\item{\code{lognormal(v1, v2)} }{If \code{X} has this distribution then \code{log(X)} has the \code{normal(v1, v2)} distribution.}
\item{\code{gamma(v1, v2)} }{Gamma distributinon with shape \code{v1} and rate \code{v2} giving mode = \code{(v1-1)/v2} and mean = \code{v1/v2}.}
\item{\code{beta(v1, v2)} }{Beta distribution with shape parameters \code{v1} and \code{v2} giving mean = \code{v1/(v1 + v2)} and the variance is \code{(v1*v2)/((v1 + v2)^2 * (v1 + v2 + 1))}}
}
}
\value{
  \item{coalmatchprob}{An object of type \code{forensicbatwing} and \code{batwing} (for now, \code{forensicbatwing} just means that there is a \code{p} column in the \code{result} matrix).}
  \item{batwing}{An object of type \code{batwing}}
  \item{batwing-object}{
    \describe{
      \item{parameters }{The parameters used}
      \item{priors }{The priors used}
      \item{result }{Matrix with information about the iterations (forensic match probability (if \code{coalmatchprob}), current value of \code{N}, \code{mu}'s, tree height \code{T}, total branch length \code{L}, population growth \code{alpha}, and loglikelihoods for times, mutations and priors).}
      \item{proposals_tree }{The number of new tree proposals.}
      \item{accepted_tree }{The number of new trees accepted.}
      \item{proposals_hyperparameters }{The number of new hyperparameters proposals.}
      \item{accepted_hyperparameters }{The number of new hyperparameters accepted.}
    }
  }
}
\source{
  \href{http://www.mas.ncl.ac.uk/~nijw/}{BATWING at Ian Wilson's homepage}
}

\references{
\href{http://www.mas.ncl.ac.uk/~nijw/batwing/batguide.pdf}{The BATWING manual}

I.J. Wilson, D.J. Balding, 
\emph{Genealogical inference from microsatellite data}, 
Genetics 150 (1998) 499-510.

I.J. Wilson, M.E. Weale, D.J. Balding, 
\emph{Inferences from DNA data: population histories, evolutionary processes and forensic match probabilities}, 
J. R. Stat. Soc. Ser. A 166 (2003) 155-201.

M.M. Andersen, A. Caliebe, A. Jochens, S. Willuweit, M. Krawczak, 
\emph{Estimating trace-suspect match probabilities for singleton Y-STR haplotypes using coalescent theory}, 
Forensic Sci. Int. Genet. (In Press, Corrected Proof \href{http://dx.doi.org/10.1016/j.fsigen.2012.11.004}{10.1016/j.fsigen.2012.11.004.}).
}
\author{
Mikkel Meyer Andersen and Ian Wilson
}

\examples{
\dontrun{
database <- matrix(c(1, 1, 2, 2, 1, 3), 3, 2)
haplotype <- c(1, 1)

# coalmatchprob:
coalmp <- coalmatchprob(database, haplotype, 
 reps = 1000, burnin = 0, treebetN = 10, Nbetsamp = 10, 
 muprior = c("normal(0.003, 0.001)",  "normal(0.005, 0.001)"), 
 Nprior = "lognormal(9, 1)", 
 alphaprior = NULL, 
 progress = TRUE, trace = FALSE)
coalmp

murange <- range(c(coalmp$result$mu1, coalmp$result$mu2))

par(mfrow = c(2, 2))
plot(coalmp)
plot(coalmp$result$N, type = "l", ylab = "N")

plot(coalmp$result$mu1, type = "l", col = "red", ylim = murange, ylab = "mu")
points(coalmp$result$mu2, type = "l", col = "blue")

hist(coalmp$result$mu1, col = "#FF000066", 
  xlim = murange, ylim = c(0, 250), main = NULL, xlab = "mu")
hist(coalmp$result$mu2, add = TRUE, col = "#0000FF66")
par(mfrow = c(1, 1))

# batwing:
bw <- batwing(database, 
 reps = 10000, burnin = 1000, treebetN = 10, Nbetsamp = 10, 
 muprior = "normal(0.003, 0.001)", 
 Nprior = "lognormal(9, 1)", 
 alphaprior = NULL, 
 progress = TRUE, trace = FALSE)
bw
}
}

\keyword{coalescent}
\keyword{coalescence}
