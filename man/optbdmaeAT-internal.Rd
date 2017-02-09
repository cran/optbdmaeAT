\name{optbdmaeAT-internal}
\alias{optbdmaeAT-internal}
\alias{Aoptbd.maeA}
\alias{Aoptbd.maeT}
\alias{MVoptbd.maeA}
\alias{MVoptbd.maeT}
\alias{Doptbd.maeA}
\alias{Doptbd.maeT}
\alias{Eoptbd.maeA}
\alias{Eoptbd.maeT}
\title{
Internal functions
}
\description{
Functions for internal usage only.
}
\usage{
 
## Computes A-optimal or near-optimal block designs
## using array exchange algorithm
Aoptbd.maeA(trt.N, blk.N, theta, nrep, itr.cvrgval)
 
## Computes A-optimal or near-optimal block designs
## using treatment exchange algorithm
Aoptbd.maeT(trt.N, blk.N, theta, nrep, itr.cvrgval)
 
## Computes MV-optimal or near-optimal block designs
## using array exchange algorithm
MVoptbd.maeA(trt.N, blk.N, theta, nrep, itr.cvrgval)
 
## Computes MV-optimal or near-optimal block designs
## using treatment exchange algorithm
MVoptbd.maeT(trt.N, blk.N, theta, nrep, itr.cvrgval)
 
## Computes D-optimal or near-optimal block designs
## using array exchange algorithm
Doptbd.maeA(trt.N, blk.N, theta, nrep, itr.cvrgval)
 
## Computes D-optimal or near-optimal block designs
## using treatment exchange algorithm
Doptbd.maeT(trt.N, blk.N, theta, nrep, itr.cvrgval)
 
## Computes E-optimal or near-optimal block designs
## using array exchange algorithm
Eoptbd.maeA(trt.N, blk.N, theta, nrep, itr.cvrgval)
 
## Computes E-optimal or near-optimal block designs
## using treatment exchange algorithm
Eoptbd.maeT(trt.N, blk.N, theta, nrep, itr.cvrgval)
}
\arguments{
  \item{trt.N}{
integer, specifying number of treatments, \code{v}. 
}
  \item{blk.N}{
integer, specifying number of arrays, \code{b}.
}
  \item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
  \item{nrep}{
integer, specifying number of replications of the optimization procedure. 
}
  \item{itr.cvrgval}{
integer, specifying number of iterations required for convergence during the exchange procedure. See \code{\link{optbdmaeAT}} documentation for details.
}
}
\details{
These functions are handled via a generic function \code{\link{optbdmaeAT}}. Please refer to the \code{\link{optbdmaeAT}} documentation for details. 
}

\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu D. B., Debusho L. K. and Haines L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
\seealso{
\code{\link{optbdmaeAT}}
}
\keyword{A-optimal block designs}
\keyword{D-optimal block designs}
\keyword{E-optimal block designs}
\keyword{MV-optimal block designs}
\keyword{Microarray experiment} 
\keyword{Treatment exchange algorithm}
\keyword{Array exchange algorithm} 
