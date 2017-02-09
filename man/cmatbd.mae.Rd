\name{cmatbd.mae}
\alias{cmatbd.mae}
\title{
Computes the treatment information matrix 
}
\description{
Computes the information matrix (C-matrix) for treatment effects under either the linear fixed effects model  or the linear mixed effects model setting for a given block design of size 2.
}
\usage{
cmatbd.mae(trt.N, blk.N, theta, des)
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
  \item{des}{
matrix, a \code{2 x b} block design with \code{b} blocks of size \code{k = 2} and \code{v} treatments.
}
}
\value{
Returns a \code{v x v} treatment information matrix (C-matrix).}
\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu D. B., Debusho L. K. and Haines L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}

\seealso{
\code{\link{optbdmaeAT}}, \code{\link{fixparbd.mae}}, \code{\link{intcbd.mae}}
}
\examples{

##Information matrix

     trt.N <- 3 
     blk.N <- 3 
     theta <- 0.2 
     dsgn <- intcbd.mae(trt.N = 3, blk.N = 3)

     cmatbd.mae(trt.N = 3, blk.N = 3, theta = 0.2, des = dsgn)
}
\keyword{Information matrix}
\keyword{C-matrix} 
