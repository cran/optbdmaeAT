\name{intcbd.mae}
\alias{intcbd.mae}
\title{
Generates initial connected block design
}
\description{
Generates a random initial connected block design for a given number of blocks \code{b} of size \code{k = 2} and the number of treatments \code{v}.}
\usage{
intcbd.mae(trt.N, blk.N)
}
\arguments{
 \item{trt.N}{
integer, specifying number of treatments, \code{v}. 
}
  \item{blk.N}{
integer, specifying number of arrays, \code{b}.
}
}
\value{
Returns a \code{2 x b} connected block design with \code{b} blocks of size \code{k = 2} and number of treatments \code{v}. 
}
\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu D. B., Debusho L. K. and Haines L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}

\seealso{
\code{\link{optbdmaeAT}}, \code{\link{cmatbd.mae}}
}
\examples{
    #Initial connected block design for 
    trt.N <- 4 #Number of treatments
    blk.N <- 4 #Number of blocks

    intcbd.mae(trt.N = 4, blk.N = 4)
}
\keyword{Initial block design}
\keyword{Connected block design}% __ONLY ONE__ keyword per line
