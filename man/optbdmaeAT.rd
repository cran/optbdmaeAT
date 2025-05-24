\name{optbdmaeAT}
\alias{optbdmaeAT}
\alias{optbdmaeAT.default}
\alias{print.optbdmaeAT}
\alias{summary.optbdmaeAT}
\alias{print.summary.optbdmaeAT}
\title{
Optimal block designs for two-colour cDNA microarray experiments
}
\description{
Used to compute A-, MV-, D- or E-optimal or near-optimal block designs for two-colour cDNA microarray experiments under either the linear fixed effects model or the linear mixed effects model settings using either the array exchange or treatment exchange algorithms of Debusho, Gemechu and Haines (2018) <doi.org/10.1080/03610918.2018.1429617>.}
\usage{
optbdmaeAT(trt.N, blk.N, theta, nrep, itr.cvrgval, Optcrit = "", Alg = "", ...)

\method{optbdmaeAT}{default}(trt.N, blk.N, theta, nrep, itr.cvrgval, Optcrit = "", Alg = "", ...)
\method{print}{optbdmaeAT}(x, ...)
\method{summary}{optbdmaeAT}(object, ...)
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
integer, specifying number of iterations required for convergence during the exchange procedure. 
}
  \item{Optcrit}{
character, specifying the optimality criteria to be used. \code{Optcrit} takes the letter \code{"A"}, \code{"MV"},\code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal block designs, respectively.
}
  \item{x}{
the object to be printed.
}
  \item{object}{
an object of class \code{"optbdmaeAT"}.
}
  \item{Alg}{
character string used to specify the algorithm to be used. Possible values of \code{Alg} are \code{Alg="trtE"} for the treatment exchange algorithm and \code{Alg="arrayE"} for the array exchange algorithm: see 'Details'.
}
  \item{\dots}{
not used.
}
}
\details{
\code{optbdmaeAT} computes optimal or near-optimal block design for the two-colour cDNA microarray experiments  
where the interest is in a comparison of all possible elementary treatment contrasts. The function computes A-, MV-, D- and E-optimal 
or near optimal block designs via calling of eight sub-functions \code{\link{Aoptbd.maeT}}, \code{\link{Aoptbd.maeA}}, 
\code{\link{MVoptbd.maeT}}, \code{\link{MVoptbd.maeA}}, \code{\link{Doptbd.maeT}}, \code{\link{Doptbd.maeA}}, 
\code{\link{Eoptbd.maeT}} and \code{\link{Eoptbd.maeA}}, respectively. Each function requires an initial connected block designs, 
generated using the function \code{\link{intcbd.mae}}.  

The minimum value of \code{trt.N} and \code{blk.N} is 3 and \code{trt.N} should be less than or equal to \code{blk.N - 1}. 
The linear fixed effects model results for given \code{trt.N} and \code{blk.N} are obtained by setting \code{theta = 0.0}.

\code{Alg} specifies the exchange algorithm of Debusho, Gemechu and Haines (2018). If \code{Alg = "trtE"}, the function 
\code{optbdmaeAT} perform  the treatment exchange procedure through deletion and addition of treatments at a time and selects a 
design with best treatment exchange with respect to the optimality criterion value. If \code{Alg = "arrayE"}, the function 
\code{optbdmaeAT} perform the array exchange procedure through deletion and addition of candidate arrays at a time and selects a 
design with best array exchange with respect to the optimality criterion value.

\code{nrep} takes a value of greater than or equal to 2. However, to ensure optimality of the resultant design, 
the \code{nrep} should be greater than or equal to 10 and in addition, as \code{trt.N} and \code{blk.N} increase, 
to ensure optimality of resultant design, it is advised to further increase the value of \code{nrep}
up to greater than or equal to 100. However, it has to be noted that as \code{trt.N} or \code{blk.N} or
 \code{nrep} or all of them increses, computer time required to generate optimal or near-optimal
block design increases.

\code{itr.cvrgval} number of iterations during exchange procedure. It takes a value between 2 and \code{blk.N}. It is used 
to speedup the computer search time by setting how long should the user should wait for the exchange process to obtain any 
different (if any) design than the one that was produced as the result of the preceding exchange of the current array in the initial 
design with candidate array. This is mainly effective if \code{blk.N} is very large. For example \code{itr.cvrgval = 2}, means the 
exchange procedure will jump to the next array test if the exchange of the two preceding arrays with candidate arrays results with the 
same efficient designs. The function  will not give error message if the users set \code{itr.cvrgval > blk.N} and it will automatically 
set \code{itr.cvrgval = blk.N}. The smaller the \code{itr.cvrgval} means the faster the exchange procedure is, but this will reduce the 
chance of getting optimal block design and users are advised to set \code{itr.cvrgval} closer to \code{blk.N}. 
}
\value{
Returns the resultant A-, MV-, D- or E-optimal or near-optimal block design with its corresponding score value and parametric combination 
saved in excel file in a working directory. In addition, the function \code{optbdmaeAT} displays the graphical layout of the resultant 
optimal or near-optimal block designs. Specifically: 

\item{call}{the method call.}         
\item{v}{number of treatments.}
\item{b}{number of blocks}
\item{theta}{theta value.}
\item{nrep}{number of replications of the optimization procedure.}  
\item{itr.cvrgval}{number of iterations required for convergence during the exchange procedure.}                          
\item{Optcrit}{optimality criteria.}                  
\item{Alg}{algorithm used.}
\item{OptdesF}{a \code{2 x blk.N} obtained optimal or near-optimal block design.}
\item{Optcrtsv}{score value of the optimality criteria \code{'Optcrit'} of the resultant optimal or near-optimal block design \code{'OptdesF'}.}
\item{file_loc, file_loc2}{location where the summary of the resultant optimal or near-optimal block design is saved in .csv format.}
\item{equireplicate}{logical value indicating whether the resultant optimal or near-optimal block design is equireplicate or not.}
\item{vtrtrep}{vector of treatment replication of the resultant optimal or near-optimal block design.}     
\item{Cmat}{the C-matrix or  treatment information matrix of the  optimal or near-optimal block design.} 

The graphical layout of the resultant optimal or near-optimal block design.   

NB: The function \code{optbdmaeAT} also saves the summary of the resultant optimal or near-optimal block design in .csv format in the working directory. 
Furthermore, the function reports only one final optimal or near-optimal block design, however, there is a possibility 
of more than one optimal or near-optimal block designs for a given parametric combination. 
The function \code{\link{graphoptbd.mae}} can be used to view and rearrange the graphical layout of the resultant 
optimal or near-optimal block design on \code{tcltk} window. Alternative to the function \code{optbdmaeAT}, a
GUI tcltk window can be used to generate optimal or near-optimal block designs, see \code{\link{mmenubd.mae}} and \code{\link{fixparbd.mae}}.   

}
\references{
Debusho, L. K., Gemechu, D. B. and Haines, L. (2018). Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Communications in Statistics - Simulation and Computation, https://doi.org/10.1080/03610918.2018.1429617}.

Gemechu D. B., Debusho L. K. and Haines L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
\seealso{
\code{\link{mmenubd.mae}}, \code{\link{fixparbd.mae}}, \code{\link{intcbd.mae}}
}
\examples{
  \donttest{
  ##To obtain the A-optimal or near-optimal block design using treatment exchange algorithm, set
  trt.N <- 3 #Number of treatments
  blk.N <- 3 #Number of blocks
  theta <- 0 #theta value
  nrep <- 5  #Number of replications
  itr.cvrgval <- 6 #Number of iterations required during the exchange procedure
  Optcrit <- "A"   #Optimality criteria
  Alg <- "trtE"    #Algorithm
  
  Aoptexample <- optbdmaeAT(trt.N = 3, blk.N = 3, theta = 0, nrep = 5, 
                            itr.cvrgval = 6, Optcrit = "A", Alg = "trtE")
  
  summary(Aoptexample)
}
}

\keyword{A-optimal block designs}
\keyword{D-optimal block designs}
\keyword{E-optimal block designs}
\keyword{MV-optimal block designs}
\keyword{Microarray experiment} 
\keyword{Treatment exchange algorithm}
\keyword{Array exchange algorithm} 
