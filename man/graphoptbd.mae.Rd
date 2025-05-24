\name{graphoptbd.mae}
\alias{graphoptbd.mae}
\title{
Creates the graphical layout of resultant optimal design
}
\description{
Creates the graphical layout of resultant A-, MV-, D- or E-optimal or near-optimal block design on a separate pop-up GUI tcltk window.
}
\usage{
graphoptbd.mae(trt.N, blk.N, theta, OptdesF, Optcrit, cbVal2)
}
\arguments{
  \item{trt.N}{
integer, specifying  number of treatments, \code{v}. 
}
  \item{blk.N}{
integer, specifying  number of arrays, \code{b}.
}
  \item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
  \item{OptdesF}{
matrix, a \code{2 x b} obtained optimal or near-optimal block design.
}
  \item{Optcrit}{
character specifying the optimality criteria to be used. Thus, \code{Optcrit} takes the letter \code{"A"}, \code{"MV"},\code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal block designs, respectively.
}
  \item{cbVal2}{
checkbox value.  It takes a value of zero or one. The default value of \code{cbVal2} is 0. Thus, if \code{cbVal2 = 0}, the function will display the graphical layout of the resultant optimal design generated using the treatment exchange algorithm. Similarly, if \code{cbVal = 1}, the function will display the graphical layout of the resultant optimal design generated using the array exchange algorithm. 
}
}
\details{
Detail discussions  concerning the constructions of a graphs can be found in \code{igraph} R package.
}
\value{
Returns the graphical layout of the resultant optimal or near-optimal block design '\code{OptdesF}' on a separate pop-up window. Furthermore, the function  \code{graphoptbd.mae} saves the graphical layout of the resultant optimal or near-optimal block design in .pdf format in a working subdirectory. 
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
\note{
When closing a pop-up window for graphical layout of the resultant designs (Graph plot),  if  the window is closed by 
clicking on the red button with "X" sign (top-right), the warning message  "Warning message: In \code{rm(list = cmd,} \code{envir = .tkplot.env):} object 'tkp \dots' not found"
will occur in R console irrespective of what command is executed next. To resolve this warning message, click 
on "\code{close}" menu that is located at the top-left of the graph plot pop-up window when closing this window. 
 }


\examples{
    ##To create the graphical layout of the D-optimal block design
    ##obtained using the treatment exchange algorithm for 
    trt.N <- 10  #Number of treatments
    blk.N <- 10  #Number of blocks
    theta <- 0.2 #theta value 
    OptdesF <- rbind(1:10, c(2:10,1)) #D-optimal design (loop design)

    graphoptbd.mae(trt.N = 10, blk.N = 10, theta = 0.2, OptdesF, Optcrit = "D", cbVal2 = 0)
}
\keyword{Graphical layout}
