\name{mmenubd.mae}
\alias{mmenubd.mae}
\title{
Creates the main menu tcltk window
}
\description{
Alternative to directly using the function \code{\link{optbdmaeAT}} on the R console, the function \code{mmenubd.mae} creates the main menu tcltk window with the list of four optimality criteria (A, MV, D and E) for search of optimal or near-optimal block designs for two-colour cDNA microarray experiments. }
\usage{
mmenubd.mae()
}
\details{
The main menu window created using the function \code{mmenubd.mae} contains the list of four optimality criteria (the A-,
MV -, D- and E) for block designs that allows user to choose an optimality criterion of interest. The function 
will then call for a function \code{fixpar.mae} for further option regarding the setup of parametric combinations and output, see \code{fixpar.mae} documentation. 
}
\author{
Dibaba B. Gemechu, Legesse K. Debusho, and Linda M. Haines
}

\seealso{
\code{\link{optbdmaeAT}}, \code{\link{fixparbd.mae}}
}
