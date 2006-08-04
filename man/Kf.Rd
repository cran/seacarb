\name{Kf}
\alias{Kf}

\title{Stability constant of hydrogen fluorid (mol/kg)}
\description{Stability constant of hydrogen fluorid (mol/kg)}
\usage{
Kf(S = 35, T = 25, P = 0, phflag = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{phflag}{The pH scale is either the total scale (0) or the free scale (1), default is 0}
}

\value{
  \item{Kf}{Stability constant of hydrogen fluorid (mol/kg)}

}
\references{Kf is given by Dickson and Riley (1979 in DOE, 1994).

Dickson, A. G. and Riley, J. P. 1979 The estimation of acid dissociation constants in seawater media from potentiometric titrations with strong base. I. The ionic product of water. \emph{Marine Chemistry} \bold{7}, 89-99.

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.}

\author{Aurelien Proye and Jean-Pierre Gattuso}


\examples{
  Kf(35,25,0,0)
}
\keyword{utilities}
