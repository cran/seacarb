\encoding{latin1}
\name{Ks}
\alias{Ks}
\title{Stability constant of hydrogen sulfate (mol/kg)}
\description{Stability constant of hydrogen sulfate (mol/kg)}
\usage{
Ks(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
  \item{Ks}{Stability constant of hydrogen sulfate (mol/kg), pHscale = free scale}
}

\references{Ks is given by Dickson (1990) in Guide to Best Practices in Ocean CO2 Measurements 2007).

Dickson, A. G. 1990 Standard potential of the reaction: AgCI(s) + 1/2H2(g) = Ag(s) + HCI(aq), and the standard acidity constant of the ion HSO4 in synthetic sea water from 273.15 to 318.15 K. \emph{Journal of Chemical Thermodynamics} \bold{22}, 113-127.

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.}

\author{
Héloïse Lavigne, Aurélien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
  Ks(S=35,T=25,P=0)
}
\keyword{utilities}
