\encoding{latin1}
\name{Kw}
\alias{Kw}
\title{Ion product of water (mol2/kg2)}
\description{Ion product of water (mol2/kg2)}
\usage{
Kw(S = 35, T = 25, P = 0, phflag = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{phflag}{The pH scale is either the total scale (0) or the free scale (1), default is 0}
}

\value{
  \item{Kw}{Ion product of water (mol2/kg2)}
}

\references{DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Millero F. J., 1995. Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica et Cosmochimica Acta} \bold{59} 661-677.}

\author{
Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
  Kw(35,25,0,0)
}
\keyword{utilities}
