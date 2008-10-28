\encoding{latin1}
\name{Kb}
\alias{Kb}
\title{Dissociation constant of boric acid (mol/kg)}

\description{Dissociation constant of boric acid (mol/kg)}
\usage{
Kb(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  }

\value{
  \item{Kb}{Dissociation constant of boric acid (mol/kg), pHscale = totale scale}

}
\references{Dickson, A. G. 1990 Thermodynamics of the dissociation of boric acid in synthetic seawater from 273.15 to 318.15 K. \emph{Deep-Sea Research} \bold{375}, 755-766.

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.}

\author{
Héloïse Lavigne, Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}


\examples{
  Kb(35,25,0)
}
\keyword{utilities}
