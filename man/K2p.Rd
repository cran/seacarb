\name{K2p}
\alias{K2p}
\title{Second dissociation constant of phosphoric acid (mol/kg)}

\description{Second dissociation constant of phosphoric acid (mol/kg)}
\usage{
K2p(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
  \item{K2p}{Second dissociation constant of phosphoric acid (mol/kg)}

}
\references{DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.}

\author{Aurelien Proye and Jean-Pierre Gattuso}

\seealso{K1p.rd, K3p.rd}

\examples{
  K2p(35,25,0)
}
\keyword{utilities}