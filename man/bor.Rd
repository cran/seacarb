\encoding{latin1}
\name{bor}
\alias{bor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Total boron concentration (mol/kg)}
\description{total boron concentration (\eqn{mol\ kg^{-1}})}
\usage{
bor(S, T, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}


\value{
  \item{bor}{total boron concentration (\eqn{mol\ kg^{-1}}))}

}
\references{DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.}

\author{
Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}


\examples{
bor(35,25,0)
}

\keyword{utilities}
