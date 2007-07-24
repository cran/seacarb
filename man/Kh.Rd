\name{Kh}
\alias{Kh}
\title{Henry's constant mol/(kg/atm)}
\description{Henry's constant mol/(kg/atm)}
\usage{
Kh(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
  \item{Kh}{Henry's constant mol/(kg/atm)}
 
}

\references{Weiss, R. F. 1974 Carbon dioxide in water and seawater: the solubility of a non-ideal gas. \emph{Marine Chemistry} \bold{2}, 203-215.}

\author{Aurelien Proye and Jean-Pierre Gattuso}


\examples{
  Kh(35,25,0)
}
\keyword{utilities}
