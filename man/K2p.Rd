\encoding{latin1}
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
\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.}

\author{
Héloïse Lavigne, Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{K1p}}, \code{\link{K3p}}.
}


\examples{
  K2p(35,25,0)
}
\keyword{utilities}
