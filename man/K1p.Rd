\encoding{latin1}
\name{K1p}
\alias{K1p}
\title{First dissociation constant of phosphoric acid (mol/kg)}

\description{First dissociation constant of phosphoric acid (mol/kg)}

\usage{
K1p(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.}

\value{
  \item{K1p}{First dissociation constant of phosphoric acid (mol/kg), pHscale = totale scale}
 
}
\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.}

\author{
Héloïse Lavigne, Aurélien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{K2p}}, \code{\link{K3p}}.
}

\examples{
  K1p(35,25,0)
}
\keyword{utilities}
