\encoding{latin1}
\name{K3p}
\alias{K3p}
\title{Third dissociation constant of phosphoric acid (mol/kg)}

\description{Third dissociation constant of phosphoric acid (mol/kg)}
\usage{
K3p(S = 35, T = 25, P = 0, pHscale = "T")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
 \item{K3p}{Third dissociation constant of phosphoric acid (mol/kg)}
}

\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.}

\author{
Héloïse Lavigne, Aurélien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{K1p}}, \code{\link{K2p}}.
}

\examples{
  K3p(35,25,0)
}
\keyword{utilities}
