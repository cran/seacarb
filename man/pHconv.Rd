\encoding{latin1}
\name{pHconv}
\alias{pHconv}
\title{Conversion factors to change pH scale}
\description{Conversion factors from the total scale to the free and seawater scales}
\usage{
pHconv(S=35,T=25, P=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
The function returns a list with 3 conversion factors :
  \item{pHtotal2SWS}{to convert from the total scale to seawater scale}
  \item{pHtotal2free}{to convert from the total scale to the free scale}
  \item{pHfree2SWS}{to convert from the free scale to the seawater scale}  
}

\references{Dickson, A.G. & F.J. Millero, 1987 A comparison of the equilibrium constants for the dissociation of carbonic acid in seawater media. \emph{Deep-Sea Research} \bold{34}:1733-1743.}
             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl}
}

\seealso{
	\code{\link{kconv}}.
}

\examples{
##To convert pH from the total scale to the seawater and free scales
## (at salinity=35, temperature=25oC and atmospheric pressure):
  pHconv(35, 25, 0)
  conv <- pHconv()
  c(pHtotal=8, pHSWS=8+conv$pHtotal2SWS, pHfree=8+conv$pHtotal2free)  
}

\keyword{utilities}
