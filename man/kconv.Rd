\encoding{latin1}
\name{kconv}
\alias{kconv}
\title{Conversion factors to change the pH scale of dissociation constants}
\description{Conversion factors from the total scale to the free and seawater scales}
\usage{
kconv(S=35,T=25, P=0, kf="pf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979), default is "pf". Note that conversion factors computation requires Kf. Here is given the choice of the formulation for Kf.} 
}

\value{
The function returns a list with 3 conversion factors :
  \item{ktotal2SWS}{to convert from the total scale to seawater scale}
  \item{ktotal2free}{to convert from the total scale to the free scale}
  \item{kfree2SWS}{to convert from the free scale to the seawater scale}  
}

\references{Dickson, A.G. & F.J. Millero, 1987 A comparison of the equilibrium constants for the dissociation of carbonic acid in seawater media. \emph{Deep-Sea Research} \bold{34}:1733-1743.}

             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl}
}

\seealso{
	\code{\link{pHconv}}.
}


\examples{
##To convert dissociation constants from the total scale to the free scale
## (at salinity=35, temperature=25oC and atmospheric pressure):
kconv(35,25,0)
conv <- kconv()
 c(K1_total=K1(),K1_SWS=K1()*conv$ktotal2SWS,K1_free=K1()*conv$ktotal2free)  
}

\keyword{utilities}
