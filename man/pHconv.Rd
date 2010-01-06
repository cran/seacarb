\encoding{latin1}
\name{pHconv}
\alias{pHconv}
\title{Conversion of pH}
\description{Converts pH from one scale to another one chosen between the total scale, the free scale and the seawater scale}
\usage{
pHconv(flag=1,pH=8.10,S=35,T=25, P=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flag}{choice of the type of conversion :
	flag=1: seawater scale to total scale
	flag=2: free scale to the total scale
	flag=3: total scale to the seawater scale
	flag=4: total scale to the free scale
	flag=5: seawater scale to the free scale
	flag=6: free scale to the seawater scale
	default is flag=1}
  \item{pH}{Enter the value of pH which need to be converted, default is 8.100}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
The function returns the values of pH converted
}

\details{The arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It is recommended to use either vectors with the same dimension or one vector for one argument and numbers for the other arguments.}

\references{
Dickson, A.G. & F.J. Millero, 1987 A comparison of the equilibrium constants for the dissociation of carbonic acid in seawater media. \emph{Deep-Sea Research} \bold{34}:1733-1743.
}
             
\author{
Heloise Lavigne and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{kconv}}.
}

\examples{
##To convert pH=8.10 from the seawater scale to the total scale
##at salinity=35, temperature=25oC and atmospheric pressure:

  pHc <- pHconv(flag=1, pH=8.10, S=35, T=25, P=0)

##note that pHc is the value of the pH converted in total scale

## By using vectors
## to convert the pH values : 8, 8.05, 8.10, 8.15, 8.20 
## from the free to the total scale

pH <- c(8, 8.05, 8.10, 8.15, 8.20)
pHc <- pHconv(flag=2, pH=pH, S=35, T=25, P=0)

## note that pHc is a vector containing the value of the pH converted 
## to the total scale
}

\keyword{utilities}
