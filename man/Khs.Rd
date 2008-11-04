\encoding{latin1}
\name{Khs}
\alias{Khs}
\title{Dissociation constant of hydrogen sulfide (mol/kg)}
\description{Dissociation constant of hydrogen sulfide (mol/kg)}
\usage{
Khs(S=35, T=25, P=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
	\item{Khs}{Dissociation constant of hydrogen sulfide on the total scale (mol/kg)}
}

\references{
Millero F. J., 1995. Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} 59: 661-677.
}
             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl}
}

%\seealso{}

\examples{
  Khs(35,25,0)
  plot(Tseq <- seq(0,30,by=0.1), Khs(T=Tseq),xlab="Temperature,dgC",ylab="Khs")
}

\keyword{utilities}
