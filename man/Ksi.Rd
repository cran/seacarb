\encoding{latin1}
\name{Ksi}
\alias{Ksi}
\title{Dissociation constant of Si(OH)4}
\description{Dissociation constant of Si(OH)4 on total scale (mol/kg)}
\usage{
Ksi(S=35, T=25, P=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.} 

\value{
	\item{Ksi}{Dissociation constant of Si(OH)4 on total scale (mol/kg)}
}

\references{
DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.
}
             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl}
}

%\seealso{}

\examples{
  Ksi(35, 25, 0)
}

\keyword{utilities}
