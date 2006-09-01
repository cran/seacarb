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
 
\value{
  \item{K1p}{First dissociation constant of phosphoric acid (mol/kg)}
 
}
\references{DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.}

\author{Aurelien Proye and Jean-Pierre Gattuso}

\seealso{
	\code{\link{K2p}}, \code{\link{K3p}}.
}

\examples{
  K1p(35,25,0)
}
\keyword{utilities}
