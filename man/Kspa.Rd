\encoding{latin1}
\name{Kspa}
\alias{Kspa}

\title{Solubility product of aragonite (mol/kg)}
\description{Solubility product of aragonite (mol/kg)}
\usage{
Kspa(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
  \item{Kspa}{Solubility product of aragonite (mol/kg)}
  
}
\references{Mucci, A. 1983 The solubility of calcite and aragonite in seawater at various salinities, temperature, and one atmosphere total pressure. \emph{American Journal of Science} \bold{283}: 780-799.}

\author{Aurelien Proye and Jean-Pierre Gattuso}

\seealso{
	\code{\link{Kspc}}.
}

\examples{
Kspa(35,25,0)
}
\keyword{utilities}
