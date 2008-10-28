\encoding{latin1}
\name{Kspc}
\alias{Kspc}

\title{Solubility product of calcite (mol/kg)}
\description{Solubility product of calcite (mol/kg)}
\usage{
Kspc(S = 35, T = 25, P = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}


\value{
  \item{Kspc}{Solubility product of calcite (mol/kg)}

}
\references{Mucci, A. 1983 The solubility of calcite and aragonite in seawater at various salinities, temperature, and one atmosphere total pressure. \emph{American Journal of Science} \bold{283}: 780-799.}

\author{
Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{Kspa}}.
}

\examples{
	Kspc(35,25,0)
}
\keyword{utilities}
