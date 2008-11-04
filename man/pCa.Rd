\encoding{latin1}
\name{pCa}
\alias{pCa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pCa}
\description{Calculates the changes in the saturation states of aragonite and calcite resulting from the manipulation of the calcium concentration}
\usage{
pCa(flag, var1, var2, Ca, S=35, T=20, P=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{flag}{select the couple of variables available. The flags which can be used are: 

flag = 1      pH and CO2 given

flag = 2      CO2 and HCO3 given

flag = 3      CO2 and CO3 given

flag = 4      CO2 and ALK given

flag = 5      CO2 and DIC given

flag = 6      pH and HCO3 given

flag = 7      pH and CO3 given

flag = 8      pH and ALK given

flag = 9      pH and DIC given

flag = 10     HCO3 and CO3 given

flag = 11     HCO3 and ALK given

flag = 12     HCO3 and DIC given

flag = 13     CO3 and ALK given

flag = 14     CO3 and DIC given

flag = 15     ALK and DIC given

flag = 21     pCO2 and pH given

flag = 22     pCO2 and HCO3 given

flag = 23     pCO2 and CO3 given

flag = 24     pCO2 and ALK given

flag = 25     pCO2 and DIC given
}
	\item{var1}{Value of the first  variable in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm}
	\item{var2}{Value of the second  variable in mol/kg, except for pH}
	\item{Ca}{Calcium concentration in mol/kg}
	\item{S}{Salinity}
	\item{T}{Temperature in degrees Celsius}
  	\item{P}{Hydrostatic pressure in bar (surface = 0)}
 }

\details{}


\value{The function returns a data frame containing the following columns:

  \item{comment}{The initial or final state water}
  \item{S}{Salinity}
  \item{T}{Temperature in degrees Celsius}
  \item{P}{Pressure in bar}
  \item{pH}{pH}
  \item{CO2}{CO2 concentration (mol/kg)}
  \item{pCO2}{pCO2, CO2 partial pressure (\eqn{\mu}atm)}
  \item{fCO2}{fCO2, CO2 fugacity (\eqn{\mu}atm)}
  \item{HCO3}{HCO3 concentration (mol/kg)}
  \item{CO3}{CO3 concentration (mol/kg)}
  \item{DIC}{DIC concentration (mol/kg)}
  \item{ALK}{ALK, total alkalinity (mol/kg)}
  \item{OmegaAragonite}{Omega aragonite, aragonite saturation state}
  \item{OmegaCalcite}{Omega calcite, calcite saturation state}

}



\references{
}

\author{
Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
pCa(flag=15, var1=2302e-6, var2=2050e-6, Ca=0.01028, S=35, T=20, P=0) # with normal Ca concentration
pCa(flag=15, var1=2302e-6, var2=2050e-6, Ca=0.01028/2, S=35, T=20, P=0) # with 0.5 * Ca concentration
}

\keyword{utilities}
