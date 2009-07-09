\encoding{latin1}
\name{psi}
\alias{psi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Molar ratio of CO2 released vs CaCO3 precipitated}
\description{Returns the molar ratio of CO2 released vs CaCO3 precipitated described by Frankignoulle et al. (1994).}
\usage{
psi(flag, var1, var2, S=35, T=20, P=0, Pt=0, Sit=0, pHscale="T", 
	kf="pf", k1k2="l")
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
	\item{var1}{enter value of the first  variable in mol/kg, except for pH and for pCO2 in \eqn{\mu}atm}
	\item{var2}{enter value of the second  variable in mol/kg, except for pH}
	\item{S}{Salinity}
	\item{T}{Temperature in degrees Celsius}
  	\item{P}{Hydrostatic pressure in bar (surface = 0)}
  	\item{Pt}{Concentration of total phosphate in mol/kg}
  	\item{Sit}{Concentration of total silicate in mol/kg}
	\item{pHscale}{choice of pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Goyet (1979), default is "pf"}
	\item{k1k2}{"l" for using K1 and K2 from Lueker et al. and "r" for using K1 and K2 from Roy et al. , default is "l"}
}


\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length. For instance, to compute parameters from one couple of variable for a range of temperatures, a vector with temperatures required can be given in enter and other arguments can be completed be one variable this variable will be used for each temperatures.}

\value{The function returns a data frame containing the following columns:
	\item{psi}{ratio of CO2 released vs CaCO3 precipitated (mol/mol)}
}

\references{

Frankignoulle, M. 1994 A complete set of buffer factors for acid/base CO2 system in seawater. \emph{Journal of Marine Systems} \bold{5}, 111-118.

Frankignoulle, M., Canon, C. & Gattuso, J.-P., 1994. Marine calcification as a source of carbon dioxide- Positive feedback of increasing atmospheric CO2. \emph{Limnology and Oceanography} \bold{2}, 458-462.
}

\author{
Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link[seacarb]{speciation}}.
}

\examples{
## Calculation using the numerical example given in Frankignoulle et al. (1994)
psi(flag=24, var1=350, var2=2400e-6, S=35, T=25, P=0, Pt=0, 
	Sit=0, pHscale="T", kf="pf", k1k2="l")
}
\keyword{utilities}

