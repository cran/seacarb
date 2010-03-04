\encoding{latin1}
\name{buffer}
\alias{buffer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Buffer parameters of the seawater carbonate system}
\description{Returns buffer parameters of the seawater carbonate system.}
\usage{
buffer(flag, var1, var2, S=35, T=25, P=0, Pt=0, Sit=0, k1k2="x", 
	kf="x", pHscale="T")
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
	\item{k1k2}{"l" for using K1 and K2 from Lueker et al. (2000), "m" from Millero et al. (2006) and "r" from Roy et al. (1993). "x" is the default flag; the default value is then "l", except if T is outside the range 2 to 35oC and/or S is outside the range 19 to 43. In these cases, the default value is "m". }
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Riley (1979). "x" is the default flag; the default value is then "pf", except if T is outside the range 9 to 33oC and/or S is outside the range 10 to 40. In these cases, the default is "dg".}
	\item{pHscale}{choice of pH scale: "T" for the total scale, "F" for the free scale and "SWS" for using the seawater scale, default is "T" (total scale)}}


\details{The Lueker et al. (2000) constant is recommended by Guide to Best Practices for Ocean CO2 Measurements (2007). It is critical to consider that each formulation is valid in specific ranges of temperature and salinity:

\emph{For K1 and K2:}
\itemize{
\item Roy et al. (1993): S ranging between 0 and 45 and T ranging between 0 and 45oC.

\item Lueker et al. (2000): S ranging between 19 and 43 and T ranging between 2 and 35oC.

\item Millero et al. (2006): S ranging between 0 and 50 and T ranging between 0 and 50oC.
}

\emph{For Kh:}
\itemize{
\item Perez and Fraga (1987): S ranging between 10 and 40 and T ranging between 9 and 33oC.

\item Dickson and Goyet (1979): S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

The arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It is recommended to use either vectors with the same dimension or one vector for one argument and numbers for the other arguments.
}

\value{The function returns a data frame containing the following columns:

\item{PhiD}{PhiD, chemical buffer factor (dpH/d[DIC]); input/output of dissolved CO2 (unit pH per mol/kg)}
  \item{BetaD}{BetaD, homogeneous or Revelle buffer factor (dln(pCO2)/dln[DIC]); input/output of dissolved CO2}
  \item{PiD}{PiD, chemical buffer factor (dpCO2/d[DIC]); input/output of dissolved CO2 (\eqn{\mu atm} per mol/kg)}
  \item{PhiB}{PhiB, chemical buffer factor (dpH/d[DIC]); from input/output of bicarbonate (unit pH per mol/kg)}
  \item{BetaB}{BetaB, homogeneous buffer factor (dln(pCO2)/dln[DIC]); input/output of bicarbonate}
  \item{PiB}{PiB, chemical buffer factor (dpCO2/d[DIC]); input/output of dissolved CO2 (\eqn{\mu atm} per mol/kg)}
  \item{PhiC}{PhiC, chemical buffer factor (dpH/d[DIC]); input/output of carbonate (unit pH per mol/kg)}
  \item{BetaC}{BetaC, homogeneous buffer factor (dln(pCO2)/dln[DIC]); input/output of carbonate}
  \item{PiC}{PiC, chemical buffer factor (dpCO2/d[DIC]); input/output of carbonate (\eqn{\mu atm} per mol/kg)}
  \item{PhiH}{PhiH, chemical buffer factor (dpH/d[ALK]); input/output of strong acid (unit pH per mol/kg)}
  \item{PiH}{PiH, chemical buffer factor (dpCO2/d[ALK]); input/output of strong acid (\eqn{\mu atm} per mol/kg)}
}

\references{
Frankignoulle, M. 1994 A complete set of buffer factors for acid/base CO2 system in seawater. \emph{Journal of Marine Systems} \bold{5}, 111-118.

Lueker, T. J. Dickson, A. G. and Keeling, C. D. 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Millero, F. J. Graham, T. B. Huang, F. Bustos-Serrano, H. and Pierrot, D. 2006 Dissociation constants of carbonic acid in seawateras a function of salinity and temperature.  \emph{Marine Chemistry} \bold{100}, 80-84.  

Perez, F. F. and Fraga, F. 1987 Association constant of fluoride and hydrogen ions in seawater. \emph{Marine Chemistry} \bold{21}, 161-168.

Roy, R. N. Roy, L. N. Vogel, K. M. Porter-Moore, C. Pearson, T. Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45ÂoC. \emph{Marine Chemistry} \bold{44}, 249-267.
}

\author{
Heloise Lavigne, Aurelien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{

## Computation with a couple of variables
buffer(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Pt=0, 
	Sit=0, pHscale="T", kf="pf", k1k2="l")

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)
kf <- c("pf", "pf", "pf")
k1k2 <- c("l", "l", "l")
pHscale <- c("T", "T", "T")
buffer(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, 
	Sit=Sit, kf=kf, k1k2=k1k2, pHscale=pHscale)


## Test for all flags 

flag <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25)

var1 <- c(8.200000, 7.477544e-06, 7.477544e-06, 7.477544e-06, 7.477544e-06, 8.2, 
	8.2, 8.2, 8.2, 0.001685024, 0.001685024, 0.001685024,  0.0002888382, 
	0.0002888382, 0.002391252, 264.2008, 264.2008, 264.2008, 264.2008, 264.2008)

var2 <- c(7.477544e-06, 0.001685024, 0.0002888382, 0.002391252, 0.001981340, 
	0.001685024, 0.0002888382, 0.002391252, 0.001981340, 0.0002888382, 0.002391252,
	0.001981340,  0.002391252, 0.001981340, 0.001981340, 8.2, 0.001685024, 
	0.0002888382, 0.002391252, 0.001981340)

buffer(flag=flag, var1=var1, var2=var2)

}
\keyword{utilities}

