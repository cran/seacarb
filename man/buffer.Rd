\encoding{latin1}
\name{buffer}
\alias{buffer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Buffer parameters of the seawater carbonate system}
\description{Returns buffer parameters of the seawater carbonate system.}
\usage{
buffer(flag, var1, var2, S=35, T=25, P=0, Pt=0, Sit=0)
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
  	\item{Pt}{enter value of the concentration of total phosphate in mol/kg}
  	\item{Sit}{enter the value of the total silicate in mol/kg}

}


\details{Vectors can be given in enter arguments. If the length of vectors entered are different the longer vector is kept and only the first value of other vectors is used. For instance, to compute parameters from one couple of variable for a range of temperatures, a vector with temperatures required can be given in enter and other arguments can be completed be one variable this variable will be used for each temperatures.}

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
}

\author{
Héloïse Lavigne, Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{

## Calcul with a couple of variables
buffer(flag=8, var1=8.2, var2=0.00234, S=35, T=25, P=0, Pt=0, Sit=0)

## Using vectors as arguments
flag <- c(8, 2, 8)
var1 <- c(8.2, 7.477544e-06, 8.2)
var2 <- c(0.002343955, 0.001649802, 2400e-6)
S <- c(35, 35, 30)
T <- c(25, 25, 30)
P <- c(0, 0, 0)
Pt <- c(0, 0, 0)
Sit <- c(0, 0, 0)

buffer(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit)

## Test for all flags 

flag <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25)

var1 <- c(8.200000, 7.477544e-06, 7.477544e-06, 7.477544e-06, 7.477544e-06, 8.2, 8.2, 8.2, 8.2, 0.001685024, 0.001685024, 0.001685024,  0.0002888382, 0.0002888382, 0.002391252, 264.2008, 264.2008, 264.2008, 264.2008, 264.2008)

var2 <- c(7.477544e-06, 0.001685024, 0.0002888382, 0.002391252, 0.001981340, 0.001685024, 0.0002888382, 0.002391252, 0.001981340, 0.0002888382, 0.002391252, 0.001981340,  0.002391252, 0.001981340, 0.001981340, 8.2, 0.001685024, 0.0002888382, 0.002391252, 0.001981340)

buffer(flag=flag, var1=var1, var2=var2)

}
\keyword{utilities}

