\encoding{latin1}
\name{pmix}
\alias{pmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pmix}
\description{Calculates the carbonate chemistry after mixing of two water samples with different pCO2}
\usage{
pmix(flag, var1, var2, pCO2s, wf, S=35, T=20, P=0, Pt=0, Sit=0, k1k2="l", kf="pf", pHscale="T")
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
	\item{var1}{Value of the first variable in mol/kg except for pH and for pCO2 in \eqn{\mu}atm}
	\item{var2}{Value of the second variable in mol/kg except for pH}
	\item{pCO2s}{Partial pressure of the high CO2 component in \eqn{\mu}atm}
	\item{wf}{Weight fraction of the high CO2 seawater vs normal seawater}
	\item{S}{Salinity}
	\item{T}{Temperature in degrees Celsius}
	\item{P}{Hydrostatic pressure in bar (surface = 0)}
	\item{Pt}{Concentration of total phosphate in mol/kg}
	\item{Sit}{Concentration of total silicate in mol/kg}
	\item{k1k2}{"l" for using K1 and K2 from Lueker et al. and "r" for using K1 and K2 from Roy and al. , default is "l"}
	\item{kf}{"pf" for using Kf from Perez and Fraga (1987) and "dg" for using Kf from Dickson and Goyet (1979), default is "pf"}
	\item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
}


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
Gattuso J.-P. and Lavigne H, 2009. Perturbation experiments to investigate the impact of ocean acidification: approaches and software tools. \emph{Biogeosciences} \bold{6}, 4413-4439.
}

\author{
Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
pmix(flag=24, var1=384, var2=2302e-6, pCO2s=1e6, wf=0.003, S=34.3, T=16, P=0, pHscale="T", kf="pf", k1k2="l")
}

\keyword{utilities}
