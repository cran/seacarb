\name{carb}
\alias{carb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameters of the seawater carbonate system}
\description{Returns parameters of the seawater carbonate system.}
\usage{
carb(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{a data frame containing the input parameters (see below)}
}


\details{The input data frame should contain the following variables in the right order:

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
}
  \item{var1}{enter value of the first variable in mol/kg except for pH}
  \item{var2}{enter value of the second variable in mol/kg except for pH}
  \item{S}{Salinity}
  \item{T}{Temperature in degrees Celsius}
  \item{P}{Hydrostatic pressure in bar (surface = 0)}
  \item{k1k2}{"r" for using K1 and K2 from Roy et al. and "m" for using K1 and K2 from Mehrbach,}
  \item{phflag}{The pH scale is either the total scale (0) or the free scale (1), default is 0}

Note that the results output is NO LONGER stored in a file.}

\value{The function returns a dataframe containing the following columns:

  \item{S}{Salinity}
  \item{T}{Temperature in degrees Celsius}
  \item{P}{Pressure in bar}
  \item{PH}{pH}
  \item{CO2}{CO2 concentration (mol/kg)}
  \item{pCO2}{pCO2, CO2 partial pressure (uatm)}
  \item{fCO2}{fCO2, CO2 fugacity (uatm)}
  \item{HCO3}{HCO3 concentration (mol/kg)}
  \item{CO3}{CO3 concentration (mol/kg)}
  \item{DIC}{DIC concentration (mol/kg)}
  \item{ALK}{ALK, total alkalinity (mol/kg)}
  \item{Oa}{Omega aragonite, aragonite saturation state}
  \item{Oc}{Omega calcite, calcite saturation state}
  \item{PhiD}{PhiD, chemical buffer factor (dpH/d[DIC]); input/output of dissolved CO2 (unit pH per mol/kg)}
  \item{BetaD}{BetaD, homogeneous buffer factor (dln(pCO2)/dln[DIC]); input/output of dissolved CO2}
  \item{PiD}{PiD, chemical buffer factor (dpCO2/d[DIC]); input/output of dissolved CO2 (uatm per mol/kg)}
  \item{PhiB}{PhiB, chemical buffer factor (dpH/d[DIC]); from input/output of bicarbonate (unit pH per mol/kg)}
  \item{BetaB}{BetaB, homogeneous buffer factor (dln(pCO2)/dln[DIC]); input/output of bicarbonate}
  \item{PiB}{PiB, chemical buffer factor (dpCO2/d[DIC]); input/output of dissolved CO2 (uatm per mol/kg)}
  \item{PhiC}{PhiC, chemical buffer factor (dpH/d[DIC]); input/output of carbonate (unit pH per mol/kg)}
  \item{BetaC}{BetaC, homogeneous buffer factor (dln(pCO2)/dln[DIC]); input/output of carbonate}
  \item{PiC}{PiC, chemical buffer factor (dpCO2/d[DIC]); input/output of carbonate (uatm per mol/kg)}
  \item{PhiH}{PhiH, chemical buffer factor (dpH/d[DIC]); input/output of strong acid (unit pH per mol/kg)}
  \item{PiH}{PiH, chemical buffer factor (dpCO2/d[DIC]); input/output of strong acid (uatm per mol/kg)}
}


\references{DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Frankignoulle, M. 1994 A complete set of buffer factors for acid/base CO2 system in seawater. \emph{Journal of Marine Systems} \bold{5}, 111-118.

Zeebe, R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}
\author{Aurelien Proye and Jean-Pierre Gattuso}

\examples{\dontrun{
x1 <- data.frame(flag=8,var1=8.2,var2=2400,S=35,T=25,P=0,k1k2='r',phflag=0)
x2 <- data.frame(flag=8,var1=8.2,var2=2400,S=30,T=30,P=0,k1k2='r',phflag=0)
x <- rbind(x1,x2)


carb(x)

\dontrun{The input parameters are pH and total alkalinity

The output is:
flag Salinity Temperature Pressure  pH          CO2     pCO2     fCO2
1     8       35          25        0 8.2 7.477544e-06 263.3691 262.5400
11    8       30          30        0 8.2 6.914443e-06 267.9528 267.1586
          HCO3          CO3         DIC    ALK OmegaAragonite OmegaCalcite
1  0.001649802 0.0003108231 0.001968102 0.0024       4.929621     7.478931
11 0.001622366 0.0003289170 0.001958198 0.0024       5.522402     8.351172
        PhiD    BetaD      PiD      PhiB    BetaB       PiB      PhiC     BetaC
1  -1369.120 8.160777 1.092065 -182.8422 1.956302 0.2617900 1003.4359 -4.248173
11 -1363.176 8.157162 1.116197 -192.5034 2.010711 0.2751386  978.1697 -4.135739
          PiC      PhiH       PiH
1  -0.5684854 -1186.278 0.8302753
11 -0.5659197 -1170.673 0.8410583


## Data entry via a file is NO LONGER possible

}
}
}
\keyword{utilities}
