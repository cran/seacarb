\name{carb}
\alias{carb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameters of the seawater carbonate system}
\description{Returns parameters of the seawater carbonate system.}
\usage{
carb(flag = 8, var1 = 8.2, var2 = 2400, S = 35, T = 25, P = 0, k1k2 = "r", 
phflag = 0, ini = "s")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flag}{select the couple of variables to be used, see below, default is 8 that is pH and total alkalinity

The flags which can be used are: 

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
  \item{var1}{enter value of the first variable in umol/kg except for pH, default is pH = 8.2}
  \item{var2}{enter value of the second variable in umol/kg except for pH, default is total alkalinity = 2400 umol/kg}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{k1k2}{"r" for using K1 and K2 from Roy et al. and "m" for using K1 and K2 from Mehrbach, default is "r"}
  \item{phflag}{The pH scale is either the total scale (0) or the free scale (1), default is 0}
  \item{ini}{"s" if the input data are given through the keyboard or "f" if there are provided in a file (see example below), default = s}
}

\value{
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
\details{Note that the results output is also stored in the file carb.out located in your working directory. It is a text file with tab delimiters. This file is created if it does not exist and the results are appended if it is already present in the directory.}


\references{DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Frankignoulle, M. 1994 A complete set of buffer factors for acid/base CO2 system in seawater. \emph{Journal of Marine Systems} \bold{5}, 111-118.

Zeebe, R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}
\author{Aurelien Proye and Jean-Pierre Gattuso}

\examples{\dontrun{##Data entry via the keyboard}
carb(8, 8.2, 2400, 35, 25, 0) 

\dontrun{The input parameters are pH=8.2, total alkalinity = 2400 umol/kg,
 Salinity = 35, temperature = 25 oC and pressure = 0 bar (that is surface water)

The output is:

Salinity:                35 
Temperature:             25 oC
Pressure:                0 bar
pH:                      8.2 
CO2:                     7.477544e-06 (mol/kg)
pCO2:                    263.3691 (uatm)
fCO2:                    262.54 (uatm)
HCO3:                    0.001649802 (mol/kg)
CO3:                     0.0003108231 (mol/kg)
DIC:                     0.001968102 (mol/kg)
ALK:                     0.0024 (mol/kg)
Omega aragonite:         4.929621 
Omega calcite:           7.478931 
PhiD:                    -1369.120 
BetaD:                   8.160777 
PiD:                     1.092065 
PhiB:                    -182.8422 
BetaB:                   1.956302 
PiB:                     0.2617900 
PhiC:                    1003.436 
BetaC:                   -4.248173 
PiC:                     -0.5684854 
PhiH:                    -1186.278 
PiH:                     0.8302753 

## Data entry via a file
carb(ini = "f")

The file must be located in your home directory and provide the column names 
(these names MUST match the names listed above) and the data must be separated
by tabs. Note that, in contrast to the keyboard input, variables used in the
file input must be in mol/kg.

For example:

flag	PH	ALK	S	T	P
8	8.2	2400e-6	35	20	0
8	8.2	2400e-6	35	20	0


}
}
\keyword{utilities}
