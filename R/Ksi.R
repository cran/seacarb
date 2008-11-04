# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

"Ksi"              <- function (S=35,T=25,P=0)

#--------------------------------------------------------------
# Dissociation constant for Si(OH)4
#--------------------------------------------------------------

{

nK <- max(length(S), length(T), length(P))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}

	
#-------Constantes----------------


tk = 273.15               # [K] (for conversion [deg C] <-> [K])
TC = T + tk               # TC [C]; T[K]

#--------------------------------------------------------------
# Dissociation constant of Si(OH)4 on total scale - DOE 1994 - correct for mol/(kg-H2O)->mol/kg soln
#--------------------------------------------------------------
Io  <- 19.924*S/(1000-1.005*S)                     # ionic strenght, mol/kg-H2O     # ionic strength
lnK <- 117.385 + 3.5913*sqrt(Io) - 1.5998*Io + 0.07871*Io*Io +
      (-8904.2 - 458.79*sqrt(Io) + 188.74*Io - 12.1652*Io*Io)/TC +
      (-19.334)*log(TC)

lnkpok0 <- rep(0, nK)

# P : applied pressure (in Bars) = (Total Pressure-1)
# Pressure Correction from Millero 1995 - same as for BOH3

for(i in (1:nK)){	
  if (P[i]> 0.0)
  {
		R = 83.131             # mol bar deg-1
		a0 <- -29.48; a1 <- 0.1622; a2 <- 2.6080e-3
    b0 <- -2.84e-3 ; b1 <- 0.0; b2 <- 0.0
    deltav  <-   a0 + a1 *T[i] + a2 *T[i]*T[i]
    deltak  <-  (b0  + b1 *T[i] + b2*T[i]*T[i])
    lnkpok0[i] <-  -(deltav/(R*TC[i]))*P[i] + (0.5*deltak /(R*TC[i]))*P[i]*P[i]
  }
}

  # from molality to molinity
  Ksi      <- exp(lnK +lnkpok0)*(1.0 - 0.001005*S)
  attr(Ksi,"unit")     = "mol/kg-soln"
  attr(Ksi,"pH scale") = "total hydrogen ion concentration"
  return(Ksi)
}  ## END Ksi
