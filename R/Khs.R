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

"Khs"              <- function (S=35,T=25,P=0)

#--------------------------------------------------------------
# Dissociation constant of H2S
#--------------------------------------------------------------

{

nK <- max(length(S), length(T), length(P))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}


tk = 273.15           # [K] (for conversion [deg C] <-> [K])
TC = T + tk           # TC [C]; T[K]

	#--------------------------------------------------------------
  # Dissociation constant of H2S on total scale - Millero 1995
 	#--------------------------------------------------------------
lnK <- 225.838 + 0.3449*sqrt(S) - 0.0274*S -13275.3/TC  -34.6435*log(TC)

lnkpok0 <- rep(0, nK)

# P : applied pressure (in Bars) = (Total Pressure-1)
#Pressure Correction from Millero 1995 - Typos corrections from Lewis and Wallace (CO2SYS)
  
	for(i in (1:nK)){		
	if (P[i]> 0.0)
  {
		R = 83.131             # mol bar deg-1
    a0 <- -14.80; a1<- 0.002; a2 <- --0.4e-3
    b0 <- 2.89e-3 ; b1 <- 0.054e-3; b2 <- 0
    deltav  <-   a0 + a1 *T[i] + a2 *T[i]*T[i]
    deltak  <-  (b0  + b1 *T[i] + b2*T[i]*T[i])
    lnkpok0[i] <-  -(deltav/(R*TC[i]))*P[i] + (0.5*deltak /(R*TC[i]))*P[i]*P[i]
  }
}
  Khs      <- exp(lnK +lnkpok0)  # K_H2S on total scale
  attr(Khs,"unit")     = "mol/kg-soln"
  attr(Khs,"pH scale") = "total hydrogen ion concentration"
  return(Khs)

}  ## END Khs
