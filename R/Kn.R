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


"Kn"              <- function (S=35, T=25, P=0)

#--------------------------------------------------------------
# Dissociation constant of ammonium 
#--------------------------------------------------------------

{

tk = 273.15           # [K] (for conversion [deg C] <-> [K])
TC = T + tk           # TC [C]; T[K]

	#--------------------------------------------------------------
  # Dissociation constant of ammonium on seawater scale - Millero 1995
 	#--------------------------------------------------------------
lnK <- -6285.33/TC+0.0001635*TC-0.25444+(0.46532-123.7184/TC)* sqrt(S)+
      (-0.01992+3.17556/TC)*S

lnkpok0 <- 0

# P : applied pressure (in Bars) = (Total Pressure-1)
# Pressure Correction from Millero 1995 - Typos corrections from Lewis and Wallace (CO2SYS)
  if (P> 0.0)
  {
		R = 83.131             # mol bar deg-1
    a0 <- -26.43; a1<- 0.0889; a2 <- -0.000905
    b0 <- -5.03e-3 ; b1 <- 0.0814e-3; b2 <- 0
    deltav  <-   a0 + a1 *T + a2 *T*T
    deltak  <-  (b0  + b1 *T + b2*T*T)
    lnkpok0 <-  -(deltav/(R*TC))*P + (0.5*deltak /(R*TC))*P*P
  }

  Kn      <- exp(lnK +lnkpok0)  # K_NH4 on SEAWATER scale
  cc  <- kconv(S,T,P)            # conversion factors from one scale to other
  Kn  <- Kn/cc$ktotal2SWS        # Kn now on total scale

  attr(Kn,"pH scale") = "total hydrogen ion concentration"
  attr(Kn,"unit")     = "mol/kg-soln"
  return(Kn)

}  ## END Kn
