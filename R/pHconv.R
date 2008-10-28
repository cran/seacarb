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


#--------------------------------------------------------------
# Conversion factors for converting pH  
# from total pH scale to free pH scale (pHtotal2free) 
# from total pH scale to seawater pH scale (pHtotal2SWS)
# and from free pH scale to seawater scale (pHfree2SWS)
# pHfree = pHtotal + pHtotal2free
# pHSWS  = pHtotal + pHtotal2SWS
# pHSWS  = pHfree  + pHfree2SWS
#--------------------------------------------------------------


"pHconv" <- function (S=35,T=25,P=0)

{
  cc <- kconv(S,T,P)     # conversion factors for dissociation constants, k 
  
  # conversion factor for pH = log10(1/conversion for k)
  return (list(pHtotal2SWS=log10(1/cc$ktotal2SWS), pHtotal2free=log10(1/cc$ktotal2free),
               pHfree2SWS=log10(1/cc$kfree2SWS)))
}

