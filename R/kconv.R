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
# Conversion factors for converting dissociation constants 
# from total pH scale to free pH scale (ktotal2free) 
# from total pH scale to seawater pH scale (ktotal2SWS)
# and from free pH scale to seawater scale (kfree2SWS)
# kfree = ktotal * ktotal2free
# kSWS  = ktotal * ktotal2SWS
# kSWS  = kfree  * kfree2SWS
#--------------------------------------------------------------


"kconv" <- function (S=35,T=25,P=0)

{
	#--------------------------------------------------------------
  # CONVERT equilibrium constants to free scale:
	#--------------------------------------------------------------
	#------------------ Ks ----------------------------------------
	#       Dickson and Goyet (1994), Chapter 5, p.13
	#       (required for total2free)
	#       Equilibrium constant for HSO4- = H+ + SO4--
	#
	#       K_S  = [H+]free [SO4--] / [HSO4-]
	#       pH-scale: free scale !!!
	#
	#       the term log(1-0.001005*S) converts from
	#       mol/kg H2O to mol/kg soln
	Ks = Ks(S=S, T=T, P=P)                 # on free pH scale
	ST  = 0.14/96.062/1.80655*S    # (mol/kg soln) total sulfate
	total2free  = 1/(1+ST/Ks)      # Kfree = Ktotal*total2free

	#---------------------------------------------------------------------
	# --------------------- Kf  ------------------------------------------
	#  Kf = [H+][F-]/[HF]
	#
	#   (Dickson and Riley, 1979 in Dickson and Goyet,
	#   1994, Chapter 5, p. 14)
	#   pH-scale: 'total'
	Kf = Kf(S=S, T=T, P=P)
	Kf  = Kf*total2free       # convert Kf from total to free pH scale


	#------- sws2free -----------------------------------------------
	#
	#       convert from pH_sws ('seawater scale`) to pH ('free`):
	#       pH_sws = pH_free - log(1+S_T/K_S(S,T)+F_T/K_F(S,T))
	FT = 7e-5*(S/35)                  # (mol/kg soln) total fluoride
	free2SWS  = 1+ST/Ks+FT/Kf         # Kfree = Ksws*sws2free
	total2SWS = total2free * free2SWS # KSWS = Ktotal*total2SWS

return (list(ktotal2SWS=total2SWS, ktotal2free=total2free,kfree2SWS=free2SWS))
}

