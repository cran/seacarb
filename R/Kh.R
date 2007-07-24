# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
"Kh" <-
function(S=35,T=25,P=0){

	
#-------Constantes----------------

#---- issues de equic----
tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
TC = T + tk;           # TC [C]; T[K]
Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
cl3 = Cl^(1/3);   
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl * Cl;   # ionic strength
iom0 = 19.924*S/(1000-1.005*S);
ST = 0.14/96.062/1.80655*S;   # (mol/kg soln) total sulfate


bor = (416.*(S/35.))* 1e-6;   # (mol/kg), DOE94   

	#---------------------------------------------------------------------
	#---------------------- Kh (K Henry) ---------------------------------
	#
	#               CO2(g) <-> CO2(aq.)
	#               Kh      = [CO2]/ p CO2
	#
	#   Weiss (1974)   [mol/kg/atm]
	#
	#                             
	
	tmp = 9345.17 / TC - 60.2409 + 23.3585 * log(TC/100);
	nKhwe74 = tmp + S*(0.023517-0.00023656*TC+0.0047036e-4*TC*TC);
	

	Kh= exp(nKhwe74);
	
	
	#attr(Kh,"unit") = "mol/kg"	
	print(Kh);
}
