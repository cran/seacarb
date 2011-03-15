# Copyright (C) 2008 Jean-Pierre Gattuso and Heloise Lavigne and Aurelien Proye
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

"Kb" <-
function(S=35,T=25,P=0,pHscale="T"){

nK <- max(length(S), length(T), length(P), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}
	
#-------Constantes----------------

#---- issues de equic----
tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
TK = T + tk;           # T [C]; TK [K]

		
	#---------------------------------------------------------------------
	# --------------------- Kb  --------------------------------------------
	#  Kbor = [H+][B(OH)4-]/[B(OH)3]
	#
	#   (Dickson, 1990 in Guide to Best Practices in Ocean CO2 Measurements 2007)
	#   pH-scale: 'total'. mol/kg-soln
	
	
	tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S^(3/2)-0.0996*S*S);
	tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
	tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S)*log(TK);
	
	lnKb = tmp1 / TK + tmp2 + tmp3 + 0.053105*sqrt(S)*TK;
	Kb <- exp(lnKb)

	# ---- Conversion from Total scale to seawater scale before pressure corrections
	factor <- kconv(S=S, T=T, P=rep(0,nK))$ktotal2SWS
	Kb <- Kb * factor

	# ------------------- Pression effect --------------------------------
	Kb <- Pcorrect(Kvalue=Kb, Ktype="Kb", T=T, S=S, P=P, pHscale="SWS")

###----------------pH scale corrections
factor <- rep(NA,nK)
pHsc <- rep(NA,nK)
for(i in (1:nK)){   
 if(pHscale[i]=="T"){factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total ; pHsc[i] <- "total scale"}
 if(pHscale[i]=="F"){factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free ; pHsc[i] <- "free scale"}
 if(pHscale[i]=="SWS"){factor[i] <- 1 ; pHsc[i] <- "seawater scale"}
Kb[i] <- Kb[i]*factor[i]
}

##------------Warnings

for(i in 1:nK){
if((T[i]>45)|(S[i]>45)|(T[i]<0)|(S[i]<5)){warning("S and/or T is outside the range of validity of the formulation available for Kb in seacarb.")}
}


	attr(Kb,"unit")     = "mol/kg-soln"
	attr(Kb,"pH scale") = pHsc
	return(Kb)
}
