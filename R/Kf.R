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
"Kf" <-
function(S=35,T=25,P=0,kf='x',pHscale="T"){

nK <- max(length(S), length(T), length(P), length(kf), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(kf)!=nK){kf <- rep(kf[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}
 
##----------Check the validity of the method regarding the T/S range

for(i in 1:nK){
if(kf[i]=='x'){
kf[i] <- 'pf'  ## Perez and Fraga by default
if((T[i]>33)|(T[i]<10)|(S[i]<10)|(S[i]>40)){kf[i] <- 'dg' }
}
}

 
#-------Constantes----------------

#---- issues de equic----
tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
TK = T + tk;           # TC [C]; T[K]
Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
cl3 = Cl^(1/3);   
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl * Cl;   # ionic strength
iom0 = 19.924*S/(1000-1.005*S);
ST = 0.14/96.062/1.80655*S;   # (mol/kg soln) total sulfate
FT = 7e-5*(S/35)    

bor = (416.*(S/35.))* 1e-6;   # (mol/kg), DOE94 

#---------------------------------------------------------------------
#---------------------- Kf Perez and Fraga ---------------------------
#  Kf = [H+][F-]/[HF]  
#
#   Perez and Fraga, 1987 in Guide to the Best Practices for Ocean CO2 Measurements
#   Dickson, Sabine and Christian, 2007, Chapter 5, p. 14)
#  
#   pH-scale: 'total'   

lnKfpf <- 874/TK - 9.68 + 0.111*S^(1/2)
Kfpf <- exp(lnKfpf)

# Conversion to free scale for pressure corrections

Ks = Ks(S=S, T=T, P=rep(0,nK))                 # on free pH scale
	ST  = 0.14/96.062/1.80655*S    # (mol/kg soln) total sulfate
	total2free  = 1/(1+ST/Ks)      # Kfree = Ktotal*total2free
	total2free <- as.numeric(total2free)	

factor <- total2free

Kfpf <- Kfpf * factor

#---------------------------------------------------------------------
# --------------------- Kf Dickson and Goyet -------------------------
#  Kf = [H+][F-]/[HF]  
#
#   (Dickson and Riley, 1979 in Dickson and Goyet, 
#   1994, Chapter 5, p. 14)
#   pH-scale: 'free' (require to convert in total scale after pressure corrections 

lnKfdg = 1590.2/TK - 12.641 + 1.525*sqrt(ION);

Kfdg <- exp(lnKfdg)

# ---------- Choice between methods (Perez and Dickson) ----------

Kf <- Kfpf

for(i in (1:nK)){
if(kf[i]=='pf'){Kf[i] <- Kfpf[i]}
if(kf[i]=='dg'){Kf[i] <- Kfdg[i]}
}

# ------------------- Pression effect --------------------------------
Kf <- Pcorrect(Kvalue=Kf, Ktype="Kf", T=T, S=S, P=P, pHscale="F") 
# All the Kf constants are on free scale, whatever the method or the pressure


###----------------pH scale corrections to retrieve the pHscale given in argument
factor <- rep(NA,nK)
pHsc <- rep(NA,nK)
for(i in (1:nK)){   
 if(pHscale[i]=="T"){factor[i] <- 1+ST[i]/Ks(S=S[i], T=T[i], P=P[i]); pHsc[i] <- "total scale"}
 if(pHscale[i]=="F"){factor[i] <- 1 ; pHsc[i] <- "free scale"}
 if(pHscale[i]=="SWS"){factor[i] <- 1 + ST[i]/Ks(S=S[i], T=T[i], P=P[i])+ FT[i]/Kf[i] ; pHsc[i] <- "seawater scale"}
Kf[i] <- Kf[i]*factor[i]
}


##------------Warnings

for(i in 1:nK){
if((kf[i]=='pf')&((T[i]>33)|(T[i]<9)|(S[i]<10)|(S[i]>40))){warning("S and/or T is outside the range of validity of the formulation chosen for Kf.")}
if((T[i]>45)|(S[i]>45)){warning("S and/or T is outside the range of validity of the formulations available for Kf in seacarb.")}
}

##---------------Attributes
method <- c()
for(i in 1:nK){
m <- "Perez and Fraga (1987)"
if(kf[i]=="dg"){m <- "Dickson and Riley (1979 in Dickson and Goyet, 1994)"}
method <- c(method, m)
}


attr(Kf,"unit")     = "mol/kg-soln"
attr(Kf,"pH scale") = pHsc
attr(Kf, "method") = method
return(Kf)
}





