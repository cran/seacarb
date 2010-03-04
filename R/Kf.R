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

for(i in (1:nK)){
if (P[i] > 0.0)
{

	RGAS = 8.314510;        # J mol-1 deg-1 (perfect Gas)  
	R = 83.14472;             # mol bar deg-1 
	                        # conversion cm3 -> m3          *1.e-6
        	                  #            bar -> Pa = N m-2  *1.e+5
	                        #                => *1.e-1 or *1/10
	# index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8,
	#        K1P 9, K2P 10, K3P 11
	
	#----- note: there is an error in Table 9 of Millero, 1995.
	#----- The coefficients -b0 and b1
	#----- have to be multiplied by 1.e-3!

	#----- there are some more errors! 
	#----- the signs (+,-) of coefficients in Millero 95 do not
	#----- agree with Millero 79
	
		
		
	a0 = c(-25.5, -15.82, -29.48, -25.60, -18.03, -9.78, -48.76, -46, -14.51, -23.12, -26.57);
	a1 = c(0.1271, -0.0219, 0.1622, 0.2324, 0.0466, -0.0090, 0.5304, 0.5304, 0.1211, 0.1758, 0.2020);
	a2 = c(0.0, 0.0, 2.608*1e-3, -3.6246*1e-3, 0.316*1e-3, -0.942*1e-3, 0.0, 0.0, -0.321*1e-3, -2.647*1e-3, -3.042*1e-3);
	b0 = c(-3.08*1e-3, 1.13*1e-3, -2.84*1e-3, -5.13*1e-3, -4.53*1e-3, -3.91*1e-3, -11.76*1e-3, -11.76*1e-3, -2.67*1e-3, -5.15*1e-3, -4.08*1e-3);
	b1 = c(0.0877*1e-3, -0.1475*1e-3, 0.0, 0.0794*1e-3, 0.09*1e-3, 0.054*1e-3, 0.3692*1e-3, 0.3692*1e-3, 0.0427*1e-3, 0.09*1e-3, 0.0714*1e-3);
	b2 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	
	deltav = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	deltak = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	lnkpok0 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	
	for (ipc in 1:length(a0))
	{
	  deltav[ipc]  =  a0[ipc] + a1[ipc] *T[i] + a2[ipc] *T[i]*T[i];
	  deltak[ipc]   = (b0[ipc]  + b1[ipc] *T[i] + b2[ipc] *T[i]*T[i]);  
	  lnkpok0[ipc]  = -(deltav[ipc] /(R*TK[i]))*P[i] + (0.5*deltak[ipc] /(R*TK[i]))*P[i]*P[i];
	}

	Kf[i] <- Kf[i]*exp(lnkpok0[6])
}
}


## pH scale correction in case of Dickson and Goyet method
## Passage from free scale to total scale
Ks <- Ks(S=S, T=T, P=P)
for(i in (1:nK)){
if(kf[i]=='dg'){
tmp2 = (1-0.001005*S[i])*(1+ST[i]/Ks[i]); #passage from "free scale" to "total scale"
Kf[i] = Kf[i]*tmp2;
}
}
 
###----------------pH scale corrections
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
if((kf[i]=='pf')&((T[i]>33)|(T[i]<9)|(S[i]<10)|(S[i]>40))){warning("The method used to compute Kf does not fit with the T/S values.")}
if((T[i]>45)|(S[i]>45)){warning("The calculation methods proposed in seacarb are not efficient to compute K1 if T exceeds 45oC and/or S exceeds 45")}
}

##---------------Attributes
method <- c()
for(i in 1:nK){
m <- "Perez and Fraga (1987)"
if(kf=="dg"){m <- "Dickson and Goyet (1979)"}
method <- c(method, m)
}


attr(Kf,"unit")     = "mol/kg-soln"
attr(Kf,"pH scale") = pHsc
attr(Kf, "method") = method
return(Kf)
}





