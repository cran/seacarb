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
"K1" <-
function(S=35,T=25,P=0,k1k2='x',pHscale="T")
{

nK <- max(length(S), length(T), length(P), length(k1k2), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(k1k2)!=nK){k1k2 <- rep(k1k2[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}

##----------Check the validity of the method regarding the T/S range

for(i in 1:nK){
if(k1k2[i]=='x'){
k1k2[i] <- 'l'  ## luecker by default
if((T[i]>35)|(T[i]<2)|(S[i]<19)|(S[i]>43)){k1k2[i] <- 'm10' }
}
}

#-------Constantes----------------

#---- issues de equic----
tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
TK = T + tk;           # TC [C]; T[K]
Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mil)
cl3 = Cl^(1/3);   
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl * Cl;   # ionic strength
iom0 = 19.924*S/(1000-1.005*S);
ST = 0.14/96.062/1.80655*S;   # (mol/kg soln) total sulfate


bor = (416.*(S/35.))* 1e-6;   # (mol/kg), DOE94   

# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#     Mehrbach et al (1973) refit by Lueker et al. (2000).
#
#(Lueker  et al., 2000 in Guide to the Best Practices for Ocean CO2 Measurements
#   Dickson, Sabine and Christian , 2007, Chapter 5, p. 13)
#
#   pH-scale: 'total'. mol/kg-soln
	
logK1lue <- (-3633.86)/TK + 61.2172 - 9.67770*log(TK) + 0.011555*S - 0.0001152*S*S
K1lue <- 10^logK1lue

# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 14)
#   pH-scale: 'total'. mol/kg-soln

tmp1 = 2.83655 - 2307.1266 / TK - 1.5529413 * log(TK);
tmp2 =         - (0.20760841 + 4.0484 / TK) * sqrt(S);
tmp3 =         + 0.08468345 * S - 0.00654208 * S * sqrt(S);   
tmp4 =         + log(1 - 0.001005 * S);

lnK1roy = tmp1 + tmp2 + tmp3 + tmp4;

        K1roy  = exp(lnK1roy);
        
# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#   Millero et al. 2006 Marine Chemistry
#   pH-scale: 'SWS scale'. mol/kg-soln

pK1o <- 6320.813/TK + 19.568224*log(TK) -126.34048
A1 <- 13.4191*S^(0.5) + 0.0331*S - (5.33e-5)*S^2
B1 <- -530.123*S^(0.5) - 6.103*S
C1 <- -2.06950*S^(0.5)

pK1 <- pK1o + A1 + B1/TK + C1*log(TK)

K1mil06 <- 10^(-pK1)

# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#   Millero 2010 Marine and Fresh water research

pK1o <- 6320.813/TK + 19.568224*log(TK) -126.34048

#   pH-scale: 'SWS scale'. mol/kg-soln

A1 <- 13.4038*S^(0.5) + 0.03206*S - (5.242e-5)*S^2
B1 <- -530.659*S^(0.5) - 5.8210*S
C1 <- -2.0664*S^(0.5)

pK1 <- pK1o + A1 + B1/TK + C1*log(TK)

K1mil10_SWS <- 10^(-pK1)

#   pH-scale: 'Total scale'. mol/kg-soln

A1 <- 13.4051*S^(0.5) + 0.03185*S - (5.218e-5)*S^2
B1 <- -531.095*S^(0.5) - 5.7789*S
C1 <- -2.0663*S^(0.5)

pK1 <- pK1o + A1 + B1/TK + C1*log(TK)

K1mil10_total <- 10^(-pK1)

#   pH-scale: 'Free scale'. mol/kg-soln

A1 <- 5.09247*S^(0.5) + 0.05574*S - (9.279e-5)*S^2
B1 <- -189.879*S^(0.5) - 11.3108*S
C1 <- -0.8080*S^(0.5)

pK1 <- pK1o + A1 + B1/TK + C1*log(TK)

K1mil10_free <- 10^(-pK1)


# ---------- Choice between methods (Lueker or Roy) ----------

K1 <- K1lue

for(i in (1:nK)){
if(k1k2[i]=='l'){K1[i] <- K1lue[i] }
if(k1k2[i]=='r'){K1[i] <- K1roy[i] }
if(k1k2[i]=='m06'){K1[i] <- K1mil06[i] }
if(k1k2[i]=='m10'){K1[i] <-  K1mil10_SWS[i]
  if((pHscale[i]=="F")&(P[i]==0)){K1[i] <-  K1mil10_free[i]}
  if((pHscale[i]=="T")&(P[i]==0)){K1[i] <-  K1mil10_total[i]}
  }
}

# ---- Conversion from Total scale to seawater scale before pressure corrections
pHsc <- rep(NA,nK)

for(i in (1:nK)){
if((k1k2[i] %in% c('l', 'r'))&(P[i]>0)){      
factor <- kconv(S=S[i], T=T[i], P=0)$ktotal2SWS
K1[i] <- K1[i] * factor
}}
## K1 is already in the sea water scale with the Millero 2006 formulation
## with the Millero 2010 formulation K1 is on SWS scale if P is greater than 0

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
	
	K1[i] = K1[i]*exp(lnkpok0[1]);

###----------------pH scale corrections  
 if(pHscale[i]=="T"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total; pHsc[i] <- "total scale"}
 if(pHscale[i]=="F"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free ; pHsc[i] <- "free scale"}
 if(pHscale[i]=="SWS"){factor <- 1 ; pHsc[i] <- "seawater scale"}
K1[i] <- K1[i]*factor

}
}

###----------------pH scale corrections in case P=0
for(i in (1:nK)){ 
  if(P[i]==0){ 
    if(k1k2[i] %in% c("l", "r")){  # whith the Luecker 2000 and Roy formulation pHscale is total scale
      if(pHscale[i]=="T"){factor <- 1; pHsc[i] <- "total scale"}
      if(pHscale[i]=="F"){factor <- kconv(S=S[i], T=T[i], P=P[i])$ktotal2free ; pHsc[i] <- "free scale"}
      if(pHscale[i]=="SWS"){factor <- kconv(S=S[i], T=T[i], P=P[i])$ktotal2SWS  ; pHsc[i] <- "seawater scale"}
    }
    if(k1k2[i]=="m06"){ # whith the Millero 2006 formulation pHscale is SWS scale
      if(pHscale[i]=="T"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total; pHsc[i] <- "total scale"}
      if(pHscale[i]=="F"){factor <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free ; pHsc[i] <- "free scale"}
      if(pHscale[i]=="SWS"){factor <- 1 ; pHsc[i] <- "seawater scale"}
    }
    if(k1k2[i]=="m10"){ # whith the Millero 2010 formulation pHscale is already adaptated to the require scale
      if(pHscale[i]=="T"){factor <- 1; pHsc[i] <- "total scale"}
      if(pHscale[i]=="F"){factor <- 1 ; pHsc[i] <- "free scale"}
      if(pHscale[i]=="SWS"){factor <- 1 ; pHsc[i] <- "seawater scale"}
    }
K1[i] <- K1[i]*factor
}
}

##---------------Attributes
method <- c()
for(i in 1:nK){
m <- "Luecker et al. (2000)"
if(k1k2[i]=="m06"){m <- "Millero et al. (2006)"}
if(k1k2[i]=="m10"){m <- "Millero (2010)"}
if(k1k2[i]=="r"){m <- "Roy et al. (1993)"}
method <- c(method, m)
}

##------------Warnings

for(i in 1:nK){
if((k1k2[i]=='l')&((T[i]>35)|(T[i]<2)|(S[i]<19)|(S[i]>43))){warning("S and/or T is outside the range of validity of the formulation chosen for K1.")}
if((k1k2[i]=='r')&((T[i]>45)|(S[i]>45))){warning("S and/or T is outside the range of validity of the formulation chosen for K1.")}
if((T[i]>50)|(S[i]>50)){warning("S and/or T is outside the range of validity of the formulations available for K1 in seacarb.")}
}


attr(K1,"unit")     = "mol/kg-soln"
attr(K1,"pH scale") = pHsc
attr(K1, "method") = method
return(K1)
}
